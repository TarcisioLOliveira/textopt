/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of textopt.
 *
 *   textopt is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   textopt is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with SolidPrep.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <SFML/Config.hpp>
#include <SFML/Graphics/Color.hpp>
#include <SFML/Graphics/Sprite.hpp>
#include <SFML/Graphics/Texture.hpp>
#include <SFML/System/Vector2.hpp>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <iostream>
#include <numeric>
#include <array>
#include <algorithm>

#include "MMASolver.hpp"
#include "param.hpp"
#include "smooth.hpp"
#include "util.hpp"
#include "render.hpp"
#include "opt_func.hpp"

void texture_map(std::vector<double>& map_z, const std::vector<double>& orig_z, double f, double ap, double vc){
    using namespace param;

    vc *= dim_scale;
    double y1 = -std::sqrt((std::pow(std::tan(alpha1)*r, 2))/(std::pow(std::tan(alpha1), 2)+1));
    double y2 =  std::sqrt((std::pow(std::tan(alpha2)*r, 2))/(std::pow(std::tan(alpha2), 2)+1));

    double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + std::tan(alpha1)*(y1+f/2);//-ap + std::tan(alpha1)*f/2;
    double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - std::tan(alpha2)*(y2+f/2); //-ap - std::tan(alpha2)*f/2;

    // Calculate whether the angles of the cutting tool intercept the surface
    // first or are stopped by the feed rate.
    //
    // Check the intersection with the feed rate edges. If they're greater than
    // zero, the cut line is narrower than the feed line, so max height is
    // zero. Otherwise, it's the height at the intersection
    double intersec1 = 0;
    double intersec2 = 0;
    if(y1 + f/2 > 0){
        // If the radius is smaller than the feed rate, check for the line
        intersec1 = smooth::min({0, line_root1});
    } else {
        // Check for the circle
        intersec1 = smooth::min({0, -std::sqrt(r*r - f*f/4) + r - ap});
    }
    if(y2 + f/2 < f){
        // If the radius is smaller than the feed rate...
        intersec2 = smooth::min({0, line_root2 + std::tan(alpha2)*f});
    } else {
        intersec2 = smooth::min({0, -std::sqrt(r*r - f*f/4) + r - ap});
    }
    // Check if one of them is greater than zero
    double max_intersec = -smooth::min({-intersec1, -intersec2});
    // If one is greater than zero, max_z must be zero. Otherwise, it's
    // the lesser one.
    max_z = smooth::min({0, max_intersec});
    min_z = 0;

    double delta_uet = vc/f_uet;
    std::vector<double> newz(tex_width);
    for(size_t Y = 0; Y < tex_height; ++Y){
        double y = static_cast<double>(Y);

        // Calculate current row
        // Apply `abs()` to prevent it from becoming less than zero
        // when close to zero
        double mult = smooth::abs(smooth::floor((y + YOFF) / f));

        // Phase differences
        double perimeter = 2*M_PI*cylinder_radius;

        double phase_diff_x = 2*M_PI*fx*mult*perimeter*dimx/vc + phix;
        double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
        double phix_cur = phase_diff_x - smooth::floor((phase_diff_x)/(2*M_PI))*2*M_PI;
        double phiz_cur = phase_diff_z - smooth::floor((phase_diff_z)/(2*M_PI))*2*M_PI;

        double xoffset_uet = mult*perimeter;

        if(y >= line_root1/std::tan(alpha1) + mult*f && y <= -line_root2/std::tan(alpha2) + mult*f){
            for(size_t X = 0; X < tex_width; ++X){
                double x = static_cast<double>(X);
                // Random oscillation
                double newx = x + Ax*std::sin(2*M_PI*fx*x*dimx/vc + phix_cur);
                double oscillation = Az*std::sin(2*M_PI*fz*newx*dimx/vc + phiz_cur);

                // Ultrasonic turning effects
                double xcirc = x + xoffset_uet;
                double mult_uet = smooth::abs(smooth::floor(xcirc / delta_uet));
                double x_uet = (xcirc - mult_uet*delta_uet)*Ax_uet/delta_uet - Ax_uet/2;
                double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)));

                newz[X] = oscillation + uet_effect;
            }

                // Tool shape
            if(y <= y1 + mult*f + f/2){
                for(size_t X = 0; X < tex_width; ++X){
                    // Get pixels
                    double& z = map_z[X + Y*tex_width];
                    double oz = orig_z[X + Y*tex_width];
                    newz[X] += -std::tan(alpha1)*(y - mult*f) + line_root1;

                    // Write results
                    z = smooth::min({oz, newz[X]});
                }
            } else if(y <= y2 + mult*f + f/2){
                for(size_t X = 0; X < tex_width; ++X){
                    // Get pixels
                    double& z = map_z[X + Y*tex_width];
                    double oz = orig_z[X + Y*tex_width];
                    newz[X] += -std::sqrt(r*r - std::pow(y - mult*f - f/2, 2)) + r - ap;

                    // Write results
                    z = smooth::min({oz, newz[X]});
                }
            } else {
                for(size_t X = 0; X < tex_width; ++X){
                    // Get pixels
                    double& z = map_z[X + Y*tex_width];
                    double oz = orig_z[X + Y*tex_width];
                    newz[X] += std::tan(alpha2)*(y - mult*f) + line_root2;

                    // Write results
                    z = smooth::min({oz, newz[X]});
                }
            }
        } else {
            for(size_t X = 0; X < tex_width; ++X){
                // Get pixels
                double& z = map_z[X + Y*tex_width];
                double oz = orig_z[X + Y*tex_width];
                z = oz;
            }
        }
    }
    min_z = *std::min_element(map_z.begin(), map_z.end());
}

void dzdvc(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdvc){
    using namespace param;

    vc *= dim_scale;
    double y1 = -std::sqrt((std::pow(std::tan(alpha1)*r, 2))/(std::pow(std::tan(alpha1), 2)+1));
    double y2 =  std::sqrt((std::pow(std::tan(alpha2)*r, 2))/(std::pow(std::tan(alpha2), 2)+1));

    double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + std::tan(alpha1)*(y1+f/2);//-ap + std::tan(alpha1)*f/2;
    double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - std::tan(alpha2)*(y2+f/2); //-ap - std::tan(alpha2)*f/2;
                                                                                      //
    double dlrdvc1 = 0;
    double dlrdvc2 = 0;

    double intersec1 = 0;
    double intersec2 = 0;
    // always zero
    double di1 = 0;
    double di2 = 0;
    if(y1 + f/2 > 0){
        intersec1 = smooth::min({0, line_root1});
    } else {
        intersec1 = smooth::min({0, -std::sqrt(r*r - f*f/4) + r - ap});
    }
    if(y2 + f/2 < f){
        intersec2 = smooth::min({0, line_root2 + std::tan(alpha2)*f});
    } else {
        intersec2 = smooth::min({0, -std::sqrt(r*r - f*f/4) + r - ap});
    }
    double max_intersec = -smooth::min({-intersec1, -intersec2});

    double dmi = 0;//-smooth::min_deriv({-intersec1, -intersec2}, -intersec1)*(-di1) - smooth::min_deriv({-intersec1, -intersec2}, -intersec2)*(-di2);
    dmax_zdvc = 0;//smooth::min_deriv({0, max_intersec}, max_intersec)*dmi;

    double delta_uet = vc/f_uet;
    double dduetdvc = 1.0/f_uet;
    std::vector<double> newz(tex_width);
    std::vector<double> dznewdvc(tex_width);
    for(size_t Y = 0; Y < tex_height; ++Y){
        double y = static_cast<double>(Y);

        double mult = smooth::abs(smooth::floor((double) (y + YOFF) / f));
        double dmult = 0;

        double perimeter = 2*M_PI*cylinder_radius;

        double phase_diff_x = 2*M_PI*fx*mult*perimeter*dimx/vc + phix;
        double dpdx = -2*M_PI*fx*mult*perimeter*dimx/(vc*vc);
        double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
        double dpdz = -2*M_PI*fz*mult*perimeter*dimx/(vc*vc);
        double phix_cur = phase_diff_x - smooth::floor((phase_diff_x)/(2*M_PI))*2*M_PI;
        double dphix_cur = dpdx - smooth::floor_deriv((phase_diff_x)/(2*M_PI))*dpdx;
        double phiz_cur = phase_diff_z - smooth::floor((phase_diff_z)/(2*M_PI))*2*M_PI;
        double dphiz_cur = dpdz - smooth::floor_deriv((phase_diff_z)/(2*M_PI))*dpdz;

        double xoffset_uet = mult*perimeter;
        double dxoff_uet = dmult*perimeter;

        if(y >= line_root1/std::tan(alpha1) + mult*f && y <= -line_root2/std::tan(alpha2) + mult*f){
            for(size_t X = 0; X < tex_width; ++X){
                double x = static_cast<double>(X);

                double newx = x + Ax*std::sin(2*M_PI*fx*x*dimx/vc + phix_cur);
                double dnewxdvc = Ax*std::cos(2*M_PI*fx*x*dimx/vc + phix_cur)*((-2)*M_PI*fx*x*dimx/(vc*vc) + dphix_cur);

                double oscillation = Az*std::sin(2*M_PI*fz*newx*dimx/vc + phiz_cur);
                double doscilldvc = Az*std::cos(2*M_PI*fz*newx*dimx/vc + phiz_cur)*(2*M_PI*fz*(dnewxdvc*dimx*vc - newx*dimx)/(vc*vc) + dphiz_cur);

                double xcirc = x + xoffset_uet;
                double dxcirc = dxoff_uet;

                double mult_uet = smooth::abs(smooth::floor(xcirc / delta_uet));
                double dmult_uet = smooth::abs_deriv(smooth::floor(xcirc / delta_uet))*smooth::floor_deriv(xcirc / delta_uet)*(dxcirc*delta_uet - dduetdvc*xcirc)/(delta_uet*delta_uet);

                double x_uet = (xcirc - mult_uet*delta_uet)*Ax_uet/delta_uet - Ax_uet/2;
                double dx_uet = Ax_uet*(dxcirc*delta_uet - dduetdvc*xcirc)/(delta_uet*delta_uet) - Ax_uet*dmult_uet;

                double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)));
                double duet = -0.5*(Az_uet/std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)))*(-2)*(x_uet/Ax_uet)*(dx_uet/Ax_uet);

                newz[X] = oscillation + uet_effect;
                dznewdvc[X] = doscilldvc + duet;
            }
            if(y <= y1 + mult*f + f/2){
                for(size_t X = 0; X < tex_width; ++X){
                    double oz = orig_z[X + Y*tex_width];
                    double& dz = dzdvc[X + Y*tex_width];
                    newz[X] += -std::tan(alpha1)*(y - mult*f) + line_root1;
                    dznewdvc[X] += dlrdvc1;
                    dz = smooth::min_deriv({oz, newz[X]}, 1)*dznewdvc[X];
                }
            } else if(y <= y2 + mult*f + f/2){
                for(size_t X = 0; X < tex_width; ++X){
                    double oz = orig_z[X + Y*tex_width];
                    double& dz = dzdvc[X + Y*tex_width];
                    newz[X] += -std::sqrt(r*r - std::pow(y - mult*f - f/2, 2)) + r - ap;
                    dznewdvc[X] += 0;
                    dz = smooth::min_deriv({oz, newz[X]}, 1)*dznewdvc[X];
                }
            } else {
                for(size_t X = 0; X < tex_width; ++X){
                    double oz = orig_z[X + Y*tex_width];
                    double& dz = dzdvc[X + Y*tex_width];
                    newz[X] +=  std::tan(alpha2)*(y - mult*f) + line_root2;
                    dznewdvc[X] += dlrdvc2;
                    dz = smooth::min_deriv({oz, newz[X]}, 1)*dznewdvc[X];
                }
            }
        } else {
            for(size_t X = 0; X < tex_width; ++X){
                double& dz = dzdvc[X + Y*tex_width];
                dz = 0;
            }
        }
    }
}

void dzdap(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdap){
    using namespace param;

    vc *= dim_scale;
    double y1 = -std::sqrt((std::pow(std::tan(alpha1)*r, 2))/(std::pow(std::tan(alpha1), 2)+1));
    double y2 =  std::sqrt((std::pow(std::tan(alpha2)*r, 2))/(std::pow(std::tan(alpha2), 2)+1));

    double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + std::tan(alpha1)*(y1+f/2);//-ap + std::tan(alpha1)*f/2;
    double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - std::tan(alpha2)*(y2+f/2); //-ap - std::tan(alpha2)*f/2;

    double dlrdap1 = -1;
    double dlrdap2 = -1;

    double intersec1 = 0;
    double intersec2 = 0;
    double di1 = 0;
    double di2 = 0;
    if(y1 + f/2 > 0){
        intersec1 = smooth::min({0, line_root1});
        di1 = smooth::min_deriv({0, line_root1}, 1)*dlrdap1;
    } else {
        intersec1 = smooth::min({0, -std::sqrt(r*r - f*f/4) + r - ap});
        di1 = smooth::min_deriv({0, -std::sqrt(r*r - f*f/4) + r - ap}, 1)*(-1);
    }
    if(y2 + f/2 < f){
        intersec2 = smooth::min({0, line_root2 + std::tan(alpha2)*f});
        di2 = smooth::min_deriv({0, line_root2 + std::tan(alpha2)*f}, 1)*dlrdap2;
    } else {
        intersec2 = smooth::min({0, -std::sqrt(r*r - f*f/4) + r - ap});
        di2 = smooth::min_deriv({0, -std::sqrt(r*r - f*f/4) + r - ap}, 1)*(-1);
    }
    double max_intersec = -smooth::min({-intersec1, -intersec2});

    double dmi = -smooth::min_deriv({-intersec1, -intersec2}, 0)*(-di1) - smooth::min_deriv({-intersec1, -intersec2}, 1)*(-di2);
    dmax_zdap = smooth::min_deriv({0, max_intersec}, 1)*dmi;

    double delta_uet = vc/f_uet;
    std::vector<double> newz(tex_width);
    for(size_t Y = 0; Y < tex_height; ++Y){
        double y = static_cast<double>(Y);

        double mult = smooth::abs(smooth::floor((double) (y + YOFF) / f));

        double perimeter = 2*M_PI*cylinder_radius;

        double phase_diff_x = 2*M_PI*fx*mult*perimeter*dimx/vc + phix;
        double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
        double phix_cur = phase_diff_x - smooth::floor((phase_diff_x)/(2*M_PI))*2*M_PI;
        double phiz_cur = phase_diff_z - smooth::floor((phase_diff_z)/(2*M_PI))*2*M_PI;

        double xoffset_uet = mult*perimeter;

        if(y >= line_root1/std::tan(alpha1) + mult*f && y <= -line_root2/std::tan(alpha2) + mult*f){
            for(size_t X = 0; X < tex_width; ++X){
                double x = static_cast<double>(X);

                double newx = x + Ax*std::sin(2*M_PI*fx*x*dimx/vc + phix_cur);
                double oscillation = Az*std::sin(2*M_PI*fz*newx*dimx/vc + phiz_cur);

                double xcirc = x + xoffset_uet;
                double mult_uet = smooth::abs(smooth::floor(xcirc / delta_uet));
                double x_uet = (xcirc - mult_uet*delta_uet)*Ax_uet/delta_uet - Ax_uet/2;
                double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)));

                newz[X] = oscillation + uet_effect;
            }
            if(y <= y1 + mult*f + f/2){
                for(size_t X = 0; X < tex_width; ++X){
                    double oz = orig_z[X + Y*tex_width];
                    double& dz = dzdap[X + Y*tex_width];
                    newz[X] += -std::tan(alpha1)*(y - mult*f) + line_root1;
                    double dznewdap = dlrdap1;
                    dz = smooth::min_deriv({oz, newz[X]}, 1)*dznewdap;
                }
            } else if(y <= y2 + mult*f + f/2){
                for(size_t X = 0; X < tex_width; ++X){
                    double oz = orig_z[X + Y*tex_width];
                    double& dz = dzdap[X + Y*tex_width];
                    newz[X] += -std::sqrt(r*r - std::pow((double)y - mult*f - f/2, 2)) + r - ap;
                    double dznewdap = -1;
                    dz = smooth::min_deriv({oz, newz[X]}, 1)*dznewdap;
                }
            } else {
                for(size_t X = 0; X < tex_width; ++X){
                    double oz = orig_z[X + Y*tex_width];
                    double& dz = dzdap[X + Y*tex_width];
                    newz[X] +=  std::tan(alpha2)*(y - mult*f) + line_root2;
                    double dznewdap = dlrdap2;
                    dz = smooth::min_deriv({oz, newz[X]}, 1)*dznewdap;
                }
            }
        } else {
            for(size_t X = 0; X < tex_width; ++X){
                double& dz = dzdap[X + Y*tex_width];
                dz = 0;
            }
        }
    }
}

void dzdf(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdf){
    using namespace param;

    vc *= dim_scale;
    double y1 = -std::sqrt((std::pow(std::tan(alpha1)*r, 2))/(std::pow(std::tan(alpha1), 2)+1));
    double y2 =  std::sqrt((std::pow(std::tan(alpha2)*r, 2))/(std::pow(std::tan(alpha2), 2)+1));

    double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + std::tan(alpha1)*(y1+f/2);//-ap + std::tan(alpha1)*f/2;
    double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - std::tan(alpha2)*(y2+f/2); //-ap - std::tan(alpha2)*f/2;

    double dlrdf1 =  std::tan(alpha1)/2;
    double dlrdf2 = -std::tan(alpha2)/2;

    double intersec1 = 0;
    double intersec2 = 0;
    double di1 = 0;
    double di2 = 0;
    if(y1 + f/2 > 0){
        intersec1 = smooth::min({0, line_root1});
        di1 = smooth::min_deriv({0, line_root1}, 1)*dlrdf1;
    } else {
        intersec1 = smooth::min({0, -std::sqrt(r*r - f*f/4) + r - ap});
        di1 = smooth::min_deriv({0, -std::sqrt(r*r - f*f/4) + r - ap}, 1)*f/(4*std::sqrt(r*r-f*f/4));
    }
    if(y2 + f/2 < f){
        intersec2 = smooth::min({0, line_root2 + std::tan(alpha2)*f});
        di2 = smooth::min_deriv({0, line_root2 + std::tan(alpha2)*f}, 1)*(dlrdf2 + std::tan(alpha2));
    } else {
        intersec2 = smooth::min({0, -std::sqrt(r*r - f*f/4) + r - ap});
        di2 = smooth::min_deriv({0, -std::sqrt(r*r - f*f/4) + r - ap}, 1)*f/(4*std::sqrt(r*r-f*f/4));
    }
    double max_intersec = -smooth::min({-intersec1, -intersec2});

    double dmi = -smooth::min_deriv({-intersec1, -intersec2}, 0)*(-di1) - smooth::min_deriv({-intersec1, -intersec2}, 1)*(-di2);
    dmax_zdf = smooth::min_deriv({0, max_intersec}, 1)*dmi;

    double delta_uet = vc/f_uet;
    double dduet = 0;
    std::vector<double> newz(tex_width);
    std::vector<double> dznewdf(tex_width);
    for(size_t Y = 0; Y < tex_height; ++Y){
        double y = static_cast<double>(Y);

        double mult = smooth::abs(smooth::floor((double) (y + YOFF) / f));
        double dmult = smooth::abs_deriv(smooth::floor((double) (y + YOFF) / f))*smooth::floor_deriv((double)(y + YOFF) / f)*(-((double)y / (f*f)));

        double perimeter = 2*M_PI*cylinder_radius;

        double phase_diff_x = 2*M_PI*fx*mult*perimeter*dimx/vc + phix;
        double dpdx = 2*M_PI*fx*dmult*perimeter*dimx/vc;
        double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
        double dpdz = 2*M_PI*fz*dmult*perimeter*dimx/vc;
        double phix_cur = phase_diff_x - smooth::floor((phase_diff_x)/(2*M_PI))*2*M_PI;
        double dphix_cur = dpdx - smooth::floor_deriv((phase_diff_x)/(2*M_PI))*dpdx;
        double phiz_cur = phase_diff_z - smooth::floor((phase_diff_z)/(2*M_PI))*2*M_PI;
        double dphiz_cur = dpdz - smooth::floor_deriv((phase_diff_z)/(2*M_PI))*dpdz;

        double xoffset_uet = mult*perimeter;
        double dxoff_uet = dmult*perimeter;

        if(y >= line_root1/std::tan(alpha1) + mult*f && y <= -line_root2/std::tan(alpha2) + mult*f){
            for(size_t X = 0; X < tex_width; ++X){
                double x = static_cast<double>(X);

                double newx = x + Ax*std::sin(2*M_PI*fx*x*dimx/vc + phix_cur);
                double dnewxdf = Ax*std::cos(2*M_PI*fx*x*dimx/vc + phix_cur)*dphix_cur;

                double oscillation = Az*std::sin(2*M_PI*fz*newx*dimx/vc + phiz_cur);
                double doscilldf = Az*std::cos(2*M_PI*fz*newx*dimx/vc + phiz_cur)*(2*M_PI*fz*dnewxdf*dimx/vc + dphiz_cur);

                double xcirc = x + xoffset_uet;
                double dxcirc = dxoff_uet;

                double mult_uet = smooth::abs(smooth::floor(xcirc / delta_uet));
                double dmult_uet = smooth::abs_deriv(smooth::floor(xcirc / delta_uet))*smooth::floor_deriv(xcirc / delta_uet)*dxcirc;

                double x_uet = (xcirc - mult_uet*delta_uet)*Ax_uet/delta_uet - Ax_uet/2;
                double dx_uet = Ax_uet*dxcirc/delta_uet - Ax_uet*dmult_uet;

                double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)));
                double duet = -0.5*(Az_uet/std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)))*(-2)*(x_uet/Ax_uet)*(dx_uet/Ax_uet);

                newz[X] = oscillation + uet_effect;
                dznewdf[X] = doscilldf + duet;
            }

            if(y <= y1 + mult*f + f/2){
                for(size_t X = 0; X < tex_width; ++X){
                    double oz = orig_z[X + Y*tex_width];
                    double& dz = dzdf[X + Y*tex_width];
                    newz[X] += -std::tan(alpha1)*(y - mult*f) + line_root1;
                    dznewdf[X] += std::tan(alpha1)*(mult + dmult*f) + dlrdf1;
                    dz = smooth::min_deriv({oz, newz[X]}, 1)*dznewdf[X];
                }
            } else if(y <= y2 + mult*f + f/2){
                for(size_t X = 0; X < tex_width; ++X){
                    double oz = orig_z[X + Y*tex_width];
                    double& dz = dzdf[X + Y*tex_width];
                    double yy = y - mult*f - f/2;
                    newz[X] += -std::sqrt(r*r - yy*yy) + r - ap;
                    dznewdf[X] += 2*yy*(mult + dmult*f + 0.5)/std::sqrt(r*r - yy*yy);
                    dz = smooth::min_deriv({oz, newz[X]}, 1)*dznewdf[X];
                }
            } else {
                for(size_t X = 0; X < tex_width; ++X){
                    double oz = orig_z[X + Y*tex_width];
                    double& dz = dzdf[X + Y*tex_width];
                    newz[X] += std::tan(alpha2)*(y - mult*f) + line_root2;
                    dznewdf[X] += -std::tan(alpha2)*(mult + dmult*f) + dlrdf2;
                    dz = smooth::min_deriv({oz, newz[X]}, 1)*dznewdf[X];
                }
            }
        } else {
            for(size_t X = 0; X < tex_width; ++X){
                double& dz = dzdf[X + Y*tex_width];
                dz = 0;
            }
        }
    }
}



int main(int argc, char* argv[]){
    using namespace param;

    size_t window_width = 600;
    size_t window_height = 600;

    double max_roughness = 60;

    auto resolution = sf::VideoMode::getDesktopMode();

    sf::RenderWindow window(sf::VideoMode(window_width, window_height), "textopt");
    window.setPosition(sf::Vector2i((resolution.width - window_width)/2, (resolution.height - window_height)/2));
    sf::Texture img;
    img.create(tex_width, tex_height);
    sf::Sprite sprite(img);
    sprite.setPosition(sf::Vector2f(window_width/2-tex_width/2, window_height/2-tex_height/2));

    std::vector<sf::Uint8> px(tex_width*tex_height*4);

    std::vector<double> map_z  (tex_width*tex_height);
    std::vector<double> orig_z (tex_width*tex_height);
    std::vector<double> df     (tex_width*tex_height);
    std::vector<double> dap    (tex_width*tex_height);
    std::vector<double> dvc    (tex_width*tex_height);

    double ch = 1;
    size_t it = 1;
    double old_surarea = 1;

    const size_t N = 3;
    double x[N] = {20, 40, 70};
    double xmax[N] = {100, 100, 100};
    double xmin[N] = {0.001, 0.001, 0.001};
    double dSa_vec[N] = {0, 0, 0};
    double dsurarea_vec[N] = {0, 0, 0};

    MMASolver mma(N, 1, 0, 1e6, 1);
    mma.SetAsymptotes(0.001, 0.7, 1.2);

    while(window.isOpen()){
        sf::Event event;
        while (window.pollEvent(event)){
            if (event.type == sf::Event::Closed){
                window.close();
            }
            if(event.type == sf::Event::Resized){
                sf::FloatRect view(0, 0, event.size.width, event.size.height);
                window.setView(sf::View(view));
            }
        }
        window.clear(sf::Color::Black);

        if(ch > 1e-8){
            double f = x[0];
            double ap = x[1];
            double vc = x[2];
            texture_map(map_z, orig_z, f, ap, vc);
            render::draw_texture(px, map_z, ap, tex_width, tex_height, render::Colorscheme::HSV);
            double surarea = -opt::surface_area(map_z);
            double roughness = opt::Sa(map_z) - max_roughness;

            dzdf(orig_z, f, ap, vc, df);
            dzdap(orig_z, f, ap, vc, dap);
            dzdvc(orig_z, f, ap, vc, dvc);

            double dsurareadf = -opt::surface_area_dz(map_z, df);
            double dSadf = opt::dSa(df, dmax_zdf);
            double dsurareadap = -opt::surface_area_dz(map_z, dap);
            double dSadap = opt::dSa(dap, dmax_zdap);
            double dsurareadvc = -opt::surface_area_dz(map_z, dvc);
            double dSadvc = opt::dSa(dvc, dmax_zdvc);

            std::cout << std::endl;
            std::cout << dsurareadf << " " << dSadf << std::endl;
            std::cout << dsurareadap << " " << dSadap << std::endl;
            std::cout << dsurareadvc << " " << dSadvc << std::endl;

            dsurarea_vec[0] = dsurareadf;
            dSa_vec[0] = dSadf;
            dsurarea_vec[1] = dsurareadap;
            dSa_vec[1] = dSadap;
            dsurarea_vec[2] = dsurareadvc;
            dSa_vec[2] = dSadvc;

            mma.Update(x, dsurarea_vec, &roughness, dSa_vec, xmin, xmax); 

            ch = std::abs(1 - surarea/old_surarea);
            old_surarea = surarea;

            std::cout << std::endl;
            std::cout << "Iteration: " << it << std::endl;
            std::cout << "Surface area: " << -surarea << std::endl;
            std::cout << "Roughness: " << roughness + max_roughness << std::endl;
            std::cout << "f: " << f << std::endl;
            std::cout << "ap: " << ap << std::endl;
            std::cout << "vc: " << vc << std::endl;
            std::cout << "Change: " << ch << std::endl;
            ++it;
        }

        img.update(px.data());
        auto wsize = window.getSize();
        window_width = wsize.x;
        window_height = wsize.y;
        sprite.setPosition(sf::Vector2f(window_width/2-tex_width/2, window_height/2-tex_height/2));

        window.draw(sprite);
        window.display();
    }

    img.copyToImage().saveToFile("result.png");

    return 0;
}
