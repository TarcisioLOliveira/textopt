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
 *   along with textopt.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <algorithm>

#include "texture.hpp"
#include "param.hpp"
#include "smooth.hpp"
#include <iostream>

namespace texture{

void map(std::vector<double>& map_z, const std::vector<double>& orig_z, double f, double ap, double vc){
    using namespace param;

    vc *= dim_scale;

    // Line parameters for tool slopes.
    const double y1 = -std::sqrt((std::pow(std::tan(alpha1)*r, 2))/(std::pow(std::tan(alpha1), 2)+1));
    const double y2 =  std::sqrt((std::pow(std::tan(alpha2)*r, 2))/(std::pow(std::tan(alpha2), 2)+1));

    const double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + std::tan(alpha1)*(y1+f/2);//-ap + std::tan(alpha1)*f/2;
    const double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - std::tan(alpha2)*(y2+f/2); //-ap - std::tan(alpha2)*f/2;

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
    const double max_intersec = -smooth::min({-intersec1, -intersec2});
    // If one is greater than zero, max_z must be zero. Otherwise, it's
    // the lesser one.
    max_z = smooth::min({0, max_intersec}) + Az;// + Az_uet;
    min_z = -(ap + Az);

    const double delta_uet = vc/f_uet;
    std::vector<double> newz(tex_width);
    for(size_t Y = 0; Y < tex_height; ++Y){
        const double y = static_cast<double>(Y);

        // Calculate current row
        // Apply `abs()` to prevent it from becoming less than zero
        // when close to zero
        const double mult = smooth::abs(smooth::floor((y + YOFF) / f));

        // Phase differences
        const double perimeter = 2*M_PI*cylinder_radius;

        const double phase_diff_x = 2*M_PI*fx*mult*perimeter*dimx/vc + phix;
        const double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
        const double phix_cur = phase_diff_x - smooth::floor((phase_diff_x)/(2*M_PI))*2*M_PI;
        const double phiz_cur = phase_diff_z - smooth::floor((phase_diff_z)/(2*M_PI))*2*M_PI;

        const double xoffset_uet = mult*perimeter;

        for(size_t X = 0; X < tex_width; ++X){
            const double x = static_cast<double>(X);
            // Random oscillation
            const double newx = x + Ax*std::sin(2*M_PI*fx*x*dimx/vc + phix_cur);
            const double oscillation = Az*std::sin(2*M_PI*fz*newx*dimx/vc + phiz_cur);

            // Ultrasonic turning effects
            const double xcirc = x + xoffset_uet;
            const double mult_uet = smooth::abs(smooth::floor(xcirc / delta_uet));
            // const double x_uet = (xcirc - mult_uet*delta_uet)*Ax_uet/delta_uet - Ax_uet/2;
            const double x_uet = Ax_uet*(xcirc/delta_uet - mult_uet - 0.5);
            const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)));

            newz[X] = oscillation + uet_effect;
        }

        // Tool shape
        double shape_z;
        if(y <= y1 + mult*f + f/2){
            shape_z = -std::tan(alpha1)*(y - mult*f) + line_root1;
        } else if(y <= y2 + mult*f + f/2){
            shape_z = -std::sqrt(r*r - std::pow(y - mult*f - f/2, 2)) + r - ap;
        } else {
            shape_z = std::tan(alpha2)*(y - mult*f) + line_root2;
        }

        for(size_t X = 0; X < tex_width; ++X){
            // Get pixels
            double& z = map_z[X + Y*tex_width];
            const double oz = orig_z[X + Y*tex_width];

            // Calculate final depth
            const double end_z = newz[X] + shape_z;

            // Write results
            z = smooth::min({oz, end_z});
        }
        // } else {
        //     for(size_t X = 0; X < tex_width; ++X){
        //         // Get pixels
        //         double& z = map_z[X + Y*tex_width];
        //         const double oz = orig_z[X + Y*tex_width];
        //         z = oz;
        //     }
        // }
    }
}

void dzdvc(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdvc){
    using namespace param;

    vc *= dim_scale;
    const double y1 = -std::sqrt((std::pow(std::tan(alpha1)*r, 2))/(std::pow(std::tan(alpha1), 2)+1));
    const double y2 =  std::sqrt((std::pow(std::tan(alpha2)*r, 2))/(std::pow(std::tan(alpha2), 2)+1));

    const double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + std::tan(alpha1)*(y1+f/2);//-ap + std::tan(alpha1)*f/2;
    const double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - std::tan(alpha2)*(y2+f/2); //-ap - std::tan(alpha2)*f/2;
    
    // const double dlrdvc1 = 0;
    // const double dlrdvc2 = 0;

    // double intersec1 = 0;
    // double intersec2 = 0;
    // // always zero
    // double di1 = 0;
    // double di2 = 0;
    // if(y1 + f/2 > 0){
    //     intersec1 = smooth::min({0, line_root1});
    // } else {
    //     intersec1 = smooth::min({0, -std::sqrt(r*r - f*f/4) + r - ap});
    // }
    // if(y2 + f/2 < f){
    //     intersec2 = smooth::min({0, line_root2 + std::tan(alpha2)*f});
    // } else {
    //     intersec2 = smooth::min({0, -std::sqrt(r*r - f*f/4) + r - ap});
    // }
    // double max_intersec = -smooth::min({-intersec1, -intersec2});

    // double dmi = 0;//-smooth::min_deriv({-intersec1, -intersec2}, -intersec1)*(-di1) - smooth::min_deriv({-intersec1, -intersec2}, -intersec2)*(-di2);
    // dmax_zdvc = 0;//smooth::min_deriv({0, max_intersec}, max_intersec)*dmi;

    const double delta_uet = vc/f_uet;
    const double dduetdvc = 1.0/f_uet;
    std::vector<double> newz(tex_width);
    std::vector<double> dznewdvc(tex_width);
    for(size_t Y = 0; Y < tex_height; ++Y){
        const double y = static_cast<double>(Y);

        const double mult = smooth::abs(smooth::floor( (y + YOFF) / f));
        // const double dmult = 0;

        const double perimeter = 2*M_PI*cylinder_radius;

        const double phase_diff_x = 2*M_PI*fx*mult*perimeter*dimx/vc + phix;
        const double dpdx = -2*M_PI*fx*mult*perimeter*dimx/(vc*vc);
        const double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
        const double dpdz = -2*M_PI*fz*mult*perimeter*dimx/(vc*vc);
        const double phix_cur = phase_diff_x - smooth::floor((phase_diff_x)/(2*M_PI))*2*M_PI;
        const double dphix_cur = dpdx - smooth::floor_deriv((phase_diff_x)/(2*M_PI))*dpdx;
        const double phiz_cur = phase_diff_z - smooth::floor((phase_diff_z)/(2*M_PI))*2*M_PI;
        const double dphiz_cur = dpdz - smooth::floor_deriv((phase_diff_z)/(2*M_PI))*dpdz;

        const double xoffset_uet = mult*perimeter;

        for(size_t X = 0; X < tex_width; ++X){
            const double x = static_cast<double>(X);

            const double newx = x + Ax*std::sin(2*M_PI*fx*x*dimx/vc + phix_cur);
            const double dnewxdvc = Ax*std::cos(2*M_PI*fx*x*dimx/vc + phix_cur)*((-2)*M_PI*fx*x*dimx/(vc*vc) + dphix_cur);

            const double oscillation = Az*std::sin(2*M_PI*fz*newx*dimx/vc + phiz_cur);
            const double doscilldvc = Az*std::cos(2*M_PI*fz*newx*dimx/vc + phiz_cur)*(2*M_PI*fz*(dnewxdvc*dimx*vc - newx*dimx)/(vc*vc) + dphiz_cur);

            const double xcirc = x + xoffset_uet;

            const double mult_uet = smooth::abs(smooth::floor(xcirc / delta_uet));
            const double dmult_uet = smooth::abs_deriv(smooth::floor(xcirc / delta_uet))*smooth::floor_deriv(xcirc / delta_uet)*((-xcirc)/(delta_uet*delta_uet))*dduetdvc;

            // const double x_uet = (xcirc - mult_uet*delta_uet)*Ax_uet/delta_uet - Ax_uet/2;
            const double x_uet = Ax_uet*(xcirc/delta_uet - mult_uet - 0.5);
            const double dx_uet = -((Ax_uet*xcirc)/(delta_uet*delta_uet))*dduetdvc - Ax_uet*dmult_uet;

            const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)));
            const double duet = -0.5*(Az_uet/std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)))*(-2)*(x_uet/Ax_uet)*(dx_uet/Ax_uet);

            newz[X] = oscillation + uet_effect;
            dznewdvc[X] = doscilldvc + duet;
        }
            
        double shape_z;
        if(y <= y1 + mult*f + f/2){
            shape_z = -std::tan(alpha1)*(y - mult*f) + line_root1;
        } else if(y <= y2 + mult*f + f/2){
            shape_z = -std::sqrt(r*r - std::pow(y - mult*f - f/2, 2)) + r - ap;
        } else {
            shape_z = std::tan(alpha2)*(y - mult*f) + line_root2;
        }

        for(size_t X = 0; X < tex_width; ++X){
            const double oz = orig_z[X + Y*tex_width];
            double& dz = dzdvc[X + Y*tex_width];

            const double end_z = newz[X] + shape_z;

            dz = smooth::min_deriv({oz, end_z}, 1)*dznewdvc[X];
        }
        // } else {
        //     for(size_t X = 0; X < tex_width; ++X){
        //         double& dz = dzdvc[X + Y*tex_width];
        //         dz = 0;
        //     }
        // }
    }
}

void dzdap(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdap){
    using namespace param;

    vc *= dim_scale;
    const double y1 = -std::sqrt((std::pow(std::tan(alpha1)*r, 2))/(std::pow(std::tan(alpha1), 2)+1));
    const double y2 =  std::sqrt((std::pow(std::tan(alpha2)*r, 2))/(std::pow(std::tan(alpha2), 2)+1));

    const double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + std::tan(alpha1)*(y1+f/2);//-ap + std::tan(alpha1)*f/2;
    const double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - std::tan(alpha2)*(y2+f/2); //-ap - std::tan(alpha2)*f/2;

    const double dlrdap1 = -1;
    const double dlrdap2 = -1;

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
    const double max_intersec = -smooth::min({-intersec1, -intersec2});

    const double dmi = -smooth::min_deriv({-intersec1, -intersec2}, 0)*(-di1) - smooth::min_deriv({-intersec1, -intersec2}, 1)*(-di2);
    dmax_zdap = smooth::min_deriv({0, max_intersec}, 1)*dmi;

    const double delta_uet = vc/f_uet;
    std::vector<double> newz(tex_width);
    for(size_t Y = 0; Y < tex_height; ++Y){
        const double y = static_cast<double>(Y);

        const double mult = smooth::abs(smooth::floor( (y + YOFF) / f));

        const double perimeter = 2*M_PI*cylinder_radius;

        const double phase_diff_x = 2*M_PI*fx*mult*perimeter*dimx/vc + phix;
        const double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
        const double phix_cur = phase_diff_x - smooth::floor((phase_diff_x)/(2*M_PI))*2*M_PI;
        const double phiz_cur = phase_diff_z - smooth::floor((phase_diff_z)/(2*M_PI))*2*M_PI;

        const double xoffset_uet = mult*perimeter;

        for(size_t X = 0; X < tex_width; ++X){
            const double x = static_cast<double>(X);

            const double newx = x + Ax*std::sin(2*M_PI*fx*x*dimx/vc + phix_cur);
            const double oscillation = Az*std::sin(2*M_PI*fz*newx*dimx/vc + phiz_cur);

            const double xcirc = x + xoffset_uet;
            const double mult_uet = smooth::abs(smooth::floor(xcirc / delta_uet));
            // const double x_uet = (xcirc - mult_uet*delta_uet)*Ax_uet/delta_uet - Ax_uet/2;
            const double x_uet = Ax_uet*(xcirc/delta_uet - mult_uet - 0.5);
            const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)));

            newz[X] = oscillation + uet_effect;
        }

        double shape_z;
        double dznewdap;
        if(y <= y1 + mult*f + f/2){
            shape_z = -std::tan(alpha1)*(y - mult*f) + line_root1;
            dznewdap = dlrdap1;
        } else if(y <= y2 + mult*f + f/2){
            shape_z = -std::sqrt(r*r - std::pow(y - mult*f - f/2, 2)) + r - ap;
            dznewdap = -1;
        } else {
            shape_z = std::tan(alpha2)*(y - mult*f) + line_root2;
            dznewdap = dlrdap2;
        }

        for(size_t X = 0; X < tex_width; ++X){
            const double oz = orig_z[X + Y*tex_width];
            double& dz = dzdap[X + Y*tex_width];

            const double end_z = newz[X] + shape_z;

            dz = smooth::min_deriv({oz, end_z}, 1)*dznewdap;
        }
        // } else {
        //     for(size_t X = 0; X < tex_width; ++X){
        //         double& dz = dzdap[X + Y*tex_width];
        //         dz = 0;
        //     }
        // }
    }
}

void dzdf(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdf){
    using namespace param;

    vc *= dim_scale;
    const double y1 = -std::sqrt((std::pow(std::tan(alpha1)*r, 2))/(std::pow(std::tan(alpha1), 2)+1));
    const double y2 =  std::sqrt((std::pow(std::tan(alpha2)*r, 2))/(std::pow(std::tan(alpha2), 2)+1));

    const double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + std::tan(alpha1)*(y1+f/2);//-ap + std::tan(alpha1)*f/2;
    const double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - std::tan(alpha2)*(y2+f/2); //-ap - std::tan(alpha2)*f/2;

    const double dlrdf1 =  std::tan(alpha1)/2;
    const double dlrdf2 = -std::tan(alpha2)/2;

    double intersec1 = 0;
    double intersec2 = 0;
    double di1 = 0;
    double di2 = 0;
    if(y1 + f/2 > 0){
        intersec1 = smooth::min({0, line_root1});
        di1 = smooth::min_deriv({0, line_root1}, 1)*dlrdf1;
    } else {
        intersec1 = smooth::min({0, -std::sqrt(r*r - f*f/4) + r - ap});
        di1 = smooth::min_deriv({0, -std::sqrt(r*r - f*f/4) + r - ap}, 1)*f/(2*std::sqrt(4*r*r-f*f));
    }
    if(y2 + f/2 < f){
        intersec2 = smooth::min({0, line_root2 + std::tan(alpha2)*f});
        di2 = smooth::min_deriv({0, line_root2 + std::tan(alpha2)*f}, 1)*(dlrdf2 + std::tan(alpha2));
    } else {
        intersec2 = smooth::min({0, -std::sqrt(r*r - f*f/4) + r - ap});
        di2 = smooth::min_deriv({0, -std::sqrt(r*r - f*f/4) + r - ap}, 1)*f/(2*std::sqrt(4*r*r-f*f));
    }
    const double max_intersec = -smooth::min({-intersec1, -intersec2});

    const double dmi = -smooth::min_deriv({-intersec1, -intersec2}, 0)*(-di1) - smooth::min_deriv({-intersec1, -intersec2}, 1)*(-di2);
    dmax_zdf = smooth::min_deriv({0, max_intersec}, 1)*dmi;

    const double delta_uet = vc/f_uet;
    //const double dduet = 0;
    std::vector<double> newz(tex_width);
    std::vector<double> dznewdf(tex_width);
    for(size_t Y = 0; Y < tex_height; ++Y){
        const double y = static_cast<double>(Y);

        const double mult = smooth::abs(smooth::floor( (y + YOFF) / f));
        const double dmult = smooth::abs_deriv(smooth::floor( (y + YOFF) / f))*smooth::floor_deriv((y + YOFF) / f)*(-(y / (f*f)));

        const double perimeter = 2*M_PI*cylinder_radius;

        const double phase_diff_x = 2*M_PI*fx*mult*perimeter*dimx/vc + phix;
        const double dpdx = 2*M_PI*fx*dmult*perimeter*dimx/vc;
        const double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
        const double dpdz = 2*M_PI*fz*dmult*perimeter*dimx/vc;
        const double phix_cur = phase_diff_x - smooth::floor((phase_diff_x)/(2*M_PI))*2*M_PI;
        const double dphix_cur = dpdx - smooth::floor_deriv((phase_diff_x)/(2*M_PI))*dpdx;
        const double phiz_cur = phase_diff_z - smooth::floor((phase_diff_z)/(2*M_PI))*2*M_PI;
        const double dphiz_cur = dpdz - smooth::floor_deriv((phase_diff_z)/(2*M_PI))*dpdz;

        const double xoffset_uet = mult*perimeter;
        const double dxoff_uet = dmult*perimeter;

        for(size_t X = 0; X < tex_width; ++X){
            const double x = static_cast<double>(X);

            const double newx = x + Ax*std::sin(2*M_PI*fx*x*dimx/vc + phix_cur);
            const double dnewxdf = Ax*std::cos(2*M_PI*fx*x*dimx/vc + phix_cur)*dphix_cur;

            const double oscillation = Az*std::sin(2*M_PI*fz*newx*dimx/vc + phiz_cur);
            const double doscilldf = Az*std::cos(2*M_PI*fz*newx*dimx/vc + phiz_cur)*(2*M_PI*fz*dnewxdf*dimx/vc + dphiz_cur);

            const double xcirc = x + xoffset_uet;
            const double dxcirc = dxoff_uet;

            const double mult_uet = smooth::abs(smooth::floor(xcirc / delta_uet));
            const double dmult_uet = smooth::abs_deriv(smooth::floor(xcirc / delta_uet))*smooth::floor_deriv(xcirc / delta_uet)*dxcirc/delta_uet;

            // const double x_uet = (xcirc - mult_uet*delta_uet)*Ax_uet/delta_uet - Ax_uet/2;
            const double x_uet = Ax_uet*(xcirc/delta_uet - mult_uet - 0.5);
            const double dx_uet = Ax_uet*dxcirc/delta_uet - Ax_uet*dmult_uet;

            const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)));
            const double duet = -0.5*(Az_uet/std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)))*(-2)*(x_uet/Ax_uet)*(dx_uet/Ax_uet);

            newz[X] = oscillation + uet_effect;
            dznewdf[X] = doscilldf + duet;
        }

        double shape_z;
        double dshape_z;
        if(y <= y1 + mult*f + f/2){
            shape_z  = -std::tan(alpha1)*(y - mult*f) + line_root1;
            dshape_z =  std::tan(alpha1)*(mult + dmult*f) + dlrdf1;
        } else if(y <= y2 + mult*f + f/2){
            const double yy = y - mult*f - f/2;
            shape_z  = -std::sqrt(r*r - yy*yy) + r - ap;
            dshape_z = -yy*(mult + dmult*f + 0.5)/std::sqrt(r*r - yy*yy);
        } else {
            shape_z  =  std::tan(alpha2)*(y - mult*f) + line_root2;
            dshape_z = -std::tan(alpha2)*(mult + dmult*f) + dlrdf2;
        }

        for(size_t X = 0; X < tex_width; ++X){
            const double oz = orig_z[X + Y*tex_width];
            double& dz = dzdf[X + Y*tex_width];

            const double end_z  = newz[X] + shape_z;
            const double dend_z = dznewdf[X] + dshape_z;

            dz = smooth::min_deriv({oz, end_z}, 1)*dend_z;
        }
        // } else {
        //     for(size_t X = 0; X < tex_width; ++X){
        //         double& dz = dzdf[X + Y*tex_width];
        //         dz = 0;
        //     }
        // }
    }
}

}
