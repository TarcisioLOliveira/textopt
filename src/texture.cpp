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
#include <cmath>
#include <iostream>
#include <vector>

namespace texture{

void map_exact(std::vector<double>& map_z, const std::vector<double>& orig_z, double f, double ap, double vc){
    using namespace param;
    using param::y1;

    vc *= 1e6/60.0;

    // Calculate maximum width
    const double ap1 = -std::sqrt(r*r - param::y1*param::y1) + r;
    const double ap2 = -std::sqrt(r*r - param::y2*param::y2) + r;
    const double ap12min = std::min(ap1, ap2);
    const double ap12max = std::max(ap1, ap2);
    if(ap < ap12min){
        // Depth fits entirely within tool radius
        w_max = 2*std::sqrt(r*r - (ap - r)*(ap - r));
    } else if(ap < ap12max){
        // Depth fits partially within tool radius
        // Can only happen for alpha1 != alpha2
        w_max = std::sqrt(r*r - (ap - r)*(ap - r));
        if(alpha1 <= alpha2){
            w_max += (ap - b1off)/tan1;
        } else {
            w_max += (ap - b2off)/tan2;
        }
    } else {
        // Depth encompasses tool's straight edges
        w_max = (ap - b1off)/tan1 + (ap - b2off)/tan2;
    }

    const size_t overlap = std::floor(w_max / f);

    const double delta_uet = vc/f_uet;

    min_z = -ap;

    const double line_root1_const = -std::sqrt(r*r - y1*y1) + r - ap + tan1*y1;
    const double line_root2_const = -std::sqrt(r*r - y2*y2) + r - ap - tan2*y2;

    const double perimeter = 2*M_PI*cylinder_radius;

    std::vector<double> z_cand(overlap + 2, 0);

    #pragma omp parallel for firstprivate(z_cand)
    for(size_t Y = 0; Y < tex_height; ++Y){
        const double y = static_cast<double>(Y);

        // Calculate current row
        const double _mult = std::floor((y - 0.5*w_max)  / f + 1);

        for(size_t X = 0; X < tex_width; ++X){
            const double x = static_cast<double>(X);

            z_cand[0] = orig_z[X + Y*tex_width];

            for(size_t m = 0; m < overlap + 1; ++m){
                // Overlap
                const double mult = _mult + m;

                // Distance travelled by tool
                const double xoffset_uet = mult*perimeter;

                // Consider distance travelled by tool
                const double xcirc = x + xoffset_uet;

                // Get distance relative to period
                const double mult_uet = std::floor(xcirc / delta_uet);
                const double x_uet = xcirc - (mult_uet + 0.5)*delta_uet;

                // Tool path
                double z_uet = 0;
                if(vc <= v_crit) {
                    // If vc <= v_crit, model as ellipses in series
                    const double H = Az_uet*std::sin(M_PI*vc/(2*v_crit));
                    const double h1 =  H + Az_uet;
                    const double h2 = -H + Az_uet;
                    const double K = delta_uet*(h1+h2)/(2*(h1-h2));
                    const double delta_1 =  K + delta_uet/2;

                    z_uet = h1*(1.0 - std::sqrt(1.0 - (4*x_uet*x_uet)/(delta_1*delta_1)));
                } else {
                    // If otherwise, model as alternating ellipses with
                    // different dimensions.
                    //
                    // Also interpolate the ellipses with a senoidal model in
                    // order to achieve better representation of the tool
                    // path, especially as Ax_uet tends to 0
                    const double H = Az_uet*std::sin(M_PI*v_crit/(2*vc));
                    const double h1 =  H + Az_uet;
                    const double h2 = -H + Az_uet;
                    const double K = delta_uet*(h1-h2)/(2*(h1+h2));
                    const double delta_1 =  K + delta_uet/2;
                    const double delta_2 = -K + delta_uet/2;

                    double z_uet1 = 0;
                    double z_uet2 = 0;
                    if(x_uet < -delta_1/2){
                        const double x_uet2 = x_uet + delta_uet/2;

                        z_uet1 = h1 + h2*std::cos(M_PI*x_uet2/delta_2);
                        z_uet2 = h1 + h2*std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                    } else if(x_uet < delta_1/2){
                        const double x_uet2 = x_uet;

                        z_uet1 = h1*(1.0 - std::cos(M_PI*x_uet2/delta_1));
                        z_uet2 = h1*(1.0 - std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_1*delta_1)));
                    } else {
                        const double x_uet2 = x_uet - delta_uet/2;

                        z_uet1 = h1 + h2*std::cos(M_PI*x_uet2/delta_2);
                        z_uet2 = h1 + h2*std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                    }
    
                    const double v_ratio = v_crit/vc;
                    z_uet = z_uet1*std::sqrt(1 - v_ratio*v_ratio) + z_uet2*(1 - std::sqrt(1 - v_ratio*v_ratio));
                }

                // Clearance angle
                const double x_uet_c = xcirc + 0.5*delta_uet - std::floor(xcirc / delta_uet + 0.5)*delta_uet;
                const double z_clear = tanc*(delta_uet - x_uet_c);

                // Final UET height
                const double z_uet_min = std::min(z_uet, z_clear);

                const double line_root1 = line_root1_const + z_uet_min;
                const double line_root2 = line_root2_const + z_uet_min;

                // Tool shape
                double shape_z;
                if(y <= y1 + mult*f){
                    shape_z = -tan1*(y - mult*f) + line_root1;
                } else if(y <= y2 + mult*f){
                    const double yy = y - mult*f;
                    shape_z = -std::sqrt(r*r - yy*yy) + r - ap + z_uet_min;
                } else {
                    shape_z = tan2*(y - mult*f) + line_root2;
                }

                z_cand[m+1] = shape_z;
            }
            // Get pixels
            double& z = map_z[X + Y*tex_width];
            // Write results
            z = *std::min_element(z_cand.begin(), z_cand.end());
        }
    }

    if(overlap == 0){
        // If it's greater than zero, max_z must be zero
        max_z = 0;
    } else {
        max_z = *std::max_element(map_z.begin(), map_z.end());
    }
}

void map(std::vector<double>& map_z, const std::vector<double>& orig_z, double f, double ap, double vc){
    using namespace param;
    using param::y1;

    vc *= 1e6/60.0;

    // Calculate maximum width
    const double ap1 = -std::sqrt(r*r - param::y1*param::y1) + r;
    const double ap2 = -std::sqrt(r*r - param::y2*param::y2) + r;
    const double ap12min = std::min(ap1, ap2);
    const double ap12max = std::max(ap1, ap2);
    if(ap < ap12min){
        // Depth fits entirely within tool radius
        w_max = 2*std::sqrt(r*r - (ap - r)*(ap - r));
    } else if(ap < ap12max){
        // Depth fits partially within tool radius
        // Can only happen for alpha1 != alpha2
        w_max = std::sqrt(r*r - (ap - r)*(ap - r));
        if(alpha1 <= alpha2){
            w_max += (ap - b1off)/tan1;
        } else {
            w_max += (ap - b2off)/tan2;
        }
    } else {
        // Depth encompasses tool's straight edges
        w_max = (ap - b1off)/tan1 + (ap - b2off)/tan2;
    }

    // This will not be differentiated, so it can be use std::floor()
    const size_t overlap = std::floor(w_max / f);

    const double delta_uet = vc/f_uet;

    min_z = -ap;

    const double line_root1_const = -std::sqrt(r*r - y1*y1) + r - ap + tan1*y1;
    const double line_root2_const = -std::sqrt(r*r - y2*y2) + r - ap - tan2*y2;

    const double perimeter = 2*M_PI*cylinder_radius;

    std::vector<double> z_cand(overlap + 2, 0);

    #pragma omp parallel for firstprivate(z_cand)
    for(size_t Y = 0; Y < tex_height; ++Y){
        const double y = static_cast<double>(Y);

        // Calculate current row
        const double _mult = smooth::floor((y - 0.5*w_max)  / f + 1);

        for(size_t X = 0; X < tex_width; ++X){
            const double x = static_cast<double>(X);

            z_cand[0] = orig_z[X + Y*tex_width];

            for(size_t m = 0; m < overlap + 1; ++m){
                // Overlap
                const double mult = _mult + m;

                // Distance travelled by tool
                const double xoffset_uet = mult*perimeter;

                // Consider distance travelled by tool
                const double xcirc = x + xoffset_uet;

                // Get distance relative to period
                const double mult_uet = smooth::floor(xcirc / delta_uet);
                const double x_uet = xcirc - (mult_uet + 0.5)*delta_uet;

                // Tool path
                double z_uet = 0;
                if(vc <= v_crit) {
                    // If vc <= v_crit, model as ellipses in series
                    const double H = Az_uet*std::sin(M_PI*vc/(2*v_crit));
                    const double h1 =  H + Az_uet;
                    const double h2 = -H + Az_uet;
                    const double K = delta_uet*(h1+h2)/(2*(h1-h2));
                    const double delta_1 =  K + delta_uet/2;

                    z_uet = h1*(1.0 - std::sqrt(1.0 - (4*x_uet*x_uet)/(delta_1*delta_1)));
                } else {
                    // If otherwise, model as alternating ellipses with
                    // different dimensions.
                    //
                    // Also interpolate the ellipses with a senoidal model in
                    // order to achieve better representation of the tool
                    // path, especially as Ax_uet tends to 0
                    const double H = Az_uet*std::sin(M_PI*v_crit/(2*vc));
                    const double h1 =  H + Az_uet;
                    const double h2 = -H + Az_uet;
                    const double K = delta_uet*(h1-h2)/(2*(h1+h2));
                    const double delta_1 =  K + delta_uet/2;
                    const double delta_2 = -K + delta_uet/2;

                    double z_uet1 = 0;
                    double z_uet2 = 0;
                    if(x_uet < -delta_1/2){
                        const double x_uet2 = x_uet + delta_uet/2;

                        z_uet1 = h1 + h2*std::cos(M_PI*x_uet2/delta_2);
                        z_uet2 = h1 + h2*std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                    } else if(x_uet < delta_1/2){
                        const double x_uet2 = x_uet;

                        z_uet1 = h1*(1.0 - std::cos(M_PI*x_uet2/delta_1));
                        z_uet2 = h1*(1.0 - std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_1*delta_1)));
                    } else {
                        const double x_uet2 = x_uet - delta_uet/2;

                        z_uet1 = h1 + h2*std::cos(M_PI*x_uet2/delta_2);
                        z_uet2 = h1 + h2*std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                    }
    
                    const double v_ratio = v_crit/vc;
                    z_uet = z_uet1*std::sqrt(1 - v_ratio*v_ratio) + z_uet2*(1 - std::sqrt(1 - v_ratio*v_ratio));
                }

                // Clearance angle
                const double x_uet_c = xcirc + 0.5*delta_uet - smooth::floor(xcirc / delta_uet + 0.5)*delta_uet;
                const double z_clear = tanc*(delta_uet - x_uet_c);

                // Final UET height
                const double z_uet_min = smooth::min({z_uet, z_clear});

                const double line_root1 = line_root1_const + z_uet_min;
                const double line_root2 = line_root2_const + z_uet_min;

                // Tool shape
                double shape_z;
                if(y <= y1 + mult*f){
                    shape_z = -tan1*(y - mult*f) + line_root1;
                } else if(y <= y2 + mult*f){
                    const double yy = y - mult*f;
                    shape_z = -std::sqrt(r*r - yy*yy) + r - ap + z_uet_min;
                } else {
                    shape_z = tan2*(y - mult*f) + line_root2;
                }

                z_cand[m+1] = shape_z;
            }
            // Get pixels
            double& z = map_z[X + Y*tex_width];
            // Write results
            z = smooth::min(z_cand);
        }
    }

    if(overlap == 0){
        // If it's greater than zero, max_z must be zero
        max_z = 0;
    } else {
        max_z = *std::max_element(map_z.begin(), map_z.end());
    }
}

void dzdvc(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdvc){
    using namespace param;
    using param::y1;

    vc *= 1e6/60.0;

    // Calculate maximum width
    const double ap1 = -std::sqrt(r*r - param::y1*param::y1) + r;
    const double ap2 = -std::sqrt(r*r - param::y2*param::y2) + r;
    const double ap12min = std::min(ap1, ap2);
    const double ap12max = std::max(ap1, ap2);
    if(ap < ap12min){
        // Depth fits entirely within tool radius
        w_max = 2*std::sqrt(r*r - (ap - r)*(ap - r));
    } else if(ap < ap12max){
        // Depth fits partially within tool radius
        // Can only happen for alpha1 != alpha2
        w_max = std::sqrt(r*r - (ap - r)*(ap - r));
        if(alpha1 <= alpha2){
            w_max += (ap - b1off)/tan1;
        } else {
            w_max += (ap - b2off)/tan2;
        }
    } else {
        // Depth encompasses tool's straight edges
        w_max = (ap - b1off)/tan1 + (ap - b2off)/tan2;
    }

    // This will not be differentiated, so it can be use std::floor()
    const size_t overlap = std::floor(w_max / f);

    const double delta_uet = vc/f_uet;
    const double ddelta_uet = 1.0/f_uet;

    min_z = -ap;

    const double line_root1_const = -std::sqrt(r*r - y1*y1) + r - ap + tan1*y1;
    const double line_root2_const = -std::sqrt(r*r - y2*y2) + r - ap - tan2*y2;

    const double perimeter = 2*M_PI*cylinder_radius;

    std::vector<double> z_cand(overlap + 2, 0);
    std::vector<double> dz_cand(overlap + 2, 0);

    #pragma omp parallel for firstprivate(z_cand,dz_cand)
    for(size_t Y = 0; Y < tex_height; ++Y){
        const double y = static_cast<double>(Y);

        // Calculate current row
        const double _mult = smooth::floor((y - 0.5*w_max)  / f + 1);

        for(size_t X = 0; X < tex_width; ++X){
            const double x = static_cast<double>(X);

            z_cand[0] = orig_z[X + Y*tex_width];

            for(size_t m = 0; m < overlap + 1; ++m){
                // Overlap
                const double mult = _mult + m;

                // Distance travelled by tool
                const double xoffset_uet = mult*perimeter;

                // Consider distance travelled by tool
                const double xcirc = x + xoffset_uet;

                // Get distance relative to period
                const double mult_uet = smooth::floor(xcirc / delta_uet);
                const double dmult_uet = smooth::floor_deriv(xcirc / delta_uet)*(-xcirc*ddelta_uet/(delta_uet*delta_uet));
                const double x_uet = xcirc - (mult_uet + 0.5)*delta_uet;
                const double dx_uet = - ((mult_uet + 0.5)*ddelta_uet + dmult_uet*delta_uet);

                // Tool path
                double z_uet = 0;
                double dz_uet = 0;
                if(vc <= v_crit) {
                    // If vc <= v_crit, model as ellipses in series
                    const double H = Az_uet*std::sin(M_PI*vc/(2*v_crit));
                    const double dH = Az_uet*std::cos(M_PI*vc/(2*v_crit))*M_PI/(2*v_crit);

                    const double h1 =  H + Az_uet;
                    const double dh1 =  dH;

                    const double h2 = -H + Az_uet;
                    const double dh2 = -dH;

                    const double K = delta_uet*(h1+h2)/(2*(h1-h2));
                    const double dK = ddelta_uet*(h1+h2)/(2*(h1-h2)) + delta_uet*((dh1+dh2)*(h1-h2) - (h1+h2)*(dh1-dh2))/(2*(h1-h2)*(h1-h2));

                    const double delta_1 =  K + delta_uet/2;
                    const double ddelta_1 = dK + ddelta_uet/2;

                    z_uet = h1*(1.0 - std::sqrt(1.0 - (4*x_uet*x_uet)/(delta_1*delta_1)));
                    const double ell = std::sqrt(1.0 - (4*x_uet*x_uet)/(delta_1*delta_1));
                    dz_uet = dh1*ell - h1*(-4*(x_uet*dx_uet*delta_1*delta_1-x_uet*x_uet*delta_1*ddelta_1)/(delta_1*delta_1*delta_1*delta_1))/ell;
                } else {
                    // If otherwise, model as alternating ellipses with
                    // different dimensions.
                    //
                    // Also interpolate the ellipses with a senoidal model in
                    // order to achieve better representation of the tool
                    // path, especially as Ax_uet tends to 0
                    const double H = Az_uet*std::sin(M_PI*v_crit/(2*vc));
                    const double dH = Az_uet*std::cos(M_PI*v_crit/(2*vc))*(-M_PI*v_crit/(2*vc*vc));

                    const double h1 =  H + Az_uet;
                    const double dh1 =  dH;

                    const double h2 = -H + Az_uet;
                    const double dh2 = -dH;

                    const double K = delta_uet*(h1-h2)/(2*(h1+h2));
                    const double dK = ddelta_uet*(h1-h2)/(2*(h1+h2)) + delta_uet*((dh1-dh2)*(h1+h2) - (h1-h2)*(dh1+dh2))/(2*(h1+h2));

                    const double delta_1 =  K + delta_uet/2;
                    const double ddelta_1 =  dK + ddelta_uet/2;

                    const double delta_2 = -K + delta_uet/2;
                    const double ddelta_2 = -dK + ddelta_uet/2;

                    double z_uet1 = 0;
                    double z_uet2 = 0;
                    double dz_uet1 = 0;
                    double dz_uet2 = 0;
                    if(x_uet < -delta_1/2){
                        const double x_uet2 = x_uet + delta_uet/2;
                        const double dx_uet2 = dx_uet + ddelta_uet/2;

                        z_uet1 = h1 + h2*std::cos(M_PI*x_uet2/delta_2);
                        dz_uet1 = dh1 + dh2*std::cos(M_PI*x_uet2/delta_2) - h2*std::sin(M_PI*x_uet2/delta_2)*M_PI*(dx_uet2*delta_2 - x_uet2*ddelta_2)/(delta_2*delta_2);

                        z_uet2 = h1 + h2*std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                        const double ell = std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                        dz_uet2 = dh1 + dh2*ell + h2*(-4*(x_uet2*dx_uet2*delta_2*delta_2-x_uet2*x_uet2*delta_2*ddelta_2)/(delta_2*delta_2*delta_2*delta_2))/ell;
                    } else if(x_uet < delta_1/2){
                        const double x_uet2 = x_uet;
                        const double dx_uet2 = dx_uet;

                        z_uet1 = h1*(1.0 - std::cos(M_PI*x_uet2/delta_1));
                        dz_uet1 = dh1 + dh1*std::cos(M_PI*x_uet2/delta_1) - h1*std::sin(M_PI*x_uet2/delta_1)*M_PI*(dx_uet2*delta_1 - x_uet2*ddelta_1)/(delta_1*delta_1);

                        z_uet2 = h1*(1.0 - std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_1*delta_1)));
                        const double ell = std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_1*delta_1));
                        dz_uet = dh1 + dh1*ell - h1*(-4*(x_uet2*dx_uet2*delta_1*delta_1-x_uet2*x_uet2*delta_1*ddelta_1)/(delta_1*delta_1*delta_1*delta_1))/ell;
                    } else {
                        const double x_uet2 = x_uet - delta_uet/2;
                        const double dx_uet2 = dx_uet - ddelta_uet/2;

                        z_uet1 = h1 + h2*std::cos(M_PI*x_uet2/delta_2);
                        dz_uet1 = dh1 + dh2*std::cos(M_PI*x_uet2/delta_2) - h2*std::sin(M_PI*x_uet2/delta_2)*M_PI*(dx_uet2*delta_2 - x_uet2*ddelta_2)/(delta_2*delta_2);

                        z_uet2 = h1 + h2*std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                        const double ell = std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                        dz_uet2 = dh1 + dh2*ell + h2*(-4*(x_uet2*dx_uet2*delta_2*delta_2-x_uet2*x_uet2*delta_2*ddelta_2)/(delta_2*delta_2*delta_2*delta_2))/ell;
                    }
    
                    const double v_ratio = v_crit/vc;
                    const double dv_ratio = -v_crit/(vc*vc);

                    const double rat = std::sqrt(1 - v_ratio*v_ratio);
                    z_uet = z_uet1*std::sqrt(1 - v_ratio*v_ratio) + z_uet2*(1 - std::sqrt(1 - v_ratio*v_ratio));
                    dz_uet = dz_uet1*rat + z_uet1*(-v_ratio*dv_ratio)/rat + dz_uet2*(1 - rat) - z_uet2*(v_ratio*dv_ratio)/rat;
                }

                // Clearance angle
                const double x_uet_c = xcirc + 0.5*delta_uet - smooth::floor(xcirc / delta_uet + 0.5)*delta_uet;
                const double dx_uet_c = 0.5*ddelta_uet - (smooth::floor(xcirc / delta_uet + 0.5)*ddelta_uet + smooth::floor_deriv(xcirc / delta_uet + 0.5)*delta_uet*(-xcirc/(delta_uet*delta_uet)));

                const double z_clear = tanc*(delta_uet - x_uet_c);
                const double dz_clear = tanc*(ddelta_uet - dx_uet_c);

                // Final UET height
                const double z_uet_min = smooth::min({z_uet, z_clear});
                const double dz_uet_min = smooth::min_deriv({z_uet, z_clear},{dz_uet, dz_clear});

                const double line_root1 = line_root1_const + z_uet_min;
                const double line_root2 = line_root2_const + z_uet_min;
                const double dline_root1 = line_root1_const + dz_uet_min;
                const double dline_root2 = line_root2_const + dz_uet_min;

                // Tool shape
                double shape_z;
                double dshape_z;
                if(y <= y1 + mult*f){
                    shape_z = -tan1*(y - mult*f) + line_root1;
                    dshape_z = dline_root1;
                } else if(y <= y2 + mult*f){
                    const double yy = y - mult*f;
                    shape_z = -std::sqrt(r*r - yy*yy) + r - ap + z_uet_min;
                    dshape_z = dz_uet_min;
                } else {
                    shape_z = tan2*(y - mult*f) + line_root2;
                    dshape_z = dline_root2;
                }

                z_cand[m+1] = shape_z;
                dz_cand[m+1] = dshape_z;
            }
            // Get pixels
            double& dz = dzdvc[X + Y*tex_width];
            // Write results
            dz = smooth::min_deriv(z_cand, dz_cand);
        }
    }
}

void dzdap(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdap){
    using namespace param;
    using param::y1;

    vc *= 1e6/60.0;

    // Calculate maximum width
    const double ap1 = -std::sqrt(r*r - param::y1*param::y1) + r;
    const double ap2 = -std::sqrt(r*r - param::y2*param::y2) + r;
    const double ap12min = std::min(ap1, ap2);
    const double ap12max = std::max(ap1, ap2);

    double dw_max;
    if(ap < ap12min){
        // Depth fits entirely within tool radius
        w_max = 2*std::sqrt(r*r - (ap - r)*(ap - r));
        dw_max = 2*(-(ap - r))/std::sqrt(r*r - (ap - r)*(ap - r));
    } else if(ap < ap12max){
        // Depth fits partially within tool radius
        // Can only happen for alpha1 != alpha2
        w_max = std::sqrt(r*r - (ap - r)*(ap - r));
        dw_max = 2*(-(ap - r))/std::sqrt(r*r - (ap - r)*(ap - r));
        if(alpha1 <= alpha2){
            w_max += (ap - b1off)/tan1;
            dw_max += 1.0/tan1;
        } else {
            w_max += (ap - b2off)/tan2;
            dw_max += 1.0/tan2;
        }
    } else {
        // Depth encompasses tool's straight edges
        w_max = (ap - b1off)/tan1 + (ap - b2off)/tan2;
        dw_max = 1.0/tan1 + 1.0/tan2;
    }

    // This will not be differentiated, so it can be use std::floor()
    const size_t overlap = std::floor(w_max / f);

    const double delta_uet = vc/f_uet;

    min_z = -ap;

    const double line_root1_const = -std::sqrt(r*r - y1*y1) + r - ap + tan1*y1;
    const double line_root2_const = -std::sqrt(r*r - y2*y2) + r - ap - tan2*y2;
    const double dline_root1_const = -1;
    const double dline_root2_const = -1;

    const double perimeter = 2*M_PI*cylinder_radius;

    std::vector<double> z_cand(overlap + 2, 0);
    std::vector<double> dz_cand(overlap + 2, 0);

    #pragma omp parallel for firstprivate(z_cand,dz_cand)
    for(size_t Y = 0; Y < tex_height; ++Y){
        const double y = static_cast<double>(Y);

        // Calculate current row
        const double _mult = smooth::floor((y - 0.5*w_max)  / f + 1);
        const double _dmult = smooth::floor_deriv((y - 0.5*w_max)  / f + 1)*(-0.5*dw_max);

        for(size_t X = 0; X < tex_width; ++X){
            const double x = static_cast<double>(X);

            z_cand[0] = orig_z[X + Y*tex_width];

            for(size_t m = 0; m < overlap + 1; ++m){
                // Overlap
                const double mult = _mult + m;
                const double dmult = _dmult;

                // Distance travelled by tool
                const double xoffset_uet = mult*perimeter;
                const double dxoffset_uet = dmult*perimeter;

                // Consider distance travelled by tool
                const double xcirc = x + xoffset_uet;
                const double dxcirc = dxoffset_uet;

                // Get distance relative to period
                const double mult_uet = smooth::floor(xcirc / delta_uet);
                const double dmult_uet = smooth::floor_deriv(xcirc / delta_uet)*dxcirc/delta_uet;

                const double x_uet = xcirc - (mult_uet + 0.5)*delta_uet;
                const double dx_uet = dxcirc - dmult_uet*delta_uet;

                // Tool path
                double z_uet = 0;
                double dz_uet = 0;
                if(vc <= v_crit) {
                    // If vc <= v_crit, model as ellipses in series
                    const double H = Az_uet*std::sin(M_PI*vc/(2*v_crit));
                    const double h1 =  H + Az_uet;
                    const double h2 = -H + Az_uet;
                    const double K = delta_uet*(h1+h2)/(2*(h1-h2));
                    const double delta_1 =  K + delta_uet/2;

                    z_uet = h1*(1.0 - std::sqrt(1.0 - (4*x_uet*x_uet)/(delta_1*delta_1)));
                    const double ell = std::sqrt(1.0 - (4*x_uet*x_uet)/(delta_1*delta_1));
                    dz_uet = h1*((4*x_uet*dx_uet)/(delta_1*delta_1))/ell;
                } else {
                    // If otherwise, model as alternating ellipses with
                    // different dimensions.
                    //
                    // Also interpolate the ellipses with a senoidal model in
                    // order to achieve better representation of the tool
                    // path, especially as Ax_uet tends to 0
                    const double H = Az_uet*std::sin(M_PI*v_crit/(2*vc));
                    const double h1 =  H + Az_uet;
                    const double h2 = -H + Az_uet;
                    const double K = delta_uet*(h1-h2)/(2*(h1+h2));
                    const double delta_1 =  K + delta_uet/2;
                    const double delta_2 = -K + delta_uet/2;

                    double z_uet1 = 0;
                    double z_uet2 = 0;
                    double dz_uet1 = 0;
                    double dz_uet2 = 0;
                    if(x_uet < -delta_1/2){
                        const double x_uet2 = x_uet + delta_uet/2;
                        const double dx_uet2 = dx_uet;

                        z_uet1 = h1 + h2*std::cos(M_PI*x_uet2/delta_2);
                        dz_uet1 = -h2*std::sin(M_PI*x_uet2/delta_2)*M_PI*dx_uet2/delta_2;

                        z_uet2 = h1 + h2*std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                        const double ell = std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                        dz_uet = h2*((4*x_uet2*dx_uet2)/(delta_2*delta_2))/ell;
                    } else if(x_uet < delta_1/2){
                        const double x_uet2 = x_uet;
                        const double dx_uet2 = dx_uet;

                        z_uet1 = h1*(1.0 - std::cos(M_PI*x_uet2/delta_1));
                        dz_uet1 = h1*std::sin(M_PI*x_uet2/delta_1)*M_PI*dx_uet2/delta_1;

                        z_uet2 = h1*(1.0 - std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_1*delta_1)));
                        const double ell = std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_1*delta_1));
                        dz_uet = -h1*((4*x_uet2*dx_uet2)/(delta_1*delta_1))/ell;
                    } else {
                        const double x_uet2 = x_uet - delta_uet/2;
                        const double dx_uet2 = dx_uet;

                        z_uet1 = h1 + h2*std::cos(M_PI*x_uet2/delta_2);
                        dz_uet1 = -h2*std::sin(M_PI*x_uet2/delta_2)*M_PI*dx_uet2/delta_2;

                        z_uet2 = h1 + h2*std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                        const double ell = std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                        dz_uet = h2*((4*x_uet2*dx_uet2)/(delta_2*delta_2))/ell;
                    }
    
                    const double v_ratio = v_crit/vc;
                    const double dv_ratio = -v_crit/(vc*vc);

                    const double rat = std::sqrt(1 - v_ratio*v_ratio);
                    z_uet = z_uet1*std::sqrt(1 - v_ratio*v_ratio) + z_uet2*(1 - std::sqrt(1 - v_ratio*v_ratio));
                    dz_uet = dz_uet1*rat + z_uet1*(-v_ratio*dv_ratio)/rat + dz_uet2*(1 - rat) - z_uet2*(v_ratio*dv_ratio)/rat;
                }

                // Clearance angle
                const double x_uet_c = xcirc + 0.5*delta_uet - smooth::floor(xcirc / delta_uet + 0.5)*delta_uet;
                const double dx_uet_c = dxcirc  - smooth::floor_deriv(xcirc / delta_uet + 0.5)*dxcirc;

                const double z_clear = tanc*(delta_uet - x_uet_c);
                const double dz_clear = tanc*(delta_uet - dx_uet_c);

                // Final UET height
                const double z_uet_min = smooth::min({z_uet, z_clear});
                const double dz_uet_min = smooth::min_deriv({z_uet, z_clear},{dz_uet, dz_clear});

                const double line_root1 = line_root1_const + z_uet_min;
                const double line_root2 = line_root2_const + z_uet_min;
                const double dline_root1 = dline_root1_const + dz_uet_min;
                const double dline_root2 = dline_root2_const + dz_uet_min;

                // Tool shape
                double shape_z;
                double dshape_z;
                if(y <= y1 + mult*f){
                    shape_z = -tan1*(y - mult*f) + line_root1;
                    dshape_z = dline_root1;
                } else if(y <= y2 + mult*f){
                    const double yy = y - mult*f;
                    shape_z = -std::sqrt(r*r - yy*yy) + r - ap + z_uet_min;
                    dshape_z = -1.0 + dz_uet_min;
                } else {
                    shape_z = tan2*(y - mult*f) + line_root2;
                    dshape_z = dline_root2;
                }

                z_cand[m+1] = shape_z;
                dz_cand[m+1] = dshape_z;
            }
            // Get pixels
            double& dz = dzdap[X + Y*tex_width];
            // Write results
            dz = smooth::min_deriv(z_cand,dz_cand);
        }
    }
}

void dzdf(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdf){
    using namespace param;
    using param::y1;

    vc *= 1e6/60.0;

    // Calculate maximum width
    const double ap1 = -std::sqrt(r*r - param::y1*param::y1) + r;
    const double ap2 = -std::sqrt(r*r - param::y2*param::y2) + r;
    const double ap12min = std::min(ap1, ap2);
    const double ap12max = std::max(ap1, ap2);

    if(ap < ap12min){
        // Depth fits entirely within tool radius
        w_max = 2*std::sqrt(r*r - (ap - r)*(ap - r));
    } else if(ap < ap12max){
        // Depth fits partially within tool radius
        // Can only happen for alpha1 != alpha2
        w_max = std::sqrt(r*r - (ap - r)*(ap - r));
        if(alpha1 <= alpha2){
            w_max += (ap - b1off)/tan1;
        } else {
            w_max += (ap - b2off)/tan2;
        }
    } else {
        // Depth encompasses tool's straight edges
        w_max = (ap - b1off)/tan1 + (ap - b2off)/tan2;
    }

    // This will not be differentiated, so it can be use std::floor()
    const size_t overlap = std::floor(w_max / f);

    const double delta_uet = vc/f_uet;

    min_z = -ap;

    const double line_root1_const = -std::sqrt(r*r - y1*y1) + r - ap + tan1*y1;
    const double line_root2_const = -std::sqrt(r*r - y2*y2) + r - ap - tan2*y2;

    const double perimeter = 2*M_PI*cylinder_radius;

    std::vector<double> z_cand(overlap + 2, 0);
    std::vector<double> dz_cand(overlap + 2, 0);

    #pragma omp parallel for firstprivate(z_cand,dz_cand)
    for(size_t Y = 0; Y < tex_height; ++Y){
        const double y = static_cast<double>(Y);

        // Calculate current row
        const double _mult = smooth::floor((y - 0.5*w_max)  / f + 1);
        const double _dmult = smooth::floor_deriv((y - 0.5*w_max)  / f + 1)*(-(y - 0.5*w_max)  / (f*f));

        for(size_t X = 0; X < tex_width; ++X){
            const double x = static_cast<double>(X);

            z_cand[0] = orig_z[X + Y*tex_width];

            for(size_t m = 0; m < overlap + 1; ++m){
                // Overlap
                const double mult = _mult + m;
                const double dmult = _dmult;

                // Distance travelled by tool
                const double xoffset_uet = mult*perimeter;
                const double dxoffset_uet = dmult*perimeter;

                // Consider distance travelled by tool
                const double xcirc = x + xoffset_uet;
                const double dxcirc = dxoffset_uet;

                // Get distance relative to period
                const double mult_uet = smooth::floor(xcirc / delta_uet);
                const double dmult_uet = smooth::floor_deriv(xcirc / delta_uet)*dxcirc/delta_uet;

                const double x_uet = xcirc - (mult_uet + 0.5)*delta_uet;
                const double dx_uet = dxcirc - dmult_uet*delta_uet;

                // Tool path
                double z_uet = 0;
                double dz_uet = 0;
                if(vc <= v_crit) {
                    // If vc <= v_crit, model as ellipses in series
                    const double H = Az_uet*std::sin(M_PI*vc/(2*v_crit));
                    const double h1 =  H + Az_uet;
                    const double h2 = -H + Az_uet;
                    const double K = delta_uet*(h1+h2)/(2*(h1-h2));
                    const double delta_1 =  K + delta_uet/2;

                    z_uet = h1*(1.0 - std::sqrt(1.0 - (4*x_uet*x_uet)/(delta_1*delta_1)));
                    const double ell = std::sqrt(1.0 - (4*x_uet*x_uet)/(delta_1*delta_1));
                    dz_uet = h1*((4*x_uet*dx_uet)/(delta_1*delta_1))/ell;
                } else {
                    // If otherwise, model as alternating ellipses with
                    // different dimensions.
                    //
                    // Also interpolate the ellipses with a senoidal model in
                    // order to achieve better representation of the tool
                    // path, especially as Ax_uet tends to 0
                    const double H = Az_uet*std::sin(M_PI*v_crit/(2*vc));
                    const double h1 =  H + Az_uet;
                    const double h2 = -H + Az_uet;
                    const double K = delta_uet*(h1-h2)/(2*(h1+h2));
                    const double delta_1 =  K + delta_uet/2;
                    const double delta_2 = -K + delta_uet/2;

                    double z_uet1 = 0;
                    double z_uet2 = 0;
                    double dz_uet1 = 0;
                    double dz_uet2 = 0;
                    if(x_uet < -delta_1/2){
                        const double x_uet2 = x_uet + delta_uet/2;
                        const double dx_uet2 = dx_uet;

                        z_uet1 = h1 + h2*std::cos(M_PI*x_uet2/delta_2);
                        dz_uet1 = -h2*std::sin(M_PI*x_uet2/delta_2)*M_PI*dx_uet2/delta_2;

                        z_uet2 = h1 + h2*std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                        const double ell = std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                        dz_uet = h2*((4*x_uet2*dx_uet2)/(delta_2*delta_2))/ell;
                    } else if(x_uet < delta_1/2){
                        const double x_uet2 = x_uet;
                        const double dx_uet2 = dx_uet;

                        z_uet1 = h1*(1.0 - std::cos(M_PI*x_uet2/delta_1));
                        dz_uet1 = h1*std::sin(M_PI*x_uet2/delta_1)*M_PI*dx_uet2/delta_1;

                        z_uet2 = h1*(1.0 - std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_1*delta_1)));
                        const double ell = std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_1*delta_1));
                        dz_uet = -h1*((4*x_uet2*dx_uet2)/(delta_1*delta_1))/ell;
                    } else {
                        const double x_uet2 = x_uet - delta_uet/2;
                        const double dx_uet2 = dx_uet;

                        z_uet1 = h1 + h2*std::cos(M_PI*x_uet2/delta_2);
                        dz_uet1 = -h2*std::sin(M_PI*x_uet2/delta_2)*M_PI*dx_uet2/delta_2;

                        z_uet2 = h1 + h2*std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                        const double ell = std::sqrt(1.0 - (4*x_uet2*x_uet2)/(delta_2*delta_2));
                        dz_uet = h2*((4*x_uet2*dx_uet2)/(delta_2*delta_2))/ell;
                    }
    
                    const double v_ratio = v_crit/vc;
                    const double dv_ratio = -v_crit/(vc*vc);

                    const double rat = std::sqrt(1 - v_ratio*v_ratio);
                    z_uet = z_uet1*std::sqrt(1 - v_ratio*v_ratio) + z_uet2*(1 - std::sqrt(1 - v_ratio*v_ratio));
                    dz_uet = dz_uet1*rat + z_uet1*(-v_ratio*dv_ratio)/rat + dz_uet2*(1 - rat) - z_uet2*(v_ratio*dv_ratio)/rat;
                }

                // Clearance angle
                const double x_uet_c = xcirc + 0.5*delta_uet - smooth::floor(xcirc / delta_uet + 0.5)*delta_uet;
                const double dx_uet_c = dxcirc  - smooth::floor_deriv(xcirc / delta_uet + 0.5)*dxcirc;

                const double z_clear = tanc*(delta_uet - x_uet_c);
                const double dz_clear = tanc*(delta_uet - dx_uet_c);

                // Final UET height
                const double z_uet_min = smooth::min({z_uet, z_clear});
                const double dz_uet_min = smooth::min_deriv({z_uet, z_clear},{dz_uet, dz_clear});

                const double line_root1 = line_root1_const + z_uet_min;
                const double line_root2 = line_root2_const + z_uet_min;
                const double dline_root1 = dz_uet_min;
                const double dline_root2 = dz_uet_min;

                // Tool shape
                double shape_z;
                double dshape_z;
                if(y <= y1 + mult*f){
                    shape_z = -tan1*(y - mult*f) + line_root1;
                    dshape_z = tan1*(dmult*f + mult) +  dline_root1;
                } else if(y <= y2 + mult*f){
                    const double yy = y - mult*f;
                    const double dyy = -(dmult*f + mult);
                    shape_z = -std::sqrt(r*r - yy*yy) + r - ap + z_uet_min;
                    dshape_z = yy*dyy/std::sqrt(r*r - yy*yy) + dz_uet_min;
                } else {
                    shape_z = tan2*(y - mult*f) + line_root2;
                    dshape_z = tan1*(dmult*f + mult) +  dline_root2;
                }

                z_cand[m+1] = shape_z;
                dz_cand[m+1] = dshape_z;
            }
            // Get pixels
            double& dz = dzdf[X + Y*tex_width];
            // Write results
            dz = smooth::min_deriv(z_cand,dz_cand);
        }
    }
}

}
