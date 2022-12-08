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

    std::copy(orig_z.begin(), orig_z.end(), map_z.begin());

    const double perimeter = 2*M_PI*(cylinder_radius - max_z);

    #pragma omp parallel for
    for(size_t Y = 0; Y < tex_height; ++Y){
        const double y = static_cast<double>(Y);

        // Calculate current row
        const double _mult = std::floor((y - 0.5*w_max)  / f + 1);

        for(size_t m = 0; m < overlap + 1; ++m){

            // Overlap
            const double mult = _mult + m;

            // Distance travelled by tool
            const double xoffset_uet = mult*perimeter;

            for(size_t X = 0; X < tex_width; ++X){
                const double x = static_cast<double>(X);

                // Consider distance travelled by tool
                const double xcirc = x + xoffset_uet;

                // Get distance relative to period
                const double mult_uet = std::floor(xcirc / delta_uet);
                const double x_uet = xcirc - (mult_uet + 0.5)*delta_uet;

                // Tool path
                double z_uet = 0;
                if(vc <= v_crit) {
                    const double H = Az_uet*std::sin(M_PI*vc/(2*v_crit));
                    const double h1 =  H + Az_uet;
                    const double h2 = -H + Az_uet;
                    const double K = delta_uet*(h1+h2)/(2*(h1-h2));
                    const double delta_1 =  K + delta_uet/2;

                    z_uet = h1*(1.0 - std::sqrt(1.0 - (4*x_uet*x_uet)/(delta_1*delta_1)));
                } else {
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

                // Get pixels
                double& z = map_z[X + Y*tex_width];

                // Write results
                z = std::min(z, shape_z);
            }
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

    // Where the slope lines intersect, so that the height at the edges is the
    // same for both
    double yc = (tan2*f+b2off-b1off)/(tan1+tan2);

    // Recalculate yc if it gets into the radius, as that skews the centering
    double s1 = -yc;
    double s2 = f-yc;
    if(alpha1 < alpha2 && s2 < y2){
        // Slope 2 must be disregarded.
        // Height of slope 1 must be equal to height of tool radius at the
        // other edge.
        const double a = tan1*tan1 + 1;
        const double b = -2*(tan1*(r-b1off)+f);
        const double c = b1off*(b1off-2*r)+f*f;
        yc = (-b - std::sqrt(b*b - 4*a*c))/(2*a);
    } else if(alpha2 < alpha1 && s1 > y1){
        // Slope 1 must be disregarded.
        // Height of slope 2 must be equal to height of tool radius at the
        // other edge.
        const double a = tan2*tan2 + 1;
        const double b = -(2*tan2*(b2off-r)+2*tan2*tan2*f);
        const double c = tan2*tan2*f*f + 2*tan2*f*(b2off-r) + b2off*(b2off-2*r);
        yc = (-b + std::sqrt(b*b - 4*a*c))/(2*a);
    }

    // Calculate whether the angles of the cutting tool intercept the surface
    // first or are stopped by the feed rate.
    //
    // Check the intersection with the feed rate edges. If they're greater than
    // zero, the cut line is narrower than the feed line, so max height is
    // zero. Otherwise, it's the height at the intersection
    //
    // Implementation of this section is based on `texture_shallow`
    double h = 0;
    s1 = -yc;
    s2 = f-yc;
    if(s1 > y1 && s2 < y2){
        // Width falls entirely within tool radius.
        yc = f/2;
        h = -std::sqrt(r*r - yc*yc) + r;
    } else if(s1 > y1){
        h = tan2*(f-yc) + b2off;
    } else if(s2 < y2){
        h = tan1*yc + b1off;
    } else {
        // Width does not intersect the radius at all.
        h = tan1*yc + b1off;
    }

    const double delta_uet = vc/f_uet;

    const double over_z = Az + Az_uet*(1.0 - std::sqrt(1.0 - (delta_uet*delta_uet)/(4*Ax_uet*Ax_uet)));

    min_z = -(ap + Az);
    // If it's greater than zero, max_z must be zero
    max_z = smooth::min({0, h+min_z+over_z+Az});
    
    // Value of z for y = 0 (global)
    const double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + tan1*(yc+y1);
    const double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - tan2*(yc+y2);

    #pragma omp parallel
    {
        std::vector<double> newz(tex_width);
        #pragma omp for
        for(size_t Y = 0; Y < tex_height; ++Y){
            const double y = static_cast<double>(Y);

            // Calculate current row
            // Apply `abs()` to prevent it from becoming less than zero
            // when close to zero
            const double mult = smooth::abs(smooth::floor((y + YOFF) / f));

            // Phase differences
            const double perimeter = 2*M_PI*(cylinder_radius - max_z);
            const double xoffset_uet = mult*perimeter;

            // It would be nice if it were possible to just cache these results
            // for a few iterations. Technically, they're constant for each row.
            // However, trying to cache them based on that made the algorithm
            // work incorrectly.
            //
            // Cheating and making `mult` actually exact _might_ be an option,
            // but last time I tried it just gave out NaNs.
            // 
            // So, for now, this will hog performance for a bit, unfortunately.
            for(size_t X = 0; X < tex_width; ++X){
                const double x = static_cast<double>(X);
                const double xcirc = x + xoffset_uet;

                // Random oscillation
                const double oscillation = Az*std::sin(2*M_PI*fz*xcirc*dimx/vc + phiz);

                // Ultrasonic turning effects
                const double mult_uet = smooth::abs(smooth::floor(xcirc / delta_uet));
                const double x_uet = xcirc - (mult_uet + 0.5)*delta_uet;
                const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)));

                newz[X] = oscillation + uet_effect;
            }

            // Tool shape
            double shape_z;
            if(y <= y1 + mult*f + yc){
                shape_z = -tan1*(y - mult*f) + line_root1;
            } else if(y <= y2 + mult*f + yc){
                const double yy = y - mult*f - yc;
                shape_z = -std::sqrt(r*r - yy*yy) + r - ap;
            } else {
                shape_z = tan2*(y - mult*f) + line_root2;
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
        }
    }
}

void dzdvc(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdvc){
    using namespace param;
    using param::y1;

    vc *= 1e6/60.0;

    double yc = (tan2*f+b2off-b1off)/(tan1+tan2);

    double s1 = -yc;
    double s2 = f-yc;
    if(alpha1 < alpha2 && s2 < y2){
        const double a = tan1*tan1 + 1;
        const double b = -2*(tan1*(r-b1off)+f);
        const double c = b1off*(b1off-2*r)+f*f;
        yc = (-b - std::sqrt(b*b - 4*a*c))/(2*a);
    } else if(alpha2 < alpha1 && s1 > y1){
        const double a = tan2*tan2 + 1;
        const double b = -(2*tan2*(b2off-r)+2*tan2*tan2*f);
        const double c = tan2*tan2*f*f + 2*tan2*f*(b2off-r) + b2off*(b2off-2*r);
        yc = (-b + std::sqrt(b*b - 4*a*c))/(2*a);
    }

    double h = 0;
    s1 = -yc;
    s2 = f-yc;
    if(s1 > y1 && s2 < y2){
        // Width falls entirely within tool radius.
        yc = f/2;
        h = -std::sqrt(r*r - yc*yc) + r;
    } else if(s1 > y1){
        h = tan2*(f-yc) + b2off;
    } else if(s2 < y2){
        h = tan1*yc + b1off;
    } else {
        // Width does not intersect the radius at all.
        h = tan1*yc + b1off;
    }

    const double delta_uet = vc/f_uet;
    const double dduetdvc = 1.0/f_uet;

    const double over_z = Az + Az_uet*(1.0 - std::sqrt(1.0 - (delta_uet*delta_uet)/(4*Ax_uet*Ax_uet)));
    const double dover_z = -0.5*(Az_uet/std::sqrt(1.0 - (delta_uet*delta_uet)/(4*Ax_uet*Ax_uet)))
                    *(-2)*(delta_uet/(2*Ax_uet))*(dduetdvc/(2*Ax_uet));

    //min_z = -(ap + Az);
    //max_z = smooth::min({0, h+min_z+over_z+Az});

    const double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + tan1*(yc+y1);
    const double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - tan2*(yc+y2);
    
    dmax_zdvc = smooth::min_deriv({0, h+min_z+over_z+Az}, 1)*dover_z;

    #pragma omp parallel
    {
        std::vector<double> newz(tex_width);
        std::vector<double> dznewdvc(tex_width);
        #pragma omp for
        for(size_t Y = 0; Y < tex_height; ++Y){
            const double y = static_cast<double>(Y);

            const double mult = smooth::abs(smooth::floor( (y + YOFF) / f));
            // const double dmult = 0;

            const double perimeter = 2*M_PI*(cylinder_radius - max_z);
            const double dperimeter = 2*M_PI*(-dmax_zdvc);

            const double xoffset_uet = mult*perimeter;
            const double dxoff_uet = mult*dperimeter;

            for(size_t X = 0; X < tex_width; ++X){
                const double x = static_cast<double>(X);
                const double xcirc = x + xoffset_uet;
                const double dxcirc = dxoff_uet;

                const double oscillation = Az*std::sin(2*M_PI*fz*xcirc*dimx/vc + phiz);
                const double doscilldvc = Az*std::cos(2*M_PI*fz*xcirc*dimx/vc + phiz)*2*M_PI*fz*dimx*(dxcirc*vc - xcirc)/(vc*vc);

                const double mult_uet = smooth::abs(smooth::floor(xcirc / delta_uet));
                const double dmult_uet = smooth::abs_deriv(smooth::floor(xcirc / delta_uet))*smooth::floor_deriv(xcirc / delta_uet)
                                                    *((dxcirc*delta_uet - xcirc*dduetdvc)/(delta_uet*delta_uet));

                const double x_uet = xcirc - (mult_uet + 0.5)*delta_uet;
                const double dx_uet = dxcirc - (dmult_uet*delta_uet + mult_uet*dduetdvc) - 0.5*dduetdvc;

                const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)));
                const double duet = -0.5*(Az_uet/std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)))*(-2)*(x_uet/Ax_uet)*(dx_uet/Ax_uet);

                newz[X] = oscillation + uet_effect;
                dznewdvc[X] = doscilldvc + duet;
            }
                
            double shape_z;
            if(y <= y1 + mult*f + yc){
                shape_z = -tan1*(y - mult*f) + line_root1;
            } else if(y <= y2 + mult*f + yc){
                const double yy = y - mult*f - yc;
                shape_z = -std::sqrt(r*r - yy*yy) + r - ap;
            } else {
                shape_z = tan2*(y - mult*f) + line_root2;
            }

            for(size_t X = 0; X < tex_width; ++X){
                const double oz = orig_z[X + Y*tex_width];
                double& dz = dzdvc[X + Y*tex_width];

                const double end_z = newz[X] + shape_z;

                dz = smooth::min_deriv({oz, end_z}, 1)*dznewdvc[X];
            }
        }
    }
}

void dzdap(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdap){
    using namespace param;
    using param::y1;

    vc *= 1e6/60.0;

    double yc = (tan2*f+b2off-b1off)/(tan1+tan2);

    double s1 = -yc;
    double s2 = f-yc;
    if(alpha1 < alpha2 && s2 < y2){
        const double a = tan1*tan1 + 1;
        const double b = -2*(tan1*(r-b1off)+f);
        const double c = b1off*(b1off-2*r)+f*f;
        yc = (-b - std::sqrt(b*b - 4*a*c))/(2*a);
    } else if(alpha2 < alpha1 && s1 > y1){
        const double a = tan2*tan2 + 1;
        const double b = -(2*tan2*(b2off-r)+2*tan2*tan2*f);
        const double c = tan2*tan2*f*f + 2*tan2*f*(b2off-r) + b2off*(b2off-2*r);
        yc = (-b + std::sqrt(b*b - 4*a*c))/(2*a);
    }

    double h = 0;
    s1 = -yc;
    s2 = f-yc;
    if(s1 > y1 && s2 < y2){
        yc = f/2;
        h = -std::sqrt(r*r - yc*yc) + r;
    } else if(s1 > y1){
        h = tan2*(f-yc) + b2off;
    } else if(s2 < y2){
        h = tan1*yc + b1off;
    } else {
        h = tan1*yc + b1off;
    }

    const double delta_uet = vc/f_uet;

    const double over_z = Az + Az_uet*(1.0 - std::sqrt(1.0 - (delta_uet*delta_uet)/(4*Ax_uet*Ax_uet)));

    // min_z = -(ap + Az);
    dmax_zdap = smooth::min_deriv({0, h+min_z+over_z+Az}, 1)*(-1);

    const double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + tan1*(yc+y1);
    const double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - tan2*(yc+y2);

    const double dlrdap1 = -1;
    const double dlrdap2 = -1;

    #pragma omp parallel
    {
        std::vector<double> newz(tex_width);
        std::vector<double> dznewdap(tex_width);
        #pragma omp for
        for(size_t Y = 0; Y < tex_height; ++Y){
            const double y = static_cast<double>(Y);

            const double mult = smooth::abs(smooth::floor( (y + YOFF) / f));

            const double perimeter = 2*M_PI*(cylinder_radius - max_z);
            const double dperimeter = -2*M_PI*dmax_zdap;

            const double xoffset_uet = mult*perimeter;
            const double dxoff_uet = mult*dperimeter;

            for(size_t X = 0; X < tex_width; ++X){
                const double x = static_cast<double>(X);

                const double xcirc = x + xoffset_uet;
                const double dxcirc = dxoff_uet;

                const double oscillation = Az*std::sin(2*M_PI*fz*xcirc*dimx/vc + phiz);
                const double doscilldap = Az*std::cos(2*M_PI*fz*xcirc*dimx/vc + phiz)*2*M_PI*fz*dxcirc*dimx/vc;

                const double mult_uet = smooth::abs(smooth::floor(xcirc / delta_uet));
                const double dmult_uet = smooth::abs_deriv(smooth::floor(xcirc / delta_uet))*smooth::floor_deriv(xcirc / delta_uet)*dxcirc/delta_uet;

                const double x_uet = xcirc - (mult_uet + 0.5)*delta_uet;
                const double dx_uet = dxcirc - dmult_uet*delta_uet;

                const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)));
                const double duet = -0.5*(Az_uet/std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)))*(-2)*(x_uet/Ax_uet)*(dx_uet/Ax_uet);

                newz[X] = oscillation + uet_effect;
                dznewdap[X] = doscilldap + duet;
            }

            double shape_z;
            double dshape_z;
            if(y <= y1 + mult*f + yc){
                shape_z = -tan1*(y - mult*f) + line_root1;
                dshape_z = dlrdap1;
            } else if(y <= y2 + mult*f + yc){
                const double yy = y - mult*f - yc;

                const double circ = std::sqrt(r*r - yy*yy);

                shape_z = -circ + r - ap;
                dshape_z = -1;
            } else {
                shape_z = tan2*(y - mult*f) + line_root2;
                dshape_z = dlrdap2;
            }

            for(size_t X = 0; X < tex_width; ++X){
                const double oz = orig_z[X + Y*tex_width];
                double& dz = dzdap[X + Y*tex_width];

                const double end_z = newz[X] + shape_z;
                const double dend_z = dznewdap[X] + dshape_z;

                dz = smooth::min_deriv({oz, end_z}, 1)*dend_z;
            }
        }
    }
}

void dzdf(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdf){
    using namespace param;
    using param::y1;

    vc *= 1e6/60;

    double yc = (tan2*f+b2off-b1off)/(tan1+tan2);
    double dyc = tan2/(tan1+tan2);

    double s1 = -yc;
    double s2 = f-yc;
    if(alpha1 < alpha2 && s2 < y2){
        const double a = tan1*tan1 + 1;
        const double b = -2*(tan1*(r-b1off)+f);
        const double db = -2;
        const double c = b1off*(b1off-2*r)+f*f;
        const double dc = 2*f;
        yc = (-b - std::sqrt(b*b - 4*a*c))/(2*a);
        dyc = (-db + 0.5*(2*b*db - 4*a*dc)/std::sqrt(b*b - 4*a*c))/(2*a);
    } else if(alpha2 < alpha1 && s1 > y1){
        const double a = tan2*tan2 + 1;
        const double b = -(2*tan2*(b2off-r)+2*tan2*tan2*f);
        const double db = -2*tan2*tan2;
        const double c = tan2*tan2*f*f + 2*tan2*f*(b2off-r) + b2off*(b2off-2*r);
        const double dc = 2*tan2*tan2*f + 2*tan2*(b2off-r);
        yc = (-b + std::sqrt(b*b - 4*a*c))/(2*a);
        dyc = (-db + 0.5*(2*b*db - 4*a*dc)/std::sqrt(b*b - 4*a*c))/(2*a);
    }

    double h = 0;
    double dh = 0;
    s1 = -yc;
    s2 = f-yc;
    if(s1 > y1 && s2 < y2){
        yc = f/2;
        dyc = 0.5;
        h = -std::sqrt(r*r - yc*yc) + r;
        dh = yc*dyc/std::sqrt(r*r - yc*yc);
    } else if(s1 > y1){
        h = tan2*(f-yc) + b2off;
        dh = tan2;
    } else if(s2 < y2){
        h = tan1*yc + b1off;
        dh = tan1;
    } else {
        h = tan1*yc + b1off;
        dh = tan1;
    }

    const double delta_uet = vc/f_uet;

    const double over_z = Az + Az_uet*(1.0 - std::sqrt(1.0 - (delta_uet*delta_uet)/(4*Ax_uet*Ax_uet)));

    dmax_zdf = smooth::min_deriv({0.0, h+min_z+over_z+Az}, 1)*dh;

    const double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + tan1*(yc+y1);
    const double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - tan2*(yc+y2);

    const double dlrdf1 =  tan1*dyc;
    const double dlrdf2 = -tan2*dyc;

    #pragma omp parallel
    {
        std::vector<double> newz(tex_width);
        std::vector<double> dznewdf(tex_width);
        #pragma omp for
        for(size_t Y = 0; Y < tex_height; ++Y){
            const double y = static_cast<double>(Y);

            const double mult = smooth::abs(smooth::floor( (y + YOFF) / f));
            const double dmult = smooth::abs_deriv(smooth::floor( (y + YOFF) / f))*smooth::floor_deriv((y + YOFF) / f)*(-(y / (f*f)));

            const double perimeter = 2*M_PI*(cylinder_radius - max_z);
            const double dperimeter = -2*M_PI*dmax_zdf;

            const double xoffset_uet = mult*perimeter;
            const double dxoff_uet = dmult*perimeter + mult*dperimeter;

            for(size_t X = 0; X < tex_width; ++X){
                const double x = static_cast<double>(X);

                const double xcirc = x + xoffset_uet;
                const double dxcirc = dxoff_uet;

                const double oscillation = Az*std::sin(2*M_PI*fz*xcirc*dimx/vc + phiz);
                const double doscilldf = Az*std::cos(2*M_PI*fz*xcirc*dimx/vc + phiz)*2*M_PI*fz*dxcirc*dimx/vc;

                const double mult_uet = smooth::abs(smooth::floor(xcirc / delta_uet));
                const double dmult_uet = smooth::abs_deriv(smooth::floor(xcirc / delta_uet))*smooth::floor_deriv(xcirc / delta_uet)*dxcirc/delta_uet;

                // const double x_uet = (xcirc - mult_uet*delta_uet)/delta_uet - 0.5;
                const double x_uet = xcirc - (mult_uet + 0.5)*delta_uet;
                const double dx_uet = dxcirc - dmult_uet*delta_uet;

                const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)));
                const double duet = -0.5*(Az_uet/std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)))*(-2)*(x_uet/Ax_uet)*(dx_uet/Ax_uet);

                newz[X] = oscillation + uet_effect;
                dznewdf[X] = doscilldf + duet;
            }

            double shape_z;
            double dshape_z;
            if(y <= y1 + mult*f + yc){
                shape_z = -tan1*(y - mult*f) + line_root1;
                dshape_z =  tan1*(mult + dmult*f) + dlrdf1;
            } else if(y <= y2 + mult*f + yc){
                const double yy = y - mult*f - yc;
                const double dyy = -mult - dmult*f - dyc;

                const double circ = std::sqrt(r*r - yy*yy);
                const double dcirc = -yy*dyy/circ;

                shape_z = -circ + r - ap;
                dshape_z = -dcirc;
            } else {
                shape_z = tan2*(y - mult*f) + line_root2;
                dshape_z = -tan2*(mult + dmult*f) + dlrdf2;
            }

            for(size_t X = 0; X < tex_width; ++X){
                const double oz = orig_z[X + Y*tex_width];
                double& dz = dzdf[X + Y*tex_width];

                const double end_z  = newz[X] + shape_z;
                const double dend_z = dznewdf[X] + dshape_z;

                dz = smooth::min_deriv({oz, end_z}, 1)*dend_z;
            }
        }
    }
}

}
