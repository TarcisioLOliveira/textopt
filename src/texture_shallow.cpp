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

#include "texture_shallow.hpp"
#include "param.hpp"
#include "smooth.hpp"
#include <iostream>

namespace texture_shallow{

void map_exact(std::vector<double>& map_z, double f, double vc){
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

    // Calculate ap taking circle center at y = 0
    //
    // Calculate whether the angles of the cutting tool intercept the surface
    // first or are stopped by the feed rate.
    //
    // Check the intersection with the feed rate edges. If they're greater than
    // zero, the cut line is narrower than the feed line, so max height is
    // zero. Otherwise, it's the height at the intersection
    s1 = -yc;
    s2 = f-yc;
    if(s1 > y1 && s2 < y2){
        // Width falls entirely within tool radius.
        yc = f/2;
        ap = -std::sqrt(r*r - yc*yc) + r;
    } else if(s1 > y1){
        ap = tan2*(f-yc) + b2off;
    } else if(s2 < y2){
        ap = tan1*yc + b1off;
    } else {
        // Width does not intersect the radius at all.
        ap = tan1*yc + b1off;
    }


    //max_z = Az;// + Az_uet;
    max_z = 0;
    min_z = -(ap + Az);

    const double delta_uet = vc/f_uet;

    const double over_z = Az + Az_uet*(1.0 - std::sqrt(1.0 - (delta_uet*delta_uet)/(4*Ax_uet*Ax_uet)));

    // Value of z for y = 0 (global)
    const double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap - over_z + tan1*(yc+y1);
    const double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - over_z - tan2*(yc+y2);
    
    #pragma omp parallel
    {
        std::vector<double> newz(tex_width);
        double prev_mult = -1;
        #pragma omp for
        for(size_t Y = 0; Y < tex_height; ++Y){
            const double y = static_cast<double>(Y);

            // Calculate current row
            // Apply `abs()` to prevent it from becoming less than zero
            // when close to zero
            const double mult = std::floor(y  / f);

            // As `mult` is always exact, we can use this optimization, meaning
            // the oscillations/ellipses along a row are calculated only once
            // per `mult`.
            if(mult > prev_mult){
                prev_mult = mult;

                // Phase differences
                const double perimeter = 2*M_PI*cylinder_radius;

                const double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
                const double phiz_cur = phase_diff_z - std::floor((phase_diff_z)/(2*M_PI))*2*M_PI;

                const double xoffset_uet = mult*perimeter;

                for(size_t X = 0; X < tex_width; ++X){
                    const double x = static_cast<double>(X);
                    // Random oscillation
                    const double oscillation = Az*std::sin(2*M_PI*fz*x*dimx/vc + phiz_cur);

                    // Ultrasonic turning effects
                    const double xcirc = x + xoffset_uet;
                    const double mult_uet = std::floor(xcirc / delta_uet);
                    const double x_uet = xcirc - (mult_uet + 0.5)*delta_uet;
                    const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)));

                    newz[X] = oscillation + uet_effect;
                }
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

                // Calculate final depth
                const double end_z = newz[X] + shape_z;

                // Write results
                z = end_z;
            }
        }
    }
}

void map(std::vector<double>& map_z, double f, double vc){
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

    // Calculate ap taking circle center at y = 0
    //
    // Calculate whether the angles of the cutting tool intercept the surface
    // first or are stopped by the feed rate.
    //
    // Check the intersection with the feed rate edges. If they're greater than
    // zero, the cut line is narrower than the feed line, so max height is
    // zero. Otherwise, it's the height at the intersection
    s1 = -yc;
    s2 = f-yc;
    if(s1 > y1 && s2 < y2){
        // Width falls entirely within tool radius.
        yc = f/2;
        ap = -std::sqrt(r*r - yc*yc) + r;
    } else if(s1 > y1){
        ap = tan2*(f-yc) + b2off;
    } else if(s2 < y2){
        ap = tan1*yc + b1off;
    } else {
        // Width does not intersect the radius at all.
        ap = tan1*yc + b1off;
    }


    // If it's greater than zero, max_z must be zero
    // max_z = Az;// + Az_uet;
    max_z = 0;
    min_z = -(ap + Az);

    const double delta_uet = vc/f_uet;

    const double over_z = Az + Az_uet*(1.0 - std::sqrt(1.0 - (delta_uet*delta_uet)/(4*Ax_uet*Ax_uet)));

    // Value of z for y = 0 (global)
    const double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap - over_z + tan1*(yc+y1);
    const double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - over_z - tan2*(yc+y2);
    
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
            const double perimeter = 2*M_PI*cylinder_radius;

            const double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
            const double phiz_cur = phase_diff_z - smooth::floor((phase_diff_z)/(2*M_PI))*2*M_PI;

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
                // Random oscillation
                const double oscillation = Az*std::sin(2*M_PI*fz*x*dimx/vc + phiz_cur);

                // Ultrasonic turning effects
                const double xcirc = x + xoffset_uet;
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

                // Calculate final depth
                const double end_z = newz[X] + shape_z;

                // Write results
                z = end_z;
            }
        }
    }
}

void dzdvc(double f, double vc, std::vector<double>& dzdvc){
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

    s1 = -yc;
    s2 = f-yc;
    if(s1 > y1 && s2 < y2){
        // Width falls entirely within tool radius.
        yc = f/2;
        ap = -std::sqrt(r*r - yc*yc) + r;
    } else if(s1 > y1){
        ap = tan2*(f-yc) + b2off;
    } else if(s2 < y2){
        ap = tan1*yc + b1off;
    } else {
        // Width does not intersect the radius at all.
        ap = tan1*yc + b1off;
    }

    dmax_zdvc = 0;
    dmin_zdvc = 0;

    const double delta_uet = vc/f_uet;
    const double dduetdvc = 1.0/f_uet;

    // const double over_z = Az + Az_uet*(1.0 - std::sqrt(1.0 - (delta_uet*delta_uet)/(4*Ax_uet_vc*Ax_uet_vc)));
    const double dover_z = -0.5*(Az_uet/std::sqrt(1.0 - (delta_uet*delta_uet)/(4*Ax_uet*Ax_uet)))
                    *(-2)*(delta_uet/(2*Ax_uet))*(dduetdvc/(2*Ax_uet));

    // const double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap - over_z + tan1*(yc+y1);
    // const double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - over_z - tan2*(yc+y2);

    const double dlrdf1 = -dover_z;
    const double dlrdf2 = -dover_z;
    #pragma omp parallel
    {
        std::vector<double> newz(tex_width);
        std::vector<double> dznewdvc(tex_width);
        #pragma omp for
        for(size_t Y = 0; Y < tex_height; ++Y){
            const double y = static_cast<double>(Y);

            const double mult = smooth::abs(smooth::floor( (y + YOFF) / f));
            // const double dmult = 0;

            const double perimeter = 2*M_PI*cylinder_radius;

            const double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
            const double dpdz = -2*M_PI*fz*mult*perimeter*dimx/(vc*vc);
            const double phiz_cur = phase_diff_z - smooth::floor((phase_diff_z)/(2*M_PI))*2*M_PI;
            const double dphiz_cur = dpdz - smooth::floor_deriv((phase_diff_z)/(2*M_PI))*dpdz;

            const double xoffset_uet = mult*perimeter;

            for(size_t X = 0; X < tex_width; ++X){
                const double x = static_cast<double>(X);

                const double oscillation = Az*std::sin(2*M_PI*fz*x*dimx/vc + phiz_cur);
                const double doscilldvc = Az*std::cos(2*M_PI*fz*x*dimx/vc + phiz_cur)*(-2*M_PI*fz*x*dimx/(vc*vc) + dphiz_cur);

                const double xcirc = x + xoffset_uet;

                const double mult_uet = smooth::abs(smooth::floor(xcirc / delta_uet));
                const double dmult_uet = smooth::abs_deriv(smooth::floor(xcirc / delta_uet))*smooth::floor_deriv(xcirc / delta_uet)*((-xcirc)/(delta_uet*delta_uet))*dduetdvc;

                const double x_uet = xcirc - (mult_uet + 0.5)*delta_uet;
                const double dx_uet = -(mult_uet*dduetdvc + dmult_uet*delta_uet) - 0.5*dduetdvc;

                const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)));
                const double duet = -0.5*(Az_uet/std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)))*(-2)*(x_uet/Ax_uet)*(dx_uet/Ax_uet);

                newz[X] = oscillation + uet_effect;
                dznewdvc[X] = doscilldvc + duet;

                // double& dz = dzdvc[X + Y*tex_width];

                // dz = dznewdvc[X];
            }

            double dshape_z;
            if(y <= y1 + mult*f + yc){
                dshape_z = dlrdf1;
            } else if(y <= y2 + mult*f + yc){
                dshape_z = 0;
            } else {
                dshape_z = dlrdf2;
            }

            for(size_t X = 0; X < tex_width; ++X){
                double& dz = dzdvc[X + Y*tex_width];

                const double dend_z = dznewdvc[X] + dshape_z;

                dz = dend_z;
            }
        }
    }
}

void dzdf(double f, double vc, std::vector<double>& dzdf){
    using namespace param;
    using param::y1;

    vc *= 1e6/60.0;

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

    double dap = 0;
    s1 = -yc;
    s2 = f-yc;
    if(s1 > y1 && s2 < y2){
        yc = f/2;
        dyc = 0.5;
        dap = yc*dyc/std::sqrt(r*r - yc*yc);
    } else if(s1 > y1){
        dap = tan2;
    } else if(s2 < y2){
        dap = tan1;
    } else {
        dap = tan1;
    }

    const double dlrdf1 = -dap + tan1*dyc;
    const double dlrdf2 = -dap - tan2*dyc;

    dmax_zdf = 0;
    dmin_zdf = -dap;

    const double delta_uet = vc/f_uet;

    //const double dduet = 0;
    #pragma omp parallel
    {
        std::vector<double> newz(tex_width);
        std::vector<double> dznewdf(tex_width);
        #pragma omp for
        for(size_t Y = 0; Y < tex_height; ++Y){
            const double y = static_cast<double>(Y);

            const double mult = smooth::abs(smooth::floor( (y + YOFF) / f));
            const double dmult = smooth::abs_deriv(smooth::floor( (y + YOFF) / f))*smooth::floor_deriv((y + YOFF) / f)*(-(y / (f*f)));

            const double perimeter = 2*M_PI*cylinder_radius;

            const double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
            const double dpdz = 2*M_PI*fz*dmult*perimeter*dimx/vc;
            const double phiz_cur = phase_diff_z - smooth::floor((phase_diff_z)/(2*M_PI))*2*M_PI;
            const double dphiz_cur = dpdz - smooth::floor_deriv((phase_diff_z)/(2*M_PI))*dpdz;

            const double xoffset_uet = mult*perimeter;
            const double dxoff_uet = dmult*perimeter;

            for(size_t X = 0; X < tex_width; ++X){
                const double x = static_cast<double>(X);

                const double oscillation = Az*std::sin(2*M_PI*fz*x*dimx/vc + phiz_cur);
                const double doscilldf = Az*std::cos(2*M_PI*fz*x*dimx/vc + phiz_cur)*dphiz_cur;

                const double xcirc = x + xoffset_uet;
                const double dxcirc = dxoff_uet;

                const double mult_uet = smooth::abs(smooth::floor(xcirc / delta_uet));
                const double dmult_uet = smooth::abs_deriv(smooth::floor(xcirc / delta_uet))*smooth::floor_deriv(xcirc / delta_uet)*dxcirc/delta_uet;

                // const double x_uet = (xcirc - mult_uet*delta_uet)/delta_uet - 0.5;
                // const double x_uet = xcirc/delta_uet - mult_uet - 0.5;
                // const double dx_uet = dxcirc/delta_uet - dmult_uet;

                // const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)));
                // const double duet = -0.5*(Az_uet/std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)))*(-2)*(x_uet/Ax_uet)*(dx_uet/Ax_uet);

                const double x_uet = xcirc - (mult_uet + 0.5)*delta_uet;
                const double dx_uet = dxcirc - dmult_uet*delta_uet;

                const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)));
                const double duet = -0.5*(Az_uet/std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)))*(-2)*(x_uet/Ax_uet)*(dx_uet/Ax_uet);

                newz[X] = oscillation + uet_effect;
                dznewdf[X] = doscilldf + duet;
            }

            double dshape_z;
            if(y <= y1 + mult*f + yc){
                dshape_z =  tan1*(mult + dmult*f) + dlrdf1;
            } else if(y <= y2 + mult*f + yc){
                const double yy = y - mult*f - yc;
                const double dyy = -mult - dmult*f - dyc;

                const double circ = std::sqrt(r*r - yy*yy);
                const double dcirc = -yy*dyy/circ;

                dshape_z = -dcirc - dap;
            } else {
                dshape_z = -tan2*(mult + dmult*f) + dlrdf2;
            }

            for(size_t X = 0; X < tex_width; ++X){
                double& dz = dzdf[X + Y*tex_width];

                const double dend_z = dznewdf[X] + dshape_z;

                dz = dend_z;
            }
        }
    }
}

}
