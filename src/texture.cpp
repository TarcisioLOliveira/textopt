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

    vc *= 1e6/60;

    const double tan1 = std::tan(alpha1);
    const double tan2 = std::tan(alpha2);

    // Line parameters for tool slopes.
    const double y1 = -std::sqrt((tan1*r*tan1*r)/(tan1*tan1+1));
    const double y2 =  std::sqrt((tan2*r*tan2*r)/(tan2*tan2+1));

    // Where the slope lines intersect, so that the height at the edges is the
    // same for both
    const double yc = tan2*f/(tan1+tan2);

    // Value of z for y = 0 (with circle center at y = 0).
    const double b1off = -std::sqrt(r*r - y1*y1) + r - ap + tan1*y1;
    const double b2off = -std::sqrt(r*r - y2*y2) + r - ap - tan2*y2;

    // Position of the intersection of the two lines relative to the circle's
    // center
    const double ycoff = (b1off - b2off)/(tan1 + tan2);

    // Value of z for y = 0 (global)
    const double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + tan1*(yc-ycoff+y1);
    const double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - tan2*(yc-ycoff+y2);

    // Calculate whether the angles of the cutting tool intercept the surface
    // first or are stopped by the feed rate.
    //
    // Check the intersection with the feed rate edges. If they're greater than
    // zero, the cut line is narrower than the feed line, so max height is
    // zero. Otherwise, it's the height at the intersection
    //
    // Both intersections always have the same height, so we only need to check
    // one of them
    double intersec = 0;
    if(yc - ycoff + y1 > 0){
        // If the intersection between line 1 and the circle is at y > 0, then
        // the intersection is equal to the line's constant.
        intersec = line_root1;
    } else {
        // Otherwise, check where it intersects the circle
        const double yi = 0 - (yc - ycoff);
        intersec = -std::sqrt(r*r - yi*yi) + r - ap;
    }
    // If it's greater than zero, max_z must be zero
    max_z = smooth::min({0, intersec}) + Az;// + Az_uet;
    min_z = -(ap + Az);

    const double delta_uet = vc/f_uet;
    
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

            const double phase_diff_x = 2*M_PI*fx*mult*perimeter*dimx/vc + phix;
            const double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
            const double phix_cur = phase_diff_x - smooth::floor((phase_diff_x)/(2*M_PI))*2*M_PI;
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
                const double newx = x + Ax*std::sin(2*M_PI*fx*x*dimx/vc + phix_cur);
                const double oscillation = Az*std::sin(2*M_PI*fz*newx*dimx/vc + phiz_cur);

                // Ultrasonic turning effects
                const double xcirc = x + xoffset_uet;
                const double mult_uet = smooth::abs(smooth::floor(xcirc / delta_uet));
                // const double x_uet = (xcirc - mult_uet*delta_uet)/delta_uet - 0.5;
                const double x_uet = xcirc/delta_uet - mult_uet - 0.5;
                const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)));

                newz[X] = oscillation + uet_effect;
            }

            // Tool shape
            double shape_z;
            if(y <= y1 + mult*f + (yc - ycoff)){
                shape_z = -tan1*(y - mult*f) + line_root1;
            } else if(y <= y2 + mult*f + (yc - ycoff)){
                const double yy = y - mult*f - (yc - ycoff);
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

    vc *= 1e6/60;

    const double tan1 = std::tan(alpha1);
    const double tan2 = std::tan(alpha2);

    const double y1 = -std::sqrt((tan1*r*tan1*r)/(tan1*tan1+1));
    const double y2 =  std::sqrt((tan2*r*tan2*r)/(tan2*tan2+1));

    const double yc = tan2*f/(tan1+tan2);

    const double b1off = -std::sqrt(r*r - y1*y1) + r - ap + tan1*y1;
    const double b2off = -std::sqrt(r*r - y2*y2) + r - ap - tan2*y2;

    const double ycoff = (b1off - b2off)/(tan1 + tan2);

    const double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + tan1*(yc-ycoff+y1);
    const double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - tan2*(yc-ycoff+y2);
    
    dmax_zdvc = 0;

    const double delta_uet = vc/f_uet;
    const double dduetdvc = 1.0/f_uet;

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

                // const double x_uet = (xcirc - mult_uet*delta_uet)/delta_uet - 0.5;
                const double x_uet = xcirc/delta_uet - mult_uet - 0.5;
                const double dx_uet = -(xcirc/(delta_uet*delta_uet))*dduetdvc - dmult_uet;

                const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)));
                const double duet = -0.5*(Az_uet/std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)))*(-2)*(x_uet/Ax_uet)*(dx_uet/Ax_uet);

                newz[X] = oscillation + uet_effect;
                dznewdvc[X] = doscilldvc + duet;
            }
                
            double shape_z;
            if(y <= y1 + mult*f + (yc - ycoff)){
                shape_z = -tan1*(y - mult*f) + line_root1;
            } else if(y <= y2 + mult*f + (yc - ycoff)){
                const double yy = y - mult*f - (yc - ycoff);
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

    vc *= 1e6/60;

    const double tan1 = std::tan(alpha1);
    const double tan2 = std::tan(alpha2);

    const double y1 = -std::sqrt((tan1*r*tan1*r)/(tan1*tan1+1));
    const double y2 =  std::sqrt((tan2*r*tan2*r)/(tan2*tan2+1));

    const double yc = tan2*f/(tan1+tan2);

    const double b1off = -std::sqrt(r*r - y1*y1) + r - ap + tan1*y1;
    const double b2off = -std::sqrt(r*r - y2*y2) + r - ap - tan2*y2;

    const double db1off = -1;
    const double db2off = -1;

    const double ycoff = (b1off - b2off)/(tan1 + tan2);
    const double dycoff = (db1off - db2off)/(tan1 + tan2);

    const double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + tan1*(yc-ycoff+y1);
    const double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - tan2*(yc-ycoff+y2);

    const double dlrdap1 = -1 + tan1*(-dycoff);
    const double dlrdap2 = -1 - tan2*(-dycoff);

    double intersec = 0;
    double di = 0;
    if(yc - ycoff + y1 > 0){
        intersec = line_root1;
        di = dlrdap1;
    } else {
        const double yi = 0 - (yc - ycoff);
        const double dyi = dycoff;

        const double circ = std::sqrt(r*r - yi*yi);
        const double dcirc = -yi*dyi/circ;

        intersec = -circ + r - ap;
        di = -dcirc - 1;
    }
    dmax_zdap = smooth::min_deriv({0, intersec}, 1)*di;

    const double delta_uet = vc/f_uet;

    #pragma omp parallel
    {
        std::vector<double> newz(tex_width);
        #pragma omp for
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
                // const double x_uet = (xcirc - mult_uet*delta_uet)/delta_uet - 0.5;
                const double x_uet = xcirc/delta_uet - mult_uet - 0.5;
                const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)));

                newz[X] = oscillation + uet_effect;
            }

            double shape_z;
            double dznewdap;
            if(y <= y1 + mult*f + (yc - ycoff)){
                shape_z = -tan1*(y - mult*f) + line_root1;
                dznewdap = dlrdap1;
            } else if(y <= y2 + mult*f + (yc - ycoff)){
                const double yy = y - mult*f - (yc - ycoff);
                const double dyy = dycoff;

                const double circ = std::sqrt(r*r - yy*yy);
                const double dcirc = -yy*dyy/circ;

                shape_z = -circ + r - ap;
                dznewdap = -dcirc - 1;
            } else {
                shape_z = tan2*(y - mult*f) + line_root2;
                dznewdap = dlrdap2;
            }

            for(size_t X = 0; X < tex_width; ++X){
                const double oz = orig_z[X + Y*tex_width];
                double& dz = dzdap[X + Y*tex_width];

                const double end_z = newz[X] + shape_z;

                dz = smooth::min_deriv({oz, end_z}, 1)*dznewdap;
            }
        }
    }
}

void dzdf(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdf){
    using namespace param;

    vc *= 1e6/60;

    const double tan1 = std::tan(alpha1);
    const double tan2 = std::tan(alpha2);

    const double y1 = -std::sqrt((tan1*r*tan1*r)/(tan1*tan1+1));
    const double y2 =  std::sqrt((tan2*r*tan2*r)/(tan2*tan2+1));

    const double yc = tan2*f/(tan1+tan2);
    const double dyc = tan2/(tan1+tan2);

    const double b1off = -std::sqrt(r*r - y1*y1) + r - ap + tan1*y1;
    const double b2off = -std::sqrt(r*r - y2*y2) + r - ap - tan2*y2;

    const double ycoff = (b1off - b2off)/(tan1 + tan2);

    const double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + tan1*(yc-ycoff+y1);
    const double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - tan2*(yc-ycoff+y2);

    const double dlrdf1 =  tan1*dyc;
    const double dlrdf2 = -tan2*dyc;

    double intersec = 0;
    double di = 0;
    if(yc - ycoff + y1 > 0){
        intersec = line_root1;
        di = dlrdf1;
    } else {
        const double yi = 0 - (yc - ycoff);
        const double dyi = -dyc;

        const double circ = std::sqrt(r*r - yi*yi);
        const double dcirc = -yi*dyi/circ;

        intersec = -circ + r - ap;
        di = -dcirc;
    }
    dmax_zdf = smooth::min_deriv({0, intersec}, 1)*di;

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

                // const double x_uet = (xcirc - mult_uet*delta_uet)/delta_uet - 0.5;
                const double x_uet = xcirc/delta_uet - mult_uet - 0.5;
                const double dx_uet = dxcirc/delta_uet - dmult_uet;

                const double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)));
                const double duet = -0.5*(Az_uet/std::sqrt(1.0 - (x_uet*x_uet)/(Ax_uet*Ax_uet)))*(-2)*(x_uet/Ax_uet)*(dx_uet/Ax_uet);

                newz[X] = oscillation + uet_effect;
                dznewdf[X] = doscilldf + duet;
            }

            double shape_z;
            double dshape_z;
            if(y <= y1 + mult*f + (yc - ycoff)){
                shape_z = -tan1*(y - mult*f) + line_root1;
                dshape_z =  tan1*(mult + dmult*f) + dlrdf1;
            } else if(y <= y2 + mult*f + (yc - ycoff)){
                const double yy = y - mult*f - (yc - ycoff);
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
