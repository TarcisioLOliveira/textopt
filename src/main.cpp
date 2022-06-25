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
#include <random>
#include <cmath>
#include <iostream>
#include <initializer_list>
#include <numeric>
#include <MMASolver.hpp>
#include <algorithm>

size_t tex_width = 300;
size_t tex_height = 300;

const double MULT = -1; // Exponent for smooth_min()
const double FLOORD = 1e-10; // smooth_floor() precision
const double ABS_EPS = 1e-15; // smooth_abs() precision
const double YOFF = 0.001; // Workaround over smooth_floor() being off by -0.5 for integers

double dim_scale = 1e6; // m to um
double dim = 1;
double dimx = dim;
double dimy = dim;
double dimz = dim;

double alpha1 = 60*M_PI/180;
double alpha2 = 60*M_PI/180;
double r = 20/dim; // um
double max_z = 0; // Used to calculate Sa, so uses smooth_min
double min_z = 0; // Used only for texture display, so uses std::min()
double dmax_zdf = 0; 
double dmax_zdap = 0; 
double dmax_zdvc = 0; 

// "Random" oscillation
// Models height variation along the tool path

// Amplitude [dim], e.g. [um]
double Ax = 0.5;
double Az = 0.1;
// Frequency [Hz]
double fx = 19*1e5;
double fz = 51*1e5;
// Phase [rad]
double phix = 20*M_PI/180;
double phiz = 0;

double cylinder_radius = 6*1e-3*dim_scale; // [um]

// Ultrasonic elliptical turning
double f_uet = 20000; // [Hz]
double Ax_uet = 1; // [um]
double Az_uet = 2; // [um]

struct Point{
    double x, y, z;
};

struct Vector{
    double x, y, z;
};

enum Colorscheme{
    GRAYSCALE,
    HSV
};

double triangle_area(std::array<Point, 3> p){
    Vector v1{p[1].x - p[0].x, p[1].y - p[0].y, p[1].z - p[0].z};
    Vector v2{p[2].x - p[0].x, p[2].y - p[0].y, p[2].z - p[0].z};

    double A = 0.5*std::sqrt(
        std::pow(v1.y*v2.z - v2.y*v1.z, 2) +
        std::pow(v1.x*v2.z - v2.x*v1.z, 2) +
        std::pow(v1.y*v2.x - v2.y*v1.x, 2)
    );

    return A;
}


double triangle_area_deriv(std::array<Point, 3> p, std::array<double, 3> dz){
    Vector v1{p[1].x - p[0].x, p[1].y - p[0].y, p[1].z - p[0].z};
    Vector v2{p[2].x - p[0].x, p[2].y - p[0].y, p[2].z - p[0].z};
    double dv1 = dz[1] - dz[0];
    double dv2 = dz[2] - dz[0];

    double A = (0.25/std::sqrt(
        std::pow(v1.y*v2.z - v2.y*v1.z, 2) +
        std::pow(v1.x*v2.z - v2.x*v1.z, 2) +
        std::pow(v1.y*v2.x - v2.y*v1.x, 2)
    ))*(2*(v1.y*v2.z - v2.y*v1.z)*(v1.y*dv2 - v2.y*dv1) +
        2*(v1.x*v2.z - v2.x*v1.z)*(v1.x*dv2 - v2.x*dv1) +
        0
    );

    return A;
}

double smooth_min(std::initializer_list<double> x){
    double frac_top = 0;
    double frac_bot = 0;
    for(double xi : x){
        double exp = std::exp(MULT*xi);
        frac_top += xi*exp;
        frac_bot += exp;
    }

    return frac_top/frac_bot;
};

double smooth_min_deriv(std::initializer_list<double> x, double xx){
    double frac_top = 0;
    double frac_bot = 0;
    for(double xi : x){
        double exp = std::exp(MULT*xi);
        frac_top += xi*exp;
        frac_bot += exp;
    }

    double sm = frac_top/frac_bot;

    return (std::exp(MULT*xx)/frac_bot)*(1+MULT*(xx-sm));
};

double smooth_floor(double v){
    // https://math.stackexchange.com/questions/2746958/smooth-floor-function
    //
    // δ = 0.01;
    // trg[x_] := 1 - 2 ArcCos[(1 - δ) Sin[2 π x]]/π;
    // sqr[x_] := 2 ArcTan[Sin[2 π x]/δ]/π;
    // swt[x_] := (1 + trg[(2 x - 1)/4] sqr[x/2])/2;
    // flr[x_] := x-swt[x]

    double p1x = (1.0 - FLOORD)*std::sin(2*M_PI*((2*v-1)/4));
    double p2x = std::sin(2*M_PI*(v/2))/FLOORD;

    double p1 = 1 - 2*std::acos(p1x)/M_PI;
    double p2 = 2 * std::atan(p2x)/M_PI;
    double sawtooth = (1 + p1 * p2)/2;

    double result = v - sawtooth;

    return result;
}

double smooth_floor_deriv(double v){
    double p1x = (1.0 - FLOORD)*std::sin(2*M_PI*((2*v-1)/4));
    double dp1x = (1.0 - FLOORD)*std::cos(2*M_PI*((2*v-1)/4))*2*M_PI*2/4;

    double p2x = std::sin(2*M_PI*(v/2))/FLOORD;
    double dp2x = (std::cos(2*M_PI*(v/2))/FLOORD)*2*M_PI/2;

    double p1 = 1 - 2*std::acos(p1x)/M_PI;
    double dp1 = (2.0/(std::sqrt(1 - p1x*p1x)*M_PI))*dp1x;
    double p2 = 2 * std::atan(p2x)/M_PI;
    double dp2 = (2.0/((p2x*p2x + 1)*M_PI))*dp2x;
    double dsawtooth = (1 + dp1 * p2 + p1 * dp2)/2;

    double result = 1.0 - dsawtooth;
    
    return result;
}

/**
 * Used mostly to fix smooth_floor() becoming negative when close to zero.
 *
 * Not really smooth, yes, but it's actually more of a workaround. It works
 * better this way, as smooth_abs_deriv() was introducing instability into
 * the optimization process.
 */
double smooth_abs(double v){
    return std::abs(v);
    // return std::sqrt(v*v+ABS_EPS);
}

double smooth_abs_deriv(double v){
    // return v/std::sqrt(v*v+ABS_EPS);
    return 1;
}

void texture_map(double*& map_z, double* orig_z, double f, double ap, double vc){
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
        intersec1 = smooth_min({0, line_root1});
    } else {
        // Check for the circle
        intersec1 = smooth_min({0, -std::sqrt(r*r - f*f/4) + r - ap});
    }
    if(y2 + f/2 < f){
        // If the radius is smaller than the feed rate...
        intersec2 = smooth_min({0, line_root2 + std::tan(alpha2)*f});
    } else {
        intersec2 = smooth_min({0, -std::sqrt(r*r - f*f/4) + r - ap});
    }
    // Check if one of them is greater than zero
    double max_intersec = -smooth_min({-intersec1, -intersec2});
    // If one is greater than zero, max_z must be zero. Otherwise, it's
    // the lesser one.
    max_z = smooth_min({0, max_intersec});
    min_z = 0;

    double delta_uet = vc/f_uet;
    for(size_t i = 0; i < tex_width*tex_height; ++i){
        double x = i % tex_width;
        double y = std::floor(((double)i) / tex_width);
            // Get pixels
            double& z = map_z[i];
            double oz = orig_z[i];

            // Calculate current row
            // Apply `abs()` to prevent it from becoming less than zero
            // when close to zero
            double mult = smooth_abs(smooth_floor((y + YOFF) / f));

            // Phase differences
            double perimeter = 2*M_PI*cylinder_radius;

            double phase_diff_x = 2*M_PI*fx*mult*perimeter*dimx/vc + phix;
            double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
            double phix_cur = phase_diff_x - smooth_floor((phase_diff_x)/(2*M_PI))*2*M_PI;
            double phiz_cur = phase_diff_z - smooth_floor((phase_diff_z)/(2*M_PI))*2*M_PI;

            double xoffset_uet = mult*perimeter;

            // Random oscillation
            double newx = x + Ax*std::sin(2*M_PI*fx*x*dimx/vc + phix_cur);
            double oscillation = Az*std::sin(2*M_PI*fz*newx*dimx/vc + phiz_cur);

            // Ultrasonic turning effects
            double xcirc = x + xoffset_uet;
            double mult_uet = smooth_abs(smooth_floor(xcirc / delta_uet));
            double x_uet = (xcirc - mult_uet*delta_uet)*Ax_uet/delta_uet - Ax_uet/2;
            double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)));

            // Tool shape
            double newz = 0;
            if(y <= y1 + mult*f + f/2){
                newz = -std::tan(alpha1)*(y - mult*f) + line_root1;
            } else if(y <= y2 + mult*f + f/2){
                newz = -std::sqrt(r*r - std::pow(y - mult*f - f/2, 2)) + r - ap;
            } else {
                newz = std::tan(alpha2)*(y - mult*f) + line_root2;
            }

            // Write results
            newz += oscillation + uet_effect;
            z = smooth_min({oz, newz});

            min_z = std::min(min_z, z);
    }
}

void dzdvc(double* orig_z, double f, double ap, double vc, double*& dzdvc){
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
        intersec1 = smooth_min({0, line_root1});
    } else {
        intersec1 = smooth_min({0, -std::sqrt(r*r - f*f/4) + r - ap});
    }
    if(y2 + f/2 < f){
        intersec2 = smooth_min({0, line_root2 + std::tan(alpha2)*f});
    } else {
        intersec2 = smooth_min({0, -std::sqrt(r*r - f*f/4) + r - ap});
    }
    double max_intersec = -smooth_min({-intersec1, -intersec2});

    double dmi = 0;//-smooth_min_deriv({-intersec1, -intersec2}, -intersec1)*(-di1) - smooth_min_deriv({-intersec1, -intersec2}, -intersec2)*(-di2);
    dmax_zdvc = 0;//smooth_min_deriv({0, max_intersec}, max_intersec)*dmi;

    double delta_uet = vc/f_uet;
    double dduetdvc = 1.0/f_uet;
    for(size_t i = 0; i < tex_width*tex_height; ++i){
        double x = i % tex_width;
        double y = std::floor(((double)i) / tex_width);
            double oz = orig_z[i];
            double& dz = dzdvc[i];

            double mult = smooth_abs(smooth_floor((double) (y + YOFF) / f));
            double dmult = 0;

            double perimeter = 2*M_PI*cylinder_radius;

            double phase_diff_x = 2*M_PI*fx*mult*perimeter*dimx/vc + phix;
            double dpdx = -2*M_PI*fx*mult*perimeter*dimx/(vc*vc);
            double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
            double dpdz = -2*M_PI*fz*mult*perimeter*dimx/(vc*vc);
            double phix_cur = phase_diff_x - smooth_floor((phase_diff_x)/(2*M_PI))*2*M_PI;
            double dphix_cur = dpdx - smooth_floor_deriv((phase_diff_x)/(2*M_PI))*dpdx;
            double phiz_cur = phase_diff_z - smooth_floor((phase_diff_z)/(2*M_PI))*2*M_PI;
            double dphiz_cur = dpdz - smooth_floor_deriv((phase_diff_z)/(2*M_PI))*dpdz;

            double xoffset_uet = mult*perimeter;
            double dxoff_uet = dmult*perimeter;

            double newx = x + Ax*std::sin(2*M_PI*fx*x*dimx/vc + phix_cur);
            double dnewxdvc = Ax*std::cos(2*M_PI*fx*x*dimx/vc + phix_cur)*((-2)*M_PI*fx*x*dimx/(vc*vc) + dphix_cur);

            double oscillation = Az*std::sin(2*M_PI*fz*newx*dimx/vc + phiz_cur);
            double doscilldvc = Az*std::cos(2*M_PI*fz*newx*dimx/vc + phiz_cur)*(2*M_PI*fz*(dnewxdvc*dimx*vc - newx*dimx)/(vc*vc) + dphiz_cur);

            double xcirc = x + xoffset_uet;
            double dxcirc = dxoff_uet;

            double mult_uet = smooth_abs(smooth_floor(xcirc / delta_uet));
            double dmult_uet = smooth_abs_deriv(smooth_floor(xcirc / delta_uet))*smooth_floor_deriv(xcirc / delta_uet)*(dxcirc*delta_uet - dduetdvc*xcirc)/(delta_uet*delta_uet);

            double x_uet = (xcirc - mult_uet*delta_uet)*Ax_uet/delta_uet - Ax_uet/2;
            double dx_uet = Ax_uet*(dxcirc*delta_uet - dduetdvc*xcirc)/(delta_uet*delta_uet) - Ax_uet*dmult_uet;

            double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)));
            double duet = -0.5*(Az_uet/std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)))*(-2)*(x_uet/Ax_uet)*(dx_uet/Ax_uet);

            double newz = 0;
            double dznewdvc = 0;
            if(y <= y1 + mult*f + f/2){
                newz = -std::tan(alpha1)*(y - mult*f) + line_root1;
                dznewdvc = dlrdvc1;
            } else if(y <= y2 + mult*f + f/2){
                newz = -std::sqrt(r*r - std::pow(y - mult*f - f/2, 2)) + r - ap;
                dznewdvc = 0;
            } else {
                newz =  std::tan(alpha2)*(y - mult*f) + line_root2;
                dznewdvc = dlrdvc2;
            }
            newz += oscillation + uet_effect;
            dznewdvc += doscilldvc + duet;
            dz = smooth_min_deriv({oz, newz}, newz)*dznewdvc;
    }
}

void dzdap(double* orig_z, double f, double ap, double vc, double*& dzdap){
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
        intersec1 = smooth_min({0, line_root1});
        di1 = smooth_min_deriv({0, line_root1}, line_root1)*dlrdap1;
    } else {
        intersec1 = smooth_min({0, -std::sqrt(r*r - f*f/4) + r - ap});
        di1 = smooth_min_deriv({0, -std::sqrt(r*r - f*f/4) + r - ap},-std::sqrt(r*r - f*f/4) + r - ap)*(-1);
    }
    if(y2 + f/2 < f){
        intersec2 = smooth_min({0, line_root2 + std::tan(alpha2)*f});
        di2 = smooth_min_deriv({0, line_root2 + std::tan(alpha2)*f}, line_root2 + std::tan(alpha2)*f)*dlrdap2;
    } else {
        intersec2 = smooth_min({0, -std::sqrt(r*r - f*f/4) + r - ap});
        di2 = smooth_min_deriv({0, -std::sqrt(r*r - f*f/4) + r - ap},-std::sqrt(r*r - f*f/4) + r - ap)*(-1);
    }
    double max_intersec = -smooth_min({-intersec1, -intersec2});

    double dmi = -smooth_min_deriv({-intersec1, -intersec2}, -intersec1)*(-di1) - smooth_min_deriv({-intersec1, -intersec2}, -intersec2)*(-di2);
    dmax_zdap = smooth_min_deriv({0, max_intersec}, max_intersec)*dmi;

    double delta_uet = vc/f_uet;
    for(size_t i = 0; i < tex_width*tex_height; ++i){
        double x = i % tex_width;
        double y = std::floor(((double)i) / tex_width);
            double oz = orig_z[i];
            double& dz = dzdap[i];

            double mult = smooth_abs(smooth_floor((double) (y + YOFF) / f));

            double perimeter = 2*M_PI*cylinder_radius;

            double phase_diff_x = 2*M_PI*fx*mult*perimeter*dimx/vc + phix;
            double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
            double phix_cur = phase_diff_x - smooth_floor((phase_diff_x)/(2*M_PI))*2*M_PI;
            double phiz_cur = phase_diff_z - smooth_floor((phase_diff_z)/(2*M_PI))*2*M_PI;

            double xoffset_uet = mult*perimeter;

            double newx = x + Ax*std::sin(2*M_PI*fx*x*dimx/vc + phix_cur);
            double oscillation = Az*std::sin(2*M_PI*fz*newx*dimx/vc + phiz_cur);

            double xcirc = x + xoffset_uet;
            double mult_uet = smooth_abs(smooth_floor(xcirc / delta_uet));
            double x_uet = (xcirc - mult_uet*delta_uet)*Ax_uet/delta_uet - Ax_uet/2;
            double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)));

            double newz = 0;
            double dznewdap = 0;
            if(y <= y1 + mult*f + f/2){
                newz = -std::tan(alpha1)*(y - mult*f) + line_root1;
                dznewdap = dlrdap1;
            } else if(y <= y2 + mult*f + f/2){
                newz = -std::sqrt(r*r - std::pow((double)y - mult*f - f/2, 2)) + r - ap;
                dznewdap = -1;
            } else {
                newz =  std::tan(alpha2)*(y - mult*f) + line_root2;
                dznewdap = dlrdap2;
            }
            newz += oscillation + uet_effect;
            dz = smooth_min_deriv({oz, newz}, newz)*dznewdap;
    }
}

void dzdf(double* orig_z, double f, double ap, double vc, double*& dzdf){
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
        intersec1 = smooth_min({0, line_root1});
        di1 = smooth_min_deriv({0, line_root1}, line_root1)*dlrdf1;
    } else {
        intersec1 = smooth_min({0, -std::sqrt(r*r - f*f/4) + r - ap});
        di1 = smooth_min_deriv({0, -std::sqrt(r*r - f*f/4) + r - ap},-std::sqrt(r*r - f*f/4) + r - ap)*f/(4*std::sqrt(r*r-f*f/4));
    }
    if(y2 + f/2 < f){
        intersec2 = smooth_min({0, line_root2 + std::tan(alpha2)*f});
        di2 = smooth_min_deriv({0, line_root2 + std::tan(alpha2)*f}, line_root2 + std::tan(alpha2)*f)*(dlrdf2 + std::tan(alpha2));
    } else {
        intersec2 = smooth_min({0, -std::sqrt(r*r - f*f/4) + r - ap});
        di2 = smooth_min_deriv({0, -std::sqrt(r*r - f*f/4) + r - ap},-std::sqrt(r*r - f*f/4) + r - ap)*f/(4*std::sqrt(r*r-f*f/4));
    }
    double max_intersec = -smooth_min({-intersec1, -intersec2});

    double dmi = -smooth_min_deriv({-intersec1, -intersec2}, -intersec1)*(-di1) - smooth_min_deriv({-intersec1, -intersec2}, -intersec2)*(-di2);
    dmax_zdf = smooth_min_deriv({0, max_intersec}, max_intersec)*dmi;

    double delta_uet = vc/f_uet;
    double dduet = 0;
    for(size_t i = 0; i < tex_width*tex_height; ++i){
        double x = i % tex_width;
        double y = std::floor(((double)i) / tex_width);
            double oz = orig_z[i];
            double& dz = dzdf[i];

            double mult = smooth_abs(smooth_floor((double) (y + YOFF) / f));
            double dmult = smooth_abs_deriv(smooth_floor((double) (y + YOFF) / f))*smooth_floor_deriv((double)(y + YOFF) / f)*(-((double)y / (f*f)));

            double perimeter = 2*M_PI*cylinder_radius;

            double phase_diff_x = 2*M_PI*fx*mult*perimeter*dimx/vc + phix;
            double dpdx = 2*M_PI*fx*dmult*perimeter*dimx/vc;
            double phase_diff_z = 2*M_PI*fz*mult*perimeter*dimx/vc + phiz;
            double dpdz = 2*M_PI*fz*dmult*perimeter*dimx/vc;
            double phix_cur = phase_diff_x - smooth_floor((phase_diff_x)/(2*M_PI))*2*M_PI;
            double dphix_cur = dpdx - smooth_floor_deriv((phase_diff_x)/(2*M_PI))*dpdx;
            double phiz_cur = phase_diff_z - smooth_floor((phase_diff_z)/(2*M_PI))*2*M_PI;
            double dphiz_cur = dpdz - smooth_floor_deriv((phase_diff_z)/(2*M_PI))*dpdz;

            double xoffset_uet = mult*perimeter;
            double dxoff_uet = dmult*perimeter;

            double newx = x + Ax*std::sin(2*M_PI*fx*x*dimx/vc + phix_cur);
            double dnewxdf = Ax*std::cos(2*M_PI*fx*x*dimx/vc + phix_cur)*dphix_cur;

            double oscillation = Az*std::sin(2*M_PI*fz*newx*dimx/vc + phiz_cur);
            double doscilldf = Az*std::cos(2*M_PI*fz*newx*dimx/vc + phiz_cur)*(2*M_PI*fz*dnewxdf*dimx/vc + dphiz_cur);

            double xcirc = x + xoffset_uet;
            double dxcirc = dxoff_uet;

            double mult_uet = smooth_abs(smooth_floor(xcirc / delta_uet));
            double dmult_uet = smooth_abs_deriv(smooth_floor(xcirc / delta_uet))*smooth_floor_deriv(xcirc / delta_uet)*dxcirc;

            double x_uet = (xcirc - mult_uet*delta_uet)*Ax_uet/delta_uet - Ax_uet/2;
            double dx_uet = Ax_uet*dxcirc/delta_uet - Ax_uet*dmult_uet;

            double uet_effect = Az_uet*(1.0 - std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)));
            double duet = -0.5*(Az_uet/std::sqrt(1.0 - std::pow(x_uet/Ax_uet, 2)))*(-2)*(x_uet/Ax_uet)*(dx_uet/Ax_uet);

            double newz = 0;
            double dznewdf = 0;
            if(y <= y1 + mult*f + f/2){
                newz = -std::tan(alpha1)*(y - mult*f) + line_root1;
                dznewdf = std::tan(alpha1)*(mult + dmult*f) + dlrdf1;
            } else if(y <= y2 + mult*f + f/2){
                double yy = y - mult*f - f/2;
                newz = -std::sqrt(r*r - yy*yy) + r - ap;
                dznewdf = 2*yy*(mult + dmult*f + 0.5)/std::sqrt(r*r - yy*yy);
            } else {
                newz = std::tan(alpha2)*(y - mult*f) + line_root2;
                dznewdf = -std::tan(alpha2)*(mult + dmult*f) + dlrdf2;
            }
            newz += oscillation + uet_effect;
            dznewdf += doscilldf + duet;
            dz = smooth_min_deriv({oz, newz}, newz)*dznewdf;
    }
}

void draw_texture(sf::Uint8*& img, double* map_z, double ap, size_t w, size_t h, Colorscheme colorscheme){
    if(colorscheme == Colorscheme::GRAYSCALE){
        for(size_t x = 0; x < tex_width; ++x){
            for(size_t y = 0; y < tex_height; ++y){
                double& z = map_z[tex_width*y + x];

                sf::Uint8 rgb = (sf::Uint8)std::round(255*(1 + (z-max_z)/(max_z - min_z)));
                img[(tex_width*y + x)*4+0] = rgb;
                img[(tex_width*y + x)*4+1] = rgb;
                img[(tex_width*y + x)*4+2] = rgb;
                img[(tex_width*y + x)*4+3] = 255;
            }
        }
    } else if(colorscheme == Colorscheme::HSV){
        for(size_t x = 0; x < tex_width; ++x){
            for(size_t y = 0; y < tex_height; ++y){
                double& z = map_z[tex_width*y + x];

                double norm = -(z-max_z)/(max_z - min_z);
                double r, g, b;
                if(norm <= 0.25){
                    r = 1;
                    g = norm*4;
                    b = 0;
                } else if(norm <= 0.5){
                    r = 1 - (norm-0.25)*4;
                    g = 1;
                    b = 0;
                } else if(norm <= 0.75){
                    r = 0;
                    g = 1;
                    b = (norm-0.5)*4;
                } else {
                    r = 0;
                    g = 1 - (norm-0.75)*4;
                    b = 1;
                }

                img[(tex_width*y + x)*4+0] = (sf::Uint8)std::round(255*r);
                img[(tex_width*y + x)*4+1] = (sf::Uint8)std::round(255*g);
                img[(tex_width*y + x)*4+2] = (sf::Uint8)std::round(255*b);
                img[(tex_width*y + x)*4+3] = 255;
            }
        }

    }
}

double Sa(double*& map_z){
    size_t area = tex_width*tex_height;
    double sum = -std::accumulate(map_z, map_z+area, 0.0, std::plus<double>());
    sum -= area*max_z; // Depth correction

    return sum/area;
}

double dSa(double*& map_z, double* dzd, double dmax){
    size_t area = tex_width*tex_height;
    double sum = -std::accumulate(dzd, dzd+area, 0.0, std::plus<double>());
    sum -= area*dmax; // Depth correction

    return sum/area;
}

double surface_area(double*& map_z){
    // For a N*M matrix of points, the actual projected surface area is (dimx*dimy)*((N-1)*(M-1))
    double A = 0;
    for(size_t x = 0; x < tex_width; x+=2){
        for(size_t y = 0; y < tex_height; y+=2){
            Point p[9];
            size_t blockw = std::min(3ul, tex_width-x);
            size_t blockh = std::min(3ul, tex_height-y);
            if(blockw <= 1 || blockh <= 1){
                continue;
            }
            for(size_t i = 0; i < blockw; ++i){
                for(size_t j = 0; j < blockh; ++j){
                    p[j*3+i] = Point{
                        dimx*(x+i),
                        dimy*(y+j),
                        dimz*map_z[tex_width*(y+j) + (x+i)]
                    };
                }
            }

            A += triangle_area({p[0], p[1], p[4]});
            A += triangle_area({p[0], p[3], p[4]});

            if(blockw == 3){
                A += triangle_area({p[1], p[2], p[4]});
                A += triangle_area({p[2], p[4], p[5]});
            }
            if(blockh == 3){
                A += triangle_area({p[3], p[4], p[6]});
                A += triangle_area({p[4], p[6], p[7]});
            }
            if(blockw == 3 && blockh == 3){
                A += triangle_area({p[4], p[5], p[8]});
                A += triangle_area({p[4], p[7], p[8]});
            }
        }
    }

    return A;
}

double surface_area_dz(double*& map_z, double* dzd){
    double A = 0;
    for(size_t x = 0; x < tex_width; x+=2){
        for(size_t y = 0; y < tex_height; y+=2){
            Point p[9];
            double dz[9];
            size_t blockw = std::min(3ul, tex_width-x);
            size_t blockh = std::min(3ul, tex_height-y);
            if(blockw <= 1 || blockh <= 1){
                continue;
            }
            for(size_t i = 0; i < blockw; ++i){
                for(size_t j = 0; j < blockh; ++j){
                    p[j*3+i] = Point{
                        dimx*(x+i),
                        dimy*(y+j),
                        dimz*map_z[tex_width*(y+j) + (x+i)]
                    };
                    dz[j*3+i] = dimz*dzd[tex_width*(y+j) + (x+i)];
                }
            }

            A += triangle_area_deriv({p[0], p[1], p[4]}, {dz[0], dz[1], dz[4]});
            A += triangle_area_deriv({p[0], p[3], p[4]}, {dz[0], dz[3], dz[4]});

            if(blockw == 3){
                A += triangle_area_deriv({p[1], p[2], p[4]}, {dz[1], dz[2], dz[4]});
                A += triangle_area_deriv({p[2], p[4], p[5]}, {dz[2], dz[4], dz[5]});
            }
            if(blockh == 3){
                A += triangle_area_deriv({p[3], p[4], p[6]}, {dz[3], dz[4], dz[6]});
                A += triangle_area_deriv({p[4], p[6], p[7]}, {dz[4], dz[6], dz[7]});
            }
            if(blockw == 3 && blockh == 3){
                A += triangle_area_deriv({p[4], p[5], p[8]}, {dz[4], dz[5], dz[8]});
                A += triangle_area_deriv({p[4], p[7], p[8]}, {dz[4], dz[7], dz[8]});
            }
        }
    }

    return A;

}

int main(int argc, char* argv[]){

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

    // std::default_random_engine generator;
    // std::uniform_int_distribution<sf::Uint8> rand(0);
    sf::Uint8* px = new sf::Uint8[tex_width*tex_height*4]();
    double* map_z = new double[tex_width*tex_height]();
    double* orig_z = new double[tex_width*tex_height]();
    double* df = new double[tex_width*tex_height]();
    double* dap = new double[tex_width*tex_height]();
    double* dvc = new double[tex_width*tex_height]();

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
    mma.SetAsymptotes(0.001, 0.01, 1.01);

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
            draw_texture(px, map_z, ap, tex_width, tex_height, Colorscheme::HSV);
            double surarea = -surface_area(map_z);
            double roughness = Sa(map_z) - max_roughness;

            dzdf(orig_z, f, ap, vc, df);
            dzdap(orig_z, f, ap, vc, dap);
            dzdvc(orig_z, f, ap, vc, dvc);

            double dsurareadf = -surface_area_dz(map_z, df);
            double dSadf = dSa(map_z, df, dmax_zdf);
            double dsurareadap = -surface_area_dz(map_z, dap);
            double dSadap = dSa(map_z, dap, dmax_zdap);
            double dsurareadvc = -surface_area_dz(map_z, dvc);
            double dSadvc = dSa(map_z, dvc, dmax_zdvc);

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

        img.update(px);
        auto wsize = window.getSize();
        window_width = wsize.x;
        window_height = wsize.y;
        sprite.setPosition(sf::Vector2f(window_width/2-tex_width/2, window_height/2-tex_height/2));

        window.draw(sprite);
        window.display();
    }

    delete[] px;
    delete[] map_z;
    delete[] orig_z;
    delete[] df;
    delete[] dap;
    delete[] dvc;

    img.copyToImage().saveToFile("result.png");

    return 0;
}
