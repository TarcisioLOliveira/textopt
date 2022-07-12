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

#ifndef SMOOTH_HPP
#define SMOOTH_HPP

#include <initializer_list>
#include <numeric>

#include "param.hpp"

namespace smooth{

inline double min(std::initializer_list<double> x){
    double frac_top = 0;
    double frac_bot = 0;
    for(double xi : x){
        double exp = std::exp(param::MULT*xi);
        frac_top += xi*exp;
        frac_bot += exp;
    }

    return frac_top/frac_bot;
};

inline double min_deriv(std::initializer_list<double> x, size_t i){
    double frac_top = 0;
    double frac_bot = 0;
    for(double xi : x){
        double exp = std::exp(param::MULT*xi);
        frac_top += xi*exp;
        frac_bot += exp;
    }

    double sm = frac_top/frac_bot;

    return (std::exp(param::MULT*(*(x.begin()+i)))/frac_bot)*(1+param::MULT*((*(x.begin()+i))-sm));
};

inline double floor(double v){
    // https://math.stackexchange.com/questions/2746958/smooth-floor-function
    //
    // δ = 0.01;
    // trg[x_] := 1 - 2 ArcCos[(1 - δ) Sin[2 π x]]/π;
    // sqr[x_] := 2 ArcTan[Sin[2 π x]/δ]/π;
    // swt[x_] := (1 + trg[(2 x - 1)/4] sqr[x/2])/2;
    // flr[x_] := x-swt[x]

    double p1x = (1.0 - param::FLOORD)*std::sin(2*M_PI*((2*v-1)/4));
    double p2x = std::sin(2*M_PI*(v/2))/param::FLOORD;

    double p1 = 1 - 2*std::acos(p1x)/M_PI;
    double p2 = 2 * std::atan(p2x)/M_PI;
    double sawtooth = (1 + p1 * p2)/2;

    double result = v - sawtooth;

    return result;
}

inline double floor_deriv(double v){
    double p1x = (1.0 - param::FLOORD)*std::sin(2*M_PI*((2*v-1)/4));
    double dp1x = (1.0 - param::FLOORD)*std::cos(2*M_PI*((2*v-1)/4))*2*M_PI*2/4;

    double p2x = std::sin(2*M_PI*(v/2))/param::FLOORD;
    double dp2x = (std::cos(2*M_PI*(v/2))/param::FLOORD)*2*M_PI/2;

    double p1 = 1 - 2*std::acos(p1x)/M_PI;
    double dp1 = (2.0/(std::sqrt(1 - p1x*p1x)*M_PI))*dp1x;
    double p2 = 2 * std::atan(p2x)/M_PI;
    double dp2 = (2.0/((p2x*p2x + 1)*M_PI))*dp2x;
    double dsawtooth = (1 + dp1 * p2 + p1 * dp2)/2;

    double result = 1.0 - dsawtooth;
    
    return result;
}

/**
 * Used mostly to fix smooth::floor() becoming negative when close to zero.
 *
 * Not really smooth, yes, but it's actually more of a workaround. It works
 * better this way, as smooth::abs_deriv() was introducing instability into
 * the optimization process.
 */
inline double abs(double v){
    return std::abs(v);
    // return std::sqrt(v*v+ABS_EPS);
}

inline double abs_deriv(double v){
    // return v/std::sqrt(v*v+ABS_EPS);
    return 1;
}

}

#endif
