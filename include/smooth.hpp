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
/**
 * @file smooth.hpp
 * @brief Smooth versions of functions.
 *
 * Smooth versions of functions and their derivatives, for use in the 
 * optimization and texture functions.
 */

#ifndef SMOOTH_HPP
#define SMOOTH_HPP

#include <initializer_list>
#include <numeric>

#include "param.hpp"

namespace smooth{

/**
 * Smooth version of `min()`, with an initializer list for argument. Uses
 * the exponential function for the calculation, as p-norm doesn't handle
 * a mix of positive and negative numbers very well. As a result, avoid
 * using it with multiple numbers.
 *
 * The exponentiation is regulated by the parameter param::MULT. Higher values
 * lead to a better approximation of `min()`, but it can very easily end up
 * beyond DOUBLE_MAX.
 *
 * @param x List of values to get the minimum value.
 *
 * @return Approximation of the minimum value of x.
 */
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

/**
 * Derivative of the smooth `min` function in relation to parameter `i`.
 *
 * @param x List of values to get the minimum value.
 * @param i Parameter for derivative.
 *
 * @return Derivative of smooth::min() in relation to parameter `i`.
 */
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

/**
 * Approximation of the `floor` function based on the sawtooth function.
 * Quality of approximation is dependent on the value of param::FLOORD. Lower
 * values increase quality of approximation. However, it seems to also increase
 * the values of the derivatives of the optimization function in relation `f`.
 * The functions seem to work correctly nonetheless, but it's important to
 * keep track of that behavior.
 *
 * It should also be noted that the function does not behave as expected when
 * the input is an integer. The result tends to be `v-0.5`. This is avoided in
 * the optimization functions by adding param::YOFF to the input.
 *
 * However, when `v == 0`, smooth::floor(v) < 0. In the functions, this is 
 * worked around by using smooth::abs(smooth::floor(v)). As the result is
 * still close to zero, it results in a good approximation.
 *
 * @param v Number to be floored.
 *
 * @return Approximation of `floor(v)`.
 */
inline double floor(double v){
    // https://math.stackexchange.com/questions/2746958/smooth-floor-function
    //
    // δ = 0.01;
    // trg[x_] := 1 - 2 ArcCos[(1 - δ) Sin[2 π x]]/π;
    // sqr[x_] := 2 ArcTan[Sin[2 π x]/δ]/π;
    // swt[x_] := (1 + trg[(2 x - 1)/4] sqr[x/2])/2;
    // flr[x_] := x-swt[x]

    const double p1x = (1.0 - param::FLOORD)*std::sin(2*M_PI*((2*v-1)/4));
    const double p2x = std::sin(2*M_PI*(v/2))/param::FLOORD;

    const double p1 = 1 - 2*std::acos(p1x)/M_PI;
    const double p2 = 2 * std::atan(p2x)/M_PI;
    const double sawtooth = (1 + p1 * p2)/2;

    const double result = v - sawtooth;

    return result;
}

/**
 * Calculates the derivative of smooth::floor() in relation to `v`.
 *
 * @param v Number to be floored.
 *
 * @return Derivative of smooth::floor(v).
 */
inline double floor_deriv(double v){
    const double p1x = (1.0 - param::FLOORD)*std::sin(2*M_PI*((2*v-1)/4));
    const double dp1x = (1.0 - param::FLOORD)*std::cos(2*M_PI*((2*v-1)/4))*2*M_PI*2/4;

    const double p2x = std::sin(2*M_PI*(v/2))/param::FLOORD;
    const double dp2x = (std::cos(2*M_PI*(v/2))/param::FLOORD)*2*M_PI/2;

    const double p1 = 1 - 2*std::acos(p1x)/M_PI;
    const double dp1 = (2.0/(std::sqrt(1 - p1x*p1x)*M_PI))*dp1x;
    const double p2 = 2 * std::atan(p2x)/M_PI;
    const double dp2 = (2.0/((p2x*p2x + 1)*M_PI))*dp2x;
    const double dsawtooth = (1 + dp1 * p2 + p1 * dp2)/2;

    const double result = 1.0 - dsawtooth;
    
    return result;
}

/**
 * Used mostly to fix smooth::floor() becoming negative when close to zero.
 *
 * Not really smooth, yes, but it's actually more of a workaround. It works
 * better this way, as smooth::abs_deriv() was introducing instability into
 * the optimization process.
 *
 * @param v Number to be taken the absolute from.
 *
 * @return Absolute of v.
 */
inline double abs(double v){
    return std::abs(v);
    // return std::sqrt(v*v+ABS_EPS);
}

/**
 * Also a workaround. Just returns 1 currently.
 *
 * @param v Number to be taken the absolute from.
 *
 * @return Derivative of smooth::abs(v).
 */
inline double abs_deriv(double v){
    // return v/std::sqrt(v*v+ABS_EPS);
    return 1;
}

}

#endif
