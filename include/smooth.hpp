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
#include <vector>

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
 * Smooth version of `min()`, with a vector for argument. Uses
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
inline double min(std::vector<double> x){
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

    const double sm = frac_top/frac_bot;
    const double xj = *(x.begin()+i);

    return (std::exp(param::MULT*xj)/frac_bot)*(1+param::MULT*(xj-sm));
};

/**
 * Derivative of the smooth `min` function with multiple derivatives.
 *
 * @param x List of values to get the minimum value.
 * @param dx List of derivatives.
 *
 * @return Derivative of smooth::min() in relation to multiple parameters.
 */
inline double min_deriv(std::initializer_list<double> x, std::initializer_list<double> dx){
    double frac_top = 0;
    double frac_bot = 0;
    for(double xi : x){
        double exp = std::exp(param::MULT*xi);
        frac_top += xi*exp;
        frac_bot += exp;
    }

    const double sm = frac_top/frac_bot;
    double result = 0;
    for(size_t i = 0; i < x.size(); ++i){
        const double xj = *(x.begin()+i);
        const double dxj = *(dx.begin()+i);

        result += (std::exp(param::MULT*xj)/frac_bot)*(1+param::MULT*(xj-sm))*dxj;
    }

    return result;
};

/**
 * Derivative of the smooth `min` function with multiple derivatives.
 *
 * @param x List of values to get the minimum value.
 * @param dx List of derivatives.
 *
 * @return Derivative of smooth::min() in relation to multiple parameters.
 */
inline double min_deriv(std::vector<double> x, std::vector<double> dx){
    double frac_top = 0;
    double frac_bot = 0;
    for(double xi : x){
        double exp = std::exp(param::MULT*xi);
        frac_top += xi*exp;
        frac_bot += exp;
    }

    const double sm = frac_top/frac_bot;
    double result = 0;
    for(size_t i = 0; i < x.size(); ++i){
        const double xj = *(x.begin()+i);
        const double dxj = *(dx.begin()+i);

        result += (std::exp(param::MULT*xj)/frac_bot)*(1+param::MULT*(xj-sm))*dxj;
    }

    return result;
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

    const double p1 = 1.0 - 2.0*std::acos(p1x)/M_PI;
    const double p2 = 2.0 * std::atan(p2x)/M_PI;
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

    const double p1 = 1.0 - 2.0*std::acos(p1x)/M_PI;
    const double dp1 = (2.0/(std::sqrt(1.0 - p1x*p1x)*M_PI))*dp1x;
    const double p2 = 2.0 * std::atan(p2x)/M_PI;
    const double dp2 = (2.0/((p2x*p2x + 1.0)*M_PI))*dp2x;
    const double dsawtooth = (dp1 * p2 + p1 * dp2)/2;

    const double result = 1.0 - dsawtooth;
    
    return result;
}

/**
 * Smooth version of the absolute function, based on the square root of squared
 * number, along with a constant, in order to smooth it at `v=0`.
 *
 * @param v Number to be taken the absolute from.
 *
 * @return Absolute of v.
 */
inline double abs(double v){
    return std::sqrt(v*v+param::ABS_EPS);
}

/**
 * Derivative of the smooth version of the absolute function.
 *
 * @param v Number to be taken the absolute from.
 *
 * @return Derivative of smooth::abs(v).
 */
inline double abs_deriv(double v){
    return v/std::sqrt(v*v+param::ABS_EPS);
}

}

#endif
