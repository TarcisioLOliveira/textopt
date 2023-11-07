/*
 *   Copyright (C) 2023 Tarcísio Ladeia de Oliveira.
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

#ifndef NEWTON_HPP
#define NEWTON_HPP

#include <cmath>
#include <iostream>
#include "param.hpp"

inline double zt(const double t){
    return param::Az_uet*std::cos(2*M_PI*param::f_uet*t);
}
inline double dzdt(const double t){
    return -param::Az_uet*2*M_PI*param::f_uet*std::sin(2*M_PI*param::f_uet*t);
}

inline double xt(const double t, const double vc){
    return -param::Ax_uet*std::sin(2*M_PI*param::f_uet*t) + vc*t;
}
inline double dxdt(const double t, const double vc){
    return -param::Ax_uet*2*M_PI*param::f_uet*std::cos(2*M_PI*param::f_uet*t) + vc;
}


inline double dtdvc(const double t, const double vc){
    return -t/dxdt(t, vc);
}

inline double newton_t(const double x0, const double t0, const double vc){
    double t = t0;
    double dt = 0;
    double x = 0, dx = 0, ti = 0;
    do{
        x = xt(t, vc) - x0;
        dx = dxdt(t, vc);
        ti = t - x/dx;

        dt = std::abs(ti - t);
        t = ti;
    } while(dt > 1e-7);

    return t;
}

#endif
