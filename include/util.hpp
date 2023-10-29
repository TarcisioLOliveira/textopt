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
 * @file util.hpp
 * @brief Utilities.
 *
 * Helper functions and classes.
 */

#ifndef UTIL_HPP
#define UTIL_HPP

#include <array>
#include <cmath>

namespace util{

/**
 * Structure representing a 3D point.
 */
struct Point{
    double x, y, z;
};

/**
 * Structure representing a 3D vector.
 */
struct Vector{
    double x, y, z;
};

/**
 * Calculates the area of a triangle.
 *
 * @param p Array of 3 Point instances.
 *
 * @return Calculated area.
 */
inline double triangle_area(const std::array<Point, 3>& p){
    const Vector v1{p[1].x - p[0].x, p[1].y - p[0].y, p[1].z - p[0].z};
    const Vector v2{p[2].x - p[0].x, p[2].y - p[0].y, p[2].z - p[0].z};

    const double expr1 =  v1.y*v2.z - v1.z*v2.y;
    const double expr2 = -v1.x*v2.z + v1.z*v2.x;
    const double expr3 =  v1.x*v2.y - v1.y*v2.x;

    const double A = 0.5*std::sqrt(
        expr1*expr1 +
        expr2*expr2 +
        expr3*expr3
    );

    return A;
}

/**
 * Calculates the derivative of the triangle area in relation to some variable.
 * For a variable 'a', does the calculation by dA/da = (dA/dz)(dz/da). The
 * variable 'dz' is used for 'dz/da'.
 *
 * @param p Array of 3 Point instances.
 * @param dz Array of derivatives of 'z',
 *
 * @return Derivative of the calculated area.
 */
inline double triangle_area_deriv(const std::array<Point, 3>& p, const std::array<double, 3>& dz){
    const Vector v1{p[1].x - p[0].x, p[1].y - p[0].y, p[1].z - p[0].z};
    const Vector v2{p[2].x - p[0].x, p[2].y - p[0].y, p[2].z - p[0].z};
    const double dv1 = dz[1] - dz[0];
    const double dv2 = dz[2] - dz[0];

    const double expr1 =  v1.y*v2.z - v1.z*v2.y;
    const double expr2 = -v1.x*v2.z + v1.z*v2.x;
    const double expr3 =  v1.x*v2.y - v1.y*v2.x;

    const double A = (0.25/std::sqrt(
        expr1*expr1 +
        expr2*expr2 +
        expr3*expr3
    ))*(2*(v1.y*v2.z - v2.y*v1.z)*( v1.y*dv2 - dv1*v2.y) +
        2*(v1.x*v2.z - v2.x*v1.z)*(-v1.x*dv2 + dv1*v2.x) +
        0
    );

    return A;
}
}

#endif
