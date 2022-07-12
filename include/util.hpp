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

#ifndef UTIL_HPP
#define UTIL_HPP

#include <array>
#include <cmath>

namespace util{

struct Point{
    double x, y, z;
};

struct Vector{
    double x, y, z;
};

inline double triangle_area(std::array<Point, 3> p){
    Vector v1{p[1].x - p[0].x, p[1].y - p[0].y, p[1].z - p[0].z};
    Vector v2{p[2].x - p[0].x, p[2].y - p[0].y, p[2].z - p[0].z};

    double A = 0.5*std::sqrt(
        std::pow(v1.y*v2.z - v2.y*v1.z, 2) +
        std::pow(v1.x*v2.z - v2.x*v1.z, 2) +
        std::pow(v1.y*v2.x - v2.y*v1.x, 2)
    );

    return A;
}


inline double triangle_area_deriv(std::array<Point, 3> p, std::array<double, 3> dz){
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
}

#endif
