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

#ifndef PARAM_HPP
#define PARAM_HPP

#include <cmath>

namespace param{

inline size_t tex_width = 300;
inline size_t tex_height = 300;

inline const double MULT = -1; // Exponent for smooth_min()
inline const double FLOORD = 1e-10; // smooth_floor() precision
inline const double ABS_EPS = 1e-15; // smooth_abs() precision
inline const double YOFF = 0.001; // Workaround over smooth_floor() being off by -0.5 for integers

inline double dim_scale = 1e6; // m to um
inline double dim = 1;
inline double dimx = dim;
inline double dimy = dim;
inline double dimz = dim;

inline double alpha1 = 60*M_PI/180;
inline double alpha2 = 60*M_PI/180;
inline double r = 20/dim; // um
inline double max_z = 0; // Used to calculate Sa, so uses smooth_min
inline double min_z = 0; // Used only for texture display, so uses std::min()
inline double dmax_zdf = 0; 
inline double dmax_zdap = 0; 
inline double dmax_zdvc = 0; 

// "Random" oscillation
// Models height variation along the tool path

// Amplitude [dim], e.g. [um]
inline double Ax = 0.5;
inline double Az = 0.1;
// Frequency [Hz]
inline double fx = 19*1e5;
inline double fz = 51*1e5;
// Phase [rad]
inline double phix = 20*M_PI/180;
inline double phiz = 0;

inline double cylinder_radius = 6*1e-3*dim_scale; // [um]

// Ultrasonic elliptical turning
inline double f_uet = 20000; // [Hz]
inline double Ax_uet = 1; // [um]
inline double Az_uet = 2; // [um]

}

#endif
