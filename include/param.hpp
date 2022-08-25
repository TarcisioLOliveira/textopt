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
 * @file param.hpp
 *
 * @brief Simulation parameters.
 *
 * Collection of global variables used as parameters for many functions.
 */

#ifndef PARAM_HPP
#define PARAM_HPP

#include <cmath>

namespace param{

// Texture dimensions
inline const size_t tex_width  = 400;
inline const size_t tex_height = 400;
inline const double base_area = (tex_width-1)*(tex_height-1);

inline const double MULT = -1;       // Exponent for smooth::min()
inline const double FLOORD = 1e-10;  // smooth::floor() precision
inline const double ABS_EPS = 1e-15; // smooth::abs() precision
inline const double YOFF = 0.001;    // Workaround over smooth::floor() being off by -0.5 for integers

// Unit scaling. Currently defaulting to um as base unit for texture.
// There is naive support for different scales in different directions,
// but it may break tool radius and ellipse shape if they're different.
// Defaults to same scale to all sides.
inline const double dim_scale = 1e3; // mm to um
inline const double dim = 1;
inline const double dimx = dim;
inline const double dimy = dim;
inline const double dimz = dim;

// Tool parameters
inline const double alpha1 = 60*M_PI/180; // Back angle (descending)
inline const double alpha2 = 60*M_PI/180; // Front angle (ascending)
inline const double r = 20/dim; // [um] Tool radius
                                
// Maximum and minimum 
inline double max_z = 0; // Used to calculate Sa, so uses smooth::min().
                         // It's actually just an estimate, and ends up (maybe)
                         // falling slightly short because of the oscillations
                         // calculated in. Using smooth::min() for the entire
                         // texture would be impossible, though.
                         
inline double min_z = 0; // Used for Sa and texture display
// Derivatives of max_z
inline double dmax_zdf = 0; 
inline double dmax_zdap = 0; 
inline double dmax_zdvc = 0; 

// "Random" oscillation
// Models height variation along the tool path

// Amplitude [dim], e.g. [um]
inline const double Ax = 0.5;
inline const double Az = 0.1;
// Frequency [Hz]
inline const double fx = 19*1e6;
inline const double fz = 51*1e6;
// Phase [rad]
inline const double phix = 20*M_PI/180;
inline const double phiz = 0;
// Cylinder radius, to simulate the phase change with the turn of the cylinder.
inline const double cylinder_radius = 6*dim_scale; // [um]

// Ultrasonic elliptical turning
inline const double f_uet = 20000; // [Hz]
inline const double Ax_uet = 1; // [um]
inline const double Az_uet = 2; // [um]

}

#endif
