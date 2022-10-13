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

inline enum class AnalysisType{
    OPT,
    PLOT,
    SINGLE
} analysis_type;

inline enum class OptMethod{
    SLP,
    MMA
} opt_method;

inline enum class PlotMethod{
    FXAP,
    VC,
} plot_method;

inline enum class SingleMethod{
    SMOOTH,
    EXACT
} single_method;

inline bool opt_ap = false;

// Window dimensions
inline size_t window_width = 600;
inline size_t window_height = 600;

// Texture dimensions
inline size_t tex_width  = 400;
inline size_t tex_height = 400;
inline double base_area = (tex_width-1)*(tex_height-1);

inline double stop = 1e-8;

inline double MULT = -1;       // Exponent for smooth::min()
inline double FLOORD = 1e-10;  // smooth::floor() precision
inline double ABS_EPS = 1e-15; // smooth::abs() precision
inline double YOFF = 0.001;    // Workaround over smooth::floor() being off by -0.5 for integers

// Unit scaling. Currently defaulting to um as base unit for texture.
// There is naive support for different scales in different directions,
// but it may break tool radius and ellipse shape if they're different.
// Defaults to same scale to all sides.
inline double dim_scale = 1e3; // mm to um
inline double dim = 1;
inline double dimx = dim;
inline double dimy = dim;
inline double dimz = dim;

// Step for plotting
inline double step = 0.1;

// Max roughness (constraint) for optimization
inline double max_roughness = 2;

// Initial values
inline double f  = 200;
inline double ap = 200;
inline double vc = 200;

// Lower bounds
inline double f_min  = 1;
inline double ap_min = 1;
inline double vc_min = 1;

// Upper bounds
inline double f_max  = 200;
inline double ap_max = 200;
inline double vc_max = 200;

// Tool parameters
inline double alpha1 = 60*M_PI/180; // Back angle (descending)
inline double alpha2 = 60*M_PI/180; // Front angle (ascending)
inline double r = 20/dim; // [um] Tool radius
                                
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

// Derivatives of min_z (texture_shallow)
inline double dmin_zdf = 0; 
inline double dmin_zdvc = 0; 


// "Random" oscillation
// Models height variation along the tool path

// Amplitude [dim], e.g. [um]
inline double Az = 0.1;
// Frequency [Hz]
inline double fz = 51*1e6;
// Phase [rad]
inline double phiz = 0;
// Cylinder radius, to simulate the phase change with the turn of the cylinder.
inline double cylinder_radius = 6*dim_scale; // [um]

// Ultrasonic elliptical turning
inline double f_uet = 20000; // [Hz]
inline double Ax_uet = 1; // [um]
inline double Az_uet = 2; // [um]
                          
inline double tan1 = std::tan(alpha1);
inline double tan2 = std::tan(alpha2);

// Line parameters for tool slopes.
inline double y1 = -tan1*r/std::sqrt(tan1*tan1+1);
inline double y2 =  tan2*r/std::sqrt(tan2*tan2+1);

// Value of z for y = 0 (with circle center at y = 0).
inline double b1off = -std::sqrt(r*r - y1*y1) + r + tan1*y1;
inline double b2off = -std::sqrt(r*r - y2*y2) + r - tan2*y2;

}

#endif
