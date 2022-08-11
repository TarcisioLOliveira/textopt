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
 * @file opt_func.hpp
 *
 * @brief Optimization functions.
 *
 * Functions used for the optimization routine (objective function and
 * restrictions) as well as their derivatives.
 */

#ifndef OPT_FUNC_HPP
#define OPT_FUNC_HPP

#include <vector>

namespace opt{

/**
 * Mean surface roughness.
 *
 * @param map_z Texture.
 *
 * @return Mean surface roughness.
 */
double Sa(const std::vector<double>& map_z);

/**
 * Derivative of mean surface roughness.
 *
 * @param dzd Derivative of texture in relation to design variable.
 * @param map_z Texture.
 * @param dmax Derivative of `max_z` in relation to design variable.
 * @param dmax Derivative of `min_z` in relation to design variable.
 *
 * @return Derivative of mean surface roughness.
 */
double dSa(const std::vector<double>& dzd, const std::vector<double>& map_z, double dmax, double dmin);

/**
 * Texture's surface area.
 * Calculates the area by diving the 3D points into multiple triangles.
 * Inspired by (JENNESS, 2004).
 *
 * @param map_z Texture.
 *
 * @return Texture's surface area.
 */
double surface_area(const std::vector<double>& map_z);

/**
 * Derivative of texture's surface area.
 *
 * @param map_z Texture.
 * @param dzd Derivative of texture in relation to design variable.
 *
 * @return Derivative of texture's surface area.
 */
double surface_area_dz(const std::vector<double>& map_z, const std::vector<double>& dzd);

}

#endif
