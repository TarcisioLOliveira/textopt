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

#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include <tuple>
#include <vector>

namespace analysis{

/**
 * Generates points from surface area and roughness within the specified range
 * for f and ap, with constant vc. Used for figuring out the resulting shape
 * and specifying the correct optimization method.
 *
 * @param f Range for feed rate
 * @param ap Range for cutting depth
 * @param vc Constant cutting velocity
 * @param step Step for getting points within ranges
 * @param map_z Memory space used for texture calculations
 * @param orig_z Original texture
 */
void plot_fxap(const std::tuple<double, double>& f, const std::tuple<double, double>& ap, const double vc, const double step, std::vector<double>& map_z, const std::vector<double>& orig_z);

/**
 * Generates points from surface area and roughness within the specified range
 * for vc, with constant f and ap. Used for figuring out the resulting shape
 * and specifying the correct optimization method. However, vc doesn't seem
 * to affect roughness and surface area too much, except with vc -> 0.
 *
 * @param f Constant feed rate
 * @param ap Constant cutting depth
 * @param vc Range for cutting velocity
 * @param step Step for getting points within ranges
 * @param map_z Memory space used for texture calculations
 * @param orig_z Original texture
 */
void plot_vc(const double f, const double ap, const std::tuple<double, double>& vc, const double step, std::vector<double>& map_z, const std::vector<double>& orig_z);

}

#endif
