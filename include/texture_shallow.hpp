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

#ifndef TEXTURE_SHALLOW_HPP
#define TEXTURE_SHALLOW_HPP

#include <vector>

namespace texture_shallow{

/**
 * Exact version of `map()`, that is, without using smooth approximations. it's
 * non-differentiable, meaning it shouldn't be used for optimization, but it's
 * faster and exact, making it useful for plotting and comparing results.
 *
 * Current tests haven't been using orig_z, so it's not currently supported.
 *
 * @param map_z Resulting texture.
 * @param orig_z Original texture.
 * @param f Feed rate.
 * @param ap Cutting depth.
 * @param vc Cutting velocity.
 */
void map_exact(std::vector<double>& map_z, double f, double vc);

/**
 * Generate texture based on cutting parameters, parameters in param.hpp, and
 * the original texture.
 *
 * Current tests haven't been using orig_z, so it's not currently supported.
 *
 * @param map_z Resulting texture.
 * @param orig_z Original texture.
 * @param f Feed rate.
 * @param vc Cutting velocity.
 */
void map(std::vector<double>& map_z, double f, double vc);

/**
 * Calculate derivative of depth (z) in relation to cutting velocity (vc).
 *
 * @param orig_z Original texture.
 * @param f Feed rate.
 * @param vc Cutting velocity.
 * @param dzdvc Resulting derivative at every point.
 */
void dzdvc(double f, double vc, std::vector<double>& dzdvc);

/**
 * Calculate derivative of depth (z) in relation to feed rate (f).
 *
 * @param orig_z Original texture.
 * @param f Feed rate.
 * @param vc Cutting velocity.
 * @param dzdf Resulting derivative at every point.
 */
void dzdf(double f, double vc, std::vector<double>& dzdf);

}

#endif
