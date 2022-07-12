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

#ifndef OPT_FUNC_HPP
#define OPT_FUNC_HPP

#include <vector>

namespace opt{

double Sa(const std::vector<double>& map_z);

double dSa(const std::vector<double>& dzd, double dmax);

double surface_area(const std::vector<double>& map_z);

double surface_area_dz(const std::vector<double>& map_z, const std::vector<double>& dzd);

}

#endif
