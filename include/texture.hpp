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

#ifndef TEXTURE_HPP
#define TEXTURE_HPP

#include <vector>

namespace texture{

void map(std::vector<double>& map_z, const std::vector<double>& orig_z, double f, double ap, double vc);

void dzdvc(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdvc);

void dzdap(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdap);

void dzdf(const std::vector<double>& orig_z, double f, double ap, double vc, std::vector<double>& dzdf);

}

#endif
