/*
 *   Copyright (C) 2024 Tarcísio Ladeia de Oliveira.
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

#ifndef NEWTON_OPT_HPP
#define NEWTON_OPT_HPP

#include <cmath>
#include <cstddef>
#include <vector>
#include <Eigen/Core>

class NewtonOpt{
    public:
    NewtonOpt(const size_t N, const size_t M);

    void update(std::vector<double>& x, const double f, const std::vector<double>& dfdx, const std::vector<double>& g, const std::vector<double>& dgdx) const;

    private:
    const size_t N, M;
};

#endif
