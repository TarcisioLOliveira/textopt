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

#include "slp.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

SLP::SLP(const size_t N, const size_t M, const std::vector<double>& xmin, const std::vector<double>& xmax):
    N(N), M(M), xmin(xmin), xmax(xmax){}

void SLP::update(std::vector<double>& x, const std::vector<double>& dfdx, const std::vector<double>& g, const std::vector<double>& dgdx) const{
    this->direct(x, dfdx, g, dgdx);
}

void SLP::direct(std::vector<double>& x, const std::vector<double>& dfdx, const std::vector<double>& g, const std::vector<double>& dgdx) const{
    // x_{i+1} = x_i + l_i*dfdx_i
    // dgdx_i*(x_{i+1}-x_i) <= g_i
    //
    // dgdx_i*(l_i*dfdx_i) <= g_i
    // l_i <= g_i/(dgdx_i*dfdx_i)
    
    // Normalize vectors (currently assuming single constraint)
    const auto norm_dfdx = this->normalized(dfdx);
    const auto norm_dgdx = this->normalized(dgdx);


    // Calculate l_i
    const double dg_df = this->dot(norm_dgdx, norm_dfdx);
    double l_i = -g[0]/dg_df;

    // Make sure l_i does not make x too low or too high.
    for(size_t i = 0; i < N; ++i){
        if(x[i] + l_i*norm_dfdx[i] < this->xmin[i]){
            l_i = 0.5*(xmin[i]-x[i])/norm_dfdx[i];
        } else if(x[i] + l_i*norm_dfdx[i] > this->xmax[i]){
            l_i = 0.5*(xmax[i]-x[i])/norm_dfdx[i];
        }
    }

    // Update x
    for(size_t i = 0; i < N; ++i){
        x[i] += l_i*norm_dfdx[i];
    }
}
