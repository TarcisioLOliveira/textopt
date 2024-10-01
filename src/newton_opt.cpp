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

#include "newton_opt.hpp"
#include <Eigen/Dense>
#include <iostream>

NewtonOpt::NewtonOpt(const size_t N, const size_t M):
    N(N), M(M){
}

void NewtonOpt::update(std::vector<double>& x, const double f, const std::vector<double>& dfdx, const std::vector<double>& g, const std::vector<double>& dgdx) const{
    Eigen::VectorXd vals;
    vals.resize(1 + M);
    Eigen::MatrixXd J;
    J.resize(N, 1 + M);
    vals[0] = f;
    for(size_t i = 0; i < M; ++i){
        vals[1+i] = g[i];
    }
    for(size_t i = 0; i < N; ++i){
        J(i, 0) = dfdx[i];
        for(size_t j = 0; j < M; ++j){
            J(i, 1+j) = dgdx[j*N + i];
        }
    }

    Eigen::VectorXd b = -J*vals;
    Eigen::MatrixXd A = J*J.transpose();

    std::cout << std::endl;
    std::cout << std::endl;
    for(size_t i = 0; i < N; ++i){
        std::cout << b[i] << " ";
    }
    std::cout << std::endl << std::endl;
    for(size_t i = 0; i < N; ++i){
        for(size_t j = 0; j < N; ++j){
            std::cout << A(i,j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    Eigen::VectorXd dx = A.fullPivLu().solve(b);
    for(size_t i = 0; i < N; ++i){
        std::cout << dx[i] << " ";
        x[i] += 1e-6*dx[i];
    }
    std::cout << std::endl << std::endl;
    for(size_t i = 0; i < N; ++i){
        std::cout << x[i] << " ";
    }
    std::cout << std::endl << std::endl;
}
