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

#ifndef SLP_HPP
#define SLP_HPP

#include <cmath>
#include <cstddef>
#include <vector>

/**
 * Successive Linear Programming implementation.
 */
class SLP{
    public:
    /**
     * Generate instance for N variables and M constraints.
     * Currently only works for M == 1.
     */
    SLP(const size_t N, const size_t M, const std::vector<double>& xmin, const std::vector<double>& xmax);

    /**
     * Update `x` based on derivatives and constraints.
     */
    void update(std::vector<double>& x, const std::vector<double>& dfdx, const std::vector<double>& g, const std::vector<double>& dgdx) const;

    private:
    const size_t N, M;
    std::vector<double> xmin;
    std::vector<double> xmax;
    const double EPS = 1e-7;

    /**
     * Update `x` based on derivatives and constraints using a "direct" method.
     */
    void direct(std::vector<double>& x, const std::vector<double>& dfdx, const std::vector<double>& g, const std::vector<double>& dgdx) const;

    /**
     * Dot product between two vectors.
     */
    inline double dot(const std::vector<double>& v1, const std::vector<double>& v2) const{
        double c = 0;
        for(size_t i = 0; i < v1.size(); ++i){
            c += v1[i]*v2[i];
        }
        return c;
    }
    /**
     * Multiplies a matrix and a vector.
     */
    inline std::vector<double> mat_dot_vec(const std::vector<double>& m, const std::vector<double>& v) const{
        const size_t l = m.size()/v.size();
        const size_t n = v.size();
        std::vector<double> r(l,0);
        for(size_t i = 0; i < l; ++i){
            for(size_t j = 0; j < n; ++j){
                r[i] += m[i*n + j]*v[j];
            }
        }
        return r;
    }
    /**
     * Scales a vector or matrix.
     */
    inline std::vector<double> scale(const double s, const std::vector<double>& v) const{
        auto sv = v;
        for(size_t i = 0; i < v.size(); ++i){
            sv[i] *= s;
        }
        return sv;
    }
    /**
     * Normalizes a vector.
     */
    inline std::vector<double> normalized(const std::vector<double>& v) const{
        double norm = 0;
        for(auto& i:v){
            norm += i*i;
        }
        norm = std::sqrt(norm);
        auto svn = v;
        for(auto& i:svn){
            i /= norm;
        }
        return svn;
    }
};

#endif
