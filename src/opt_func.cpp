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

#include <numeric>

#include "opt_func.hpp"
#include "param.hpp"
#include "util.hpp"
#include "smooth.hpp"

namespace opt{

double Sa(const std::vector<double>& map_z){
    using namespace param;

    const double avg = (max_z + min_z)/2;
    const size_t area = tex_width*tex_height;
    double Sa = 0;
    for(double z:map_z){
        Sa += smooth::abs(z - avg);
    }
    Sa /= area;

    return Sa;
}

double dSa(const std::vector<double>& dzd, const std::vector<double>& map_z, double dmax, double dmin){
    using namespace param;

    const double avg = (max_z + min_z)/2;
    const double davg = (dmax + dmin)/2;
    const size_t area = tex_width*tex_height;
    double dSa = 0;
    for(size_t i = 0; i < area; ++i){
        dSa += smooth::abs_deriv(map_z[i] - avg)*(dzd[i] - davg);
    }
    dSa /= area;

    return dSa;
}

double surface_area(const std::vector<double>& map_z){
    using namespace param;

    // For a N*M matrix of points, the actual projected surface area is (dimx*dimy)*((N-1)*(M-1))
    double A = 0;
    for(size_t x = 0; x < tex_width; x+=2){
        for(size_t y = 0; y < tex_height; y+=2){
            util::Point p[9];
            const size_t blockw = std::min(3ul, tex_width-x);
            const size_t blockh = std::min(3ul, tex_height-y);
            if(blockw <= 1 || blockh <= 1){
                continue;
            }
            for(size_t i = 0; i < blockw; ++i){
                for(size_t j = 0; j < blockh; ++j){
                    p[j*3+i] = util::Point{
                        dimx*(x+i),
                        dimy*(y+j),
                        dimz*map_z[tex_width*(y+j) + (x+i)]
                    };
                }
            }

            A += util::triangle_area({p[0], p[1], p[4]});
            A += util::triangle_area({p[0], p[3], p[4]});

            if(blockw == 3){
                A += util::triangle_area({p[1], p[2], p[4]});
                A += util::triangle_area({p[2], p[4], p[5]});
            }
            if(blockh == 3){
                A += util::triangle_area({p[3], p[4], p[6]});
                A += util::triangle_area({p[4], p[6], p[7]});
            }
            if(blockw == 3 && blockh == 3){
                A += util::triangle_area({p[4], p[5], p[8]});
                A += util::triangle_area({p[4], p[7], p[8]});
            }
        }
    }

    return A;
}

double surface_area_dz(const std::vector<double>& map_z, const std::vector<double>& dzd){
    using namespace param;

    double A = 0;
    for(size_t x = 0; x < tex_width; x+=2){
        for(size_t y = 0; y < tex_height; y+=2){
            util::Point p[9];
            double dz[9];
            const size_t blockw = std::min(3ul, tex_width-x);
            const size_t blockh = std::min(3ul, tex_height-y);
            if(blockw <= 1 || blockh <= 1){
                continue;
            }
            for(size_t i = 0; i < blockw; ++i){
                for(size_t j = 0; j < blockh; ++j){
                    p[j*3+i] = util::Point{
                        dimx*(x+i),
                        dimy*(y+j),
                        dimz*map_z[tex_width*(y+j) + (x+i)]
                    };
                    dz[j*3+i] = dimz*dzd[tex_width*(y+j) + (x+i)];
                }
            }

            A += util::triangle_area_deriv({p[0], p[1], p[4]}, {dz[0], dz[1], dz[4]});
            A += util::triangle_area_deriv({p[0], p[3], p[4]}, {dz[0], dz[3], dz[4]});

            if(blockw == 3){
                A += util::triangle_area_deriv({p[1], p[2], p[4]}, {dz[1], dz[2], dz[4]});
                A += util::triangle_area_deriv({p[2], p[4], p[5]}, {dz[2], dz[4], dz[5]});
            }
            if(blockh == 3){
                A += util::triangle_area_deriv({p[3], p[4], p[6]}, {dz[3], dz[4], dz[6]});
                A += util::triangle_area_deriv({p[4], p[6], p[7]}, {dz[4], dz[6], dz[7]});
            }
            if(blockw == 3 && blockh == 3){
                A += util::triangle_area_deriv({p[4], p[5], p[8]}, {dz[4], dz[5], dz[8]});
                A += util::triangle_area_deriv({p[4], p[7], p[8]}, {dz[4], dz[7], dz[8]});
            }
        }
    }

    return A;

}

}
