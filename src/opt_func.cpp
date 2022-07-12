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

#include <numeric>

#include "opt_func.hpp"
#include "param.hpp"
#include "util.hpp"

namespace opt{

double Sa(const std::vector<double>& map_z){
    using namespace param;

    size_t area = tex_width*tex_height;
    double sum = -std::accumulate(map_z.begin(), map_z.end(), 0.0, std::plus<double>());
    sum -= area*max_z; // Depth correction

    return sum/area;
}

double dSa(const std::vector<double>& dzd, double dmax){
    using namespace param;

    size_t area = tex_width*tex_height;
    double sum = -std::accumulate(dzd.begin(), dzd.end(), 0.0, std::plus<double>());
    sum -= area*dmax; // Depth correction

    return sum/area;
}

double surface_area(const std::vector<double>& map_z){
    using namespace param;

    // For a N*M matrix of points, the actual projected surface area is (dimx*dimy)*((N-1)*(M-1))
    double A = 0;
    for(size_t x = 0; x < tex_width; x+=2){
        for(size_t y = 0; y < tex_height; y+=2){
            util::Point p[9];
            size_t blockw = std::min(3ul, tex_width-x);
            size_t blockh = std::min(3ul, tex_height-y);
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
            size_t blockw = std::min(3ul, tex_width-x);
            size_t blockh = std::min(3ul, tex_height-y);
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
