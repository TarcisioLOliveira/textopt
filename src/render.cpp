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

#include "render.hpp"
#include "param.hpp"
#include <exception>

namespace render{

void draw_texture(std::vector<sf::Uint8>& img, const std::vector<double>& map_z, Colorscheme colorscheme){
    using namespace param;

    if(colorscheme == Colorscheme::GRAYSCALE){
        for(size_t i = 0; i < tex_width*tex_height; ++i){
            const double z = map_z[i];

            const sf::Uint8 rgb = static_cast<sf::Uint8>(255*(1 - (max_z - z)/(max_z - min_z)));
            img[i*4+0] = rgb;
            img[i*4+1] = rgb;
            img[i*4+2] = rgb;
        }
    } else if(colorscheme == Colorscheme::HSV){
        // Uses the HSV colorscheme, going from red to blue passing through
        // green, imitating a common colorscheme used in surface topography.
        const std::vector<double> R{1, 1, 0, 0, 0};
        const std::vector<double> G{0, 1, 1, 1, 0};
        const std::vector<double> B{0, 0, 0, 1, 1};
        const size_t L = R.size()-1;
        for(size_t i = 0; i < tex_width*tex_height; ++i){
            const double z = map_z[i];

            double norm = (max_z - z)/(max_z - min_z);
            const double low = std::floor(norm*L);
            const double rem = norm*L - low;
            const double r = R[low] + rem*(R[low+1]-R[low]);
            const double g = G[low] + rem*(G[low+1]-G[low]);
            const double b = B[low] + rem*(B[low+1]-B[low]);

            img[i*4+0] = static_cast<sf::Uint8>(255*r);
            img[i*4+1] = static_cast<sf::Uint8>(255*g);
            img[i*4+2] = static_cast<sf::Uint8>(255*b);
        }
    }
}

}
