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

#include "render.hpp"
#include "param.hpp"

namespace render{

void draw_texture(std::vector<sf::Uint8>& img, const std::vector<double>& map_z, double ap, size_t w, size_t h, Colorscheme colorscheme){
    using namespace param;

    if(colorscheme == Colorscheme::GRAYSCALE){
        for(size_t x = 0; x < tex_width; ++x){
            for(size_t y = 0; y < tex_height; ++y){
                double z = map_z[tex_width*y + x];

                sf::Uint8 rgb = (sf::Uint8)std::round(255*(1 + (z-max_z)/(max_z - min_z)));
                img[(tex_width*y + x)*4+0] = rgb;
                img[(tex_width*y + x)*4+1] = rgb;
                img[(tex_width*y + x)*4+2] = rgb;
                img[(tex_width*y + x)*4+3] = 255;
            }
        }
    } else if(colorscheme == Colorscheme::HSV){
        for(size_t x = 0; x < tex_width; ++x){
            for(size_t y = 0; y < tex_height; ++y){
                double z = map_z[tex_width*y + x];

                double norm = -(z-max_z)/(max_z - min_z);
                norm = std::max(0.0, std::min(1.0, norm));
                double r, g, b;
                if(norm <= 0.25){
                    r = 1;
                    g = norm*4;
                    b = 0;
                } else if(norm <= 0.5){
                    r = 1 - (norm-0.25)*4;
                    g = 1;
                    b = 0;
                } else if(norm <= 0.75){
                    r = 0;
                    g = 1;
                    b = (norm-0.5)*4;
                } else {
                    r = 0;
                    g = 1 - (norm-0.75)*4;
                    b = 1;
                }

                img[(tex_width*y + x)*4+0] = (sf::Uint8)std::round(255*r);
                img[(tex_width*y + x)*4+1] = (sf::Uint8)std::round(255*g);
                img[(tex_width*y + x)*4+2] = (sf::Uint8)std::round(255*b);
                img[(tex_width*y + x)*4+3] = 255;
            }
        }

    }
}

}
