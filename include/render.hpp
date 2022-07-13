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
/**
 * @file render.hpp
 *
 * @brief Texture rendering.
 *
 * Functions and utilities to render the calculated texture.
 */

#ifndef RENDER_HPP
#define RENDER_HPP

#include <SFML/Graphics.hpp>

namespace render{

/**
 * Available colorschemes.
 */
enum Colorscheme{
    GRAYSCALE,
    HSV
};

/**
 * Draws the texture using the desired colorscheme.
 * 
 * @param img Image vector.
 * @param map_z Texture.
 * @param colorscheme Colorscheme to be used.
 */
void draw_texture(std::vector<sf::Uint8>& img, const std::vector<double>& map_z, Colorscheme colorscheme);

}

#endif
