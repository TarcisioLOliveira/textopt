/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of textopt.
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with SolidPrep.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <SFML/Config.hpp>
#include <SFML/Graphics/Color.hpp>
#include <SFML/Graphics/Sprite.hpp>
#include <SFML/Graphics/Texture.hpp>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <random>
#include <cmath>
#include <iostream>
#include <initializer_list>

size_t tex_width = 300;
size_t tex_height = 300;

double alpha = 60*M_PI/180;
double r = 2; // um

double smooth_min(std::initializer_list<double> x){
    double mult = -8;
    double frac_top = 0;
    double frac_bot = 0;
    for(double xi : x){
        double exp = std::exp(mult*xi);
        frac_top += xi*exp;
        frac_bot += exp;
    }

    return frac_top/frac_bot;
};

void texture_map(double*& map_z, double f, double ap, double vc){
    // double dp1 = ap/std::tan(alpha);
    // double p1 = tex_height/2 - dp1;
    // double p4 = tex_height/2 + dp1;

    // std::cout << p1 << " " << p4 << std::endl;
    for(size_t x = 0; x < tex_width; ++x){
        size_t mult = 1;
        for(size_t y = 0; y < tex_height; ++y){
            double& z = map_z[tex_width*y + x];
            // if(y >= p1 && y <= p4){
            if(y < mult*f){
                if(y <= (mult-1)*f + f/2){
                    z = smooth_min({z, -std::tan(alpha)*(y - (mult-1)*f)});
                } else {
                    z = smooth_min({z, std::tan(alpha)*(y - ((mult-1)*f + f))});
                }
            } else {
                ++mult;
                --y;
                continue;
            }
            //}
        }
    }
}

void draw_texture(sf::Uint8*& img, double* map_z, double ap, size_t w, size_t h){
    for(size_t x = 0; x < tex_width; ++x){
        for(size_t y = 0; y < tex_height; ++y){
            double& z = map_z[tex_width*y + x];

            sf::Uint8 rgb = (sf::Uint8)std::round(255*(1 + z/ap));
            img[(tex_width*y + x)*4+0] = rgb;
            img[(tex_width*y + x)*4+1] = rgb;
            img[(tex_width*y + x)*4+2] = rgb;
            img[(tex_width*y + x)*4+3] = 255;
        }
    }
}

int main(int argc, char* argv[]){

    size_t window_width = 800;
    size_t window_height = 600;
    sf::RenderWindow window(sf::VideoMode(window_width, window_height), "textopt");
    sf::Texture img;
    img.create(tex_width, tex_height);
    sf::Sprite sprite(img);
    sprite.setPosition(sf::Vector2f(window_width/2-tex_width/2, window_height/2-tex_height/2));

    // std::default_random_engine generator;
    // std::uniform_int_distribution<sf::Uint8> rand(0);
    sf::Uint8* px = new sf::Uint8[tex_width*tex_height*4]();
    double* map_z = new double[tex_width*tex_height]();

    double f = 30;
    double ap = std::tan(alpha)*f/2;
    texture_map(map_z, f, ap, 10);
    draw_texture(px, map_z, ap, tex_width, tex_height);

    while(window.isOpen()){
        sf::Event event;
        while (window.pollEvent(event)){
            if (event.type == sf::Event::Closed){
                window.close();
            }
            if(event.type == sf::Event::Resized){
                sf::FloatRect view(0, 0, event.size.width, event.size.height);
                window.setView(sf::View(view));
            }
        }
        window.clear(sf::Color::Black);

        img.update(px);
        auto wsize = window.getSize();
        window_width = wsize.x;
        window_height = wsize.y;
        sprite.setPosition(sf::Vector2f(window_width/2-tex_width/2, window_height/2-tex_height/2));

        window.draw(sprite);
        window.display();
    }

    img.copyToImage().saveToFile("result.png");

    return 0;
}
