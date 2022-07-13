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

#include <SFML/Config.hpp>
#include <SFML/Graphics/Color.hpp>
#include <SFML/Graphics/Sprite.hpp>
#include <SFML/Graphics/Texture.hpp>
#include <SFML/System/Vector2.hpp>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <iostream>

#include "MMASolver.hpp"
#include "param.hpp"
#include "render.hpp"
#include "opt_func.hpp"
#include "texture.hpp"

int main(int argc, char* argv[]){
    using namespace param;

    size_t window_width = 600;
    size_t window_height = 600;

    double max_roughness = 60;

    auto resolution = sf::VideoMode::getDesktopMode();

    sf::RenderWindow window(sf::VideoMode(window_width, window_height), "textopt");
    window.setPosition(sf::Vector2i((resolution.width - window_width)/2, (resolution.height - window_height)/2));
    sf::Texture img;
    img.create(tex_width, tex_height);
    sf::Sprite sprite(img);
    sprite.setPosition(sf::Vector2f(window_width/2-tex_width/2, window_height/2-tex_height/2));

    std::vector<sf::Uint8> px(tex_width*tex_height*4);

    std::vector<double> map_z  (tex_width*tex_height);
    std::vector<double> orig_z (tex_width*tex_height);
    std::vector<double> df     (tex_width*tex_height);
    std::vector<double> dap    (tex_width*tex_height);
    std::vector<double> dvc    (tex_width*tex_height);

    double ch = 1;
    size_t it = 1;
    double old_surarea = 1;

    const size_t N = 3;
    double x[N] = {20, 40, 70};
    double xmax[N] = {100, 100, 100};
    double xmin[N] = {0.001, 0.001, 0.001};
    double dSa_vec[N] = {0, 0, 0};
    double dsurarea_vec[N] = {0, 0, 0};

    MMASolver mma(N, 1, 0, 1e6, 1);
    mma.SetAsymptotes(0.001, 0.7, 1.2);

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

        if(ch > 1e-8){
            const double f = x[0];
            const double ap = x[1];
            const double vc = x[2];
            texture::map(map_z, orig_z, f, ap, vc);
            render::draw_texture(px, map_z, ap, tex_width, tex_height, render::Colorscheme::HSV);
            const double surarea = -opt::surface_area(map_z);
            const double roughness = opt::Sa(map_z) - max_roughness;

            texture::dzdf(orig_z, f, ap, vc, df);
            texture::dzdap(orig_z, f, ap, vc, dap);
            texture::dzdvc(orig_z, f, ap, vc, dvc);

            const double dsurareadf = -opt::surface_area_dz(map_z, df);
            const double dSadf = opt::dSa(df, dmax_zdf);
            const double dsurareadap = -opt::surface_area_dz(map_z, dap);
            const double dSadap = opt::dSa(dap, dmax_zdap);
            const double dsurareadvc = -opt::surface_area_dz(map_z, dvc);
            const double dSadvc = opt::dSa(dvc, dmax_zdvc);

            std::cout << std::endl;
            std::cout << dsurareadf << " " << dSadf << std::endl;
            std::cout << dsurareadap << " " << dSadap << std::endl;
            std::cout << dsurareadvc << " " << dSadvc << std::endl;

            dsurarea_vec[0] = dsurareadf;
            dSa_vec[0] = dSadf;
            dsurarea_vec[1] = dsurareadap;
            dSa_vec[1] = dSadap;
            dsurarea_vec[2] = dsurareadvc;
            dSa_vec[2] = dSadvc;

            mma.Update(x, dsurarea_vec, &roughness, dSa_vec, xmin, xmax); 

            ch = std::abs(1 - surarea/old_surarea);
            old_surarea = surarea;

            std::cout << std::endl;
            std::cout << "Iteration: " << it << std::endl;
            std::cout << "Surface area: " << -surarea << std::endl;
            std::cout << "Roughness: " << roughness + max_roughness << std::endl;
            std::cout << "f: " << f << std::endl;
            std::cout << "ap: " << ap << std::endl;
            std::cout << "vc: " << vc << std::endl;
            std::cout << "Change: " << ch << std::endl;
            ++it;
        }

        img.update(px.data());
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
