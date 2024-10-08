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
#include "analysis.hpp"
#include "slp.hpp"
#include "config.hpp"
#include "texture_shallow.hpp"
#include "newton_opt.hpp"

int main(int argc, char* argv[]){
    using namespace param;

    if(argc > 1){
        config::load(argv[1]);
    } else {
        std::cout << "Error: missing path to configuration file." << std::endl;
        throw;
    }

    std::vector<double> map_z  (tex_width*tex_height);
    std::vector<double> orig_z (tex_width*tex_height);

    if(analysis_type == AnalysisType::PLOT){
        if(plot_method == PlotMethod::FXAP){
            if(opt_f){
                analysis::plot_fxap({f_min, f_max}, {ap_min, ap_max}, vc, step, map_z, orig_z);
            } else {
                analysis::plot_ap_shallow({ap_min, ap_max}, vc, step, map_z, orig_z);
            }
        } else if(plot_method == PlotMethod::VC){
            if(opt_f){
                analysis::plot_vc(f, ap, {vc_min, vc_max}, step, map_z, orig_z);
            } else {
                analysis::plot_vc_shallow(ap, {vc_min, vc_max}, step, map_z, orig_z);
            }
        }
    } else {
        std::vector<double> df (tex_width*tex_height);
        std::vector<double> dap(tex_width*tex_height);
        std::vector<double> dvc(tex_width*tex_height);

        auto resolution = sf::VideoMode::getDesktopMode();

        sf::Texture img;
        img.create(tex_width, tex_height);
        sf::Sprite sprite(img);
        sprite.setPosition(sf::Vector2f(window_width/2-tex_width/2, window_height/2-tex_height/2));

        std::vector<sf::Uint8> px(tex_width*tex_height*4, 255);

        sf::RenderWindow window(sf::VideoMode(window_width, window_height), "textopt");
        window.setPosition(sf::Vector2i((resolution.width - window_width)/2, (resolution.height - window_height)/2));
        
        if(analysis_type == AnalysisType::SINGLE){
            std::vector<double> dSa_vec{0, 0, 0};
            std::vector<double> dsurarea_vec{0, 0, 0};

            std::cout << std::endl;
            std::cout << "===========================" << std::endl;
            std::cout << "========   EXACT    =======" << std::endl;
            std::cout << "===========================" << std::endl;

            if(param::opt_f){
                texture::map_exact(map_z, orig_z, f, ap, vc);
            } else {
                texture_shallow::map_exact(map_z, orig_z, ap, vc);
            }

            if(param::single_method == param::SingleMethod::EXACT){
                std::cout << "== (CURRENTLY RENDERED)  ==" << std::endl;
                std::cout << "===========================" << std::endl;
                render::draw_texture(px, map_z, render::Colorscheme::HSV);
            }
            const double surarea_ex = -opt::surface_area(map_z);
            const double roughness_ex = opt::Sa(map_z) - max_roughness;

            std::cout << std::endl;
            std::cout << "Surface area: " << -surarea_ex << std::endl;
            std::cout << "Roughness: " << roughness_ex + max_roughness << std::endl;
            std::cout << std::endl;
            std::cout << "f: " << f << std::endl;
            std::cout << "ap: " << ap << std::endl;
            std::cout << "vc: " << vc << std::endl;
            std::cout << std::endl;
            std::cout << "max z: " << param::max_z << " (red)" << std::endl;
            std::cout << "mid z: " << (param::max_z+param::min_z)/2.0 << " (green)" << std::endl;
            std::cout << "min z: " << param::min_z << " (blue)" << std::endl;
            std::cout << std::endl;
            std::cout << "avg z: " << param::z_avg << std::endl;

            std::cout << std::endl;
            std::cout << "===========================" << std::endl;
            std::cout << "========   SMOOTH   =======" << std::endl;
            std::cout << "===========================" << std::endl;

            if(param::opt_f){
                texture::map(map_z, orig_z, f, ap, vc);
            } else {
                texture_shallow::map(map_z, orig_z, ap, vc);
            }
            if(param::single_method == param::SingleMethod::SMOOTH){
                std::cout << "== (CURRENTLY RENDERED)  ==" << std::endl;
                std::cout << "===========================" << std::endl;
                render::draw_texture(px, map_z, render::Colorscheme::HSV);
            }
            const double surarea = -opt::surface_area(map_z);
            const double roughness = opt::Sa(map_z) - max_roughness;

            if(param::opt_f){
                texture::dzdf(orig_z, f, ap, vc, df);
                texture::dzdap(orig_z, f, ap, vc, dap);
                texture::dzdvc(orig_z, f, ap, vc, dvc);

                const double dsurareadf = -opt::surface_area_dz(map_z, df);
                const double dSadf = opt::dSa(df, map_z);
                const double dsurareadap = -opt::surface_area_dz(map_z, dap);
                const double dSadap = opt::dSa(dap, map_z);
                const double dsurareadvc = -opt::surface_area_dz(map_z, dvc);
                const double dSadvc = opt::dSa(dvc, map_z);

                std::cout << std::endl;
                std::cout << "dAdf:   " << dsurareadf << std::endl;
                std::cout << "dSadf:  " << dSadf << std::endl;
                std::cout << std::endl;
                std::cout << "dAdap:  " << dsurareadap << std::endl;
                std::cout << "dSadap: " << dSadap << std::endl;
                std::cout << std::endl;
                std::cout << "dAdvc:  " << dsurareadvc << std::endl;
                std::cout << "dSadvc: " << dSadvc << std::endl;
            } else {
                texture_shallow::dzdap(orig_z, ap, vc, dap);
                texture_shallow::dzdvc(orig_z, ap, vc, dvc);

                const double dsurareadap = -opt::surface_area_dz(map_z, dap);
                const double dSadap = opt::dSa(dap, map_z);
                const double dsurareadvc = -opt::surface_area_dz(map_z, dvc);
                const double dSadvc = opt::dSa(dvc, map_z);

                std::cout << std::endl;
                std::cout << "dAdap:   " << dsurareadap << std::endl;
                std::cout << "dSadap:  " << dSadap << std::endl;
                std::cout << std::endl;
                std::cout << "dAdvc:  " << dsurareadvc << std::endl;
                std::cout << "dSadvc: " << dSadvc << std::endl;
            }

            std::cout << std::endl;
            std::cout << "Surface area: " << -surarea << std::endl;
            std::cout << "Roughness: " << roughness + max_roughness << std::endl;
            std::cout << std::endl;
            std::cout << "f: " << f << std::endl;
            std::cout << "ap: " << ap << std::endl;
            std::cout << "vc: " << vc << std::endl;
            std::cout << std::endl;
            std::cout << "max z: " << param::max_z << " (red)" << std::endl;
            std::cout << "mid z: " << (param::max_z+param::min_z)/2.0 << " (green)" << std::endl;
            std::cout << "min z: " << param::min_z << " (blue)" << std::endl;
            std::cout << std::endl;
            std::cout << "avg z: " << param::z_avg << std::endl;

            std::cout << std::endl;
            std::cout << "===========================" << std::endl;
            std::cout << "=====  DIFFERENCE   =======" << std::endl;
            std::cout << "===========================" << std::endl;

            std::cout << std::endl;
            std::cout << "Surface area: " << -surarea+surarea_ex << std::endl;
            std::cout << "Roughness: " << roughness - roughness_ex << std::endl;

            std::cout << std::endl;

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

                img.update(px.data());
                auto wsize = window.getSize();
                window_width = wsize.x;
                window_height = wsize.y;
                sprite.setPosition(sf::Vector2f(window_width/2-tex_width/2, window_height/2-tex_height/2));

                window.draw(sprite);
                window.display();
            }

            img.copyToImage().saveToFile("result.png");

        } else if(analysis_type == AnalysisType::OPT){
            if(param::opt_f){
                double ch = 1;
                size_t it = 1;
                double old_surarea = 1;

                const size_t N = 3;
                std::vector<double> x{f, ap, vc};
                std::vector<double> xmin{f_min, ap_min, vc_min};
                std::vector<double> xmax{f_max, ap_max, vc_max};
                // Workaround. When the angles are different, you can't assume that both
                // edges have the same height when at least one side is within the tool
                // radius zone. Not sure if it's worth it to adapt to this case, as it
                // may involve a larger refactor, but this will do for now.
                if(alpha1 != alpha2){
                    xmin[0] = std::max(2*r, xmin[0]);
                }
                std::vector<double> dSa_vec{0, 0, 0};
                std::vector<double> dsurarea_vec{0, 0, 0};

                // Only works if starting from an exterior point, for some reason
                SLP slp(N, 1, xmin, xmax);
                MMASolver mma(N, 1, 0, 1e7, 1);
                mma.SetAsymptotes(0.1, 0.7, 1.2);
                NewtonOpt nopt(N, 1);

                bool printed_finished = false;

                double surarea = 0;
                double roughness = 0;
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

                    if(ch > stop){
                        f = x[0];
                        ap = x[1];
                        vc = x[2];
                        texture::map(map_z, orig_z, f, ap, vc);
                        render::draw_texture(px, map_z, render::Colorscheme::HSV);
                        surarea = -opt::surface_area(map_z);
                        roughness = opt::Sa(map_z) - max_roughness;

                        texture::dzdf(orig_z, f, ap, vc, df);
                        texture::dzdap(orig_z, f, ap, vc, dap);
                        texture::dzdvc(orig_z, f, ap, vc, dvc);

                        const double dsurareadf = -opt::surface_area_dz(map_z, df);
                        const double dSadf = opt::dSa(df, map_z);
                        const double dsurareadap = -opt::surface_area_dz(map_z, dap);
                        const double dSadap = opt::dSa(dap, map_z);
                        const double dsurareadvc = -opt::surface_area_dz(map_z, dvc);
                        const double dSadvc = opt::dSa(dvc, map_z);

                        std::cout << std::endl;
                        std::cout << "===========================" << std::endl;
                        std::cout << "Iteration: " << it << std::endl;
                        std::cout << "===========================" << std::endl;

                        std::cout << std::endl;
                        std::cout << "dAdf:   " << dsurareadf << std::endl;
                        std::cout << "dSadf:  " << dSadf << std::endl;
                        std::cout << std::endl;
                        std::cout << "dAdap:  " << dsurareadap << std::endl;
                        std::cout << "dSadap: " << dSadap << std::endl;
                        std::cout << std::endl;
                        std::cout << "dAdvc:  " << dsurareadvc << std::endl;
                        std::cout << "dSadvc: " << dSadvc << std::endl;

                        dsurarea_vec[0] = dsurareadf;
                        dSa_vec[0] = dSadf;
                        dsurarea_vec[1] = dsurareadap;
                        dSa_vec[1] = dSadap;
                        dsurarea_vec[2] = dsurareadvc;
                        dSa_vec[2] = dSadvc;

                        if(opt_method == OptMethod::SLP){
                            slp.update(x, dsurarea_vec, {roughness}, dSa_vec);
                        } else if(opt_method == OptMethod::MMA){
                            mma.Update(x.data(), dsurarea_vec.data(), &roughness, dSa_vec.data(), xmin.data(), xmax.data());
                        } else if(opt_method == OptMethod::NEWTON){
                            nopt.update(x, surarea, dsurarea_vec, {roughness}, dSa_vec);
                        }

                        ch = std::abs(1 - surarea/old_surarea);
                        old_surarea = surarea;

                        std::cout << std::endl;
                        std::cout << "Surface area: " << -surarea << std::endl;
                        std::cout << "Roughness: " << roughness + max_roughness << std::endl;
                        std::cout << std::endl;
                        std::cout << "f: " << f << std::endl;
                        std::cout << "ap: " << ap << std::endl;
                        std::cout << "vc: " << vc << std::endl;
                        std::cout << std::endl;
                        std::cout << "max z: " << param::max_z << " (red)" << std::endl;
                        std::cout << "mid z: " << (param::max_z+param::min_z)/2.0 << " (green)" << std::endl;
                        std::cout << "min z: " << param::min_z << " (blue)" << std::endl;
                        std::cout << std::endl;
                        std::cout << "avg z: " << param::z_avg << std::endl;
                        std::cout << std::endl;
                        std::cout << "Change: " << ch << std::endl;
                        ++it;
                    } else if(!printed_finished){
                        printed_finished = true;
                        std::vector<double> dSa_vec{0, 0, 0};
                        std::vector<double> dsurarea_vec{0, 0, 0};

                        std::cout << std::endl;
                        std::cout << "===========================" << std::endl;
                        std::cout << "======   FINISHED   =======" << std::endl;
                        std::cout << "===========================" << std::endl;
                        std::cout << "===   EXACT RESULTS:   ====" << std::endl;
                        std::cout << "===========================" << std::endl;

                        texture::map_exact(map_z, orig_z, f, ap, vc);

                        const double surarea_ex = -opt::surface_area(map_z);
                        const double roughness_ex = opt::Sa(map_z) - max_roughness;

                        std::cout << std::endl;
                        std::cout << "Surface area: " << -surarea_ex << std::endl;
                        std::cout << "Roughness: " << roughness_ex + max_roughness << std::endl;
                        std::cout << std::endl;
                        std::cout << "f: " << f << std::endl;
                        std::cout << "ap: " << ap << std::endl;
                        std::cout << "vc: " << vc << std::endl;
                        std::cout << std::endl;
                        std::cout << "max z: " << param::max_z << " (red)" << std::endl;
                        std::cout << "mid z: " << (param::max_z+param::min_z)/2.0 << " (green)" << std::endl;
                        std::cout << "min z: " << param::min_z << " (blue)" << std::endl;
                        std::cout << std::endl;
                        std::cout << "avg z: " << param::z_avg << std::endl;

                        std::cout << std::endl;
                        std::cout << "===========================" << std::endl;
                        std::cout << "=====  DIFFERENCE   =======" << std::endl;
                        std::cout << "===========================" << std::endl;

                        std::cout << std::endl;
                        std::cout << "Surface area: " << -surarea+surarea_ex << std::endl;
                        std::cout << "Roughness: " << roughness - roughness_ex << std::endl;

                        std::cout << std::endl;
                    }

                    img.update(px.data());
                    auto wsize = window.getSize();
                    window_width = wsize.x;
                    window_height = wsize.y;
                    sprite.setPosition(sf::Vector2f(window_width/2-tex_width/2, window_height/2-tex_height/2));

                    window.draw(sprite);
                    window.display();
                }
            } else {
                double ch = 1;
                size_t it = 1;
                double old_surarea = 1;

                const size_t N = 2;
                std::vector<double> x{ap, vc};
                std::vector<double> xmin{ap_min, vc_min};
                std::vector<double> xmax{ap_max, vc_max};

                std::vector<double> dSa_vec{0, 0};
                std::vector<double> dsurarea_vec{0, 0};

                // Only works if starting from an exterior point, for some reason
                SLP slp(N, 1, xmin, xmax);
                MMASolver mma(N, 1, 0, 1e7, 1);
                mma.SetAsymptotes(0.1, 0.7, 1.2);

                bool printed_finished = false;

                double surarea = 0;
                double roughness = 0;
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

                    if(ch > stop){
                        ap = x[0];
                        vc = x[1];
                        texture_shallow::map(map_z, orig_z, ap, vc);
                        render::draw_texture(px, map_z, render::Colorscheme::HSV);
                        surarea = -opt::surface_area(map_z);
                        roughness = opt::Sa(map_z) - max_roughness;

                        texture_shallow::dzdap(orig_z, ap, vc, dap);
                        texture_shallow::dzdvc(orig_z, ap, vc, dvc);

                        const double dsurareadap = -opt::surface_area_dz(map_z, dap);
                        const double dSadap = opt::dSa(dap, map_z);
                        const double dsurareadvc = -opt::surface_area_dz(map_z, dvc);
                        const double dSadvc = opt::dSa(dvc, map_z);

                        std::cout << std::endl;
                        std::cout << "===========================" << std::endl;
                        std::cout << "Iteration: " << it << std::endl;
                        std::cout << "===========================" << std::endl;

                        std::cout << std::endl;
                        std::cout << "dAdap:   " << dsurareadap << std::endl;
                        std::cout << "dSadap:  " << dSadap << std::endl;
                        std::cout << std::endl;
                        std::cout << "dAdvc:  " << dsurareadvc << std::endl;
                        std::cout << "dSadvc: " << dSadvc << std::endl;

                        dsurarea_vec[0] = dsurareadap;
                        dSa_vec[0] = dSadap;
                        dsurarea_vec[1] = dsurareadvc;
                        dSa_vec[1] = dSadvc;

                        if(opt_method == OptMethod::SLP){
                            slp.update(x, dsurarea_vec, {roughness}, dSa_vec);
                        } else if(opt_method == OptMethod::MMA){
                            mma.Update(x.data(), dsurarea_vec.data(), &roughness, dSa_vec.data(), xmin.data(), xmax.data());
                        }

                        ch = std::abs(1 - surarea/old_surarea);
                        old_surarea = surarea;

                        std::cout << std::endl;
                        std::cout << "Surface area: " << -surarea << std::endl;
                        std::cout << "Roughness: " << roughness + max_roughness << std::endl;
                        std::cout << std::endl;
                        std::cout << "f: " << f << std::endl;
                        std::cout << "ap: " << ap << std::endl;
                        std::cout << "vc: " << vc << std::endl;
                        std::cout << std::endl;
                        std::cout << "max z: " << param::max_z << " (red)" << std::endl;
                        std::cout << "mid z: " << (param::max_z+param::min_z)/2.0 << " (green)" << std::endl;
                        std::cout << "min z: " << param::min_z << " (blue)" << std::endl;
                        std::cout << std::endl;
                        std::cout << "avg z: " << param::z_avg << std::endl;
                        std::cout << std::endl;
                        std::cout << "Change: " << ch << std::endl;
                        ++it;
                    } else if(!printed_finished){
                        printed_finished = true;
                        std::vector<double> dSa_vec{0, 0, 0};
                        std::vector<double> dsurarea_vec{0, 0, 0};

                        std::cout << std::endl;
                        std::cout << "===========================" << std::endl;
                        std::cout << "======   FINISHED   =======" << std::endl;
                        std::cout << "===========================" << std::endl;
                        std::cout << "===   EXACT RESULTS:   ====" << std::endl;
                        std::cout << "===========================" << std::endl;

                        texture_shallow::map_exact(map_z, orig_z, ap, vc);

                        const double surarea_ex = -opt::surface_area(map_z);
                        const double roughness_ex = opt::Sa(map_z) - max_roughness;

                        std::cout << std::endl;
                        std::cout << "Surface area: " << -surarea_ex << std::endl;
                        std::cout << "Roughness: " << roughness_ex + max_roughness << std::endl;
                        std::cout << std::endl;
                        std::cout << "f: " << f << std::endl;
                        std::cout << "ap: " << ap << std::endl;
                        std::cout << "vc: " << vc << std::endl;
                        std::cout << std::endl;
                        std::cout << "max z: " << param::max_z << " (red)" << std::endl;
                        std::cout << "mid z: " << (param::max_z+param::min_z)/2.0 << " (green)" << std::endl;
                        std::cout << "min z: " << param::min_z << " (blue)" << std::endl;
                        std::cout << std::endl;
                        std::cout << "avg z: " << param::z_avg << std::endl;

                        std::cout << std::endl;
                        std::cout << "===========================" << std::endl;
                        std::cout << "=====  DIFFERENCE   =======" << std::endl;
                        std::cout << "===========================" << std::endl;

                        std::cout << std::endl;
                        std::cout << "Surface area: " << -surarea+surarea_ex << std::endl;
                        std::cout << "Roughness: " << roughness - roughness_ex << std::endl;

                        std::cout << std::endl;
                    }

                    img.update(px.data());
                    auto wsize = window.getSize();
                    window_width = wsize.x;
                    window_height = wsize.y;
                    sprite.setPosition(sf::Vector2f(window_width/2-tex_width/2, window_height/2-tex_height/2));

                    window.draw(sprite);
                    window.display();
                }
            }

            img.copyToImage().saveToFile("result.png");
        }
    }

    return 0;
}
