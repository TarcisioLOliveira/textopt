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
#include <numeric>
#include <MMASolver.hpp>
#include <algorithm>

size_t tex_width = 300;
size_t tex_height = 300;

const double EXP = 100;

double dimx = 1;
double dimy = 1;
double dimz = 1;

double alpha = 60*M_PI/180;
double r = 2; // um
double max_z = 0; // Used to calculate Sa, so uses smooth_min
double min_z = 0; // Used only for texture display, so uses std::min()
double dmax_zdf = 0; 
                  
struct Point{
    double x, y, z;
};

struct Vector{
    double x, y, z;
};

enum Colorscheme{
    GRAYSCALE,
    HSV
};

double triangle_area(std::array<Point, 3> p){
    Vector v1{p[1].x - p[0].x, p[1].y - p[0].y, p[1].z - p[0].z};
    Vector v2{p[2].x - p[0].x, p[2].y - p[0].y, p[2].z - p[0].z};

    double A = 0.5*std::sqrt(
        std::pow(v1.y*v2.z - v2.y*v1.z, 2) +
        std::pow(v1.x*v2.z - v2.x*v1.z, 2) +
        std::pow(v1.y*v2.x - v2.y*v1.x, 2)
    );

    return A;
}


double triangle_area_deriv(std::array<Point, 3> p, std::array<double, 3> dz){
    Vector v1{p[1].x - p[0].x, p[1].y - p[0].y, p[1].z - p[0].z};
    Vector v2{p[2].x - p[0].x, p[2].y - p[0].y, p[2].z - p[0].z};
    double dv1 = dz[1] - dz[0];
    double dv2 = dz[2] - dz[0];

    double A = (0.25/std::sqrt(
        std::pow(v1.y*v2.z - v2.y*v1.z, 2) +
        std::pow(v1.x*v2.z - v2.x*v1.z, 2) +
        std::pow(v1.y*v2.x - v2.y*v1.x, 2)
    ))*(2*(v1.y*v2.z - v2.y*v1.z)*(v1.y*dv2 - v2.y*dv1) +
        2*(v1.x*v2.z - v2.x*v1.z)*(v1.x*dv2 - v2.x*dv1) +
        0
    );

    return A;
}

double smooth_min(std::initializer_list<double> x){
    double sum = 0;
    for(double xi : x){
        sum += std::pow(-xi, EXP);
    }

    return -std::pow(sum, 1.0/EXP);
};

double smooth_min_deriv(std::initializer_list<double> x, double xx){
    double sum = 0;
    for(double xi : x){
        sum += std::pow(-xi, EXP);
    }
    double result = -std::pow(sum, 1.0/EXP - 1);
    result *= -std::pow(-xx, EXP-1);

    return result;
};

void texture_map(double*& map_z, double* orig_z, double f, double ap, double vc){
    // double dp1 = ap/std::tan(alpha);
    // double p1 = tex_height/2 - dp1;
    // double p4 = tex_height/2 + dp1;

    // std::cout << p1 << " " << p4 << std::endl;
    double line_root = ap - std::tan(alpha)*f/2;
    max_z = smooth_min({0, -line_root});
    min_z = 0;
    for(size_t x = 0; x < tex_width; ++x){
        size_t mult = 1;
        for(size_t y = 0; y < tex_height; ++y){
            double& z = map_z[tex_width*y + x];
            double oz = orig_z[tex_width*y + x];
            if(y < mult*f){
                if(y <= (mult-1)*f + f/2){
                    z = smooth_min({oz, -std::tan(alpha)*(y - (mult-1)*f) - line_root});
                } else {
                    z = smooth_min({oz, std::tan(alpha)*(y - ((mult-1)*f + f)) - line_root});
                }
                min_z = std::min(min_z, z);
            } else {
                ++mult;
                --y;
            }
        }
    }
}


void dzdf(double*& map_z, double* orig_z, double f, double ap, double vc, double*& dzdf){
    double line_root = ap - std::tan(alpha)*f/2;
    double dlrdf = -std::tan(alpha)/2;
    max_z = smooth_min({0, -line_root});
    dmax_zdf = smooth_min_deriv({0, -line_root}, -line_root)*dlrdf;
    min_z = 0;
    for(size_t x = 0; x < tex_width; ++x){
        size_t mult = 1;
        for(size_t y = 0; y < tex_height; ++y){
            double& z = map_z[tex_width*y + x];
            double oz = orig_z[tex_width*y + x];
            double& dz = dzdf[tex_width*y + x];
            if(y < mult*f){
                if(y <= (mult-1)*f + f/2){
                    // z = smooth_min({z, -std::tan(alpha)*(y - (mult-1)*f) - line_root});
                    double znew = -std::tan(alpha)*(y - (mult-1)*f) - line_root;
                    double dznewdf = std::tan(alpha)*(mult-1) - dlrdf;
                    dz = smooth_min_deriv({oz, znew}, znew)*dznewdf;
                } else {
                    // z = smooth_min({z, std::tan(alpha)*(y - ((mult-1)*f + f)) - line_root});
                    double znew = -std::tan(alpha)*(y - (mult-1)*f + f) - line_root;
                    double dznewdf = std::tan(alpha)*mult - dlrdf;
                    dz = smooth_min_deriv({oz, znew}, znew)*dznewdf;
                    if(std::isnan(dz)){
                        std::cout << "here" << std::endl;
                        dz = smooth_min_deriv({oz, znew}, znew)*dznewdf;
                    }
                }
            } else {
                ++mult;
                --y;
            }
        }
    }
}

void draw_texture(sf::Uint8*& img, double* map_z, double ap, size_t w, size_t h, Colorscheme colorscheme){
    if(colorscheme == Colorscheme::GRAYSCALE){
        for(size_t x = 0; x < tex_width; ++x){
            for(size_t y = 0; y < tex_height; ++y){
                double& z = map_z[tex_width*y + x];

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
                double& z = map_z[tex_width*y + x];

                double norm = -(z-max_z)/(max_z - min_z);
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

double Sa(double*& map_z){
    size_t area = tex_width*tex_height;
    double sum = std::abs(std::accumulate(map_z, map_z+area, 0.0, std::plus<double>()));
    sum += area*max_z; // Depth correction

    return sum/area;
}

double dSa(double*& map_z, double* dzd){
    size_t area = tex_width*tex_height;
    double sum = std::abs(std::accumulate(dzd, dzd+area, 0.0, std::plus<double>()));
    sum += area*dmax_zdf; // Depth correction

    return sum/area;
}

double surface_area(double*& map_z){
    // For a N*M matrix of points, the actual projected surface area is (dimx*dimy)*((N-1)*(M-1))
    double A = 0;
    for(size_t x = 0; x < tex_width; x+=2){
        for(size_t y = 0; y < tex_height; y+=2){
            Point p[9];
            size_t blockw = std::min(3ul, tex_width-x);
            size_t blockh = std::min(3ul, tex_height-y);
            if(blockw <= 1 || blockh <= 1){
                continue;
            }
            for(size_t i = 0; i < blockw; ++i){
                for(size_t j = 0; j < blockh; ++j){
                    p[i*3+j] = Point{
                        dimx*(x+i),
                        dimy*(y+j),
                        dimz*map_z[tex_width*(y+j) + (x+i)]
                    };
                }
            }

            A += triangle_area({p[0], p[1], p[4]});
            A += triangle_area({p[0], p[3], p[4]});

            if(blockw == 3){
                A += triangle_area({p[1], p[2], p[4]});
                A += triangle_area({p[2], p[4], p[5]});
            }
            if(blockh == 3){
                A += triangle_area({p[3], p[4], p[6]});
                A += triangle_area({p[4], p[6], p[7]});
            }
            if(blockw == 3 && blockh == 3){
                A += triangle_area({p[4], p[5], p[8]});
                A += triangle_area({p[4], p[7], p[8]});
            }
        }
    }

    return A;
}

double surface_area_dz(double*& map_z, double* dzd){
    double A = 0;
    for(size_t x = 0; x < tex_width; x+=2){
        for(size_t y = 0; y < tex_height; y+=2){
            Point p[9];
            double dz[9];
            size_t blockw = std::min(3ul, tex_width-x);
            size_t blockh = std::min(3ul, tex_height-y);
            if(blockw <= 1 || blockh <= 1){
                continue;
            }
            for(size_t i = 0; i < blockw; ++i){
                for(size_t j = 0; j < blockh; ++j){
                    p[i*3+j] = Point{
                        dimx*(x+i),
                        dimy*(y+j),
                        dimz*map_z[tex_width*(y+j) + (x+i)]
                    };
                    dz[i*3+j] = dimz*dzd[tex_width*(y+j) + (x+i)];
                }
            }

            A += triangle_area_deriv({p[0], p[1], p[4]}, {dz[0], dz[1], dz[4]});
            A += triangle_area_deriv({p[0], p[3], p[4]}, {dz[0], dz[3], dz[4]});

            if(blockw == 3){
                A += triangle_area_deriv({p[1], p[2], p[4]}, {dz[1], dz[2], dz[4]});
                A += triangle_area_deriv({p[2], p[4], p[5]}, {dz[2], dz[4], dz[5]});
            }
            if(blockh == 3){
                A += triangle_area_deriv({p[3], p[4], p[6]}, {dz[3], dz[4], dz[6]});
                A += triangle_area_deriv({p[4], p[6], p[7]}, {dz[4], dz[6], dz[7]});
            }
            if(blockw == 3 && blockh == 3){
                A += triangle_area_deriv({p[4], p[5], p[8]}, {dz[4], dz[5], dz[8]});
                A += triangle_area_deriv({p[4], p[7], p[8]}, {dz[4], dz[7], dz[8]});
            }
        }
    }

    return A;

}

int main(int argc, char* argv[]){

    size_t window_width = 800;
    size_t window_height = 600;

    double max_roughness = 6;

    sf::RenderWindow window(sf::VideoMode(window_width, window_height), "textopt");
    sf::Texture img;
    img.create(tex_width, tex_height);
    sf::Sprite sprite(img);
    sprite.setPosition(sf::Vector2f(window_width/2-tex_width/2, window_height/2-tex_height/2));

    // std::default_random_engine generator;
    // std::uniform_int_distribution<sf::Uint8> rand(0);
    sf::Uint8* px = new sf::Uint8[tex_width*tex_height*4]();
    double* map_z = new double[tex_width*tex_height]();
    double* orig_z = new double[tex_width*tex_height]();
    double* df = new double[tex_width*tex_height]();

    double ch = 1;
    size_t it = 1;
    double old_surarea = 1;

    const size_t N = 3;
    double x[N] = {30, 40, 10};
    double xmax[N] = {100, 100, 100};
    double xmin[N] = {0.001, 0.001, 0.001};
    double dSa_vec[N] = {0, 0, 0};
    double dsurarea_vec[N] = {0, 0, 0};

    MMASolver mma(1, 1, 0, 1e4, 1);
    mma.SetAsymptotes(0.001, 0.7, 1.001);

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
            double f = x[0];
            double ap = x[1];
            double vc = x[2];
            texture_map(map_z, orig_z, f, ap, vc);
            draw_texture(px, map_z, ap, tex_width, tex_height, Colorscheme::HSV);
            double surarea = -surface_area(map_z);
            double roughness = Sa(map_z) - max_roughness;

            dzdf(map_z, orig_z, f, ap, vc, df);
            double dsurareadf = -surface_area_dz(map_z, df);
            double dSadf = dSa(map_z, df);
            std::cout << dsurareadf << " " << dSadf << std::endl;

            dsurarea_vec[0] = dsurareadf;
            dSa_vec[0] = dSadf;

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

        img.update(px);
        auto wsize = window.getSize();
        window_width = wsize.x;
        window_height = wsize.y;
        sprite.setPosition(sf::Vector2f(window_width/2-tex_width/2, window_height/2-tex_height/2));

        window.draw(sprite);
        window.display();
    }

    delete[] px;
    delete[] map_z;
    delete[] orig_z;
    delete[] df;

    img.copyToImage().saveToFile("result.png");

    return 0;
}
