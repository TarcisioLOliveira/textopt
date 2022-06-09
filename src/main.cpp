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

const double MULT = -1;

double dim_scale = 1e6; // m to um
double dim = 1;
double dimx = dim;
double dimy = dim;
double dimz = dim;

double alpha1 = 60*M_PI/180;
double alpha2 = 60*M_PI/180;
double r = 20/dim; // um
double max_z = 0; // Used to calculate Sa, so uses smooth_min
double min_z = 0; // Used only for texture display, so uses std::min()
double dmax_zdf = 0; 

// "Random" oscillation
// Models height variation along the tool path

// Amplitude [dim], e.g. [um]
double Ax = 10;
double Az = 5;
// Frequency [Hz]
double fx = 15*dim_scale/10;
double fz = 20*dim_scale/10;
// Phase [rad]
double phix = 20*M_PI/180;
double phiz = 0;
                  
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
    double frac_top = 0;
    double frac_bot = 0;
    for(double xi : x){
        double exp = std::exp(MULT*xi);
        frac_top += xi*exp;
        frac_bot += exp;
    }

    return frac_top/frac_bot;
};

double smooth_min_deriv(std::initializer_list<double> x, double xx){
    double frac_top = 0;
    double frac_bot = 0;
    for(double xi : x){
        double exp = std::exp(MULT*xi);
        frac_top += xi*exp;
        frac_bot += exp;
    }

    double sm = frac_top/frac_bot;

    return (std::exp(MULT*xx)/frac_bot)*(1+MULT*(xx-sm));
};

void texture_map(double*& map_z, double* orig_z, double f, double ap, double vc){
    double y1 = -std::sqrt((std::pow(std::tan(alpha1)*r, 2))/(std::pow(std::tan(alpha1), 2)+1));
    double y2 =  std::sqrt((std::pow(std::tan(alpha2)*r, 2))/(std::pow(std::tan(alpha2), 2)+1));

    double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + std::tan(alpha1)*(y1+f/2);//-ap + std::tan(alpha1)*f/2;
    double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - std::tan(alpha2)*(y2+f/2); //-ap - std::tan(alpha2)*f/2;

    // Calculate whether the angles of the cutting tool intercept the surface
    // first or are stopped by the feed rate.
    //
    // Check the intersection with the feed rate edges. If they're greater than
    // zero, the cut line is narrower than the feed line, so max height is
    // zero. Otherwise, it's the height at the intersection
    double intersec1 = 0;
    double intersec2 = 0;
    if(y1 + f/2 > 0){
        // If the radius is smaller than the feed rate, check for the line
        intersec1 = smooth_min({0, line_root1});
    } else {
        // Check for the circle
        intersec1 = smooth_min({0, -std::sqrt(r*r - f*f/4) + r - ap});
    }
    if(y2 + f/2 < f){
        // If the radius is smaller than the feed rate...
        intersec2 = smooth_min({0, line_root2 + std::tan(alpha2)*f});
    } else {
        intersec2 = smooth_min({0, -std::sqrt(r*r - f*f/4) + r - ap});
    }
    // Check if one of them is greater than zero
    double max_intersec = -smooth_min({-intersec1, -intersec2});
    // If one is greater than zero, max_z must be zero. Otherwise, it's
    // the lesser one.
    max_z = smooth_min({0, max_intersec});
    min_z = 0;
    for(size_t x = 0; x < tex_width; ++x){
        size_t mult = 1;
        for(size_t y = 0; y < tex_height; ++y){
            double& z = map_z[tex_width*y + x];
            double oz = orig_z[tex_width*y + x];
            double newx = x + Ax*std::sin(2*M_PI*fx*x*dimx/(vc*dim_scale) + phix);
            z = oz;
            if(y <= mult*f){
                // If it's within the area that's actually cut
                if(y >=  line_root1/std::tan(alpha1) + (mult-1)*f &&
                   y <= -line_root2/std::tan(alpha2) + (mult-1)*f){
                    if(y <= y1 + (mult-1)*f + f/2){
                        z = smooth_min({oz, -std::tan(alpha1)*((double)y - (mult-1)*f) + line_root1 + Az*std::sin(2*M_PI*fz*newx*dimx/(vc*dim_scale) + phiz)});
                    } else if(y <= y2 + (mult-1)*f + f/2){
                        z = smooth_min({oz, -std::sqrt(r*r - std::pow((double)y - (mult-1)*f - f/2, 2)) + r - ap + Az*std::sin(2*M_PI*fz*newx*dimx/(vc*dim_scale) + phiz)});
                    } else {
                        z = smooth_min({oz,  std::tan(alpha2)*((double)y - (mult-1)*f) + line_root2 + Az*std::sin(2*M_PI*fz*newx*dimx/(vc*dim_scale) + phiz)});
                    }
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
    double y1 = -std::sqrt((std::pow(std::tan(alpha1)*r, 2))/(std::pow(std::tan(alpha1), 2)+1));
    double y2 =  std::sqrt((std::pow(std::tan(alpha2)*r, 2))/(std::pow(std::tan(alpha2), 2)+1));

    double line_root1 = -std::sqrt(r*r - y1*y1) + r - ap + std::tan(alpha1)*(y1+f/2);//-ap + std::tan(alpha1)*f/2;
    double line_root2 = -std::sqrt(r*r - y2*y2) + r - ap - std::tan(alpha2)*(y2+f/2); //-ap - std::tan(alpha2)*f/2;

    double dlrdf1 =  std::tan(alpha1)/2;
    double dlrdf2 = -std::tan(alpha2)/2;

    double intersec1 = 0;
    double intersec2 = 0;
    double di1 = 0;
    double di2 = 0;
    if(y1 + f/2 > 0){
        intersec1 = smooth_min({0, line_root1});
        di1 = smooth_min_deriv({0, line_root1}, line_root1)*dlrdf1;
    } else {
        intersec1 = smooth_min({0, -std::sqrt(r*r - f*f/4) + r - ap});
        di1 = smooth_min_deriv({0, -std::sqrt(r*r - f*f/4) + r - ap},-std::sqrt(r*r - f*f/4) + r - ap)*f/(4*std::sqrt(r*r-f*f/4));
    }
    if(y2 + f/2 < f){
        intersec2 = smooth_min({0, line_root2 + std::tan(alpha2)*f});
        di2 = smooth_min_deriv({0, line_root2 + std::tan(alpha2)*f}, line_root2 + std::tan(alpha2)*f)*(dlrdf2 + std::tan(alpha2));
    } else {
        intersec2 = smooth_min({0, -std::sqrt(r*r - f*f/4) + r - ap});
        di2 = smooth_min_deriv({0, -std::sqrt(r*r - f*f/4) + r - ap},-std::sqrt(r*r - f*f/4) + r - ap)*f/(4*std::sqrt(r*r-f*f/4));
    }
    double max_intersec = -smooth_min({-intersec1, -intersec2});
    max_z = smooth_min({0, max_intersec});

    double dmi = -smooth_min_deriv({-intersec1, -intersec2}, -intersec1)*(-di1) - smooth_min_deriv({-intersec1, -intersec2}, -intersec2)*(-di2);
    dmax_zdf = smooth_min_deriv({0, max_intersec}, max_intersec)*dmi;
    min_z = 0;
    for(size_t x = 0; x < tex_width; ++x){
        size_t mult = 1;
        for(size_t y = 0; y < tex_height; ++y){
            double& z = map_z[tex_width*y + x];
            double oz = orig_z[tex_width*y + x];
            double& dz = dzdf[tex_width*y + x];
            double newx = x + Ax*std::sin(2*M_PI*fx*x*dimx/(vc*dim_scale) + phix);
            dz = 0;
            if(y <= mult*f){
                if(y >=  line_root1/std::tan(alpha1) + (mult-1)*f &&
                   y <= -line_root2/std::tan(alpha2) + (mult-1)*f){
                    if(y <= y1 + (mult-1)*f + f/2){
                        double znew = -std::tan(alpha1)*(y - (mult-1)*f) + line_root1 + Az*std::sin(2*M_PI*fz*newx*dimx/(vc*dim_scale) + phiz);
                        double dznewdf = -std::tan(alpha1)*(mult-1) + dlrdf1;
                        dz = smooth_min_deriv({oz, znew}, znew)*dznewdf;
                    } else if(y <= y2 + (mult-1)*f + f/2){
                        double znew = -std::sqrt(r*r - std::pow((double)y - (mult-1)*f - f/2, 2)) + r - ap + Az*std::sin(2*M_PI*fz*newx*dimx/(vc*dim_scale) + phiz);
                        double dznewdf = ((mult-1)+0.5)/std::sqrt(r*r - std::pow((double)y - (mult-1)*f - f/2, 2));
                        dz = smooth_min_deriv({oz, znew}, znew)*dznewdf;
                    } else {
                        double znew =  std::tan(alpha2)*(y - (mult-1)*f) + line_root2 + Az*std::sin(2*M_PI*fz*newx*dimx/(vc*dim_scale) + phiz);
                        double dznewdf =  std::tan(alpha2)*mult + dlrdf2;
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
    double sum = -std::accumulate(map_z, map_z+area, 0.0, std::plus<double>());
    sum -= area*max_z; // Depth correction

    return sum/area;
}

double dSa(double*& map_z, double* dzd){
    size_t area = tex_width*tex_height;
    double sum = -std::accumulate(dzd, dzd+area, 0.0, std::plus<double>());
    sum -= area*dmax_zdf; // Depth correction

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
                    p[j*3+i] = Point{
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
                    p[j*3+i] = Point{
                        dimx*(x+i),
                        dimy*(y+j),
                        dimz*map_z[tex_width*(y+j) + (x+i)]
                    };
                    dz[j*3+i] = dimz*dzd[tex_width*(y+j) + (x+i)];
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

    double max_roughness = 60;

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
    double x[N] = {20, 40, 10};
    double xmax[N] = {100, 100, 100};
    double xmin[N] = {0.001, 0.001, 0.001};
    double dSa_vec[N] = {0, 0, 0};
    double dsurarea_vec[N] = {0, 0, 0};

    MMASolver mma(1, 1, 0, 1e6, 1);
    mma.SetAsymptotes(0.001, 0.1, 1.001);

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
