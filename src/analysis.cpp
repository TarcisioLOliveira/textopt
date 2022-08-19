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

#include <fstream>
#include <sstream>
#include <iostream>
#include "analysis.hpp"
#include "texture.hpp"
#include "opt_func.hpp"

namespace analysis{

void plot_fxap(const std::tuple<double, double>& f, const std::tuple<double, double>& ap, const double vc, const double step, std::vector<double>& map_z, const std::vector<double>& orig_z){
    std::stringstream result;
    const double diff = std::get<1>(f) - std::get<0>(f);
    std::cout << "Calculating plot points..." << std::endl;
    for(double fi = std::get<0>(f); fi <= std::get<1>(f); fi += step){
        std::cout << "\r" << (fi - std::get<0>(f))*100/diff << "%      " << std::flush;
        for(double api = std::get<0>(ap); api <= std::get<1>(ap); api += step){
            texture::map(map_z, orig_z, fi, api, vc);
            const double surarea = opt::surface_area(map_z);
            const double roughness = opt::Sa(map_z);
            result << surarea << " " << roughness << std::endl;
        }
    }
    std::cout << "\r100%      " << std::flush;
    std::ofstream file;
    file.open("plot_fxap.txt");
    file << std::get<0>(f) << " " << std::get<1>(f) << " " << step << std::endl;
    file << std::get<0>(ap) << " " << std::get<1>(ap) << " " << step << std::endl;
    file << result.str();
    file.close();
}

}
