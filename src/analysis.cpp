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
    // Cache resulting string into memory, as it's faster than continually
    // writing to disk

    const double diff_f = std::get<1>(f) - std::get<0>(f);
    const double diff_ap = std::get<1>(ap) - std::get<0>(ap);
    const size_t size = (diff_f/step + 1)*(diff_ap/step + 1);
    std::vector<double> results(size*2);
    auto r = results.begin();

    std::cout << "Calculating plot points..." << std::endl;
    for(double fi = std::get<0>(f); fi <= std::get<1>(f)+1e-7; fi += step){
        std::cout << "\r" << (fi - std::get<0>(f))*100/diff_f << "%      " << std::flush;
        for(double api = std::get<0>(ap); api <= std::get<1>(ap)+1e-7; api += step){
            texture::map(map_z, orig_z, fi, api, vc);
            const double surarea = opt::surface_area(map_z);
            const double roughness = opt::Sa(map_z);
            *r = surarea;
            ++r;
            *r = roughness;
            ++r;
        }
    }
    std::cout << "\r100%      " << std::endl;
    std::string resstr;
    resstr.reserve(results.size()*25);
    std::stringstream result(resstr);
    r = results.begin();
    while(r < results.end()){
        result << *r << " " << *(r+1) << std::endl;
        r += 2;
    }

    std::ofstream file;
    file.open("plot_fxap.txt");
    // Write ranges
    file << std::get<0>(f) << " " << std::get<1>(f) << " " << step << std::endl;
    file << std::get<0>(ap) << " " << std::get<1>(ap) << " " << step << std::endl;

    // Write values for surface area and roughness
    file << result.str();
    file.close();
}

void plot_vc(const double f, const double ap, const std::tuple<double, double>& vc, const double step, std::vector<double>& map_z, const std::vector<double>& orig_z){
    // Cache resulting string into memory, as it's faster than continually
    // writing to disk
    const double diff = std::get<1>(vc) - std::get<0>(vc);
    const size_t size = (diff/step + 1);
    std::vector<double> results(size*2);
    auto r = results.begin();

    std::cout << "Calculating plot points..." << std::endl;
    for(double vci = std::get<0>(vc); vci <= std::get<1>(vc)+1e-7; vci += step){
        std::cout << "\r" << (vci - std::get<0>(vc))*100/diff << "%      " << std::flush;
        texture::map(map_z, orig_z, f, ap, vci);
        const double surarea = opt::surface_area(map_z);
        const double roughness = opt::Sa(map_z);
        *r = surarea;
        ++r;
        *r = roughness;
        ++r;
    }
    std::string resstr;
    resstr.reserve(results.size()*25);
    std::stringstream result(resstr);
    r = results.begin();
    while(r < results.end()){
        result << *r << " " << *(r+1) << std::endl;
        r += 2;
    }

    std::cout << "\r100%      " << std::endl;
    std::ofstream file;
    file.open("plot_vc.txt");
    // Write ranges
    file << std::get<0>(vc) << " " << std::get<1>(vc) << " " << step << std::endl;

    // Write values for surface area and roughness
    file << result.str();
    file.close();

}

}
