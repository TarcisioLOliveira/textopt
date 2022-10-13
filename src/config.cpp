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

#include "config.hpp"
#include "param.hpp"
#include <yaml-cpp/exceptions.h>
#include <yaml-cpp/node/impl.h>
#include <yaml-cpp/node/node.h>
#include <yaml-cpp/node/type.h>

namespace config{

void load(const std::string& path){
 
    YAML ::Node config;
    try{ 
        config = YAML::LoadFile(path);
    } catch(YAML::BadFile& e){
        std::cout << std::endl << "Error: inserted file path leads to a non existent file." << std::endl << std::endl;
        throw;
    }

    required(config, "analysis", YAML::NodeType::Map);
    auto at = get_scalar<std::string>(config["analysis"], "type");
    auto am = get_scalar<std::string>(config["analysis"], "method");
    param::opt_ap = get_scalar<bool>(config["analysis"], "opt_ap");
    if(at == "opt"){
        param::analysis_type = param::AnalysisType::OPT;
        if(am == "slp"){
            param::opt_method = param::OptMethod::SLP;
        } else if(am == "mma"){
            param::opt_method = param::OptMethod::MMA;
        } else {
            std::cout << std::endl << "Unknown method for " << at << ": " << am << std::endl << std::endl;
            throw;
        }
        param::max_roughness = get_scalar<double>(config["analysis"], "max_roughness");
        param::stop = get_scalar<double>(config["analysis"], "stop");
    } else if(at == "plot"){
        param::analysis_type = param::AnalysisType::PLOT;
        if(am == "fxap"){
            param::plot_method = param::PlotMethod::FXAP;
        } else if(am == "vc"){
            param::plot_method = param::PlotMethod::VC;
        } else {
            std::cout << std::endl << "Unknown method for " << at << ": " << am << std::endl << std::endl;
            throw;
        }
        param::step = get_scalar<double>(config["analysis"], "step");
    } else if(at == "single"){
        param::analysis_type = param::AnalysisType::SINGLE;
        if(am == "exact"){
            param::single_method = param::SingleMethod::EXACT;
        } else if(am == "smooth"){
            param::single_method = param::SingleMethod::SMOOTH;
        } else {
            std::cout << std::endl << "Unknown method for " << at << ": " << am << std::endl << std::endl;
            throw;
        }
    } else {
        std::cout << std::endl << "Unknown analysis type: " << at << std::endl << std::endl;
        throw;
    }

    required(config, "window", YAML::NodeType::Map);
    param::window_width = get_scalar<size_t>(config["window"], "window_width");
    param::window_height = get_scalar<size_t>(config["window"], "window_height");

    required(config, "texture", YAML::NodeType::Map);
    param::tex_width = get_scalar<size_t>(config["texture"], "tex_width");
    param::tex_height = get_scalar<size_t>(config["texture"], "tex_height");
    param::base_area = (param::tex_width-1)*(param::tex_height-1);
    param::dim = get_scalar<double>(config["texture"], "dim");
    param::dimx = param::dimy = param::dimz = param::dim;

    required(config, "cutting", YAML::NodeType::Map);
    if(param::analysis_type != param::AnalysisType::PLOT){
        param::f  = get_scalar<double>(config["cutting"], "f" );
        param::ap = get_scalar<double>(config["cutting"], "ap");
        param::vc = get_scalar<double>(config["cutting"], "vc");
    }
    if(param::analysis_type != param::AnalysisType::SINGLE){
        required(config["cutting"], "f_range", YAML::NodeType::Sequence);
        required(config["cutting"], "ap_range", YAML::NodeType::Sequence);
        required(config["cutting"], "vc_range", YAML::NodeType::Sequence);

        bool full_range = true;
        if(param::analysis_type == param::AnalysisType::OPT &&
           param::opt_method == param::OptMethod::SLP){
            required_size(config["cutting"], "f_range", 1);
            required_size(config["cutting"], "ap_range", 1);
            required_size(config["cutting"], "vc_range", 1);
            full_range = false;
        } else {
            required_size(config["cutting"], "f_range", 2);
            required_size(config["cutting"], "ap_range", 2);
            required_size(config["cutting"], "vc_range", 2);
        }

        auto fit = config["cutting"]["f_range"].begin();
        auto apit = config["cutting"]["ap_range"].begin();
        auto vcit = config["cutting"]["vc_range"].begin();
        param::f_min  = fit->as<double>();
        param::ap_min = apit->as<double>();
        param::vc_min = vcit->as<double>();
        if(full_range){
            ++fit;
            ++apit;
            ++vcit;
            param::f_max  = fit->as<double>();
            param::ap_max = apit->as<double>();
            param::vc_max = vcit->as<double>();
        }
    }

    required(config, "approx_constants", YAML::NodeType::Map);
    param::MULT = get_scalar<double>(config["approx_constants"], "MULT");
    param::FLOORD = get_scalar<double>(config["approx_constants"], "FLOORD");
    param::ABS_EPS = get_scalar<double>(config["approx_constants"], "ABS_EPS");
    param::YOFF = get_scalar<double>(config["approx_constants"], "YOFF");

    required(config, "tool", YAML::NodeType::Map);
    param::r = get_scalar<double>(config["tool"], "r") / param::dim;
    param::alpha1 = get_scalar<double>(config["tool"], "alpha1") * M_PI / 180.0;
    param::alpha2 = get_scalar<double>(config["tool"], "alpha2") * M_PI / 180.0;

    {
        using namespace param;
        using param::y1;

        tan1 = std::tan(param::alpha1);
        tan2 = std::tan(param::alpha2);

        y1 = -tan1*r/std::sqrt(tan1*tan1+1);
        y2 =  tan2*r/std::sqrt(tan2*tan2+1);

        b1off = -std::sqrt(r*r - y1*y1) + r + tan1*y1;
        b2off = -std::sqrt(r*r - y2*y2) + r - tan2*y2;
    }

    required(config, "oscillation", YAML::NodeType::Map);
    param::Az = get_scalar<double>(config["oscillation"], "Az");
    param::fz = get_scalar<double>(config["oscillation"], "fz");
    param::phiz = get_scalar<double>(config["oscillation"], "phiz") * M_PI / 180.0;
    param::cylinder_radius = 1000*get_scalar<double>(config["oscillation"], "cylinder_radius");

    required(config, "ellipse", YAML::NodeType::Map);
    param::f_uet = get_scalar<double>(config["ellipse"], "f_uet");
    param::Ax_uet = get_scalar<double>(config["ellipse"], "Ax_uet");
    param::Az_uet = get_scalar<double>(config["ellipse"], "Az_uet");
}

}
