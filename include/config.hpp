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
 * @file config.hpp
 *
 * @brief Loads configuration file.
 *
 * Populates the global variables in param.hpp with parameters present in the
 * specified YAML configuration file.
 */

#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <iostream>
#include <string>
#include <yaml-cpp/exceptions.h>
#include <yaml-cpp/node/type.h>
#include <yaml-cpp/yaml.h>

namespace config{

/**
 * Loads optimization parameters from a file.
 *
 * @param path Path to YAML configuration file.
 */
void load(const std::string& path);

/**
 * Initializes param::dV()
 */
void init_dV();

inline bool valid(const YAML::Node& n, YAML::NodeType::value type){
    return n.Type() == type;
}

inline bool exists(const YAML::Node& n, const std::string& obj){
    try{
        return !n[obj].IsNull();
    } catch(YAML::InvalidNode& e){
        return false;
    }
    return true;
}

inline bool required(const YAML::Node& n, const std::string obj, YAML::NodeType::value type){
    if(!exists(n, obj)){
        std::cout << std::endl << "Missing parameter: " << obj << std::endl << std::endl;
        throw;
    } else if(!valid(n[obj], type)){
        std::cout << std::endl << "Parameter \"" << obj << "\" has invalid type" << std::endl << std::endl;
        throw;
    }
    return true;
}

inline bool required_size(const YAML::Node& n, const std::string obj, size_t size){
    required(n, obj, YAML::NodeType::Sequence);
    if(n[obj].size() < size){
        std::cout << std::endl << "Incorrect number of elements for sequence " << obj << ": has " << n[obj].size() << ", must be at least " << size << std::endl << std::endl;
        throw;
    }
    return true;
}

template<class T>
inline T get_scalar(const YAML::Node& n, const std::string obj){
    try{
        return n[obj].as<T>();
    } catch(YAML::InvalidNode& e){
        std::cout << std::endl << "Error: missing parameter: " << obj << std::endl << std::endl;
        throw;
    } catch(YAML::TypedBadConversion<T>& e){
        std::cout << std::endl << "Error: parameter has invalid type: " << obj << std::endl << std::endl;
        throw;
    }
    return 0;
}

}

#endif
