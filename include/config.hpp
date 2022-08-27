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
#include <yaml-cpp/node/type.h>
#include <yaml-cpp/yaml.h>

namespace config{

void load(const std::string& path);

inline bool valid(const YAML::Node& n, YAML::NodeType::value type){
    return n.Type() == type;
}

inline bool null(const YAML::Node& n, YAML::NodeType::value type){
    return n.Type() == YAML::NodeType::Null;
}

inline bool required(const YAML::Node& n, const std::string obj, YAML::NodeType::value type){
    if(null(n[obj], type)){
        std::cout << "Missing parameter: " << obj << std::endl;
        throw;
    } else if(!valid(n[obj], type)){
        std::cout << "Parameter \"" << obj << "\" has invalid type" << std::endl;
        throw;
    }
    return true;
}

inline bool required_size(const YAML::Node& n, const std::string obj, size_t size){
    required(n, obj, YAML::NodeType::Sequence);
    if(n[obj].size() < size){
        std::cout << "Incorrect number of elements for sequence " << obj << ": has " << n[obj].size() << ", must be at least " << size << std::endl;
        throw;
    }
    return true;
}

template<class T>
inline T get_scalar(const YAML::Node& n, const std::string obj){
    if(required(n, obj, YAML::NodeType::Scalar)){
        return n[obj].as<T>();
    }
    return 0;
}

template<class T>
inline T get_array(const YAML::Node& n, const std::string obj){
    if(required(n, obj, YAML::NodeType::Sequence)){
        return n[obj].as<T>();
    }
}

}

#endif
