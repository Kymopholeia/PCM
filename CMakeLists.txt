cmake_minimum_required(VERSION 3.26)
project(Lebedev)

set(CMAKE_CXX_STANDARD 17)

add_executable(Lebedev main.cpp
        Lebedev_grid.cpp
        Lebedev_grid.h
        BasisSets.cpp
        BasisSets.h
        json.hpp
        real_water.json
        Point_class.cpp
        Point_class.h
        Pcm.cpp
        Pcm.h)
