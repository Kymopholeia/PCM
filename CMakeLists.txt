cmake_minimum_required(VERSION 3.26)
project(Lebedev)

set(CMAKE_CXX_STANDARD 17)
find_package (Eigen3 5.0.1 REQUIRED NO_MODULE)

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
        Pcm.h
        ObaraSaika.cpp
        ObaraSaika.h
        BasicFunctions.cpp
        BasicFunctions.h
        ClassesStructs.cpp
        ClassesStructs.h
        Integrals.cpp
        Integrals.h
        Overlap.cpp
        Overlap.h
        Kinetic.cpp
        Kinetic.h
        Nuclear.cpp
        Nuclear.h
        hidrogen.json
        metane.json
        acetylene.json
        benzene_test.json
        DensityMatrix.cpp
        DensityMatrix.h
        FockMatrix.cpp
        FockMatrix.h
        SCF.cpp
        SCF.h)
target_link_libraries (Lebedev Eigen3::Eigen)