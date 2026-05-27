//
// Created by Domokos Kyra on 2026. 05. 26.
//

#ifndef LEBEDEV_LEBEDEV_GRID_H
#define LEBEDEV_LEBEDEV_GRID_H
#include <iostream>
#include <cmath>
#include <vector>
#include "BasisSets.h"
#include "Point_class.h"
using namespace std;

class Tessera{
public:
    Point r;
    double area;
    void cal_area();
};

class Lebedev_grid{
    const double A2=(1/sqrt(2));
    const double A3=(1/sqrt(3));
    const double B1=0.1851156353447362;
    const double B2=0.6904210483822922;
    const double B3=0.3956894730559419;
    const double C1=0.4783690288121502;

    const double a1w=0.003828270494937162;
    const double a3w=0.009793737512487512;
    const double b1w=0.008211737283191111;
    const double b2w=0.009942814891178103;
    const double b3w=0.009595471336070963;
    const double c1w=0.009694996361663028;

public:
    vector<Point> sphere_points;
    vector<Point> grid_points;
//    void write_out_grid_points();
    void Lebedev_grid_();
    void build_sphere(BasisSet &bs, Lebedev_grid &g);
//    double get_weight();

    vector<Tessera> tesserae;
    void cal_tessera(Lebedev_grid &g, BasisSet bs);
};



//build the tesserae

#endif //LEBEDEV_LEBEDEV_GRID_H
