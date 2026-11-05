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
//    const double A1=(1/sqrt(2));
//    const double A2=(1/sqrt(2));
//    const double A3=(1/sqrt(3));
//    const double B1=0.1851156353447362;
//    const double B2=0.6904210483822922;
//    const double B3=0.3956894730559419;
//    const double C1=0.4783690288121502;

    const double A1=8.54591172878* pow(10,-4);
    const double A2=(1/sqrt(2));
    const double A3=3.59911928502*pow(10,-3);
    const double B1=3.65004580768*pow(10,-3);
    const double B2=3.060482260142*pow(10,-3);
    const double B3=3.57672966173*pow(10,-3);
    const double B4=3.44978842429*pow(10,-3);
    const double B5=3.10895312238*pow(10,-3);
    const double B6=2.35210141366*pow(10,-3);
    const double C1=3.60082093222*pow(10,-3);
    const double C2=2.98234496317*pow(10,-3);
    const double D1=3.57154055427*pow(10,-3);
    const double D2=3.39231220501*pow(10,-3);

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
