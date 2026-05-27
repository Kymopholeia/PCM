//
// Created by Domokos Kyra on 2026. 05. 27..
//

#ifndef LEBEDEV_POINT_CLASS_H
#define LEBEDEV_POINT_CLASS_H

class Point{
public:
    double x;
    double y;
    double z;
    double weight;
    double Area;
    int no_of_atom;
    void write_point();
    void cal_area();
};

#endif //LEBEDEV_POINT_CLASS_H
