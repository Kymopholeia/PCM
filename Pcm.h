//
// Created by Domokos Kyra on 2026. 05. 28.
//

#ifndef LEBEDEV_PCM_H
#define LEBEDEV_PCM_H

#include "Lebedev_grid.h"
#include "BasisSets.h"
#include "Point_class.h"
#include "Eigen/Dense"

class Pcm{
public:
    int s_mx_size;
    Eigen::MatrixXd S_matrix;
    vector<double> pot_vect;
    vector<double>v_fixed;
    vector<double>v_elec;
    void cal_smatrix(Lebedev_grid g);
    void build_pot_vector(BasisSet bs,Lebedev_grid g);
    double dist_bet_po(Point a, Point b);
    void cal
};


#endif //LEBEDEV_PCM_H
