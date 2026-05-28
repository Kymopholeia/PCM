//
// Created by Domokos Kyra on 10/11/2024.
//

#ifndef ONLAB_NEWHEADER_DENSITYMATRIX_H
#define ONLAB_NEWHEADER_DENSITYMATRIX_H


#include "ClassesStructs.h"
#include "BasicFunctions.h"
#include "Eigen/dense"
#include "BasisSets.h"


class Density_Matrix{
public:


    Eigen::MatrixXd p;
    int basis_size;
    int mx_size;
    int n_elect;

    Eigen::MatrixXd constructFromCoeff(Eigen::MatrixXd& C,Eigen::MatrixXd S);


};

#endif //ONLAB_NEWHEADER_DENSITYMATRIX_H
