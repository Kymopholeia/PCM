//
// Created by Domokos Kyra on 10/11/2024.
//

#include "DensityMatrix.h"

Eigen::MatrixXd Density_Matrix::constructFromCoeff( Eigen::MatrixXd& C,Eigen::MatrixXd S){
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(mx_size, mx_size);
    int n_occ=n_elect/2;

    for(int mu = 0; mu < mx_size; mu++) {
        for(int nu = 0; nu < mx_size; nu++) {
            double sum = 0.0;
            for(int i = 0; i < n_occ; i++) {
                sum += C(mu,i) * C(nu,i);
            }
            P(mu,nu) = 2.0 * sum;
        }
    }
    return P;
}
