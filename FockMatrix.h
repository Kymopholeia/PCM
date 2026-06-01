//
// Created by Domokos Kyra on 11/11/2024.
//

#ifndef ONLAB_NEWHEADER_FOCKMATRIX_H
#define ONLAB_NEWHEADER_FOCKMATRIX_H

#include "ClassesStructs.h"
#include "BasicFunctions.h"
#include "Eigen/dense"
#include "BasisSets.h"

class Fock_mx{
public:
    Eigen::MatrixXd F;
    Eigen::MatrixXd H_core;
    int basis_size_;
    int mx_size;
    double CoulombE,ExchangeE;
//    double energy;

    Eigen::MatrixXd calculateCoreHamiltonian(Eigen::MatrixXd kin, Eigen::MatrixXd nuc);
    Eigen::MatrixXd calculateFockMatrix(Eigen::MatrixXd p, vector<El_rep_indexing> eri);
    Eigen::MatrixXd calculateFockMatrix(Eigen::MatrixXd p, vector<El_rep_indexing> eri, vector<double> V_PCM);
    double calculateEnergy(Eigen::MatrixXd p,vector<El_rep_indexing> eri);
};

#endif //ONLAB_NEWHEADER_FOCKMATRIX_H
