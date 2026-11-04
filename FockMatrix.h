//
// Created by Domokos Kyra on 11/11/2024.
//

#ifndef ONLAB_NEWHEADER_FOCKMATRIX_H
#define ONLAB_NEWHEADER_FOCKMATRIX_H

#include "ClassesStructs.h"
#include "BasicFunctions.h"
#include "Eigen/dense"
#include "BasisSets.h"
#include "Pcm.h"

class Fock_mx{
public:
    Eigen::MatrixXd F;
    Eigen::MatrixXd H_core;
    Eigen::MatrixXd V_pcm;
    int basis_size_;
    int mx_size;


    Eigen::MatrixXd calculateCoreHamiltonian(Eigen::MatrixXd& kin, Eigen::MatrixXd& nuc);
    Eigen::MatrixXd calculateFockMatrix(Eigen::MatrixXd& p, vector<vector<vector<vector<double>>>>& eri);
    Eigen::MatrixXd calculateFockMatrix(Eigen::MatrixXd& p, vector<vector<vector<vector<double>>>>& eri, Eigen::VectorXd& q, vector<Eigen::MatrixXd>& Nuclear_pcm_matrix);
    double calculateEnergy(Eigen::MatrixXd& p,vector<vector<vector<vector<double>>>>& eri);
    double calculateEnergy(Eigen::MatrixXd& p,vector<vector<vector<vector<double>>>>& eri, double& E_solve);
};

#endif //ONLAB_NEWHEADER_FOCKMATRIX_H
