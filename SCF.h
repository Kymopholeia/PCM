//
// Created by Domokos Kyra on 11/11/2024.
//

#ifndef ONLAB_NEWHEADER_SCF_H
#define ONLAB_NEWHEADER_SCF_H

#include "ClassesStructs.h"
#include "BasicFunctions.h"
#include "Eigen/dense"
#include "BasisSets.h"
//#include "Coulomb_Exchange.h"
#include "DensityMatrix.h"
#include "FockMatrix.h"
#include "Overlap.h"
#include "Kinetic.h"
#include "Nuclear.h"
#include "ObaraSaika.h"


class SCF{
public:
    int max_iter;
    double conv_tresh;
    double total_energy;
    double nuc_rep_en;
    Eigen::MatrixXd S,T,V;  //overlap,kinetic,nuclear matrix
    vector<El_rep_indexing> ERI;
    int mx_size; //electron number
    int basis_size;

    SCF(double conv_threshold, int max_iter);
    void runSCF(BasisSet& Water,Eigen::MatrixXd& P);
    double calculateNuclearRepulsionEnergy(BasisSet& Water) const;
};
#endif //ONLAB_NEWHEADER_SCF_H
