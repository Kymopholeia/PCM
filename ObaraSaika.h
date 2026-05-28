//
// Created by Domokos Kyra on 17/09/2024.
//

#ifndef ONLAB_NEWHEADER_OBARASAIKA_H
#define ONLAB_NEWHEADER_OBARASAIKA_H

#include "Integrals.h"
class ObaraSaika : public integral
{

public:
    vector<El_rep_indexing> ERI;
    ObaraSaika(BasisSet Test_set);
    vector<vector<vector<vector<double>>>> Eri2;
    virtual double calculate_PrimInt(double Al, double Be, double Ga, double De, int mm, double Tt);

};

#endif //ONLAB_NEWHEADER_OBARASAIKA_H
