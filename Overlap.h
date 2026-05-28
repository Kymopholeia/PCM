//
// Created by Domokos Kyra on 17/09/2024.
//

#ifndef ONLAB_NEWHEADER_OVERLAP_H
#define ONLAB_NEWHEADER_OVERLAP_H

#include "Integrals.h"

class Overlap : public integral
{
public:
    Overlap(BasisSet Test_set);
    Overlap(BasisFunction A, BasisFunction C);
    virtual double calculate_PrimInt(double Al, double Be, double Ga, double De, int mm, double Tt);
};

#endif //ONLAB_NEWHEADER_OVERLAP_H
