//
// Created by Domokos Kyra on 17/09/2024.
//

#ifndef ONLAB_NEWHEADER_NUCLEAR_H
#define ONLAB_NEWHEADER_NUCLEAR_H

#include "Integrals.h"
#include "Lebedev_grid.h"


class Nuclear : public integral
{

public:
    Nuclear(BasisSet Test_set);
    Nuclear(BasisSet bs, Lebedev_grid g );
    virtual double calculate_PrimInt(double Al, double Be, double Ga, double De, int mm, double Tt);
};

#endif //ONLAB_NEWHEADER_NUCLEAR_H
