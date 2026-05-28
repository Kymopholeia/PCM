//
// Created by Domokos Kyra on 17/09/2024.
//

#ifndef ONLAB_NEWHEADER_KINETIC_H
#define ONLAB_NEWHEADER_KINETIC_H

#include "Integrals.h"
class Kinetic : public integral
{

public:
    Kinetic (BasisSet Test_set);
//    Kinetic(double Ax,double Ay, double Az, int asumA, double al,double Cx,double Cy, double Cz, int asumC, double ga);
    virtual double calculate_PrimInt(double Al, double Be, double Ga, double De, int mm, double Tt);
};

#endif //ONLAB_NEWHEADER_KINETIC_H
