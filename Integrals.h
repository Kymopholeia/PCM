//
// Created by Domokos Kyra on 17/09/2024.
//

#ifndef ONLAB_NEWHEADER_INTEGRALS_H
#define ONLAB_NEWHEADER_INTEGRALS_H



#include "ClassesStructs.h"
#include "BasicFunctions.h"
#include "Eigen/dense"
#include "BasisSets.h"

class integral{
public:
    const double e =2.7182818284590452353602874713526624977;
    const double pi = 3.14159265;
    int type;
    point A, B, C, D,A_pcm_nuc;
    int angsumA, angsumB, angsumC, angsumD,AZ,BZ,CZ,DZ;
    double alpha, beta, gamma, delta;
    Eigen::MatrixXd integral_matrix;
    Eigen::MatrixXd norm_mx;
    vector<storage> Stored_matrixes;
    vector<Matrix> needed;
    vector<Matrix> made_matrixe;
    vector<equation> equations;
    vector <pair<int,double> > cal_T;
    void StarterPack(int AZn, double Ax,double Ay, double Az, int asumA, double al,int BZn,double Bx,double By, double Bz, int asumB, double be,int CZn,double Cx,double Cy, double Cz, int asumC, double ga,int DZn,double Dx,double Dy, double Dz, int asumD, double de);
    void set_type(int t);
    bool stored_it(Matrix M);
    void Make_matrixes(int a,int b,int c,int d, int m);
    virtual double calculate_PrimInt(double Al, double Be, double Ga, double De, int mm, double Tt)=0;
    double calculation(int m,double alpha, double beta,double gamma, double delta,  double sigma_p, double sigma_q, int type, point A, point B, point C, point D,double QmT);
    void make_equations(Matrix Ma,double alpha, double beta, double gamma, double delta,point A, point B, point C, point D,int type);
    int what_row(int col, Matrix M);
    void finder(Matrix M,double data);
    void download(double &f, Matrix M);
    void what_2(storage s,double alpha,double beta,double gamma,double delta,int type, point A, point B, point C, point D);
    void nom_mx_cal();
    double cal_norm(double alp);



};

#endif //ONLAB_NEWHEADER_INTEGRALS_H
