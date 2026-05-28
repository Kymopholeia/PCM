//
// Created by Domokos Kyra on 17/09/2024.
//

#ifndef ONLAB_NEWHEADER_CLASSESSTRUCTS_H
#define ONLAB_NEWHEADER_CLASSESSTRUCTS_H

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
//#include "Eigen/dense"
//#include "BasisSets.h"

using namespace std;
//extern const double e =2.7182818284590452353602874713526624977;
//extern const double pi = 3.14159265;


class point{
public:
    double x, y, z;
};

class Matrix {
    double MatrixArray[3][4]={0};
public:
    int m=0;
    void change_data (int row, int column, double data);
    void write_out_Matrix();
    double get_data(int row, int col);
    double col_sum(int col);
    bool same_matrix(Matrix M);
};

struct Pair{
    Matrix M;
    bool solved=false;
};

struct storage{
    Matrix M;
    vector<Pair> need_to;
    bool solved;
    double result_mx;
public:
    bool good();
    void check_good();
    void ki_store();
};

class equation {
public:
    point PQ, AD;
    Matrix f_0, f_1, a_0, a_1, b_0, b_1, c_1, d_1, orig;
    double xi_eta, aa, bb, cc, dd;
    point RR;

    void ki_equ_matrix();
    storage make_store();

};

class El_rep_indexing{
public:
    int Fia,Fib,Fic,Fid;
    int shell_size_A, shell_size_B, shell_size_C;
    double value;
    void write_out_el_rep_vector();
};

#endif //ONLAB_NEWHEADER_CLASSESSTRUCTS_H
