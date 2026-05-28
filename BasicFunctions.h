//
// Created by Domokos Kyra on 17/09/2024.
//

#ifndef ONLAB_NEWHEADER_BASICFUNCTIONS_H
#define ONLAB_NEWHEADER_BASICFUNCTIONS_H

#include "ClassesStructs.h"
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;


double Pitagoras(double x, double y, double z);
double distance_between_points(point aA, point bB);
double fv_in_point(double u, double T, int m);
double normal_simpson(double from_s, double to_s, double T, int m);
double better_simpson_with_step(double step, double T, int m);
point multiply_vector(point Al, double alpha_l);
point adding_of_vectors_and_multiply(point Aa, point Bb, double alpha_a, double beta_b);
double calculate_R(double sigma_p, double sigma_q,double alpha, double beta, double gamma, double delta, int type, point A, point B, point C, point D);
point Cal_R(double alpha, double beta, double gamma, double delta,double sigma_p, double sigma_q, int type, point A, point B, point C, point D);
double calculate_U(double sigma_p, double sigma_q, double alpha, double beta,double gamma, double delta, point A, point B, point C, point D);
double calculate_v(double sigma_p, double sigma_q);
double calculate_T(double Rq,double sigma_p, double sigma_q);
int what_row(int col, Matrix M);
int return_type_num(int type);

#endif //ONLAB_NEWHEADER_BASICFUNCTIONS_H
