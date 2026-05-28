//
// Created by Domokos Kyra on 17/09/2024.
//

#include "BasicFunctions.h"
const double e =2.7182818284590452353602874713526624977;
const double pi = 3.14159265;

double Pitagoras(double x, double y, double z)
{
    double distance;
    distance = sqrt(x*x + y*y+ z*z);
    return distance;
}

double distance_between_points(point aA, point bB){
    point new_vector;
    new_vector.x=(bB.x-aA.x);
    new_vector.y=(bB.y-aA.y);
    new_vector.z=(bB.z-aA.z);
    double distance =Pitagoras(new_vector.x, new_vector.y, new_vector.z);
    return distance;
}

double fv_in_point(double u, double T, int m){
    double value;
    double pow_1=2*m;
    value=pow(u,pow_1);
    double pow_2=-T*pow(u,2);
    value=value*pow(e,pow_2);
    return value;}

double normal_simpson(double from_s, double to_s, double T, int m){
    double a=from_s;
    double b=to_s;
    double ab_2=(from_s+to_s)/2;
    double fv_a=fv_in_point(from_s, T,m);
    double fv_ab_2=fv_in_point(ab_2, T,m);
    double fv_b=fv_in_point(to_s, T,m);
    double c_0=(b-a)/6;
    double c_1=(4*(b-a))/6;
    double value=(c_0*fv_a)+(c_1*fv_ab_2)+(c_0*fv_b);
    return value;
}

double better_simpson_with_step(double step, double T, int m){
    double value=0;
    for(double i=0;i<(1-(1/step));i=i+(1/step)){
        value+=normal_simpson(i,i+(1/step),T,m);
    }
    return value;
}

point multiply_vector(point Al, double alpha_l) {
    point newA;
    newA.x=Al.x*alpha_l;
    newA.y=Al.y*alpha_l;
    newA.z=Al.z*alpha_l;
    return newA;
}

point adding_of_vectors_and_multiply(point Aa, point Bb, double alpha_a, double beta_b){
    point newA,newB,new_A_plus_B;
    newA=multiply_vector(Aa,alpha_a);
    newB=multiply_vector(Bb,beta_b);

    new_A_plus_B.x=newA.x+newB.x;
    new_A_plus_B.y=newA.y+newB.y;
    new_A_plus_B.z=newA.z+newB.z;

    return new_A_plus_B;
}

double calculate_R(double sigma_p, double sigma_q,double alpha, double beta, double gamma, double delta, int type, point A, point B, point C, point D){
    double eta=gamma+delta;
    double EtaRecip = 1.0 / (eta + eta);
    point PP=adding_of_vectors_and_multiply(A,B,alpha,beta);
    PP=multiply_vector(PP,sigma_p);
    point QQ=adding_of_vectors_and_multiply(C,D,gamma,delta);
    QQ=multiply_vector(QQ,sigma_q);
    if(type==3){
        PP = A;
        QQ = (adding_of_vectors_and_multiply(C,D,gamma,delta));
        QQ= multiply_vector(QQ,2*EtaRecip);
    }
    double RR=pow(distance_between_points(PP,QQ),2);
    return RR;
}

point Cal_R(double alpha, double beta, double gamma, double delta,double sigma_p, double sigma_q, int type, point A, point B, point C, point D)
{
    if(type==3) {
        double eta = gamma + delta;
        sigma_p = 0;
        sigma_q = 1.0 / (eta);
    }
    point P=adding_of_vectors_and_multiply(A,B,alpha,beta);
    P=multiply_vector(P,sigma_p);
    if(type==3){P=A;}
    point Q=adding_of_vectors_and_multiply(C,D,gamma,delta);
    Q=multiply_vector(Q,sigma_q);

    point new_vector;
    new_vector.x=(P.x-Q.x);
    new_vector.y=(P.y-Q.y);
    new_vector.z=(P.z-Q.z);

    return new_vector;
}

double calculate_U(double sigma_p, double sigma_q, double alpha, double beta,double gamma, double delta, point A, point B, point C, point D){
    int D_a=1;
    int D_b=1;
    int D_c=1;
    int D_d=1;
    double pow_e_p=-alpha*beta*sigma_p*pow(distance_between_points(A,B),2);
    double U_p=sqrt(pow((pi*sigma_p),3))*D_a*D_b*pow(e,pow_e_p);
    double pow_e_q=-gamma*delta*sigma_q*pow(distance_between_points(C,D),2);
    double U_q=sqrt(pow((pi*sigma_q),3))*D_c*D_d*pow(e,pow_e_q);
    double U=U_p*U_q;
    return U;
}
double calculate_v(double sigma_p, double sigma_q){
    double v=1/(sigma_q+sigma_p);
    return v;
}
double calculate_T(double Rq,double sigma_p, double sigma_q){
    double v=calculate_v(sigma_p, sigma_q);
    double T=v*Rq;
    return T;
}
int what_row(int col, Matrix M){

    for(int i=0;i<3;i++){
        if(M.get_data(i,col)>0)
            return i;
    }
}

int return_type_num(int type){
    switch (type) {
        case 0:
            return 1;
            break;
        case 1:
            return 3;
            break;
        case 2:
            return 6;
            break;
        case 3:
            return 10;
            break;
    }
    return -1;
}