//
// Created by Domokos Kyra on 2026. 05. 28..
//

#include "Pcm.h"
#include "Nuclear.h"

double Pcm::dist_bet_po(Point a, Point b) {
    double xx=(a.x-b.x)*(a.x-b.x);
    double yy=(a.y-b.y)*(a.y-b.y);
    double zz=(a.z-b.z)*(a.z-b.z);
    return sqrt(xx+yy+zz);
}

void Pcm::cal_smatrix(Lebedev_grid g) {
    S_matrix=Eigen::MatrixXd::Zero(s_mx_size,s_mx_size);
    for(int i=0;i<g.tesserae.size();i++){
        for(int j=0;j<g.tesserae.size();j++){
            if(i!=j){
                S_matrix(i,j)=1/(dist_bet_po(g.tesserae[i].r,g.tesserae[j].r));
            }
            if(i==j){
                S_matrix(i,j)=1.0694*sqrt((4*M_PI)/(g.tesserae[i].area));
            }
        }
    }
}

void Pcm::build_pot_vector(BasisSet bs,Lebedev_grid g) {
    double fix_element;
        if(v_fixed.size()!=0){
            for(int i=0;i<g.tesserae.size();i++){
                fix_element=bs.atoms[g.tesserae[i].r.no_of_atom].Element/
                        fabs(dist_bet_po(g.tesserae[i].r,bs.atoms[g.tesserae[i].r.no_of_atom].centrum));
                v_fixed.push_back(fix_element);
            }
        }
        double vi_elec;



}


