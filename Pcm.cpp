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

void Pcm::cal_smatrix(Lebedev_grid& g) {
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

void Pcm::build_pot_vector(BasisSet& bs,Lebedev_grid& g, Eigen::MatrixXd& P) {
    v_elec.clear();
    V_electric.clear();
    v_fixed.clear();
    pot_vect.clear();
    double fix_element;
        if(v_fixed.size()==0){
            for(int i=0;i<g.tesserae.size();i++){
                fix_element=0;
                for(int j=0;j<bs.atoms.size();j++){
                    fix_element += bs.atoms[j].Element /
                                  fabs(dist_bet_po(g.tesserae[i].r, bs.atoms[j].centrum));

                }
                v_fixed.push_back(fix_element);
            }
        }


        if(Nuclear_pcm_matrix.empty()){
            Nuclear N(bs,g);
            Nuclear_pcm_matrix=N.Nuclear_PCM;
        }
        for(int i=0;i<Nuclear_pcm_matrix.size();i++){
            v_elec.push_back(0);
            for(int j=0;j<P.cols();j++){
                for(int k=0;k<P.cols();k++){
                    v_elec[i]+=(P(j,k)*Nuclear_pcm_matrix[i](j,k));
                }
            }
             v_elec[i]*=-1;
        }

        for(int i=0;i<v_fixed.size();i++){

            V_electric.push_back(v_fixed[i]+v_elec[i]);
        }

    Eigen::VectorXd V_eigen = Eigen::Map<Eigen::VectorXd>(V_electric.data(), V_electric.size());
    Eigen::VectorXd b = -0.987243 * V_eigen;

    Eigen::LLT<Eigen::MatrixXd> llt(S_matrix);
    q = llt.solve(b);
    cout<<"q. size: "<<q.size()<<endl;

// Need v as Eigen vector for dot product:
    Eigen::VectorXd v_eigen = Eigen::Map<Eigen::VectorXd>(V_electric.data(), V_electric.size());
    E_solv = 0.5 * q.dot(v_eigen);
}


