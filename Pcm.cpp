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

        Nuclear N(bs,g);
        Nuclear_pcm_matrix=N.Nuclear_PCM;
        cout<<"Here we are inside the build pot vect: Nuclear pcm.size: "<<Nuclear_pcm_matrix.size()<<endl;

        for(int i=0;i<Nuclear_pcm_matrix.size();i++){
            v_elec.push_back(0);
            for(int j=0;j<P.cols();j++){
                for(int k=0;k<P.cols();k++){
                    v_elec[i]+=P(j,k)*Nuclear_pcm_matrix[i](j,k);
                }
            }
             v_elec[i]*=-1;
        }

        for(int i=0;i<v_fixed.size();i++){

            V_electric.push_back(v_fixed[i]+v_elec[i]);
        }
//        cout<<"V_electric: "<<endl;
//        for(int i=0;i<V_electric.size();i++){
//            cout<<V_electric[i]<<endl;
//        }
    Eigen::VectorXd V_eigen = Eigen::Map<Eigen::VectorXd>(V_electric.data(), V_electric.size());
    Eigen::VectorXd b = -0.987 * V_eigen;

    Eigen::LLT<Eigen::MatrixXd> llt(S_matrix);
    q = llt.solve(b);
//    cout<<q<<endl;
    cout<<"q. size: "<<q.size()<<endl;
    double q_sum = q.sum();
    double q_min = q.minCoeff();
    double q_max = q.maxCoeff();

// Need v as Eigen vector for dot product:
    Eigen::VectorXd v_eigen = Eigen::Map<Eigen::VectorXd>(V_electric.data(), V_electric.size());
    E_solv = 0.5 * q.dot(v_eigen);

//    std::cout << "Sum of charges (Gauss's law): " << q_sum << std::endl;
//    std::cout << "Min charge: " << q_min << std::endl;
//    std::cout << "Max charge: " << q_max << std::endl;
//    std::cout << "Solvation energy estimate: " << E_solv << " Hartree" << std::endl;
//    std::cout << "                         = " << E_solv * 627.509 << " kcal/mol" << std::endl;
//        cout<<N.integral_matrix<<endl;


//        Nuclear_pcm_matrix=N.Nuclear_PCM;
//    //    for(int i=0;i<nuclear.Nuclear_PCM.size();i++){
//        for(int j=0;j<Nuclear_pcm_matrix[259].rows();j++){
//            for(int k=0;k<Nuclear_pcm_matrix[259].cols();k++){
//    //                if(j==k){cout<<nuclear.Nuclear_PCM[i](j,k)<<' ';}
//                cout<<Nuclear_pcm_matrix[259](j,k)<<" ";
//            }
//            cout<<endl;
//        }


//    cout << "min of v_nuc: " << *min_element(v_fixed.begin(), v_fixed.end())
//         << " max of v_nuc: " << *max_element(v_fixed.begin(), v_fixed.end()) << endl;
//
//    double min_d = 1e30, max_d = 0;
//    for (auto& t : g.tesserae) {
//        double d_min_to_any_atom = 1e30;
//        for (auto& a : bs.atoms) {
//            double dx = t.r.x - a.Centrum.x;
//            double dy = t.r.y - a.Centrum.y;
//            double dz = t.r.z - a.Centrum.z;
//            double d = std::sqrt(dx*dx + dy*dy + dz*dz);
//            if (d < d_min_to_any_atom) d_min_to_any_atom = d;
//        }
//        if (d_min_to_any_atom < min_d) min_d = d_min_to_any_atom;
//        if (d_min_to_any_atom > max_d) max_d = d_min_to_any_atom;
//    }
//    std::cout << "Closest tessera-to-nucleus distance: " << min_d << "\n";
//    std::cout << "Farthest tessera-to-nucleus distance: " << max_d << "\n";
//
//    auto min_it = std::min_element(v_fixed.begin(), v_fixed.end());
//    int min_idx = std::distance(v_fixed.begin(), min_it);
//    int i=min_idx;
//    cout<<"Element where the min index values is: "<<bs.atoms[g.tesserae[i].r.no_of_atom].Element<<endl;
//    cout<<"atom point position: "<<bs.atoms[g.tesserae[i].r.no_of_atom].Centrum.x<<','<<bs.atoms[g.tesserae[i].r.no_of_atom].Centrum.y<<','<<bs.atoms[g.tesserae[i].r.no_of_atom].Centrum.z<<endl;
//    cout<<"tessera point position: "<<g.tesserae[i].r.x<<','<<g.tesserae[i].r.y<<','<<g.tesserae[i].r.z<<endl;
//    cout<<"absolute of distance of the points "<<fabs(dist_bet_po(g.tesserae[i].r,bs.atoms[g.tesserae[i].r.no_of_atom].centrum))<<endl;
//    cout<<"min index: "<<min_idx<<endl;
//    cout<<"v_fixed at min index: "<<v_fixed[min_idx]<<endl;


}


