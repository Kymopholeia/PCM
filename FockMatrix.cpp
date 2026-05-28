//
// Created by Domokos Kyra on 11/11/2024.
//

#include "FockMatrix.h"

Eigen::MatrixXd Fock_mx::calculateCoreHamiltonian(Eigen::MatrixXd kin, Eigen::MatrixXd nuc) {
    H_core=kin+nuc;
    cout<<"h_core:"<<endl;
    cout<<H_core<<endl;
    return H_core;
}

Eigen::MatrixXd Fock_mx::calculateFockMatrix(Eigen::MatrixXd p,vector<El_rep_indexing> eri) {
    auto start = std::chrono::high_resolution_clock::now();
cout<<"Fock matrix started:"<<endl;
//    cout<<"here"<<endl;
    F=H_core;
    cout<<"eri. size:  "<<eri.size()<<endl;
    for(int mu = 0; mu < mx_size; mu++) {
        for(int nu = 0; nu < mx_size; nu++) {
            double sum=0.0;
            for(int lambda = 0; lambda < mx_size; lambda++) {
                for(int sigma = 0; sigma < mx_size; sigma++) {
                    double val1, val2;

                    for(const auto& e : eri) {

                        if(e.Fia == mu && e.Fib == nu &&
                        e.Fic == lambda && e.Fid == sigma) {
                            val1= e.value;
                        }

                        if(e.Fia == mu && e.Fic == nu &&
                        e.Fib == lambda && e.Fid == sigma) {
                            val2=e.value;
                        }
                    }

//                        val1=eri[mu][nu][lambda][sigma];
//                        val2=eri[mu][lambda][nu][sigma];

                    double P_ls = p(lambda,sigma);
                    sum +=0.5* P_ls * (2.0 * val1 - val2);
                }

                auto now = std::chrono::steady_clock::now();
                auto elapsed = duration_cast<std::chrono::seconds>(now - start).count();
                std::cout << "\rElapsed time: " << elapsed << " seconds" << std::flush;
            }
            F(mu,nu) += sum;
        }
    }



    cout<<endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " << duration.count() << " seconds" << endl;
    cout<<"Fock matrix done."<<endl;
    cout<<endl;

//    cout<<"Here is the FOCK matrix: "<<endl;
//    cout<<F<<endl;
//   double energy= calculateEnergy(p,eri);
//    cout<<"Energy now is: "<<energy<<endl;
    return F;
}



double Fock_mx::calculateEnergy(Eigen::MatrixXd P,vector<El_rep_indexing> eri) {
//  TOTOD basis size helyett n_electron
   Eigen::MatrixXd Co(mx_size,mx_size);
   Co=Eigen::MatrixXd::Zero(mx_size,mx_size);
   Eigen::MatrixXd Ex(mx_size,mx_size);
   Ex=Eigen::MatrixXd::Zero(mx_size,mx_size);

    for(auto e : eri) {
        int mu = e.Fia;
        int nu= e.Fib;
        int lambda = e.Fic;
        int sigma= e.Fid;
        double value = e.value;

        Co(mu,nu) += P(lambda,sigma) * value;

    }
    Co=Co*0.5;

    for(auto e : eri) {
        int mu = e.Fia;
        int nu= e.Fib;
        int lambda = e.Fic;
        int sigma= e.Fid;
        double value = e.value;

        Ex(mu,lambda) += value * P(nu,sigma);

    }
    Ex=Ex*-0.25;

//    p*coulom p*hcore p*exchange
//     1.23      -2.4    -0.614
//-1.7945
//with hidrogen
//      30.68   -77.27    -5..06
//-34
//    cout<<"Total Coulomb energy:  "<<(P*Co).trace()<<endl;
//    cout<<"One electron energy :  "<<(P*H_core).trace()<<endl;
//    cout<<"HF Exchange energy :  "<<(P*Ex).trace()<<endl;
//    cout<<"Energy: "<<(P*Co).trace()+(P*H_core).trace()+(P*Ex).trace()<<endl;

    return (P*Co).trace()+(P*H_core).trace()+(P*Ex).trace();
}