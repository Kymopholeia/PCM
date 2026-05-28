//
// Created by Domokos Kyra on 17/09/2024.
//

#include "BasicFunctions.h"
#include "ObaraSaika.h"

ObaraSaika::ObaraSaika(BasisSet Test_set) {
    cout<<"Electron Repulsion Started"<<endl;
    auto start = std::chrono::high_resolution_clock::now();
    vector<El_rep_indexing> El_rep_vec;
    BasisFunction basisA, basisB,basisC, basisD;
    Eri2.resize(Test_set.matrix_size, vector<vector<vector<double>>>(Test_set.matrix_size, vector<vector<double>>(Test_set.matrix_size, vector<double>(Test_set.matrix_size, 0.0))));
    cout<<"eri2 size:  "<<Eri2.size()<<endl;

//    std::ofstream outputFile("output.txt");
    vector<double> rests;
//    Ebben a vektorban vannak az egyes shellquartettek eredmény tömbjei
    for (int a = 0; a < Test_set.basis.size(); a++){
        for(int b = 0; b < Test_set.basis.size(); b++){
            for (int c = 0; c < Test_set.basis.size(); c++) {
                for (int d = 0; d < Test_set.basis.size(); d++) {
                    cal_T.clear();
                    basisA = Test_set.basis[a];
                    basisB = Test_set.basis[b];
                    basisC = Test_set.basis[c];
                    basisD = Test_set.basis[d];


//                cout<<"basis function A:  ";basisA.write_out_bF();cout<<endl;
//                cout<<"basis function B:  ";basisB.write_out_bF();cout<<endl;
//                cout<<"basis function C:  ";basisC.write_out_bF();cout<<endl;
//                cout<<"basis function D:  ";basisD.write_out_bF();

//                cout<<"basis function A fine index:  "<<basisA.FI<<"   and type:   "<<basisA.type<<endl;
//                cout<<"basis function B fine index:  "<<basisB.FI<<"   and type:   "<<basisB.type<<endl;
//                cout<<"basis function C fine index:  "<<basisC.FI<<"   and type:   "<<basisC.type<<endl;
//                cout<<"basis function D fine index:  "<<basisD.FI<<"   and type:   "<<basisD.type<<endl;


                    set_type(0);
                    StarterPack(basisA.element, basisA.Centrum.x, basisA.Centrum.y, basisA.Centrum.z, basisA.type, basisA.exponent,
                                basisB.element, basisB.Centrum.x, basisB.Centrum.y, basisB.Centrum.z, basisB.type, basisB.exponent,
                                basisC.element, basisC.Centrum.x, basisC.Centrum.y, basisC.Centrum.z,basisC.type, basisC.exponent,
                                basisD.element, basisD.Centrum.x, basisD.Centrum.y,basisD.Centrum.z, basisD.type, basisD.exponent);
                    int Mm = 0;
                    made_matrixe.clear();
                    rests.clear();

//
//                  IDE kell egy vektor ahova el tudom menteni az eredményeket,
//                  ez mehet azonnal mint integral matrix helyére
//                  mert itt nem kell belehelyezni
//                  indexelni kell majd jól!!!!
//                  esetleg rests törölni?
//
//

                    Make_matrixes(basisA.type, basisB.type, basisC.type, basisD.type, Mm);
//                    cout<<"made matrixes: "<<made_matrixe.size()<<endl;

                    for (size_t m = 0; m < made_matrixe.size(); m++) {
                        made_matrixe[m];
                        equations.clear();
                        needed.clear();
                        Stored_matrixes.clear();
                        make_equations(made_matrixe[m], alpha, beta, gamma, delta, A, B, C, D, type);

                        for (size_t i = 0; i < needed.size(); i++) {
                            make_equations(needed[i], alpha, beta, gamma, delta, A, B, C, D, type);
                        }


                        for (size_t i = 0; i < equations.size(); i++) {
                            Stored_matrixes.push_back(equations[i].make_store());

                        }

                        for (size_t i = Stored_matrixes.size(); i > 0; i--) {
                            what_2(Stored_matrixes[i - 1], alpha, beta, gamma, delta, type, A, B, C, D);
                            for (size_t j = 0; j < Stored_matrixes.size(); j++) {
                                Stored_matrixes[j].check_good();
                            }
                        }

                        if (Stored_matrixes[0].solved) {
//                    cout<<"Stored matrix result mx: "<<endl;
//                    cout << Stored_matrixes[0].result_mx << endl;
//                    cout<<"-----------------------------------~~~---------------------------------"<<endl;
                            rests.push_back(Stored_matrixes[0].result_mx);

                        }

                    }



                         El_rep_indexing tmp;
                    int t=0;
                        for(int i=0;i< return_type_dec(basisA.type);i++){
                            for(int j=0;j< return_type_dec(basisB.type);j++){
                                for(int k=0;k< return_type_dec(basisC.type);k++){
                                    for(int l=0;l< return_type_dec(basisD.type);l++){
//
                                        tmp.Fia = basisA.FI + i;
                                        tmp.Fib = basisB.FI + j;
                                        tmp.Fic = basisC.FI + k;
                                        tmp.Fid = basisD.FI + l;

                                        tmp.value=rests[t];
                                        ERI.push_back(tmp);
                                        Eri2[tmp.Fia][tmp.Fib][tmp.Fic][tmp.Fid]=tmp.value;
                                        t++;

                                    }
                                }
                            }
                        }

//cout<<"-----------------------------------~~~---------------------------------"<<endl;
                    auto now = std::chrono::steady_clock::now();
                    auto elapsed = duration_cast<std::chrono::seconds>(now - start).count();
                    std::cout << "\rElapsed time: " << elapsed << " seconds" << std::flush;
                }
            }
        }
    }

    cout<<endl;
    cout<<"Eri2. size: "<<Eri2.size()<<endl;

    cout<<endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " << duration.count() << " seconds" << endl;
    cout<<"Electron repulsion done."<<endl;
    cout<<endl;
}

double ObaraSaika::calculate_PrimInt(double Al, double Be, double Ga, double De, int mm, double Tt) {
    mm=0;
    double Omt;
    double aux=1/(Al+Ga);
    double theta=(Al*Ga)/(Al+Ga);
    double theta2m=pow(theta*2,mm);
    Omt=theta2m*aux *sqrt(aux)*exp(-Tt)* sqrt(pow(pi,3));
    return Omt;
}

