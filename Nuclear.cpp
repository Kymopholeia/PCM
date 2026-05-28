//
// Created by Domokos Kyra on 17/09/2024.
//

#include "BasicFunctions.h"
#include "Nuclear.h"

Nuclear::Nuclear(BasisSet Test_set) {
    cout<<"Nuclear started"<<endl;
    auto start = std::chrono::high_resolution_clock::now();
    Eigen::Index I = Test_set.matrix_size;
    integral_matrix.resize(I, I);

    BasisFunction basisC, basisD;
    Atom Atom_A;
//    std::ofstream outputFile("output.txt");
    vector<double> rests;
//    cout<<"Atoms vector size:  "<<Test_set.atoms.size()<<endl;
    for (int a = 0; a < Test_set.atoms.size(); a++){
        for (int c = 0; c < Test_set.basis.size(); c++){
            for(int d = 0; d < Test_set.basis.size(); d++){
            cal_T.clear();
            Atom_A = Test_set.atoms[a];
            basisC = Test_set.basis[c];
            basisD = Test_set.basis[d];
            set_type(3);

            StarterPack(Atom_A.Element, Atom_A.Centrum.x, Atom_A.Centrum.y, Atom_A.Centrum.z, 0, 0, 1,0.0,0.0, 0.0, 0, 0.0, basisC.element, basisC.Centrum.x, basisC.Centrum.y, basisC.Centrum.z, basisC.type,basisC.exponent, basisD.element, basisD.Centrum.x, basisD.Centrum.y, basisD.Centrum.z, basisD.type,basisD.exponent);
            int Mm = 0;
            made_matrixe.clear();
            rests.clear();
            Eigen::MatrixXd Result_matrix(return_type_num(basisC.type), return_type_num(basisD.type));
//            cout<<"Res_mx size:  "<<Result_matrix.rows()<<"x"<<Result_matrix.cols()<<endl;

            Make_matrixes(0, 0, basisC.type, basisD.type, Mm);

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
//                cout<<"rest size:"<<rests.size()<<endl;
            int t = 0;
            for (int r = 0; r < Result_matrix.cols(); r++) {
                for (int co = 0; co < Result_matrix.rows(); co++) {
//                    cout<<"here"<<endl;
                    Result_matrix(co, r) = rests[t];
                    t += 1;

                }
            }

//here
            for (int row = 0; row < Result_matrix.rows(); row++) {
                for (int col = 0; col < Result_matrix.cols(); col++) {
                    integral_matrix(basisC.FI + row, basisD.FI + col) += Result_matrix(row, col);
                }
            }

                auto now = std::chrono::steady_clock::now();
                auto elapsed = duration_cast<std::chrono::seconds>(now - start).count();
                std::cout << "\rElapsed time: " << elapsed << " seconds" << std::flush;
            }
        }
    }

    cout<<endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " << duration.count() << " seconds" << endl;
    cout<<"Nuclear done."<<endl;
    cout<<endl;
//    cout<<"here is the NUCLEAR matrix: "<<endl<<integral_matrix<<endl;

}

Nuclear::Nuclear(BasisSet bs, Lebedev_grid g) {
    cout<<"Nuclear started"<<endl;
    auto start = std::chrono::high_resolution_clock::now();
    Eigen::Index I = bs.matrix_size;
    integral_matrix.resize(I, I);

    BasisFunction basisC, basisD;
    Atom Atom_A;
//    std::ofstream outputFile("output.txt");
    vector<double> rests;
//    cout<<"Atoms vector size:  "<<Test_set.atoms.size()<<endl;
    for (int a = 0; a < bs.atoms.size(); a++){
        for (int c = 0; c < bs.basis.size(); c++){
            for(int d = 0; d < bs.basis.size(); d++){
                cal_T.clear();
                Atom_A = bs.atoms[a];
                basisC = bs.basis[c];
                basisD = bs.basis[d];
                set_type(4);

                StarterPack(Atom_A.Element, Atom_A.Centrum.x, Atom_A.Centrum.y, Atom_A.Centrum.z, 0, 0, 1,0.0,0.0, 0.0, 0, 0.0, basisC.element, basisC.Centrum.x, basisC.Centrum.y, basisC.Centrum.z, basisC.type,basisC.exponent, basisD.element, basisD.Centrum.x, basisD.Centrum.y, basisD.Centrum.z, basisD.type,basisD.exponent);
                int Mm = 0;
                made_matrixe.clear();
                rests.clear();
                Eigen::MatrixXd Result_matrix(return_type_num(basisC.type), return_type_num(basisD.type));
//            cout<<"Res_mx size:  "<<Result_matrix.rows()<<"x"<<Result_matrix.cols()<<endl;

                Make_matrixes(0, 0, basisC.type, basisD.type, Mm);

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
//                cout<<"rest size:"<<rests.size()<<endl;
                int t = 0;
                for (int r = 0; r < Result_matrix.cols(); r++) {
                    for (int co = 0; co < Result_matrix.rows(); co++) {
//                    cout<<"here"<<endl;
                        Result_matrix(co, r) = rests[t];
                        t += 1;

                    }
                }

//here
                for (int row = 0; row < Result_matrix.rows(); row++) {
                    for (int col = 0; col < Result_matrix.cols(); col++) {
                        integral_matrix(basisC.FI + row, basisD.FI + col) += Result_matrix(row, col);
                    }
                }

                auto now = std::chrono::steady_clock::now();
                auto elapsed = duration_cast<std::chrono::seconds>(now - start).count();
                std::cout << "\rElapsed time: " << elapsed << " seconds" << std::flush;
            }
        }
    }

    cout<<endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " << duration.count() << " seconds" << endl;
    cout<<"Nuclear done."<<endl;
    cout<<endl;
//    cout<<"here is the NUCLEAR matrix: "<<endl<<integral_matrix<<endl;

}

double Nuclear::calculate_PrimInt(double Al, double Be, double Ga, double De, int mm, double Tt) {

}
