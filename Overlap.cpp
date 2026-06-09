//
// Created by Domokos Kyra on 17/09/2024.
//


#include "BasicFunctions.h"
#include "Overlap.h"

Overlap::Overlap(BasisSet Test_set) {
    cout<<"Overlap started"<<endl;
    auto start = std::chrono::high_resolution_clock::now();

    Eigen::Index I=Test_set.matrix_size;
    integral_matrix.resize(I,I);

    BasisFunction basisA,basisC;

    vector<double> rests;
    for(int a=0;a<Test_set.basis.size();a++){
        for (int c = 0; c < Test_set.basis.size(); c++) {
            cal_T.clear();
            basisA=Test_set.basis[a];
            basisC=Test_set.basis[c];
//            cout<<endl;
//            basisA.write_out_bF();
//            basisC.write_out_bF();

            set_type(1);
            StarterPack(basisA.element, basisA.Centrum.x, basisA.Centrum.y, basisA.Centrum.z, basisA.type, basisA.exponent, 1, 0.0,0.0,0.0, 0, 0.0, basisC.element, basisC.Centrum.x, basisC.Centrum.y, basisC.Centrum.z, basisC.type,basisC.exponent,1, 0.0, 0.0, 0.0, 0, 0.0);
            int Mm = 0;
            made_matrixe.clear();
            Make_matrixes(basisA.type, 0, basisC.type, 0, Mm);
            rests.clear();
            Eigen::MatrixXd Result_matrix(return_type_num(basisA.type), return_type_num(basisC.type));
//            cout<<"result mx size:  "<<Result_matrix.rows()<<"x"<<Result_matrix.cols()<<endl;

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
//                    cout << Stored_matrixes[0].result_mx << endl;
                    rests.push_back(Stored_matrixes[0].result_mx);


            }



            }
            int t=0;
            for(int r=0;r<Result_matrix.cols();r++){
                for(int co=0;co<Result_matrix.rows();co++){
                    Result_matrix(co,r)=rests[t];
                    t+=1;
                }
            }

//here
//            cout<<"INTEGRAL MX SIZE:  "<<integral_matrix.size()<<endl;
            for(int row=0;row<Result_matrix.rows();row++){
                for(int col=0;col<Result_matrix.cols();col++){
//                    cout<<"row: "<<row<<endl;
//                    cout<<"col: "<<col<<endl;
//                    cout<<"basisA.FI+row: "<<basisA.FI+row<<endl;
//                    cout<<"basisC.FI+col: "<<basisC.FI+col<<endl;
                    integral_matrix(basisA.FI+row,basisC.FI+col)=Result_matrix(row,col);
                }
            }

            auto now = std::chrono::steady_clock::now();
            auto elapsed = duration_cast<std::chrono::seconds>(now - start).count();
            std::cout << "\rElapsed time: " << elapsed << " seconds" << std::flush;
        }
    }


    cout<<endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " << duration.count() << " seconds" << endl;
    cout<<"Overlap done."<<endl;
    cout<<endl;
////integral_matrix[basisA.RI,basisC.RI]=Stored_matrixes[0].result_mx;
    cout<<"here is the OVERLAP matrix: "<<endl<<integral_matrix<<endl;

}

double Overlap::calculate_PrimInt(double Al, double Be, double Ga, double De, int mm, double Tt) {
//    mm=0;
//    double Omt;
//    double aux=1/(Al+Ga);
//    double theta=(Al*Ga)/(Al+Ga);
//    double theta2m=pow(theta*2,mm);
//    Omt=theta2m*aux *sqrt(aux)*exp(-Tt)* sqrt(pow(pi,3));
//    return Omt;
    double N1 = cal_norm(Al);
    double N2 = cal_norm(Ga);
    mm = 0;
    double aux = 1.0/(Al + Ga);
    double theta = (Al * Ga)/(Al + Ga);
    double theta2m = pow(theta * 2.0, mm);

    // Átfedési integrál
    double overlap = theta2m * aux * sqrt(aux) * exp(-Tt) * sqrt(pow(M_PI, 3));
    return N1 * N2 * overlap;

}