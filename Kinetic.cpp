//
// Created by Domokos Kyra on 17/09/2024.
//


#include "BasicFunctions.h"
#include "Kinetic.h"

//Kinetic::Kinetic(double Ax, double Ay, double Az, int asumA, double al, double Cx, double Cy, double Cz, int asumC,double ga) {
//    set_type(2);
//    StarterPack(1,Ax,Ay,Az,asumA,al,1,0.0,0.0,0.0,0,0.0,1,Cx,Cy,Cz,asumC,ga,1,0.0,0.0,0.0,0,0.0);
//    int Mm = 0;
//
//    Make_matrixes(angsumA, angsumB, angsumC, angsumD, Mm);
//    cout<<endl;
//    cout << "Kinetic:" << endl;
//    cout<<"python3 kinetic_reference.py "<<angsumA<<' '<<alpha<<' '<<A.x<<".0"<<' '<<A.y<<".0"<<' '<<A.z<<".0"<<' '<<angsumC<<' '<<gamma<<' '<<C.x<<".0"<<' '<<C.y<<".0"<<' '<<C.z<<".0"<<' '<<'0'<<endl;
//
//    cout<<endl;
//    cout << "Made matrixes: " << made_matrixe.size() << endl;
//
//    for (size_t m = 0; m < made_matrixe.size(); m++) {
////      cout<<"m. "<<m<<endl;
//
//        made_matrixe[m];
//        equations.clear();
//        needed.clear();
//        Stored_matrixes.clear();
//        make_equations(made_matrixe[m], alpha, beta, gamma, delta, A, B, C, D,type);
//
//        for (size_t i = 0; i < needed.size(); i++) {
//            make_equations(needed[i], alpha, beta, gamma, delta, A, B, C, D, type);
//        }
//
//
//        for (size_t i = 0; i < equations.size(); i++) {
//            Stored_matrixes.push_back(equations[i].make_store());
//
//        }
//
//        for (size_t i = Stored_matrixes.size(); i > 0; i--) {
//            what_2(Stored_matrixes[i - 1], alpha, beta, gamma, delta, type, A, B, C, D);
//            for (size_t j = 0; j < Stored_matrixes.size(); j++) {
//                Stored_matrixes[j].check_good();
//            }
//        }
//
//        if (Stored_matrixes[0].solved) {
//            cout << Stored_matrixes[0].result_mx << endl;
//        }
//
//    }
//}
//

Kinetic::Kinetic(BasisSet Test_set) {
    cout<<"Kinetic started"<<endl;
    auto start = std::chrono::high_resolution_clock::now();

    Eigen::Index I=Test_set.matrix_size;
    integral_matrix.resize(I,I);

    BasisFunction basisA,basisC;

    vector<double> rests;
    for(int a=0;a<Test_set.basis.size();a++){
        for (int c = 0; c < Test_set.basis.size(); c++) {
            cal_T.clear();
//            cout<<"a: "<<a<<"   "<<"c:  "<<c<<endl;

            basisA=Test_set.basis[a];
            basisC=Test_set.basis[c];

            set_type(2);
            StarterPack(basisA.element, basisA.Centrum.x, basisA.Centrum.y, basisA.Centrum.z, basisA.type, basisA.exponent, 1, 0.0,0.0,0.0, 0, 0.0, basisC.element, basisC.Centrum.x, basisC.Centrum.y, basisC.Centrum.z, basisC.type,basisC.exponent,1, 0.0, 0.0, 0.0, 0, 0.0);
            int Mm = 0;
            made_matrixe.clear();
            Make_matrixes(basisA.type, 0, basisC.type, 0, Mm);
//            cout << endl;

            rests.clear();
            Eigen::MatrixXd Result_matrix(return_type_num(basisA.type), return_type_num(basisC.type));

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
            for(int row=0;row<Result_matrix.rows();row++){
                for(int col=0;col<Result_matrix.cols();col++){
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
    cout<<"Kinetic done."<<endl;
    cout<<endl;
//    cout<<endl<<endl;
//    cout<<"here is the KINETIC matrix: "<<endl<<endl<<integral_matrix<<endl;

}

double Kinetic::calculate_PrimInt(double Al, double Be, double Ga, double De, int mm, double Tt) {
//    double KmT;
//    double aux=1/(Al+Ga);
//    double theta=(Al*Ga)/(Al+Ga);
//    double theta2m=pow(theta*2,mm+1);
//    KmT=-0.5*aux *sqrt(aux)*exp(-Tt)* sqrt(pow(pi,3))*(theta2m)*(2* Tt-(2 *(mm)+3));
//    return KmT;

    double N1 = cal_norm(Al);
    double N2 = cal_norm(Ga);

    // Az eredeti számításaid
    double aux = 1.0/(Al + Ga);
    double theta = (Al * Ga)/(Al + Ga);
    double theta2m = pow(theta * 2.0, mm + 1);

    // Kinetikus energia integrál
    double KmT = -0.5 * aux * sqrt(aux) * exp(-Tt) * sqrt(pow(pi, 3))
                 * theta2m * (2.0 * Tt - (2.0 * mm + 3));

    // Normalizált eredmény
    return N1 * N2 * KmT;

}
