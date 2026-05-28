//
// Created by Domokos Kyra on 17/09/2024.
//

#include "ClassesStructs.h"
#include <iomanip>

using namespace std;

void Matrix::change_data(int row, int column, double data) {
    MatrixArray[row][column]=data;
}

void Matrix::write_out_Matrix() {
    cout<<"m= "<<m<<endl;
    for(int i=0;i<3;i++){
        for(int j=0;j<4;j++){
            cout<< MatrixArray[i][j]<<" ";
        }
        cout<<endl;
    }
}

double Matrix::get_data(int row, int col) {
    return MatrixArray[row][col];
}

double Matrix::col_sum(int col) {
    double res=0;
    for(int i=0;i<3;i++){
        res+=MatrixArray[i][col];
    }
    return res;
}

bool Matrix::same_matrix(Matrix M) {
    bool same=true;
    for(int i=0;i<3;i++){
        for(int j=0;j<4;j++){
            if(MatrixArray[i][j]!=M.get_data(i,j)||m!=M.m){
                same=false;
            }
        }
    }
    return same;
}

bool storage::good() {
    int a=need_to.size();
    for(size_t i=0;i<need_to.size();i++){
        if(!need_to[i].solved) {return false;}
    }
    return true;
}

void storage::check_good() {
    if(good()) { solved = true; }
    else{solved=false;}
}

void storage::ki_store() {
    M.write_out_Matrix();
    cout<<endl;
    for(size_t i=0;i<need_to.size();i++){
        cout<<"need matrix: "<<endl;
        need_to[i].M.write_out_Matrix();
        cout<<"Solved:  "<<need_to[i].solved<<endl;
        cout<<endl;
    }
}

void equation::ki_equ_matrix() {
    orig.write_out_Matrix();
    cout<<endl;
    f_0.write_out_Matrix();
    cout<<endl;
    f_1.write_out_Matrix();
    cout<<endl;
    a_0.write_out_Matrix();
    cout<<endl;
    a_1.write_out_Matrix();
    cout<<endl;
    b_0.write_out_Matrix();
    cout<<endl;
    b_1.write_out_Matrix();
    cout<<endl;
    c_1.write_out_Matrix();
    cout<<endl;
    d_1.write_out_Matrix();
    cout<<endl;
}

storage equation::make_store() {
    storage s;
    s.M = orig;
    Pair p;

    p.M = f_0;
    p.solved = false;
    s.need_to.push_back(p);

    p.M = f_1;
    p.solved = false;
    s.need_to.push_back(p);
//        cout<<"aa:  "<<aa<<endl;

    int Y = -1;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            if (orig.get_data(i, j) - f_1.get_data(i, j) != 0) {
                Y = j;
            }
        }
    }

    if (Y == 0 || Y == 1) {
        if (aa > 0) {
            p.M = a_0;
//            a_0.write_out_Matrix();
            p.solved = false;
            s.need_to.push_back(p);

            p.M = a_1;
            p.solved = false;
            s.need_to.push_back(p);
        }

        if (bb > 0) {
            p.M = b_0;
            p.solved = false;
            s.need_to.push_back(p);

            p.M = b_1;
            p.solved = false;
            s.need_to.push_back(p);
        }

        if (cc > 0) {
            p.M = c_1;
            p.solved = false;
            s.need_to.push_back(p);
        }

        if (dd > 0) {
            p.M = d_1;
            p.solved = false;
            s.need_to.push_back(p);
        }
    }
    else {
        if (cc > 0) {
            p.M = a_0;
//            a_0.write_out_Matrix();
            p.solved = false;
            s.need_to.push_back(p);

            p.M = a_1;
            p.solved = false;
            s.need_to.push_back(p);
        }

        if (dd > 0) {
            p.M = b_0;
            p.solved = false;
            s.need_to.push_back(p);

            p.M = b_1;
            p.solved = false;
            s.need_to.push_back(p);
        }

        if (aa > 0) {
            p.M = c_1;
            p.solved = false;
            s.need_to.push_back(p);
        }

        if (bb > 0) {
            p.M = d_1;
            p.solved = false;
            s.need_to.push_back(p);
        }
    }

    s.solved = s.good();
    if (s.M.col_sum(0) + s.M.col_sum(1) + s.M.col_sum(2) + s.M.col_sum(3) == 0) {
        storage s0;
        s0.M = orig;
        s0.solved = s0.good();
        return s0;
    }
    return s;
}

void El_rep_indexing::write_out_el_rep_vector() {
    cout<<"Fia: "<<Fia<<" Fib: "<<Fib<<" Fic: "<<Fic<<" Fid: "<<Fid;
    cout<<"  value:  "<<value<<endl;
    cout<<endl;
}