//
// Created by Domokos Kyra on 17/09/2024.
//

#include "Integrals.h"

void integral::StarterPack(int AZn, double Ax,double Ay, double Az, int asumA, double al,int BZn,double Bx,double By, double Bz, int asumB, double be,int CZn,double Cx,double Cy, double Cz, int asumC, double ga,int DZn,double Dx,double Dy, double Dz, int asumD, double de) {
    angsumA = asumA;
    A.x = Ax;
    A.y =Ay;
    A.z = Az;
    alpha = al;
    AZ=AZn;
    angsumB = asumB;
    B.x = Bx;
    B.y = By;
    B.z = Bz;
    beta = be;
    BZ=BZn;
    angsumC = asumC;
    C.x = Cx;
    C.y = Cy;
    C.z = Cz;
    gamma = ga;
    CZ=CZn;
    angsumD = asumD;
    D.x = Dx;
    D.y = Dy;
    D.z = Dz;
    delta = de;
    DZ=DZn;
}
void integral::set_type(int t){
    type=t;
}
bool integral::stored_it(Matrix M) {
    for(int i=0;i<needed.size();i++){
        if(M.same_matrix(needed[i])){
            return true;
        }
    }
    return false;
}

void integral::Make_matrixes(int a, int b, int c, int d, int m) {
    Matrix tmp;
    tmp.m=m;
    if(a+b+c+d==0){
        made_matrixe.push_back(tmp);
    }
    else{
        vector<point> one_a, one_b, one_c, one_d;
        for (int l = 0; l < 4; l++) {
            int am;
            if (l == 0) am = a;
            if (l == 1) am = b;
            if (l == 2) am = c;
            if (l == 3) am = d;
            for(int i = am; i >= 0; --i) {
                for(int j = am - i; j >= 0; --j) {
                    int k = am - i - j;
                    point tmp_p;
                    tmp_p.x = i;
                    tmp_p.y = j;
                    tmp_p.z = k;
                    if (l == 0) one_a.push_back(tmp_p);
                    if (l == 1) one_b.push_back(tmp_p);
                    if (l == 2) one_c.push_back(tmp_p);
                    if (l == 3) one_d.push_back(tmp_p);
                }
            }
        }

        for (int i = 0; i < one_a.size(); i++) {
            for (int j = 0; j < one_b.size(); j++) {
                for (int k = 0; k < one_c.size(); k++) {
                    for (int l = 0; l < one_d.size(); l++) {
                        Matrix tmp_m;
                        tmp_m.m=m;
                        tmp_m.change_data(0, 0, one_a[i].x);
                        tmp_m.change_data(1, 0, one_a[i].y);
                        tmp_m.change_data(2, 0, one_a[i].z);

                        tmp_m.change_data(0, 1, one_b[j].x);
                        tmp_m.change_data(1, 1, one_b[j].y);
                        tmp_m.change_data(2, 1, one_b[j].z);

                        tmp_m.change_data(0, 2, one_c[k].x);
                        tmp_m.change_data(1, 2, one_c[k].y);
                        tmp_m.change_data(2, 2, one_c[k].z);

                        tmp_m.change_data(0, 3, one_d[l].x);
                        tmp_m.change_data(1, 3, one_d[l].y);
                        tmp_m.change_data(2, 3, one_d[l].z);

//                        tmp_m.write_out_Matrix();
                        made_matrixe.push_back(tmp_m);
                    }
                }
            }
        }
    }
}

double integral::calculation(int m, double alpha, double beta, double gamma, double delta, double sigma_p,double sigma_q, int type, point A, point B, point C, point D, double QmT) {

    for(int i=0;i<cal_T.size();i++){
        if(cal_T[i].first==m){
            return cal_T[i].second;
        }
    }
    if(type==3){
        double Eta=gamma+delta;
        sigma_p=0;
        sigma_q=1/(Eta+Eta);
    }

    double Rq=calculate_R(sigma_p,sigma_q,alpha,beta,gamma, delta,type,A,B,C,D);
//     cout<<"R  "<<Rq<<endl;
    double U=calculate_U(sigma_p,sigma_q,alpha,beta,gamma, delta,A,B,C,D);
//     cout<<"U  "<<U<<endl;
    double T=calculate_T(Rq,sigma_p,sigma_q);
    if(type==3){
        T=(gamma+delta)*Rq;
    }
//    cout<<"T   "<<T<<endl;
    double v=calculate_v(sigma_p,sigma_q);
    if(type==3)v=(gamma+delta);
//    cout<<"v   "<<v<<endl;
    double FmT=better_simpson_with_step(30,T,m);
//    cout<<"Fmt: "<<FmT<<endl;
    double resoult=0;
    double U_q=sqrt(pow((pi*(1/alpha+gamma)),3))*pow(e,((-alpha*gamma*pow(distance_between_points(C,A),2))*(1/alpha+gamma)));
    if(type==0){
        double N1 = cal_norm(alpha);
        double N2 = cal_norm(beta);
        double N3 = cal_norm(gamma);
        double N4 = cal_norm(delta);

        // Az eredeti számításod
        resoult = U * sqrt(2.0 * v) * sqrt(2.0 / pi) * FmT;

        // Normalizált eredmény
        resoult= N1 * N2 * N3 * N4 * resoult;}
    if(type==1) {
        resoult = QmT;
    }
    if(type==2) {
        resoult=QmT;
    }
    if(type==3) {
        double N1 = cal_norm(gamma);
        double N2 = cal_norm(delta);
        double eta=gamma+delta;
        double sigma_qq = 1 / (2 * (gamma+delta));
        double pow_e_q=-gamma*delta*sigma_qq*2*pow(distance_between_points(C,D),2);
        double U_qq=pow(e,pow_e_q);
        double Uu=1*U_qq;
//        cout<<"AZ: "<<AZ<<"CZ:  "<<CZ<<"DZ:  "<<DZ<<endl;
        resoult=-Uu*(1/eta)*2*pi*FmT*pow(2*eta,m)*AZ;
        resoult=resoult*N1*N2;
    }

    cal_T.push_back(make_pair(m,resoult));
    return resoult;
}

void integral::make_equations(Matrix Ma, double alpha, double beta, double gamma, double delta, point A, point B, point C, point D, int type) {
    equation equ;
    equ.orig=Ma;
    double eta = gamma + delta;
    double xi = alpha + beta;
    double eta_xi = (eta / (xi + eta));
    double xi_eta = (xi / (xi + eta));
    if(type==2||type==3){
        eta_xi=1/(2*eta);
        xi_eta=1/(2*xi);
    }
    equ.xi_eta=eta_xi;
    double sigmap= 1 / (alpha + beta);
    double sigmaq = 1 / (gamma + delta);
    equ.RR= Cal_R(alpha,beta,gamma,delta, sigmap,sigmaq,type,A,B,C,D);
    Matrix new_matrix=Ma;
    point P,Q;
    if(type==0){
        P=adding_of_vectors_and_multiply(A,B,alpha,beta);
        P=multiply_vector(P,sigmap);
        Q=adding_of_vectors_and_multiply(C,D,gamma,delta);
        Q=multiply_vector(Q,sigmaq);
    }
    if(type>0){
        P=A;Q=C;
    }
    if(type==3){
        P=A;
        Q= adding_of_vectors_and_multiply(C,D,gamma,delta);
        Q=multiply_vector(Q,(1/(eta))); //nem kell egy *2 az eta elé?
//        sigmap=0;
//        sigmaq=(1/(2*eta));
    }

    int w_col=4;
    int row=0;

    if(Ma.col_sum(3)>0){w_col=3;equ.PQ = Q; equ.AD = D;
        row=what_row(w_col,Ma);
        equ.cc=((Ma.get_data(row,2))/(2*eta));
        equ.dd=((Ma.get_data(row,3)-1)/(2*eta));
        equ.aa=((Ma.get_data(row,0))/(2*(xi+eta)));
        equ.bb=(Ma.get_data(row,1)/(2*(xi+eta)));
        if(type==2||type==3){
            equ.aa=((Ma.get_data(row,0))/(4*(xi*eta)));
            equ.bb=(Ma.get_data(row,1)/(4*(xi*eta)));
        }

    }
    if(Ma.col_sum(2)>0){w_col=2;equ.PQ = Q; equ.AD = C;
        row=what_row(w_col,Ma);
        equ.cc=((Ma.get_data(row,2)-1)/(2*eta));
        equ.dd=((Ma.get_data(row,3))/(2*eta));
        equ.aa=((Ma.get_data(row,0))/(2*(xi+eta)));
        equ.bb=(Ma.get_data(row,1)/(2*(xi+eta)));
        if(type==2||type==3){
            equ.aa=((Ma.get_data(row,0))/(4*(xi*eta)));
            equ.bb=(Ma.get_data(row,1)/(4*(xi*eta)));
        }
    }
    if(Ma.col_sum(1)>0){w_col=1;equ.PQ = P; equ.AD = B;
        row=what_row(w_col,Ma);
        equ.aa=((Ma.get_data(row,0))/(2*xi));
        equ.bb=((Ma.get_data(row,1)-1)/(2*xi));
        equ.cc=((Ma.get_data(row,2))/(2*(xi+eta)));
        equ.dd=((Ma.get_data(row,3))/(2*(xi+eta)));
        if(type==2||type==3){
            equ.cc=((Ma.get_data(row,2))/(4*(xi*eta)));
            equ.dd=(Ma.get_data(row,3)/(4*(xi*eta)));
        }
    }
    if(Ma.col_sum(0)>0){w_col=0;equ.PQ = P; equ.AD = A;
        row=what_row(w_col,Ma);
        equ.aa=((Ma.get_data(row,0)-1)/(2*xi));
        equ.bb=((Ma.get_data(row,1))/(2*xi));
        equ.cc=((Ma.get_data(row,2))/(2*(xi+eta)));
        equ.dd=((Ma.get_data(row,3))/(2*(xi+eta)));
        if(type==2||type==3){
            equ.cc=((Ma.get_data(row,2))/(4*(xi*eta)));
            equ.dd=(Ma.get_data(row,3)/(4*(xi*eta)));
        }
    }

    if(w_col!=4){
        new_matrix = Ma;
        new_matrix.change_data(row, w_col, Ma.get_data(row, w_col) - 1);
        new_matrix.m = Ma.m;
        equ.f_0 = new_matrix;

        new_matrix = Ma;
        new_matrix.change_data(row, w_col, Ma.get_data(row, w_col) - 1);
        new_matrix.m = Ma.m + 1;
        equ.f_1 = new_matrix;
    }

    if(w_col==0||w_col==1){
        if(equ.aa>0)
        {
            new_matrix=Ma;
            if(w_col==0){
                new_matrix.change_data(row, w_col,Ma.get_data(row,w_col)-2);
            }
            if(w_col==1){
                new_matrix.change_data(row,0,Ma.get_data(row,0)-1);
                new_matrix.change_data(row,1,Ma.get_data(row,1)-1);
            }

            new_matrix.m=Ma.m;
            equ.a_0=new_matrix;

            new_matrix=Ma;
            if(w_col==0){new_matrix.change_data(row, w_col,Ma.get_data(row,w_col)-2);}
            if(w_col==1){
                new_matrix.change_data(row,0,Ma.get_data(row,0)-1);
                new_matrix.change_data(row,1,Ma.get_data(row,1)-1);
            }
            new_matrix.m=Ma.m+1;
            equ.a_1=new_matrix;
        }
        if(equ.bb>0)
        {
            new_matrix=Ma;
            if(w_col==1){
                new_matrix.change_data(row, w_col,Ma.get_data(row,w_col)-2);
            }
            if(w_col==0){
                new_matrix.change_data(row,w_col,Ma.get_data(row,w_col)-1);
                new_matrix.change_data(row,1,Ma.get_data(row,1)-1);
            }
            new_matrix.m=Ma.m;
            equ.b_0=new_matrix;

            new_matrix=Ma;
            if(w_col==1){
                new_matrix.change_data(row, 1,Ma.get_data(row,1)-2);
            }
            if(w_col==0){
                new_matrix.change_data(row,0,Ma.get_data(row,0)-1);
                new_matrix.change_data(row,1,Ma.get_data(row,1)-1);
            }
            new_matrix.m=Ma.m+1;
            equ.b_1=new_matrix;
        }
        if(equ.cc>0){
            new_matrix=Ma;
            new_matrix.change_data(row,w_col,Ma.get_data(row,w_col)-1);
            new_matrix.change_data(row,2,Ma.get_data(row,2)-1);
            new_matrix.m=Ma.m+1;
            equ.c_1=new_matrix;
        }
        if(equ.dd>0){

            new_matrix=Ma;
            new_matrix.change_data(row,w_col,Ma.get_data(row,w_col)-1);
            new_matrix.change_data(row,3,Ma.get_data(row,3)-1);
            new_matrix.m=Ma.m+1;
            equ.d_1=new_matrix;
        }


    }
    if(w_col==2||w_col==3){
        if(equ.cc>0)
        {
            new_matrix=Ma;
            if(w_col==2){
                new_matrix.change_data(row, w_col,Ma.get_data(row,w_col)-2);
            }
            else{
                new_matrix.change_data(row,2,Ma.get_data(row,2)-1);
                new_matrix.change_data(row,3,Ma.get_data(row,3)-1);
            }

            new_matrix.m=Ma.m;
            equ.a_0=new_matrix;

            new_matrix=Ma;
            if(w_col==2){new_matrix.change_data(row, w_col,Ma.get_data(row,w_col)-2);}
            else{
                new_matrix.change_data(row,2,Ma.get_data(row,2)-1);
                new_matrix.change_data(row,3,Ma.get_data(row,3)-1);
            }
            new_matrix.m=Ma.m+1;
            equ.a_1=new_matrix;
        }
        if(equ.dd>0)
        {
            new_matrix=Ma;
            if(w_col==3){
                new_matrix.change_data(row, w_col,Ma.get_data(row,w_col)-2);
            }
            else{
                new_matrix.change_data(row,2,Ma.get_data(row,2)-1);
                new_matrix.change_data(row,3,Ma.get_data(row,3)-1);
            }
            new_matrix.m=Ma.m;
            equ.b_0=new_matrix;

            new_matrix=Ma;
            if(w_col==3){
                new_matrix.change_data(row, w_col,Ma.get_data(row,w_col)-2);
            }
            else{
                new_matrix.change_data(row,2,Ma.get_data(row,2)-1);
                new_matrix.change_data(row,3,Ma.get_data(row,3)-1);
            }
            new_matrix.m=Ma.m+1;
            equ.b_1=new_matrix;
        }
        if(equ.aa>0){
            new_matrix=Ma;
            new_matrix.change_data(row,w_col,Ma.get_data(row,w_col)-1);
            new_matrix.change_data(row,0,Ma.get_data(row,0)-1);
            new_matrix.m=Ma.m+1;
            equ.c_1=new_matrix;
        }

        if(equ.bb>0){

            new_matrix=Ma;
            new_matrix.change_data(row,w_col,Ma.get_data(row,w_col)-1);
            new_matrix.change_data(row,1,Ma.get_data(row,1)-1);
            new_matrix.m=Ma.m+1;
            equ.d_1=new_matrix;
        }

    }
//cout<<"here"<<endl;
    equ.orig.m=Ma.m;
//    equ.ki_equ_matrix();
    equations.push_back(equ);
    if(!stored_it(equ.f_0)){needed.push_back(equ.f_0); }
    if(!stored_it(equ.f_1)){needed.push_back(equ.f_1); }
    if(!stored_it(equ.a_0)){needed.push_back(equ.a_0); }
    if(!stored_it(equ.a_1)){needed.push_back(equ.a_1); }
    if(!stored_it(equ.b_0)){needed.push_back(equ.b_0); }
    if(!stored_it(equ.b_1)){needed.push_back(equ.b_1); }
    if(!stored_it(equ.c_1)){needed.push_back(equ.c_1); }
    if(!stored_it(equ.d_1)){needed.push_back(equ.d_1); }
}

int integral::what_row(int col, Matrix M) {
    for(int i=0;i<3;i++){
        if(M.get_data(i,col)>0)
            return i;
    }
}

void integral::finder(Matrix M, double data) {
    for(size_t i=0;i<Stored_matrixes.size();i++){
        if(M.same_matrix(Stored_matrixes[i].M)){
            Stored_matrixes[i].solved=true;
            Stored_matrixes[i].result_mx=data;
        }
    }
    for(size_t i=0;i<Stored_matrixes.size();i++){
        for(size_t j=0;j<Stored_matrixes[i].need_to.size();j++){
            if(M.same_matrix(Stored_matrixes[i].need_to[j].M)){
                Stored_matrixes[i].need_to[j].solved=true;
            }
        }
    }
    for(size_t i=0;i<Stored_matrixes.size();i++){
        Stored_matrixes[i].check_good();
    }
}

void integral::download(double &f, Matrix M) {
    for(size_t i=0;i<Stored_matrixes.size();i++){
        if(M.same_matrix(Stored_matrixes[i].M)){
            f=Stored_matrixes[i].result_mx;
        }
    }
}

void integral::what_2(storage s, double alpha, double beta, double gamma, double delta, int type, point A, point B,point C, point D) {
    double sigmap= 1 / (alpha + beta);
    double sigmaq = 1 / (gamma + delta);
    double eta,xi,eta_xi,xi_eta;
    eta = gamma + delta;
    xi = alpha + beta;
    eta_xi = (eta / (xi + eta));
    xi_eta = (xi / (xi + eta));
    if(type==3){
        sigmap=0;
        sigmaq=(1/(eta+eta));
        xi=0;
    }
//    s.M.write_out_Matrix();

    point P,Q,R;
    if(type==0){
        P=adding_of_vectors_and_multiply(A,B,alpha,beta);
        P=multiply_vector(P,sigmap);
        Q=adding_of_vectors_and_multiply(C,D,gamma,delta);
        Q=multiply_vector(Q,sigmaq);
        R= Cal_R(alpha,beta,gamma,delta,sigmap,sigmaq,type,A,B,C,D);
    }
    if(type>0){
        P=A;
        Q=C;
        R= Cal_R(alpha,beta,gamma,delta,sigmap,sigmaq,type,A,B,C,D);
    }
    if(type==3){
        P=A;
        Q= adding_of_vectors_and_multiply(C,D,gamma,delta);
        Q=multiply_vector(Q,(1/(eta+eta))*2);
        R= Cal_R(alpha,beta,gamma,delta,sigmap,sigmaq,type,A,B,C,D);
    }
    if(!s.good()) {

        for (size_t i = 0; i < s.need_to.size(); i++){
            if(!s.need_to[i].solved){
                for(size_t j=0;j<Stored_matrixes.size();j++){
                    if(s.need_to[i].M.same_matrix(Stored_matrixes[j].M)){
                        what_2(Stored_matrixes[j],alpha,beta,gamma,delta,type, A,B,C,D);

                    }
                }
            }
        }

    }

    Matrix M=s.M;
    double value=0;
    double aa,bb,cc,dd;
    double f_0,f_1,a_0,a_1,b_0,b_1,c_1,d_1;
    double pa,r;
    if(type==2||type==3){
        eta_xi=1/(2*xi);
        xi_eta=1/(2*eta);
    }



    if (M.col_sum(0) + M.col_sum(1) + M.col_sum(2) + M.col_sum(3) == 0) {
//        M.write_out_Matrix();
        if (type==0){
            value = calculation(M.m, alpha, beta, gamma, delta, sigmap, sigmaq, type, A, B, C, D,0.0);
        }
        if (type==1){
            value = calculation(M.m, alpha, beta, gamma, delta, sigmap, sigmaq, type, A, B, C, D,
                                calculate_PrimInt(alpha, beta, gamma, delta, M.m, calculate_T(
                                        calculate_R(sigmap, sigmaq, alpha, beta, gamma, delta, type, A, B, C, D),
                                        sigmap, sigmaq)));
        }
        if (type==2){
            value = calculation(M.m, alpha, beta, gamma, delta, sigmap, sigmaq, type, A, B, C, D,
                                calculate_PrimInt(alpha, beta, gamma, delta, M.m, calculate_T(
                                        calculate_R(sigmap, sigmaq, alpha, beta, gamma, delta, type, A, B, C, D),
                                        sigmap, sigmaq)));
        }
        if (type==3){
            value = calculation(M.m, alpha, beta, gamma, delta, sigmap, sigmaq, type, A, B, C, D,0.0);
        }
        finder(M, value);

    }
    if (M.col_sum(0) + M.col_sum(1) + M.col_sum(2) + M.col_sum(3) == 1) {
//            cout<<"01"<<endl;
        int x, y, z;
        double f0, f1;
        Matrix M_0, M_1;


        for (size_t i = 0; i < Stored_matrixes.size(); i++) {
            if (M.same_matrix(Stored_matrixes[i].M)) {
                M_0 = Stored_matrixes[i].need_to[0].M;
//                cout<<"M_0:  ";M_0.write_out_Matrix();
                M_1 = Stored_matrixes[i].need_to[1].M;
//                cout<<"M_1:  ";M_1.write_out_Matrix();
            }
        }

        if (M.col_sum(0) == 1) {
            if (M.get_data(0, 0) == 1) {
                pa = P.x - A.x;
                r = -R.x;
            }
            if (M.get_data(1, 0) == 1) {
                pa = P.y - A.y;
                r = -R.y;

            }
            if (M.get_data(2, 0) == 1) {
                pa = P.z - A.z;
                r = -R.z;
//                cout<<"here!"<<endl;
            }
        }
        if (M.col_sum(1) == 1) {
            if (M.get_data(0, 1) == 1) {
                pa = P.x - B.x;
                r = -R.x;
            }
            if (M.get_data(1, 1) == 1) {
                pa = P.y - B.y;
                r = -R.y;
            }
            if (M.get_data(2, 1) == 1) {
                pa = P.z - B.z;
                r = -R.z;
            }
        }
        if (M.col_sum(2) == 1) {
            if (M.get_data(0, 2) == 1) {
                pa = Q.x - C.x;
                r = R.x;
                eta_xi=xi_eta;
            }
            if (M.get_data(1, 2) == 1) {
                pa = Q.y - C.y;
                r = R.y;
                eta_xi=xi_eta;
            }
            if (M.get_data(2, 2) == 1) {
                pa = Q.z - C.z;
                r = R.z;
                eta_xi=xi_eta;
            }
        }
        if (M.col_sum(3) == 1) {
            if (M.get_data(0, 3) == 1) {
                pa = Q.x - D.x;
                r = R.x;
                eta_xi=xi_eta;
            }
            if (M.get_data(1, 3) == 1) {
                pa = Q.y - D.y;
                r = R.y;
                eta_xi=xi_eta;
            }
            if (M.get_data(2, 3) == 1) {
                pa = Q.z - D.z;
                r = R.z;
                eta_xi=xi_eta;
            }
        }
//        cout<<"pa  "<<pa<<endl;
//        cout<<"r  "<<r<<endl;
//        cout<<"1/2eta:  "<<eta_xi<<endl;
//        cout<<"M_0:  ";M_0.write_out_Matrix();
        if (type==0) {
            f0 = calculation(M_0.m, alpha, beta, gamma, delta, sigmap, sigmaq, type, A, B, C, D,
                             calculate_PrimInt(alpha, beta, gamma, delta, M_0.m, calculate_T(
                                     calculate_R(sigmap, sigmaq, alpha, beta, gamma, delta, type, A, B, C, D), sigmap,
                                     sigmaq)));
//        cout<<"M_1:  ";M_1.write_out_Matrix();
            f1 = calculation(M_1.m, alpha, beta, gamma, delta, sigmap, sigmaq, type, A, B, C, D,
                             calculate_PrimInt(alpha, beta, gamma, delta, M_1.m, calculate_T(
                                     calculate_R(sigmap, sigmaq, alpha, beta, gamma, delta, type, A, B, C, D), sigmap,
                                     sigmaq)));
        }
        if (type==1) {
            f0 = calculation(M_0.m, alpha, beta, gamma, delta, sigmap, sigmaq, type, A, B, C, D,
                             calculate_PrimInt(alpha, beta, gamma, delta, M_0.m, calculate_T(
                                     calculate_R(sigmap, sigmaq, alpha, beta, gamma, delta, type, A, B, C, D), sigmap,
                                     sigmaq)));
//        cout<<"M_1:  ";M_1.write_out_Matrix();
            f1 = calculation(M_1.m, alpha, beta, gamma, delta, sigmap, sigmaq, type, A, B, C, D,
                             calculate_PrimInt(alpha, beta, gamma, delta, M_1.m, calculate_T(
                                     calculate_R(sigmap, sigmaq, alpha, beta, gamma, delta, type, A, B, C, D), sigmap,
                                     sigmaq)));
        }
        if (type==2) {
            f0 = calculation(M_0.m, alpha, beta, gamma, delta, sigmap, sigmaq, type, A, B, C, D,
                             calculate_PrimInt(alpha, beta, gamma, delta, M_0.m, calculate_T(
                                     calculate_R(sigmap, sigmaq, alpha, beta, gamma, delta, type, A, B, C, D), sigmap,
                                     sigmaq)));
//        cout<<"M_1:  ";M_1.write_out_Matrix();
            f1 = calculation(M_1.m, alpha, beta, gamma, delta, sigmap, sigmaq, type, A, B, C, D,
                             calculate_PrimInt(alpha, beta, gamma, delta, M_1.m, calculate_T(
                                     calculate_R(sigmap, sigmaq, alpha, beta, gamma, delta, type, A, B, C, D), sigmap,
                                     sigmaq)));
        }
        if (type==3) {
            f0 = calculation(M_0.m, alpha, beta, gamma, delta, sigmap, sigmaq, type, A, B, C, D,0.0);
//        cout<<"M_1:  ";M_1.write_out_Matrix();
            f1 = calculation(M_1.m, alpha, beta, gamma, delta, sigmap, sigmaq, type, A, B, C, D,0.0);
        }

        value = (pa * f0) + ( eta_xi* r * f1);

        finder(M, value);
    }
    if (M.col_sum(0) + M.col_sum(1) + M.col_sum(2) + M.col_sum(3) > 1) {

        Matrix N = s.need_to[0].M;
//        M.write_out_Matrix();

        int X = -1;
        int Y = -1;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 4; j++) {
                if (M.get_data(i, j) - N.get_data(i, j) != 0) {
                    X = i;
                    Y = j;
                }
            }
        }
        if (X == -1 || Y == -1)cout << "ERROR";
//        cout<<"X= "<<X<<"Y=  "<<Y<<endl;

        if (Y == 0) { //a
            if (X == 0) {
                pa = P.x - A.x;
                r = -R.x;
            }
            if (X == 1) {
                pa = P.y - A.y;
                r = -R.y;
            }
            if (X == 2) {
                pa = P.z - A.z;
                r = -R.z;
            }
        }

        if (Y == 1) { //b
            if (X == 0) {
                pa = P.x - B.x;
                r = -R.x;
            }
            if (X == 1) {
                pa = P.y - B.y;
                r = -R.y;
            }
            if (X == 2) {
                pa = P.z - B.z;
                r = -R.z;
            }
        }

        if (Y == 2) { //c
            if (X == 0) {
                pa = Q.x - C.x;
                r = R.x;
            }
            if (X == 1) {
                pa = Q.y - C.y;
                r = R.y;
            }
            if (X == 2) {
                pa = Q.z - C.z;
                r = R.z;
            }
        }

        if (Y == 3) { //d
            if (X == 0) {
                pa = Q.x - D.x;
                r = R.x;
            }
            if (X == 1) {
                pa = Q.y - D.y;
                r = R.y;
            }
            if (X == 2) {
                pa = Q.z - D.z;
                r = R.z;
            }
        }

        if (Y == 0) {
            aa = ((M.get_data(X, 0) - 1) / (2 * xi));
            bb = (M.get_data(X, 1) / (2 * xi));
            cc = (M.get_data(X, 2) / (2 * (xi + eta)));
            dd = (M.get_data(X, 3) / (2 * (xi + eta)));
            if(type==2||type==3){
                cc = (M.get_data(X, 2) / (4 * (xi * eta)));
                dd = (M.get_data(X, 3) / (4 * (xi * eta)));
            }

        }
        if (Y == 1) {
            aa = (M.get_data(X, 0) / (2 * xi));
            bb = ((M.get_data(X, 1) - 1) / (2 * xi));
            cc = (M.get_data(X, 2) / (2 * (xi + eta)));
            dd = (M.get_data(X, 3) / (2 * (xi + eta)));
            if(type==2||type==3){
                cc = (M.get_data(X, 2) / (4 * (xi * eta)));
                dd = (M.get_data(X, 3) / (4 * (xi * eta)));
            }
        }
        if (Y == 2) {
            cc = ((M.get_data(X, 2) - 1) / (2 * eta));
            dd = ((M.get_data(X, 3)) / (2 * eta));
            aa = ((M.get_data(X, 0)) / (2 * (xi + eta)));
            bb = (M.get_data(X, 1) / (2 * (xi + eta)));
            if(type==2||type==3){
                aa = (M.get_data(X, 0) / (4 * (xi * eta)));
                bb = (M.get_data(X, 1) / (4 * (xi * eta)));
            }
            eta_xi=xi_eta;
        }
        if (Y == 3) {
            cc = ((M.get_data(X, 2)) / (2 * eta));
            dd = ((M.get_data(X, 3) - 1) / (2 * eta));
            aa = ((M.get_data(X, 0)) / (2 * (xi + eta)));
            bb = (M.get_data(X, 1) / (2 * (xi + eta)));
            if(type==2||type==3){
                aa = (M.get_data(X, 0) / (4 * (xi * eta)));
                bb = (M.get_data(X, 1) / (4 * (xi * eta)));
            }
            eta_xi=xi_eta;
        }
        download(f_0, s.need_to[0].M);
        download(f_1, s.need_to[1].M);

        a_0 = 0;
        a_1 = 0;
        b_0 = 0;
        b_1 = 0;
        c_1 = 0;
        d_1 = 0;

        if (Y == 0 || Y == 1) {
            if (aa > 0) {
                if (bb > 0) {
                    if (cc > 0) {
                        if (dd > 0) {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                            download(b_0, s.need_to[4].M);
                            download(b_1, s.need_to[5].M);
                            download(c_1, s.need_to[6].M);
                            download(d_1, s.need_to[7].M);
                        } else {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                            download(b_0, s.need_to[4].M);
                            download(b_1, s.need_to[5].M);
                            download(c_1, s.need_to[6].M);
                        }
                    } else {
                        if (dd > 0) {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                            download(b_0, s.need_to[4].M);
                            download(b_1, s.need_to[5].M);
                            download(d_1, s.need_to[6].M);
                        } else {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                            download(b_0, s.need_to[4].M);
                            download(b_1, s.need_to[5].M);
                        }
                    }
                } else {
                    if (cc > 0) {
                        if (dd > 0) {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                            download(c_1, s.need_to[4].M);
                            download(d_1, s.need_to[5].M);
                        } else {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                            download(c_1, s.need_to[4].M);
                        }
                    } else {
                        if (dd > 0) {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                            download(d_1, s.need_to[4].M);
                        } else {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                        }
                    }
                }

            } else {
                if (bb > 0) {
                    if (cc > 0) {
                        if (dd > 0) {
                            download(b_0, s.need_to[2].M);
                            download(b_1, s.need_to[3].M);
                            download(c_1, s.need_to[4].M);
                            download(d_1, s.need_to[5].M);
                        } else {
                            download(b_0, s.need_to[2].M);
                            download(b_1, s.need_to[3].M);
                            download(c_1, s.need_to[4].M);
                        }
                    } else {
                        if (dd > 0) {
                            download(b_0, s.need_to[2].M);
                            download(b_1, s.need_to[3].M);
                            download(d_1, s.need_to[4].M);
                        } else {
                            download(b_0, s.need_to[2].M);
                            download(b_1, s.need_to[3].M);
                        }
                    }
                } else {
                    if (cc > 0) {
                        if (dd > 0) {
                            download(c_1, s.need_to[2].M);
                            download(d_1, s.need_to[3].M);
                        } else {
                            download(c_1, s.need_to[2].M);
//                            s.need_to[2].M.write_out_Matrix();
                        }
                    } else {
                        if (dd > 0) {
                            download(d_1, s.need_to[2].M);
                        }
                    }
                }
            }
        }
        else {
            if (cc > 0) {
                if (dd > 0) {
                    if (aa > 0) {
                        if (bb > 0) {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                            download(b_0, s.need_to[4].M);
                            download(b_1, s.need_to[5].M);
                            download(c_1, s.need_to[6].M);
                            download(d_1, s.need_to[7].M);
                        } else {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                            download(b_0, s.need_to[4].M);
                            download(b_1, s.need_to[5].M);
                            download(c_1, s.need_to[6].M);
                        }

                    } else {
                        if (bb > 0) {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                            download(b_0, s.need_to[4].M);
                            download(b_1, s.need_to[5].M);
                            download(d_1, s.need_to[6].M);
                        } else {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                            download(b_0, s.need_to[4].M);
                            download(b_1, s.need_to[5].M);
                        }

                    }

                } else {
                    if (aa > 0) {
                        if (bb > 0) {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                            download(c_1, s.need_to[4].M);
                            download(d_1, s.need_to[5].M);
                        } else {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                            download(c_1, s.need_to[4].M);
                        }

                    } else {
                        if (bb > 0) {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                            download(d_1, s.need_to[4].M);
                        } else {
                            download(a_0, s.need_to[2].M);
                            download(a_1, s.need_to[3].M);
                        }
                    }
                }

            } else {
                if (dd > 0) {
                    if (aa > 0) {
                        if (bb > 0) {
                            download(b_0, s.need_to[2].M);
                            download(b_1, s.need_to[3].M);
                            download(c_1, s.need_to[4].M);
                            download(d_1, s.need_to[5].M);
                        } else {
                            download(b_0, s.need_to[2].M);
                            download(b_1, s.need_to[3].M);
                            download(c_1, s.need_to[4].M);
                        }

                    } else {
                        if (bb > 0) {
                            download(b_0, s.need_to[2].M);
                            download(b_1, s.need_to[3].M);
                            download(d_1, s.need_to[4].M);
                        } else {
                            download(b_0, s.need_to[2].M);
                            download(b_1, s.need_to[3].M);
                        }

                    }

                } else {
                    if (aa > 0) {
                        if (bb > 0) {
                            download(c_1, s.need_to[2].M);
                            download(d_1, s.need_to[3].M);
                        } else {
                            download(c_1, s.need_to[2].M);
                        }

                    } else {
                        if (bb > 0) {
                            download(d_1, s.need_to[2].M);
                        }
                    }
                }
            }
        }
        if (Y == 0 || Y == 1) {
//            cout<<"here"<<endl;

            value = (pa * f_0) + (eta_xi * r * f_1) + (aa * (a_0 - (eta_xi * a_1))) +
                    (bb * (b_0 - (eta_xi * b_1))) + (cc * c_1) + (dd * d_1);
        }
        else {
            if(type==3){
                value = (pa * f_0) + (xi_eta * r * f_1) + (cc * (a_0 - (xi_eta * a_1))) +
                        (dd * (b_0 - (xi_eta * b_1)));
            }
            else{
                value = (pa * f_0) + (xi_eta * r * f_1) + (cc * (a_0 - (xi_eta * a_1))) +
                        (dd * (b_0 - (xi_eta * b_1))) + (aa * c_1) + (bb * d_1);
            }
        }
        finder(M, value);
    }
}

void integral::nom_mx_cal() {

}

double integral::cal_norm(double alp) {
    return pow(2.0 * alp / M_PI, 0.75);

}