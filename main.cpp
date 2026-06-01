#include <iostream>
#include <vector>
#include "cmath"
#include "Lebedev_grid.h"
#include "BasisSets.h"
#include "Pcm.h"
#include <iomanip>
//#include "ObaraSaika.h"
//#include "Overlap.h"
//#include "Kinetic.h"
//#include "Nuclear.h"
#include "BasisSets.h"
//#include "Eigen/dense"
//#include "Coulomb_Exchange.h"
//#include "FockMatrix.h"
#include "SCF.h"


using namespace std;
using json = nlohmann::json;


int main() {
    int TYPE=0;
    if(TYPE==0){
    cout<<"Normal SCF"<<endl;
    cout << setprecision(12);
    BasisSet Water;
//    Water.filename="real_water.json";
//    Water.filename="ammonia.json";
//    Water.filename="ammonia.json";
    Water.filename="hydrogen_peroxide.json";
    Water.Make_Basis_Functions();
    cout<<"mx size:  "<<Water.matrix_size<<endl;
    Eigen::MatrixXd P=Eigen::MatrixXd::Zero(Water.matrix_size,Water.matrix_size);
    SCF SCf(1e-8,100);
    SCf.runSCF(Water,P);
    }
    if(TYPE==1){
        cout<<"PCM SCF"<<endl;
        cout << setprecision(12);
        BasisSet Water;
//    Water.filename="real_water.json";
//    Water.filename="metane.json";
//    Water.filename="ammonia.json";
        Water.filename="hydrogen_peroxide.json";
        Water.Make_Basis_Functions();
        Lebedev_grid Grid;
        Grid.Lebedev_grid_();
        Grid.build_sphere(Water,Grid);
        Grid.cal_tessera(Grid,Water);
        cout<<"mx size:  "<<Water.matrix_size<<endl;
        Eigen::MatrixXd P=Eigen::MatrixXd::Zero(Water.matrix_size,Water.matrix_size);
        SCF SCf(1e-8,100);
        SCf.runSCF(Water,P,Grid);
    }

//



//    Lebedev_grid Grid;
//    Grid.Lebedev_grid_();
////    Grid.write_out_grid_points();
//    BasisSet Bs;
//    Bs.filename="real_water.json";
////    Bs.filename="Oxygen.json";
//    Bs.Make_Basis_Functions();
////    Bs.write_out_BasisSet();
//    Grid.build_sphere(Bs,Grid);
//    Grid.cal_tessera(Grid,Bs);
//    double kk=0;
//    for(Tessera k:Grid.tesserae){cout<<k.area<<endl;kk=kk+k.area;}
//    cout<<endl<<kk<<endl;



//    cout << setprecision(12);
//    BasisSet Water;
//    Water.filename="real_water.json";
//    Water.Make_Basis_Functions();
//    Lebedev_grid Grid;
//    Grid.Lebedev_grid_();
//    Grid.build_sphere(Water,Grid);
//    Grid.cal_tessera(Grid,Water);
//    cout<<"mx size:  "<<Water.matrix_size<<endl;
//    Eigen::MatrixXd P=Eigen::MatrixXd::Zero(Water.matrix_size,Water.matrix_size);
//    SCF SCf(1e-8,100);
//    SCf.runSCF(Water,P);
    return 0;
}





//std::vector<int> per_atom_count(Water.atoms.size(), 0);
//for (const auto& t : Grid.tesserae) {
//per_atom_count[t.r.no_of_atom]++;
//}
//for (size_t i = 0; i < per_atom_count.size(); ++i) {
//std::cout << "Atom " << i << " (" << Water.atoms[i].element_name
//<< "): " << per_atom_count[i] << " tesserae" << std::endl;
//}
//double total_area = 0.0;
//for (const auto& t : Grid.tesserae) {
//total_area += t.area;
//}
//std::cout << "Total cavity surface area: " << total_area << " Bohr²" << std::endl;
//std::cout << "Total tesserae: " << Grid.tesserae.size() << std::endl;


//S MATRIX TESTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//    cout<<pcm.S_matrix<<endl;
//    bool good=true;
//    for(int i=0;i<pcm.s_mx_size;i++){
//        for(int j=0;j<pcm.s_mx_size;j++){
////            cout<<abs(pcm.S_matrix(i,j)-pcm.S_matrix(j,i))<<endl;
////            if(abs(pcm.S_matrix(i,j)-pcm.S_matrix(j,i))>1e-12){
////                good=false;
//                if(i==j){
//                    cout<<pcm.S_matrix(i,j)<<endl;
//            }
//        }
//    }
//    cout<<good<<endl;
//Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(pcm.S_matrix);
//double lambda_min=es.eigenvalues()(0);
//double lambda_max=es.eigenvalues()(es.eigenvalues().size()-1);
//double cond=lambda_max/lambda_min;
//cout<<"Condition number: "<<cond<<"\n";

//double min_dist=1e30;
//int imin=-1, jmin=-1;
//for(int i=0;i<110;i++){
//    for(int j=i+1;j<110;j++){
//        double alma=pcm.dist_bet_po(Grid.tesserae[i].r,Grid.tesserae[j].r);
//        if(alma<min_dist){min_dist=alma;imin=i;jmin=j;}
//    }
//}
//cout<<"Min tess dist: "<<min_dist<<" between "<<imin<<" and "<<jmin<<endl;