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
////    double kk=0;
////    for(Tessera k:Grid.tesserae){cout<<k.area<<endl;kk=kk+k.area;}
////    cout<<endl<<kk<<endl;

    cout << setprecision(12);
    BasisSet Water;
    Water.filename="real_water.json";
    Water.Make_Basis_Functions();
    cout<<"mx size:  "<<Water.matrix_size<<endl;
    Eigen::MatrixXd P=Eigen::MatrixXd::Zero(Water.matrix_size,Water.matrix_size);
    SCF SCf(1e-8,100);
    SCf.runSCF(Water,P);
    return 0;
}
