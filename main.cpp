#include <iostream>
#include <vector>
#include "cmath"
#include "Lebedev_grid.h"
#include "BasisSets.h"
#include "Pcm.h"

using namespace std;
using json = nlohmann::json;



int main() {
    Lebedev_grid Grid;
    Grid.Lebedev_grid_();
//    Grid.write_out_grid_points();
    BasisSet Bs;
    Bs.filename="real_water.json";
//    Bs.filename="Oxygen.json";
    Bs.Make_Basis_Functions();
//    Bs.write_out_BasisSet();
    Grid.build_sphere(Bs,Grid);
    Grid.cal_tessera(Grid,Bs);
//    double kk=0;
//    for(Tessera k:Grid.tesserae){cout<<k.area<<endl;kk=kk+k.area;}
//    cout<<endl<<kk<<endl;
    return 0;
}
