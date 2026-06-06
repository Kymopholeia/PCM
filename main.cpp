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

    auto start = std::chrono::high_resolution_clock::now();
    if(TYPE==0){
    cout<<"Normal SCF"<<endl;
    cout << setprecision(12);
    BasisSet Water;
//    Water.filename="real_water.json";
//    Water.filename="hcn_h_polarized.json";
//    Water.filename="ammonia.json";
    Water.filename="methanol.json";
    Water.Make_Basis_Functions();
    cout<<"mx size:  "<<Water.matrix_size<<endl;
    Eigen::MatrixXd P=Eigen::MatrixXd::Zero(Water.matrix_size,Water.matrix_size);
    SCF SCf(1e-8,100);
    SCf.runSCF(Water,P);
        // === One-shot PCM energy test ===
// Run AFTER vacuum SCF has converged on methanol.
// This computes E_pol = (1/2) * q . v from the vacuum-converged P,
// without coupling PCM into the Fock matrix. If this gives a reasonable
// answer (~ few kcal/mol), the SCF coupling is the problem. If it also
// blows up, then the cavity/basis is the cause.
//        Lebedev_grid Grid;
//        Grid.Lebedev_grid_();
//        Grid.build_sphere(Water,Grid);
//        Grid.cal_tessera(Grid,Water);
//        Pcm Pc;
//        Pc.s_mx_size=Grid.tesserae.size();
//        cout<<Pc.s_mx_size<<endl;
//        Pc.cal_smatrix(Grid);
//        P=SCf.density_P;
//        std::cout << "=== P DIAGNOSTIC BEFORE ONE-SHOT ===" << std::endl;
//        std::cout << "P diagonal: ";
//        for (int i = 0; i < P.rows(); i++) {
//            std::cout << P(i,i) << " ";
//        }
//        std::cout << std::endl;
//        std::cout << "P(0,0..5) row: ";
//        for (int j = 0; j < 5; j++) {
//            std::cout << P(0,j) << " ";
//        }
//        Overlap O(Water);
//
//        std::cout << std::endl;
//        std::cout << "(P*S).trace() = " << (P * O.integral_matrix).trace() << std::endl;
//        std::cout << "===================================" << std::endl;
//
//        std::cout << std::endl;
//        std::cout << "=== ONE-SHOT PCM TEST (vacuum density) ===" << std::endl;
//        std::cout << "Vacuum SCF has converged; computing PCM energy from P_vacuum" << std::endl;
//        std::cout << "without iterating, to isolate basis/cavity vs. coupling effects." << std::endl;
//        std::cout << std::endl;
//
//// Step 1: rebuild v_elec from the vacuum-converged P (and v_nuc which is already there)
////         This uses your existing pcm.build_pot_vector code path.
//        Pc.build_pot_vector(Water, Grid, P);
//
//// Step 2: solve S q = -0.987 v ONCE, using the vacuum-converged v
//        Eigen::VectorXd V_eigen = Eigen::Map<Eigen::VectorXd>(
//                Pc.V_electric.data(), Pc.V_electric.size());
//        Eigen::VectorXd b_oneshot = -0.987 * V_eigen;
//        Eigen::LLT<Eigen::MatrixXd> llt_oneshot(Pc.S_matrix);
//        Eigen::VectorXd q_oneshot = llt_oneshot.solve(b_oneshot);
//
//// Step 3: report diagnostics
//        double sum_q_oneshot = q_oneshot.sum();
//        double q_min_oneshot = q_oneshot.minCoeff();
//        double q_max_oneshot = q_oneshot.maxCoeff();
//
//// Step 4: compute the one-shot polarization energy
//        double E_pol_oneshot = 0.5 * q_oneshot.dot(V_eigen);
//
//        std::cout << "Sum of charges (Gauss's law): " << sum_q_oneshot << std::endl;
//        std::cout << "  (should be near 0 for neutral solute)" << std::endl;
//        std::cout << "Min charge: " << q_min_oneshot << std::endl;
//        std::cout << "Max charge: " << q_max_oneshot << std::endl;
//        std::cout << std::endl;
//        std::cout << "One-shot polarization energy:" << std::endl;
//        std::cout << "  E_pol = 0.5 * q . v = " << E_pol_oneshot << " Hartree" << std::endl;
//        std::cout << "                     = " << E_pol_oneshot * 627.509 << " kcal/mol" << std::endl;
//        std::cout << std::endl;
//        std::cout << "Compare to full PCM-SCF E_solv:" << std::endl;
//        std::cout << "  (your previously reported value, e.g. for methanol: -1527 kcal/mol)" << std::endl;
//        std::cout << std::endl;
//        std::cout << "Interpretation:" << std::endl;
//        std::cout << "  - If one-shot ~ few kcal/mol AND full-SCF blows up:" << std::endl;
//        std::cout << "    -> the iterative coupling is the bug" << std::endl;
//        std::cout << "  - If one-shot also blows up (large magnitude):" << std::endl;
//        std::cout << "    -> basis-set / cavity issue (response is correct, basis is wrong)" << std::endl;
//        std::cout << "==========================================" << std::endl;
//
//        // Right after vacuum SCF converges on methanol
    }
    if(TYPE==1){
        cout<<"PCM SCF"<<endl;
        cout << setprecision(12);
        BasisSet Water;
//    Water.filename="ion_test.json";
//    Water.filename="real_water.json";
//    Water.filename="ammonia.json";
//        Water.filename="methanol.json";
        Water.Make_Basis_Functions();
        Lebedev_grid Grid;
        Grid.Lebedev_grid_();
        Grid.build_sphere(Water,Grid);
        Grid.cal_tessera(Grid,Water);
        cout<<"mx size:  "<<Water.matrix_size<<endl;
        std::vector<int> per_atom_count(Water.atoms.size(), 0);
        for (const auto& t : Grid.tesserae) {
            per_atom_count[t.r.no_of_atom]++;
        }
        int atoms_with_tesserae = 0;
        for (size_t i = 0; i < per_atom_count.size(); ++i) {
            std::cout << "  Atom " << i << " (" << Water.atoms[i].element_name
                      << "): " << per_atom_count[i] << " tesserae";
            if (per_atom_count[i] == 0) std::cout << "  <-- SWALLOWED";
            std::cout << std::endl;
            if (per_atom_count[i] > 0) atoms_with_tesserae++;
        }
        std::cout << "Atoms with tesserae: " << atoms_with_tesserae
                  << " / " << Water.atoms.size() << std::endl;
        std::cout << "Total tesserae: " << Grid.tesserae.size() << std::endl;
        std::cout << "============================================" << std::endl;

        // Compute cavity volume via divergence theorem
        double V_cavity = 0.0;
        Point r_center;  // pick any interior point - e.g., center of mass
        r_center.x = 0; r_center.y = 0; r_center.z = 0;
        for (const auto& a : Water.atoms) {
            r_center.x += a.centrum.x / Water.atoms.size();
            r_center.y += a.centrum.y / Water.atoms.size();
            r_center.z += a.centrum.z / Water.atoms.size();
        }

        double A_tot=0;
        for (size_t i = 0; i < Grid.tesserae.size(); ++i) {
            A_tot+=Grid.tesserae[i].area;
            Point r_t = Grid.tesserae[i].r;
            int parent = r_t.no_of_atom;
            Point r_parent = Water.atoms[parent].centrum;
            double R = Water.atoms[parent].get_VdW_radii();

            // Outward normal at tessera = (r_t - r_parent) / R
            double nx = (r_t.x - r_parent.x) / R;
            double ny = (r_t.y - r_parent.y) / R;
            double nz = (r_t.z - r_parent.z) / R;

            // (r_t - r_center) . n_hat
            double dx = r_t.x - r_center.x;
            double dy = r_t.y - r_center.y;
            double dz = r_t.z - r_center.z;
            double dot = dx*nx + dy*ny + dz*nz;

            V_cavity += dot * Grid.tesserae[i].area / 3.0;
        }

        std::cout << "Cavity volume: " << V_cavity << " Bohr^3" << std::endl;

// Compute the volume of one atom-sphere alone
        double V_C_alone = 4.0/3.0 * M_PI * pow(Water.atoms[0].get_VdW_radii(), 3);
        std::cout << "Volume of C sphere alone: " << V_C_alone << " Bohr^3" << std::endl;
        std::cout << "Ratio V_cavity / V_C_alone: " << V_cavity / V_C_alone << std::endl;

        double sphericity = pow(M_PI, 1.0/3.0) * pow(6.0 * V_cavity, 2.0/3.0) / A_tot;
        std::cout << "Sphericity: " << sphericity << std::endl;
        std::cout << "  1.0 = perfect sphere (cavity is ball-shaped)" << std::endl;
        std::cout << "  <0.95 = clearly non-spherical" << std::endl;
        cout<<endl;
        cout<<"==============================================================================="<<endl;
        cout<<endl;
        cout<<endl;
        Eigen::MatrixXd P=Eigen::MatrixXd::Zero(Water.matrix_size,Water.matrix_size);
        SCF SCf(1e-8,100);
        SCf.runSCF(Water,P,Grid);
    }
    cout<<endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    cout<<"==============================================================================="<<endl;
    std::cout << "Total Execution time: " << duration.count() << " seconds" << endl;
    cout<<"Calculation done."<<endl;
    cout<<endl;
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