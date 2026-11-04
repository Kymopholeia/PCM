#include "SCF.h"
#include "BasicFunctions.h"
#include "Lebedev_grid.h"

SCF::SCF(double conv_threshold, int max_iter) :
        conv_tresh(conv_threshold), max_iter(max_iter) {}


void SCF::runSCF(BasisSet& Water, Eigen::MatrixXd& P_initial) {
    int electorn_num=0;
    for(auto a: Water.atoms){
    electorn_num+=a.Element;
    }
    nuc_rep_en= calculateNuclearRepulsionEnergy(Water);
    cout<<"nuclear rep energy: "<<nuc_rep_en<<endl;
    mx_size = Water.matrix_size;
    std::cout << "Matrix size: " << mx_size << std::endl;
    basis_size=Water.basis.size();
    cout<<"basis size: "<<basis_size<<endl;

    // Integrálok számítása
    Overlap overlap(Water);
//    cout<< "overlap mx: "<<endl<<overlap.integral_matrix<<endl;
    Kinetic kinetic(Water);
//    cout<< "kinetic mx: "<<endl<<kinetic.integral_matrix<<endl;

    Nuclear nuclear(Water);
//    cout<< "nuclear mx: "<<endl<<nuclear.integral_matrix<<endl;
    ObaraSaika electron_rep(Water);
    cout<<"ERI2 size is: "<<electron_rep.Eri2.size()<<endl;

    Fock_mx fock;
    fock.basis_size_ = Water.basis.size();
    fock.mx_size=mx_size;
    fock.calculateCoreHamiltonian(kinetic.integral_matrix, nuclear.integral_matrix);

    Eigen::MatrixXd P = P_initial;
    double energy_prev = 0.0;

    // SCF iterációk
    for(int iter = 0; iter < max_iter; iter++) {
        Eigen::MatrixXd F = fock.calculateFockMatrix(P, electron_rep.Eri2);
//        cout<<"Fock matrix: "<<endl<<F<<endl;

        double energy = fock.calculateEnergy(P,electron_rep.Eri2);
        std::cout << "Iteration " << iter << ", Energy: " << energy+nuc_rep_en << std::endl;

        // Konvergencia ellenőrzés
        if(iter > 0 && std::abs(energy - energy_prev) < conv_tresh) {
            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(F, overlap.integral_matrix);
            std::cout << "SCF converged!" << std::endl;
            std::cout << "Final energy: " << energy << std::endl;
            std::cout << "Iterations: " << iter + 1 << std::endl;

            total_energy=nuc_rep_en+energy;
            cout<<endl;
            cout<<endl;
            cout<<"The total energy is: "<<total_energy<<endl;
            cout<<endl;
            cout<<endl;
//            cout<<"Eltérés (sample molekula esetén):  "<<(abs(-34.9401447199-total_energy))/34.9401447199<<" ( "<<conv_tresh<<" kéne lennie )"<<endl;
            return;
        }
        energy_prev = energy;

        // Új sajátértékek és sajátvektorok számítása
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(F, overlap.integral_matrix);
        Eigen::MatrixXd C = es.eigenvectors();
//        cout<<"eigen vektorok:  "<<endl<<es.eigenvectors()<<endl;
//        cout<<"C matrix: "<<endl<<C<<endl;
//        cout<<"eigen values:  "<<endl<<es.eigenvalues()<<endl;
        std::cout << "Methanol orbital energies (vacuum, occupied):" << std::endl;
        for (int i = 0; i < 9; i++) {  // 9 occupied orbitals for 18 electrons
            std::cout << "  " << es.eigenvalues()[i] << std::endl;
        }

        // Új sűrűségmátrix építése
        Density_Matrix density;
        density.mx_size=mx_size;
        density.basis_size=basis_size;
        density.n_elect=electorn_num;
        P = density.constructFromCoeff(C,overlap.integral_matrix);
//        cout<<"Density matrix: "<<endl<<P<<endl;
        std::cout << "Iteration completed" << std::endl;
    }

    // Ha nem konvergált
    std::cout << "SCF did not converge within " << max_iter << " iterations" << std::endl;
}

void SCF::runSCF(BasisSet& Water, Eigen::MatrixXd& P_initial, Lebedev_grid& G) {
    int electorn_num=0;
    for(auto a: Water.atoms){
        electorn_num+=a.Element;
    }
    nuc_rep_en= calculateNuclearRepulsionEnergy(Water);
    cout<<"nuclear rep energy: "<<nuc_rep_en<<endl;
    mx_size = Water.matrix_size;
    std::cout << "Matrix size: " << mx_size << std::endl;
    basis_size=Water.basis.size();
    cout<<"basis size: "<<basis_size<<endl;

    // Integrálok számítása
    Overlap overlap(Water);
//    cout<< "overlap mx: "<<endl<<overlap.integral_matrix<<endl;
    Kinetic kinetic(Water);
//    cout<< "kinetic mx: "<<endl<<kinetic.integral_matrix<<endl;
    Nuclear nuclear(Water);

//    cout<< "nuclear mx: "<<endl<<nuclear.integral_matrix<<endl;
    ObaraSaika electron_rep(Water);
//    cout<<"ERI2 size is: "<<electron_rep.Eri2.size()<<endl;

    Fock_mx fock;
    fock.basis_size_ = Water.basis.size();
    fock.mx_size=mx_size;
    fock.calculateCoreHamiltonian(kinetic.integral_matrix, nuclear.integral_matrix);

    Eigen::MatrixXd P = P_initial;
    double energy_prev = 0.0;

    Pcm pcm;
    pcm.s_mx_size=G.tesserae.size();
    cout<<pcm.s_mx_size<<endl;
    pcm.cal_smatrix(G);
//    smatrix test can be placed here;
    pcm.build_pot_vector(Water,G,P);


    // SCF iterációk
    for(int iter = 0; iter < max_iter; iter++) {
        Eigen::MatrixXd F = fock.calculateFockMatrix(P, electron_rep.Eri2,pcm.q,pcm.Nuclear_pcm_matrix);
//        cout<<"Fock matrix: "<<endl<<F<<endl;

        double energy = fock.calculateEnergy(P,electron_rep.Eri2,pcm.E_solv);
        std::cout << "Iteration " << iter << ", Energy: " << energy << std::endl;

        // Konvergencia ellenőrzés
        if(iter > 0 && std::abs(energy - energy_prev) < conv_tresh) {
            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(F, overlap.integral_matrix);
            std::cout << "SCF converged!" << std::endl;
            std::cout << "Final energy: " << energy << std::endl;
            std::cout << "Iterations: " << iter + 1 << std::endl;

            total_energy=nuc_rep_en+energy;
            cout<<endl;
            cout<<endl;
            cout<<"The total energy is: "<<total_energy<<endl;
            cout<<"The PCM solvation energy is: "<<pcm.E_solv*627.5095<<endl;
            cout<<endl;
            cout<<endl;
//            cout<<"Eltérés (sample molekula esetén):  "<<(abs(-34.9401447199-total_energy))/34.9401447199<<" ( "<<conv_tresh<<" kéne lennie )"<<endl;
            return;
        }
        energy_prev = energy;

        // Új sajátértékek és sajátvektorok számítása
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(F, overlap.integral_matrix);
        Eigen::MatrixXd C = es.eigenvectors();
//        cout<<"eigen vektorok:  "<<endl<<es.eigenvectors()<<endl;
//        cout<<"C matrix: "<<endl<<C<<endl;
//        cout<<"eigen values:  "<<endl<<es.eigenvalues()<<endl;

        // Új sűrűségmátrix építése
        Density_Matrix density;
        density.mx_size=mx_size;
        density.basis_size=basis_size;
        density.n_elect=electorn_num;
        P = density.constructFromCoeff(C,overlap.integral_matrix);
//        cout<<"Density matrix: "<<endl<<P<<endl;
        pcm.build_pot_vector(Water,G,P);
        std::cout << "Iteration completed" << std::endl;
        double E_pol = pcm.E_solv;
        std::cout << "Iter " << iter
                  << "  E (returned) = " << energy
                  << "  Σq = " << pcm.q.sum()
                  << "  E_pol = " << E_pol
                  << "  Tr(P·V_pcm) = " << (P * fock.V_pcm).trace()
                  << std::endl;
    }

    // Ha nem konvergált
    std::cout << "SCF did not converge within " << max_iter << " iterations" << std::endl;
}

double SCF::calculateNuclearRepulsionEnergy(BasisSet& Water) const
{
    double nuclearRepulsionEnergy = 0;
    // we iterate through every pair of atoms
    for (int i = 0; i < Water.atoms.size(); i++)
    {
        for (int j = 0; j < i; j++)
        {
            nuclearRepulsionEnergy += (double)(Water.atoms[i].Element) * (double)(Water.atoms[j].Element) /
                    distance_between_points(Water.atoms[i].Centrum,Water.atoms[j].Centrum);
        }
    }
    return nuclearRepulsionEnergy;
}