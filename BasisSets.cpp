//
// Created by Domokos Kyra on 26/09/2024.
//

#include "BasisSets.h"
#include "json.hpp"
#include "fstream"


using namespace std;
using json = nlohmann::json;

int BasisFunction::decode_type(std::string Type) {
    if(Type=="s")return 0;
    if(Type=="p")return 1;
    if(Type=="d")return 2;
    if(Type=="f")return 3;

    return 0;
}

void BasisFunction::write_out_bF() {
    std::cout << "RI: " << RI  << std::endl;
    std::cout << "FI: " << FI  << std::endl;
    std::cout << "Element name: " << Element_name << std::endl;
    std::cout << "Element: " <<element << std::endl;
    std::cout << "Position: (" << Centrum.x << ", " << Centrum.y << ", " <<Centrum.z << ")" << std::endl;
    std::cout << "  Type: " <<type << ", Exponent: " << exponent << std::endl;
    cout<<endl;
}

bool compare_by_type(const BasisFunction &a, const BasisFunction &b) {
    return a.type < b.type;  // Sort in ascending order of 'type'
}

// Method to sort the basis vector by the type of basis function
void BasisSet::sort_basis_by_type() {
    std::sort(basis.begin(), basis.end(), compare_by_type);
}

int Element_DataBase::get_element() {
    if(name=="O") { return 8; }
    if(name=="H") { return 1; }
    if(name=="C") { return 6; }
    if(name=="N") { return 7; }

    return -1;
}

int Atom::Get_element(string elem) {
    Element_DataBase EDB;
    EDB.name=elem;
    if(EDB.get_element()==-1) { cout << "New element type, please update the Element_Database.cpp file" << endl; }
    else{return EDB.get_element(); }
    return -1;
}

double Atom::get_VdW_radii() {
    double tuner=1.2;
    switch (Element) {
        case 1: return 1.2*1.8897259886*tuner;//angstrom H
        case 6: return 1.7*1.8897259886*tuner; //C
        case 7: return 1.55*1.8897259886*tuner; //N
        case 8: return 1.52*1.8897259886*tuner; //O
        case -1: return -1;
    }
    return-1;
}

void BasisSet::add_basis_set() {
    Atom atom;
    atom.Element = 8;
    atom.Centrum.x = 0.0;
    atom.Centrum.y =0.0;
    atom.Centrum.z =0.0;

    atom.centrum.x = 0.0;
    atom.centrum.y =0.0;
    atom.centrum.z =0.0;

            BasisFunction basisFunc;
    basisFunc.element = atom.Element;
    basisFunc.Centrum = atom.Centrum;
    basisFunc.centrum = atom.centrum;
    basisFunc.type = 1;
    basisFunc.exponent = 2.0;

    // Add the basis function to the basis list
    basis.push_back(basisFunc);

}

int BasisSet::cal_mat_size() {

//    for ciklus+swich case a type ra
//    matrix_size = basis[basis.size()-1].FI+ return_type_num(basis[basis.size()-1].type);
int size=0;
    for (int i=0;i<basis.size();i++){
        size+= return_type_dec(basis[i].type);
    }
    matrix_size=size;
//    cout<<"Here is the matrix size: "<<matrix_size<<'x'<<matrix_size<<endl;

}

void BasisSet::Make_Basis_Functions() {
    // Open the JSON file
//    cout << fs::current_path();
    ifstream file("../"+filename);
//
    // Parse the JSON file
    json jsonData;
    file >> jsonData;

    // Extract the molecule name
    Name = jsonData["molecule"]["name"];
    // Extract atom and basis function data
    for (const auto& atom_json : jsonData["molecule"]["atoms"]) {
        // Create an Atom instance
        Atom atom;
        atom.element_name=atom_json["element"];
        atom.Element = atom.Get_element(atom_json["element"]);
        atom.Centrum.x = atom_json["position"]["x"];
        atom.Centrum.y = atom_json["position"]["y"];
        atom.Centrum.z = atom_json["position"]["z"];

//        code mergeing problem hardcode solution
        atom.centrum.x = atom_json["position"]["x"];
        atom.centrum.y = atom_json["position"]["y"];
        atom.centrum.z = atom_json["position"]["z"];

        atom.Centrum.x= atom.Centrum.x*1.8897259886;
        atom.Centrum.y= atom.Centrum.y*1.8897259886;
        atom.Centrum.z= atom.Centrum.z*1.8897259886;
//        code mergeing problem hardcode solution
        atom.centrum.x= atom.centrum.x*1.8897259886;
        atom.centrum.y= atom.centrum.y*1.8897259886;
        atom.centrum.z= atom.centrum.z*1.8897259886;

        // Add the Atom to the atom list
        atoms.push_back(atom);

        // Create BasisFunction instances for each basis in the atom
        for (const auto& basis_json : atom_json["basis"]) {
            BasisFunction basisFunc;
            basisFunc.Element_name=atom.element_name;
            basisFunc.element = atom.Element;
            basisFunc.Centrum = atom.Centrum;
            basisFunc.type = basisFunc.decode_type(basis_json["type"]);
            basisFunc.exponent = basis_json["exponent"];

            // Add the basis function to the basis list
            basis.push_back(basisFunc);
        }
    }
//    add_basis_set();
    std::reverse(basis.begin(), basis.end());
    sort_basis_by_type();
    cal_mat_size();

    for(int i=0;i<basis.size();i++){
        basis[i].RI=i;
        if(i>0){
            if(basis[i-1].type>0) {
                basis[i].FI = basis[i - 1].FI + (basis[i - 1].type * 3);
            }
            else{
                basis[i].FI = basis[i - 1].FI +1;
            }
        }
        else{
            basis[i].FI=i;
        }
    }
}



void BasisSet::write_out_BasisSet() {

    for (const auto& basiss : basis){
        std::cout << "RI: " << basiss.RI  << std::endl;
        std::cout << "FI: " << basiss.FI  << std::endl;
        std::cout << "Element name: " << basiss.Element_name << std::endl;
        std::cout << "Element: " << basiss.element << std::endl;
        std::cout << "Position: (" << basiss.Centrum.x << ", " << basiss.Centrum.y << ", " << basiss.Centrum.z << ")" << std::endl;
        std::cout << "  Type: " << basiss.type << ", Exponent: " << basiss.exponent << std::endl;
        cout<<endl;
    }

}

void BasisSet::write_out_JsonFile() {
    std::cout << "Molecule Name: " << Name << std::endl;


    for (const auto& atom : atoms) {
        std::cout << "Element: " << atom.Element << std::endl;
        std::cout << "Position: (" << atom.Centrum.x << ", " << atom.Centrum.y << ", " << atom.Centrum.z << ")" << std::endl;
        std::cout << "Basis Set: " << std::endl;
        for (const auto& basiss : basis) {
            if(basiss.element==atom.Element) {
                std::cout << "  Type: " << basiss.type << ", Exponent: " << basiss.exponent << std::endl;
            }
        }

        std::cout << std::endl;
    }
}

int return_type_dec(int kyra_type){
    switch (kyra_type) {
        case 0: return 1;
        case 1: return 3;
        case 2: return 6;
        case 3: return 8;
    }
    return 0;
}