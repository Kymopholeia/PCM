//
// Created by Domokos Kyra on 26/09/2024.
//

#ifndef ONLAB_NEWHEADER_BASISSETS_H
#define ONLAB_NEWHEADER_BASISSETS_H

//#include "Lebedev_grid.h"
#include "Point_class.h"
#include "iostream"
#include "json.hpp"
#include "fstream"


using json = nlohmann::json;
using namespace std;

class BasisFunction{
public:
    int RI,FI;
    string Element_name;
    int element;
    Point Centrum;
    int type;
    int decode_type(string type);
    double exponent;
    void write_out_bF();

};

class Atom{
public:
    string element_name;
    int Element;
    Point Centrum;
    int Get_element(string elem);
    double get_VdW_radii();
};

class BasisSet{
public:
    string Name;
    string filename;
    int matrix_size;
    vector<Atom> atoms;
    vector<BasisFunction> basis;
    int cal_mat_size();
    void add_basis_set();
    void sort_basis_by_type();
    void Make_Basis_Functions();
    void write_out_BasisSet();
    void write_out_JsonFile();
};
class Element_DataBase{
public:
    string name;
    int get_element();
};

int return_type_dec(int kyra_type);
#endif //ONLAB_NEWHEADER_BASISSETS_H
