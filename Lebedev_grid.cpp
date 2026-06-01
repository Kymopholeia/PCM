//
// Created by Domokos Kyra on 2026. 05. 26..
//

#include "Lebedev_grid.h"

void Lebedev_grid::Lebedev_grid_() {
    Point P;
    vector<double> point;
    point.push_back(0);
    point.push_back(0);
    point.push_back(0);
//    a1 grid points. 6 total [(1,0,0) and all permutations]
    for (int i=0;i<3;i++){
        P.x=0;P.y=0;P.z=0;
        point[0]=0;point[1]=0;point[2]=0;
        point[i]=1;
        P.x=point[0];P.y=point[1];P.z=point[2];P.weight=a1w;
        grid_points.push_back(P);
        point[i]=-1;
        P.x=point[0];P.y=point[1];P.z=point[2];P.weight=a1w;
        grid_points.push_back(P);
    }
//    a2 grid points 12 in total, for 110 its not needed.
//    for (int i = 0; i < 3; i++) {
//        for (int j = i + 1; j < 3; j++) {
//            vector<double> point;
//            P.x=0;P.y=0;P.z=0;
//            point[i] = A2*1;
//            point[j] = A2*1;
//            grid_points.push_back(P);
//            point[j] = A2*-1;
//            grid_points.push_back(P);
//            point[i] = A2*-1;
//            grid_points.push_back(P);
//            point[j] = A2*1;
//            grid_points.push_back(P);
//        }
//    }
    point[0] = A3;   // = 1/sqrt(3)
    point[1] = A3;
    point[2] = A3;
    for (int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            for (int k=0;k<2;k++){
                point[2]*=-1;
                P.x=point[0];P.y=point[1];P.z=point[2];P.weight=a3w;
                grid_points.push_back(P);
            }
            point[1]*=-1;
        }
        point[0]*=-1;
    }
//b1 b2 and b3 for three time 24
    double B1l;
    double B1m;
    for(int z=0;z<3;z++) {
        if (z == 0) {B1l = B1;}
        if (z == 1) {B1l = B2;}
        if (z == 2) {B1l = B3;}
        B1m = sqrt(1 - (2 * pow(B1l, 2)));
        point[0]=B1l;
        point[1]=B1l;
        point[2]=B1m;
        P.x=point[0];P.y=point[1];P.z=point[2];
        for (int f = 0; f < 3; f++) {
            for (int g = f + 1; g < 3; g++) {
                point[f] = B1l;
                point[g] = B1l;
                int h = 3 - f - g;
                point[h] = B1m;
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        for (int k = 0; k < 2; k++) {
                            point[2] *= -1;
                            if (z==0) { P.x=point[0];P.y=point[1];P.z=point[2]; P.weight=b1w;
                                        grid_points.push_back(P);}
                            if (z==1) { P.x=point[0];P.y=point[1];P.z=point[2]; P.weight=b2w;
                                        grid_points.push_back(P);}
                            if (z==2) { P.x=point[0];P.y=point[1];P.z=point[2];P.weight=b3w;
                                        grid_points.push_back(P);}
                        }
                        point[1] *= -1;
                    }
                    point[0] *= -1;
                }
            }
        }
    }
    double Cp1=C1;
    double Cq1=sqrt(1-pow(Cp1,2));
    point[0]=0;
    point[1]=0;
    point[2]=0;
    for(int k=0;k<2;k++){
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 3; j++) {
                int h=3-i-j;
                point[i] = Cp1;
                point[j] = Cq1;
                point[h]=0;
                P.x=point[0];P.y=point[1];P.z=point[2];P.weight=c1w;
                grid_points.push_back(P);
                point[j] *= -1;
                point[h]=0;
                P.x=point[0];P.y=point[1];P.z=point[2]; P.weight=c1w;
                grid_points.push_back(P);
                point[i] *= -1;
                point[h]=0;
                P.x=point[0];P.y=point[1];P.z=point[2];P.weight=c1w;
                grid_points.push_back(P);
                point[j] = Cq1;
                point[h]=0;
                P.x=point[0];P.y=point[1];P.z=point[2];P.weight=c1w;
                grid_points.push_back(P);
            }
        }
        Cp1=Cq1;
        Cq1=C1;
    }

//write_out_grid_points();
}

//void Lebedev_grid::write_out_grid_points() {
//    for(size_t i=0;i<grid_points.size();i++){
//        cout<<i+1<<' '<<'('<<grid_points[i].x<<","<<grid_points[i].y<<","<<grid_points[i].z<<")"<<endl;
//    }
//}

void Lebedev_grid::build_sphere(BasisSet &bs, Lebedev_grid &g) {
    vector<Point> centers;
    vector<double> radiis;
    Point P;
    for (Atom a: bs.atoms){
        centers.push_back(a.centrum);
        radiis.push_back(a.get_VdW_radii());
    }
//    for(double k:radiis){cout<<k<<endl;}
    int i=0;
    for(Point &pc:centers){
        for(Point &pg:g.grid_points){
            P.x=pc.x+(pg.x*radiis[i]);
            P.y=pc.y+(pg.y*radiis[i]);
            P.z=pc.z+(pg.z*radiis[i]);
            P.weight=pg.weight;
            P.no_of_atom=i;
            bool keep=true;
            for(int neighbor_idx = 0; neighbor_idx < centers.size(); ++neighbor_idx) {
                if (neighbor_idx == i) continue;  // skip self
                Point& tmp = centers[neighbor_idx];
                double dx = tmp.x - P.x;
                double dy = tmp.y - P.y;
                double dz = tmp.z - P.z;
                double dist2 = dx*dx + dy*dy + dz*dz;
                if (dist2 < radiis[neighbor_idx]*radiis[neighbor_idx]) {  // ← use neighbor's radius
                    keep = false;
                    break;
                }
            }
            if(keep){g.sphere_points.push_back(P);}
        }
        i+=1;
    }
//    int j=1;
//    for(Point tmp:g.sphere_points){
//        cout<<j;
//        tmp.write_point();
//        j+=1;
//    }

}

void Lebedev_grid::cal_tessera(Lebedev_grid &g, BasisSet bs) {
    Tessera t;
    for(Point p:g.sphere_points){
        t.r=p;
        t.area=4*M_PI*(bs.atoms[p.no_of_atom].get_VdW_radii()*bs.atoms[p.no_of_atom].get_VdW_radii())*p.weight;
        tesserae.push_back(t);
    }
    cout<<"Tessera points: "<<tesserae.size()<<endl;

}


