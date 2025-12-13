#include "../headers/Particule.h"


// Constructeur
Particule::Particule(double xc, double yc, double zc) : x(xc), y(yc), z(zc){}

// Printeur pour debug
void Particule::afficheParticule(){
    std::cout << "x,y,z : " << x << "," << y << "," << z << "\n";
}

// Getters
double Particule::coorx() const{return x;}
double Particule::coory() const{return y;}
double Particule::coorz() const{return z;}

// Setter
void Particule::update_coor(double xx, double yy, double zz){
    x = xx;
    y = yy;
    z = zz;
}
