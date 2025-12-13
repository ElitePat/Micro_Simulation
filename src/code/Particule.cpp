#include "../headers/Particule.h"


// Constructeur
Particule::Particule(double xc, double yc, double zc) : x(xc), y(yc), z(zc){}

// Printeur pour debug
void Particule::afficheParticule(){
    std::cout << "x,y,z : " << x << "," << y << "," << z << "\n";
}

// Getters
double Particule::coorx(){return x;}
double Particule::coory(){return y;}
double Particule::coorz(){return z;}

// Setter
void Particule::update_coor(double xx, double yy, double zz){
    x = xx;
    y = yy;
    z = zz;
}

// lecture de N particules dans un fichier (contenant N particules ou plus ...)