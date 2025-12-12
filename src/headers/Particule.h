#ifndef PARTICLE_H
#define PARTICLE_H

#include <string>
#include <iostream>

#define N_particules_total 1000


// définition d'une particule
class Particule{
private:
    // Position de la particule
    double x;
    double y;
    double z;

public:
    // Constructeur
    Particule(double xc, double yc, double zc);
    // Destructeur
    ~Particule() = default;

    // Printeur pour debug
    void afficheParticule();

    // Getters
    double coorx();
    double coory();
    double coorz();

    // Setter
    void update_coor(double xx, double yy, double zz);

    // lecture de N particules dans un fichier (contenant N particules ou plus ...)


};


#endif // PARTICLE_H