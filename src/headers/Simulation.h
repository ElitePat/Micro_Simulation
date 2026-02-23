#ifndef SIMULATION_H
#define SIMULATION_H


#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>

#include "Particule.h"


#define N_particules_total 1000
#define N_sym 27                        // vecteurs de translation possibles

// les constantes suivantes sont des doubles
#define L 31.09   
#define EPSILON 0.53
#define R 2.850
#define RCUT 10.0                       // rayon de coupure

// Pour Velocity Verlet
#define CONVERSION_FORCE 0.0001*4.186
#define M 18                            // masse d'une particule
#define CONSTANTE_R 8.31e-7
#define N 3*N_particules_total-3        // Nombre de degrées de liberté du système

//
#define T 10 // pas de temps maximal


class Simulation{
private:
    /* Variables */
    // nombre local de particules
    int N_particules_local;

    // valeur de l'energie (selon le terme de Leonard-Jones)
    double ulj = 0;
    
    // valeur de l'energie cinétique du système
    double ec = 0;

    // valeur de la température cinétique du système
    double tc = 0;

    // temps actuel de la simulation
    int t;

    // Liste des particules dans la simulation
    std::vector<Particule> *list_particules;
    // Liste des positions précdentes des particules
    std::vector<Particule> *list_particules_prec;

    // forces pour chaque particule
    std::vector<std::vector<double>> *list_forces;

    // vitesses de chaque particule dans les 3 axes
    std::vector<std::vector<double>> *list_v;
    // vitesses precedentes de chaque particule
    std::vector<std::vector<double>> *list_v_prec;

    // moment cinétique pour chaque particule
    std::vector<std::vector<double>> *list_mc;

    // vecteurs de translation
    std::vector<std::vector<double>> *trans_vect;


    /* Méthodes */
    // initialisation des vecteurs de translation
    void trans_vect_init();

    // lecture de N particules dans un fichier (contenant N particules ou plus ...)
    int lireP(const std::string filepath);

    // lecture de N moments cinétiques dans un fichier (contenant N particules ou plus ...)
    int lireM(const std::string filepath);

    // carré de distances entre 2 partcules
    double carre_dist(Particule const& p1, Particule const& p2, std::vector<double> const& vec);

    // calcul des forces agissant sur chacune des particules pour potentiel de Lennard Jones
    void energieLJ();

    // distances entre 2 partcules
    double distX(Particule const& p1, Particule const& p2);
    double distY(Particule const& p1, Particule const& p2);
    double distZ(Particule const& p1, Particule const& p2);

    // algorithme de velocity Verlet sans contrôle de la température
    void vverlet();

    // calcul de l'energie cinétique et de la température cinétique du système
    void cinetic_ET();

    // nombre moyen de particules situées à une distance inférieure à Rc
    double pproches(int const& rcut);

public:
    // Constructeur
    Simulation();
    // Destructeur
    ~Simulation();

    // Affiche les information/état de la simulation au momment de l'appel
    void printInfo();

    // Lancer la simulation
    int run(std::string const& filepath_xyz, std::string const& filepath_mc);
};


#endif // SIMULATION_H