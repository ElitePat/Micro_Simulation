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
#define N_sym 27 // vecteurs de translation possibles
#define L 42.0 // c'est un double !


class Simulation{
private:
    /* Variables */
    // nombre local de particules
    int N_particules_local;

    // Liste des particules dans la simulation
    std::vector<Particule> *list_particules;

    // forces pour chaque particule
    std::vector<std::vector<double>> *list_forces;

    // vecteurs de translation
    std::vector<std::vector<double>> *trans_vect;


    /* Méthodes */
    // initialisation des vecteurs de translation
    void trans_vect_init();

    // lecture de N particules dans un fichier (contenant N particules ou plus ...)
    int lireP(const std::string filepath);

    // carré de distances entre 2 partcules
    double carre_dist(Particule const& p1, Particule const& p2);

    // calcul des forces agissant sur chacune des particules pour potentiel de Lennard Jones
    void energieLJ();

public:
    // Constructeur
    Simulation();
    // Destructeur
    ~Simulation();

    // Affiche les information/état de la simulation au momment de l'appel
    void printInfo();

    // Lancer la simulation
    int run(std::string const& filepath);
};


#endif // SIMULATION_H