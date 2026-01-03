#include <iostream>
#include <string>

#include "../headers/Simulation.h"


const std::string file = "../data/particule";

int main(int argc, char** argv){

    std::cout << "Projet de Simulation Microscopique\n";

    // Création de la simulation
    Simulation* simu = new Simulation{};

    // Lancement de la simulation
    int res = simu->run(file);

    // Terminaison
    delete simu;
    std::cout << "Fin du programme" << std::endl;
    return res;
}

/* Continuer debogage! */