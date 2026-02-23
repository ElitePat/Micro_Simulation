#include <iostream>
#include <string>

#include "../headers/Simulation.h"


const std::string file_xyz = "../data/examen_270226_xyz";
const std::string file_mc = "../data/examen_270226_mci";

int main(int argc, char** argv){

    std::cout << "Projet de Simulation Microscopique\n";

    // Création de la simulation
    Simulation* simu = new Simulation{};

    // Lancement de la simulation
    int res = simu->run(file_xyz,file_mc);

    // Terminaison
    delete simu;
    std::cout << "Fin du programme" << std::endl;
    return res;
}