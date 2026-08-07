#include <iostream>
#include <string>

#include "../headers/Simulation.h"


const std::string file_xyz = "../data/examen_270226_xyz";
const std::string file_mc = "../data/examen_270226_mci";

int main(int argc, char** argv){

    // affectation de la graine au générateur aléatoire (une seule fois pour tout le code)
    std::random_device rand;

    std::cout << "Programme du projet de Simulation Microscopique\n";

    // Création de la simulation
    Simulation* simu = new Simulation{};

    // Génération des données
    int res = simu->alea_gen(file_xyz);
    //int res = simu->file_gen(file_xyz,file_mc);

    if(!res){ // Si la génération s'est bien faite

        // Lancement de la simulation
        res = simu->run();

    }
    
    // Terminaison
    delete simu;
    std::cout << "Fin du programme" << std::endl;
    return res;
}