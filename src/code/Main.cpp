#include <iostream>
#include <string>

#include "../headers/Simulation.h"


const std::string file_xyz = "../data/examen_270226_xyz"; // positions
const std::string file_mc = "../data/examen_270226_mci"; // moments cinétiques

int main(int argc, char** argv){

    // affectation de la graine au générateur aléatoire (une seule fois pour tout le code)
    std::random_device rand;

    std::cout << "Programme du projet de Simulation Microscopique\n";

    // Création de la simulation
    Simulation* simu = new Simulation{};

    // Génération des données
    /* En entrée soit le fichier des positions soit les deux fichiers ! */
    int res = simu->data_gen(file_xyz);

    if(!res){ // Si la génération s'est bien faite

        // Lancement de la simulation
        res = simu->run();

    }
    
    // Terminaison
    delete simu;
    std::cout << "Fin du programme" << std::endl;
    return res;
}