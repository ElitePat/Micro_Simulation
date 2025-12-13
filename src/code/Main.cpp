#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "../headers/Particule.h"


#define N_particules_total 1000


const std::string file = "../data/particule";
int N_particules_local;
int sort_prog = 0;


// à mettre dans la classe qui contient les Particules:
// lecture de N particules dans un fichier (contenant N particules ou plus ...)
void lireP(const std::string filepath, std::vector<Particule>& listP){

    // ouvrir un fichier en mode input
    std::ifstream fs(filepath);
    if(!fs.is_open()){ // rapport de l'erreur
        std::cerr << "ERREUR: le fichier n'a pas pu être ouvert\n";
        sort_prog = 1;
        return;
    }

    // varaibles temporaires
    std::string ligne, cell;
    double tmpx, tmpy, tmpz;
    int i = 0;

    // Pour ignorer la première ligne
    std::getline(fs,ligne);
    // Tant que le fichier n'est pas à la fin
    while(std::getline(fs,ligne)){

        std::stringstream sligne(ligne);
        
        // Ce n'est peut-être pas comme ça qu'on fait mais au moins ça marche :)
        i = 0;
        while(std::getline(sligne, cell, ' ')){
            try{
                if(stoi(cell) == 2){
                    continue;
                }
                switch(i){
                    case 0:
                        tmpx = stod(cell);
                        break;
                    case 1:
                        tmpy = stod(cell);
                        break;
                    case 2:
                        tmpz = stod(cell);
                        break;
                    default:
                        std::cerr << "Attention: cas inconnu abordé !";
                }
            }catch(std::invalid_argument const& ex){
                continue;
            }
            ++i;
        }        

        // creation nouvelle particule
        listP.push_back(Particule{tmpx,tmpy,tmpz});
    }

    // on ferme le fichier en sortie
    fs.close();
}


int main(int argc, char** argv){

    std::cout << "Projet de Simulation Microscopique\n";

    std::vector<Particule> list_particules;
    list_particules.reserve(N_particules_total);

    // Lecture des particules ...
    lireP(file,list_particules);

    // Affichage de ce qu'on a lu (debug)
    ///*
    for(Particule p : list_particules){
        p.afficheParticule();
    }
    //*/

    std::cout << "Fin du programme" << std::endl;
    return sort_prog;
}