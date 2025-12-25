#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>

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


// carré de distances entre 2 partcules (r_ij)
double carre_dist(Particule const& p1, Particule const& p2){
    return std::pow(p2.coorx()-p1.coorx(),2) + std::pow(p2.coory()-p1.coory(),2) + std::pow(p2.coorz()-p1.coorz(),2);
}

// calcul des forces agissant sur chacune des particules pour potentiel de Lennard Jones
void energieLJ(std::vector<Particule> const& lp, std::vector<std::vector<double>>& lf){
    double fx,fy,fz,r_ij,for_all;
    for(int i=0; i < N_particules_total; i++){
        fx=0, fy=0, fz=0;
        for(int j=0; j < N_particules_total; j++){
            if(i == j)
                continue;
            r_ij = carre_dist(lp.at(i),lp.at(j));
            //std::cout << "valeur distance etre " << i << " et " << j << " = " << r_ij << "\n";
            
            // on sait que 3.0^2 = 9.0
            //for_all = -48 * 0.2 * (std::pow((9.0/r_ij),6) - std::pow((9.0/r_ij),3));
            for_all = -48 * 0.2 * ((std::pow(3.0,12)/std::pow(r_ij,6)) - (std::pow(3.0,6)/std::pow(r_ij,3)));

            fx += for_all * ((lp.at(i).coorx()-lp.at(j).coorx()) / r_ij);
            fy += for_all * ((lp.at(i).coory()-lp.at(j).coory()) / r_ij);
            fz += for_all * ((lp.at(i).coorz()-lp.at(j).coorz()) / r_ij);
            //std::cout << "valeur force fx " << fx << "\n"; // debug line
        }
        lf.at(i).at(0) = fx;
        lf.at(i).at(1) = fy;
        lf.at(i).at(2) = fz;
    }
}


int main(int argc, char** argv){

    std::cout << "Projet de Simulation Microscopique\n";

    std::vector<Particule> list_particules;
    list_particules.reserve(N_particules_total);

    std::vector<std::vector<double>> list_forces(N_particules_total, std::vector<double>(3, 0));
    /* Debug code for initialisation of list_forces
    for(auto ok : list_forces)
        std::cout << ok[0] << ok[1]<< ok[2] << "\n";
    //*/
    //list_forces.reserve(N_particules_total);

    // Lecture des particules ...
    lireP(file,list_particules);
    //list_particules.at(0).afficheParticule(); // debug line
    //list_particules.at(1).afficheParticule(); // debug line

    // Affichage de ce qu'on a lu (debug)
    /*
    for(Particule p : list_particules){
        p.afficheParticule();
    }
    //*/

    // Calcul de l'energie du système
    energieLJ(list_particules,list_forces);
    //std::cout << "Energie pour chaque point dy système\n";
    /* Debug code
    std::cout << "Forces: fx fy fz\n";
    for(auto ok : list_forces)
        std::cout << ok[0] << " " << ok[1] << " "  << ok[2] << "\n";
    std::cout << "==================\n";
    ///*/

    // On vérifie que la somme des forces agissant sur toutes les particules est nulle
    double sumfx=0, sumfy=0, sumfz=0;
    for(auto force : list_forces){
        sumfx += force.at(0);
        sumfy += force.at(1);
        sumfz += force.at(2);
    }

    std::cout << "Somme des forces pour x, y, z: " << sumfx << ", " << sumfy << ", "  << sumfz << "\n";
    //std::cout << "Soit en somme -> " << (sumfx+sumfy) << "\n"; // debug line

    std::cout << "Fin du programme" << std::endl;
    return sort_prog;
}

/* Reprendre les calculs à la main */