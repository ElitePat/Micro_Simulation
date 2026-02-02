#include "../headers/Simulation.h"

// Constructeur
Simulation::Simulation(){

    list_particules = new std::vector<Particule>(N_particules_total);
    list_forces = new std::vector<std::vector<double>>(N_particules_total, std::vector<double>(3)); // soit 3 valeurs pour chaque particule
    trans_vect = new std::vector<std::vector<double>>(N_sym, std::vector<double>(3));

}

// Destructeur
Simulation::~Simulation(){
    delete list_particules;
    delete list_forces;
    std::cout << "Mémoire de la Simulation liberé !\n";
}

// Affiche les information/état de la simulation au momment de l'appel
void Simulation::printInfo(){}


// initialisation des vecteurs de translation
void Simulation::trans_vect_init(){
    int n1=0, n0=0;
    double deb=-1*L; // debut des itérations
    double c1=deb, c2=deb, c0=deb;
    for(int i=0; i<N_sym; ++i){
        if(n0 > N_sym/3-1){
            c0 += L;
            n0 = 0;
        }
        trans_vect->at(i).at(0) = c0;
        ++n0;

        if(n1 > 2){
            c1 += L;
            if(c1 > L){
                c1 = -1*L;
            }
            n1 = 0;
        }
        trans_vect->at(i).at(1) = c1;
        ++n1;

        trans_vect->at(i).at(2) = c2;
        c2 += L;
        if(c2 > L){
            c2 = -1*L;
        }
    }
}


// lecture de N particules dans un fichier (contenant N particules ou plus ...)
int Simulation::lireP(const std::string filepath){

    // ouvrir un fichier en mode input
    std::ifstream fs(filepath);
    if(!fs.is_open()){ // rapport de l'erreur
        std::cerr << "ERREUR: le fichier n'a pas pu être ouvert\n";
        return 1;
    }

    // varaibles temporaires
    std::string ligne, cell;
    double tmpx, tmpy, tmpz;
    int i=0,p=0;

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
        list_particules->at(p).update_coor(tmpx,tmpy,tmpz);
        ++p;
    }

    // on ferme le fichier en sortie
    fs.close();
    return 0;
}



// carré de distances entre 2 partcules (r_ij)
double Simulation::carre_dist(Particule const& p1, Particule const& p2){
    double xx = p2.coorx()-p1.coorx();
    double yy = p2.coory()-p1.coory();
    double zz = p2.coorz()-p1.coorz();
    return xx*xx + yy*yy + zz*zz;
}

// calcul des forces agissant sur chacune des particules pour potentiel de Lennard Jones
void Simulation::energieLJ(){
    double fx, fy, fz, r_ij, for_all, temp1, temp2, tulj, r2, r6, r12, r_ij3, r_ij6;
    int i=0, j;
    //int n = 0; // debug line

    // Calcul de r^6 et de r^12
    r2 = r*r; // r^2
    r6 = r2*r2; // r^4
    r6 *= r2; // r^6
    r12 = r6*r6; // r^12

    for(std::vector<double> i_sym : *trans_vect){

        //std::cout << "Vecteur <" << i_sym.at(0) << ";"  << i_sym.at(1) << ";"  << i_sym.at(2) << ">\n"; // debug line

        for(Particule p1 : *list_particules){
            //std::cout << "\n======================: i = " << i << "\n"; // debug line

            fx=0, fy=0, fz=0, j=0;
            for(Particule p2: *list_particules){
                if(j != i){
                    r_ij = carre_dist(p1,p2);
                    //std::cout << "valeur distance etre " << i << " et " << j << " = " << r_ij << "\n";

                    // Calcul de r_ij^3 et de r_ij^6
                    r_ij3 = r_ij * r_ij;
                    r_ij3 *= r_ij;
                    r_ij6 = r_ij3 * r_ij3;

                    temp1 = r12 / r_ij6;
                    temp2 = r6 / r_ij3;
                    
                    // on sait que 3.0^2 = 9.0
                    //for_all = -48 * 0.2 * (std::pow((9.0/r_ij),6) - std::pow((9.0/r_ij),3));
                    for_all = -48 * epsilon * (temp1 - temp2);


                    fx += for_all * ((p1.coorx()-p2.coorx()) / r_ij);
                    fy += for_all * ((p1.coory()-p2.coory()) / r_ij);
                    fz += for_all * ((p1.coorz()-p2.coorz()) / r_ij);
                    //std::cout << "valeur force fx " << fx << "\n"; // debug line
                    //std::cout << " " << j << " "; // debug line
                
                    if(j > i){

                        tulj += temp1 - (2 * temp2);
                        tulj *= epsilon;
                    }
                }
                ++j;
            }

            /* Pour chaque vecteur on va reécrire la valeur de list_forces. En gros, on ne garde que
            la valeur des list_forces calculés pour le dernier vecteur */
            list_forces->at(i).at(0) = fx;
            list_forces->at(i).at(1) = fy;
            list_forces->at(i).at(2) = fz;
            
            // reset de j
            //std::cout << "j = "; // debug line
            j = 0;
            ++i;

        }
        
        tulj *= 4;
        ulj += tulj;

        // reset des variables
        tulj = 0;
        i = 0;
    }

    // moyenne sur les vecteurs ... BOF
    ulj /= N_sym;

    /* Non pertinent
    for(auto f: *list_forces){
        f.at(0) /= N_sym;
        f.at(1) /= N_sym;
        f.at(2) /= N_sym;
    }*/

}


// Lance la simulation
int Simulation::run(std::string const& filepath){

    std::cout << "==========Début d'execution de la Simulation==========\n\n";

    // initialisation des vecteurs de translation
    trans_vect_init();
    /* Debug code 
    for(auto t : *trans_vect){
        std::cout << t[0] << " " << t[1] << " "  << t[2] << "\n"; // debug line
    }
    //*/

    // Lecture des particules
    int test = lireP(filepath);
    if(test){
        std::cout << "Échec de la simulation !\n";
        return 1;
    }
    /* Debug code
    // Affichage de ce qu'on a lu (debug)
    for(Particule p : list_particules){
        p.afficheParticule();
    }
    //*/

    // Calcul de l'energie du système
    energieLJ();
    // On vérifie que la somme des forces agissant sur toutes les particules est nulle
    double sumfx=0, sumfy=0, sumfz=0;
    //std::cout << "Energie pour chaque point du système\nForces: fx fy fz\n";
    for(auto force : *list_forces){
        //std::cout << force[0] << " " << force[1] << " "  << force[2] << "\n"; // debug line
        sumfx += force[0];
        sumfy += force[1];
        sumfz += force[2];
    }
    //std::cout << "==================\n";
    std::cout << "Somme des forces pour x, y, z: " << sumfx << ", " << sumfy << ", "  << sumfz << "\n";
    std::cout << "ULJ = " << ulj << "\n";
    //std::cout << "Soit en somme -> " << (sumfx+sumfy) << "\n"; // debug line

    std::cout << "\n==========Fin de l'execution de la Simulation==========\n";
    return 0;
}