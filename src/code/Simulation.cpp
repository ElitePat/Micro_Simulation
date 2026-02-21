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
    #ifndef NDEBUG // Debug code"
    std::cout << "Mémoire de la Simulation liberé !\n";
    #endif
}

// Affiche les information/état de la simulation au momment de l'appel
void Simulation::printInfo(){}



// initialisation des vecteurs de translation
void Simulation::trans_vect_init(){
/*  Modifier cette fonction si j'ai le temps:
    - n'accepter N_sym que si ça racine cubique est entière !
    - repenser alors comment doit être réalisé le parcours
    et l'affectation des valeurs dans trans_vect !
        - par exemple: n'affecter que des 0, 1 et -1
        puis tout multiplier par L !
*/
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
double Simulation::carre_dist(Particule const& p1, Particule const& p2, std::vector<double> const& vec){
    double xx = p2.coorx()-(p1.coorx()+vec[0]);
    double yy = p2.coory()-(p1.coory()+vec[1]);
    double zz = p2.coorz()-(p1.coorz()+vec[2]);
    return xx*xx + yy*yy + zz*zz;
}

// calcul des forces agissant sur chacune des particules pour potentiel de Lennard Jones
void Simulation::energieLJ(){

    // Déclaration des variables
    double r_ij, for_all, temp1, temp2, r_ij3, r_ij6;
    int i=0, j;
    //int n = 0; // debug line

    // Calcul de r^6 et de r^12
    double r2 = R*R;        // r^2
    double r6 = r2*r2;      // r^4
    r6 *= r2;               // r^6
    double r12 = r6*r6;     // r^12

    // Itération sur les conditions périodiques
    for(std::vector<double> i_sym : *trans_vect){
        //std::cout << "Vecteur <" << i_sym.at(0) << ";"  << i_sym.at(1) << ";"  << i_sym.at(2) << ">\n"; // debug line

        i = 0;
        for(Particule p1 : *list_particules){
            //std::cout << "\n======================: i = " << i << "\n"; // debug line

            j=0;
            for(Particule p2: *list_particules){
                if(j != i){
                    r_ij = carre_dist(p1,p2,i_sym);
                    //std::cout << "valeur distance etre " << i << " et " << j << " = " << r_ij << "\n"; // debug line

                    // Calcul de r_ij^3 et de r_ij^6
                    r_ij3 = r_ij * r_ij;    // r_ij^2
                    r_ij3 *= r_ij;          // r_ij^3
                    r_ij6 = r_ij3 * r_ij3;  // r_ij^6

                    temp1 = r12 / r_ij6;
                    temp2 = r6 / r_ij3;
                    
                    for_all = -48 * EPSILON * (temp1 - temp2);

                    list_forces->at(i).at(0) += for_all * ((p1.coorx()-p2.coorx()) / r_ij);
                    list_forces->at(i).at(1) += for_all * ((p1.coory()-p2.coory()) / r_ij);
                    list_forces->at(i).at(2) += for_all * ((p1.coorz()-p2.coorz()) / r_ij);
                    //std::cout << "valeur force fx " << fx << "\n"; // debug line
                    //std::cout << " " << j << " "; // debug line
                
                    // application de la troncature dans le calcul du terme de Leonard-Jones
                    if((r_ij < RCUT) && (j > i)){
                        ulj += temp1 - (2 * temp2);
                    }
                }
                ++j;
            }
            ++i;
        }
    }

    // valeur finale de ULJ, version tronqué
    ulj = (EPSILON * 2) * ulj;
}


// Lance la simulation
int Simulation::run(std::string const& filepath){

    std::cout << "==========Début d'execution de la Simulation==========\n\n";

    // initialisation des vecteurs de translation
    trans_vect_init();
    #ifndef NDEBUG // Debug code
    std::cout << "Vecteur -> conditions périodiques:\n";
    for(auto t : *trans_vect){
        std::cout << t[0] << " " << t[1] << " "  << t[2] << "\n"; // debug line
    }
    #endif

    // Lecture des particules
    int test = lireP(filepath);
    if(test){
        std::cout << "Échec de la simulation !\n";
        return 1;
    }
    #ifndef NDEBUG // Debug code
    std::cout << "\nAffichage des points lus dans le fichier: \npoints en  x,y,z\n";
    for(Particule p : *list_particules){
        p.afficheParticule();
    }
    #endif

    // Calcul de l'energie du système
    energieLJ();
    // On vérifie que la somme des forces agissant sur toutes les particules est nulle
    double sumfx=0, sumfy=0, sumfz=0;
    #ifndef NDEBUG // Debug line
    std::cout << "\nEnergie pour chaque point du système, forces fx fy fz:\n";
    #endif
    for(auto force : *list_forces){
        #ifndef NDEBUG // for Debug 
        std::cout << force[0] << "," << force[1] << ","  << force[2] << "\n"; // debug line
        #endif
        sumfx += force[0];
        sumfy += force[1];
        sumfz += force[2];
    }
    //std::cout << "==================\n";
    std::cout << "\nSomme des forces pour x, y, z: " << sumfx << ", " << sumfy << ", "  << sumfz << "\n";
    std::cout << "ULJ = " << ulj << "\n";
    //std::cout << "Soit en somme -> " << (sumfx+sumfy) << "\n"; // debug line

    std::cout << "\n==========Fin de l'execution de la Simulation==========\n";
    return 0;
}