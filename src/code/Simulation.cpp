#include "../headers/Simulation.h"

// Constructeur
Simulation::Simulation(){

    /* L'initialisation des listes ci-dessous ne fait que remplir les vecteurs avec des zéros! */

    // initialisation des particules du système
    list_particules = new std::vector<Particule>(N_particules_total);
    // initialisation des positions antérieures des particules du système
    list_particules_prec = new std::vector<Particule>(N_particules_total);

    // pour chaque particule on represente les 3 forces sur les axes x, y, et z soit 3 valeurs
    list_forces = new std::vector<std::vector<double>>(N_particules_total, std::vector<double>(3));
    // même logique pour les vitesses de chaque particule
    list_v = new std::vector<std::vector<double>>(N_particules_total, std::vector<double>(3));
    list_v_prec = new std::vector<std::vector<double>>(N_particules_total, std::vector<double>(3));
    // même logique pour les moments cinétiques de chaque particule
    list_mc = new std::vector<std::vector<double>>(N_particules_total, std::vector<double>(3));

    // initialisation des vecteurs de translation, pour les conditions périodiques
    trans_vect = new std::vector<std::vector<double>>(N_sym, std::vector<double>(3));

}

// Destructeur
Simulation::~Simulation(){
    delete list_particules;
    delete list_forces;
    delete list_particules_prec;
    delete list_v;
    delete list_v_prec;
    delete list_mc;
    delete trans_vect;
    #ifndef NDEBUG // Debug code"
    std::cout << "Mémoire de la Simulation liberé !\n";
    #endif
}

// Affiche les information/état de la simulation au momment de l'appel
void Simulation::printInfo(){ // debug function
    /*
    #ifndef NDEBUG
    
    std::cout << "Vecteur -> conditions périodiques:\n";
    for(auto t : *trans_vect){
        //std::cout << t[0] << " " << t[1] << " "  << t[2] << "\n"; // debug line
    }

    std::cout << "\nAffichage des points lus dans le fichier: \npoints en  x,y,z\n";
    for(Particule p : *list_particules){
        //p.afficheParticule();
    }

    std::cout << "\nList Particules\n";
    for(auto force : *list_particules){
        //std::cout << force.coorx() << "," << force.coory() << ","  << force.coorx() << "\n";
    }

    std::cout << "\nEnergie pour chaque Particule, forces fx fy fz:\n";
    for(auto force : *list_forces){
        //std::cout << force[0] << "," << force[1] << ","  << force[2] << "\n";
    }
    
    #endif
    */

    // Variables locales
    double sumfx=0, sumfy=0, sumfz=0;

    // On vérifie que la somme des forces agissant sur toutes les particules est nulle
    for(auto force : *list_forces){
        sumfx += force[0];
        sumfy += force[1];
        sumfz += force[2];
    }
    std::cout << "\nSomme des forces pour x, y, z: " << sumfx << ", " << sumfy << ", "  << sumfz << "\n";
    std::cout << "ULJ = " << ulj << "\n";
    //std::cout << "Soit en somme -> " << (sumfx+sumfy) << "\n"; // debug line
    std::cout << "ENERGIE_CINETIQUE = " << ec << " et Température cinétique = " << tc << "\n";

}



// initialisation des vecteurs de translation
void Simulation::trans_vect_init(){
    /* Modifier cette fonction si j'ai le temps:
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
    /* On introduit certaines variables intermédiares pour attenuer l'accumulation d'erreurs
    dues à la sommation de termes d'ordre de grandeurs très differents! */
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

            j = 0;
            for(Particule p2: *list_particules){
                if(j != i){
                    r_ij = carre_dist(p1,p2,i_sym);
                    //std::cout << "valeur distance etre " << i << " et " << j << " = " << r_ij << "\n"; // debug line

                    if(r_ij < RCUT){ // application du rayon de troncature

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
                        //std::cout << "valeur force fy " << list_forces->at(i).at(1) << "\n"; // debug line
                        //std::cout << " " << j << " "; // debug line
                    
                        // calcul du terme de Leonard-Jones
                        if(j > i){
                            ulj += temp1 - (2 * temp2);
                        }
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

// distances entre 2 partcules
double Simulation::distX(Particule const& p1, Particule const& p2){
    double xx = p2.coorx()-p1.coorx();
    return xx*xx;
}
// distances entre 2 partcules
double Simulation::distY(Particule const& p1, Particule const& p2){
    double yy = p2.coory()-p1.coory();
    return yy*yy;
}
// distances entre 2 partcules
double Simulation::distZ(Particule const& p1, Particule const& p2){
    double zz = p2.coorz()-p1.coorz();
    return zz*zz;
}

// algorithme de velicty Verlet sans contrôle de la température
void Simulation::vverlet(){
    /* Pour chaque pas de temps nous permet de calculer la vitesse de chaque particule et d'actualiser sa position
    pour le pas de temps suivant! */

    // actualisation pour chaque particule du système
    for(int i=0; i<N_particules_total; ++i){
        // nouvelle positions
        list_particules->at(i).update_coor((-1 * list_v->at(i).at(0) * list_particules_prec->at(i).coorx()) / (1 - list_v->at(i).at(0)),
                                        (-1 * list_v->at(i).at(1) * list_particules_prec->at(i).coory()) / (1 - list_v->at(i).at(1)),
                                        (-1 * list_v->at(i).at(2) * list_particules_prec->at(i).coorz()) / (1 - list_v->at(i).at(2)));
    
        // nouvelles vitesses vx = sqrt(rx^2 + r_(-1)x^2)
        list_v->at(i).at(0) = sqrt(distX(list_particules->at(i),list_particules_prec->at(i)));
        list_v->at(i).at(1) = sqrt(distY(list_particules->at(i),list_particules_prec->at(i)));
        list_v->at(i).at(2) = sqrt(distZ(list_particules->at(i),list_particules_prec->at(i)));
    }
}

// calcul de ec et de tc
void Simulation::cinetic_ET(){
    /* La valeur de la masse étant constante on peut la sortir de la boucle de sommation pour grandement réduire
    le nombre de divisions */
    double temp = 0;
    for(auto p : *list_mc){
        temp += p.at(0) * p.at(0) + p.at(1) * p.at(1) + p.at(2) * p.at(2);
    }
    ec = temp / (2 * CONVERSION_FORCE * M);
    tc = ec / (N * CONSTANTE_R);
}


// Lance la simulation
int Simulation::run(std::string const& filepath){

    std::cout << "==========Début d'execution de la Simulation==========\n\n";

    // initialisation des vecteurs de translation
    trans_vect_init();

    // Lecture des particules dans le fichier particule
    int test = lireP(filepath);
    if(test){
        std::cout << "Échec de la simulation !\n";
        return 1;
    }  

    // pour chaque pas de temps
    for(t=0; t<T; ++t){

        #ifndef NDEBUG // Debug line
        std::cout << "\n-----------------Itération n°" << t << "-----------------\n";
        #endif

        // Remise a zéro de forces sur chaque particule
        for(int i=0;i<N_particules_total;++i){
            list_forces->at(i).at(0) = 0;
            list_forces->at(i).at(1) = 0;
            list_forces->at(i).at(2) = 0;
            //std::cout << list_forces->at(i).at(0) << "," << list_forces->at(i).at(1) << ","  << list_forces->at(i).at(2) << "\n"; // debug line
        }
        // Calcul de l'energie du système
        energieLJ();

        // sauvegarde de l'itéation précdente pour chaque particule du système
        for(int i=0; i<N_particules_total; ++i){
            // position des particules
            list_particules_prec->at(i).update_coor(list_particules->at(i).coorx(),
                                                    list_particules->at(i).coory(),
                                                    list_particules->at(i).coorz());
            // anciennes vitesses
            list_v_prec->at(i).at(0) = list_v->at(i).at(0);
            list_v_prec->at(i).at(1) = list_v->at(i).at(1);
            list_v_prec->at(i).at(2) = list_v->at(i).at(2);
        }

        // calcul des vitesses
        vverlet();

        // calcul du moment cinétique
        for(int i=0; i<N_particules_total; ++i){
            list_mc->at(i).at(0) = list_v->at(i).at(0) * M * CONVERSION_FORCE;
            list_mc->at(i).at(1) = list_v->at(i).at(1) * M * CONVERSION_FORCE;
            list_mc->at(i).at(2) = list_v->at(i).at(2) * M * CONVERSION_FORCE;
        }

        // calcul de l'energie et de la temperature cinetique
        cinetic_ET();

        // suite ...

        printInfo();
    }

    std::cout << "\n==========Fin de l'execution de la Simulation==========\n";
    return 0;
}