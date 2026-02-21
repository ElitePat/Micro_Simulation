# CODE DE SIMULATION MICROSCOPIQUE

Projet académique de simulation d'un système microscopique en C++ du cours Introduction à la Simulation Microscopique du M2 CHPS de l'UVSQ (Paris-Saclay).

## Objectif

Développer un code de simulation microscopique basé sur un potentiel de Lennard-Jones pour étudier un fluide de particules homogène (i.e. toutes les particules sont identiques).
Code écrit en C++, pour un processeur Intel, dans WSL.

## Structure du projet

- Le dossier [data](data/) contient les differents jeux de données.
- Le dossier [src](src/) contient les fichiers source du programme.

## Procédure de compilation et d'exection

Pour le mode DEBUG on se place dans le dossier:
```bash
mkdir build-debug # si le dossier n'existe pas on le crée avant
cd build-debug
```
Pour construire le projet:
```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

Pour le mode RELEASE on se place dans le dossier:
```bash
mkdir build-release # si le dossier n'existe pas on le crée avant
cd build-release
```
Pour construire le projet:
```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
```

Les prochaines étapes sont les mêmes pour les deux modes

Compiler le projet:
```bash
cmake --build .
```
(Si on veut voir le temps d'execution du programme et se rassurer que tout se passe bien)
```bash
ctest
```
Executer le projet
```bash
./Main.exe
```
Pour le mode DEBUG il est préférable de rediriger l'output vers un fichier de sortie à part pour eviter de remplir le terminal. Ainsi mieux vaut executer ceci
```bash
./Main.exe > ../data/log.txt
```
Pour nettoyer le dossier build (correctement) faire:
```bash
cd build
cmake --build . --target my_clean && cd ../
```