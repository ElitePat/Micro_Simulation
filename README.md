# CODE DE SIMULATION MICROSCOPIQUE

Projet académique de simulation d'un système microscopique en C++ du cours Introduction à la Simulation Microscopique du M2 CHPS de l'UVSQ (Paris-Saclay).

## Objectif

Développer un code de simulation microscopique basé sur un potentiel de Lennard-Jones pour étudier un fluide de particules homogène (i.e. toutes les particules sont identiques).
Code écrit en C++, pour un processeur Intel, dans WSL.

## Structure du projet

- Le dossier [data](data/) contient les differents jeux de données.
- Le dossier [src](src/) contient les fichiers source du programme.

## Procédure de compilation et d'exection

On crée un dossier `build` et on se place dedans:
```bash
mkdir build
cd build
```
Construire le projet:
```bash
cmake ..
```
Compiler le projet:
```bash
cmake --build .
```
(Si on veut voir la longueur de l'execution du programme et se rassurer que tout se passe bien)
```bash
ctest
```
Executer le projet
```bash
./Main.exe
```