# Fichier CMakeLists.txt à la racine du projet

# Spécification de la version minimale de CMake requise:
cmake_minimum_required(VERSION 3.28)

# enable modern C++ (en a-t-on besoin ?)
# set(CMAKE_CXX_STANDARD 20)

# Nom du projet
project(MICRO_SIMULATION CXX) # Indique que c'est un projet C++


# --- Définition des répertoires sources et headers ---

# Répertoire contenant les fichiers .cpp
# Basé sur la structure : HP_Mandelbrot/src/code/
set(CODE_DIR ${CMAKE_SOURCE_DIR}/src/code)

# Répertoire contenant les fichiers .h
# Basé sur la structure : HP_Mandelbrot/src/headers/
set(HEADERS_DIR ${CMAKE_SOURCE_DIR}/src/headers)


# --- Liste de TOUS les fichiers sources (.cpp) pour l'exécutable "Main.exe" ---

# Doit inclure TOUS les .cpp dont le code à besoin pour compiler et les lier.
set(MAIN_SRC
    ${CODE_DIR}/Main.cpp
    ${CODE_DIR}/Particule.cpp
    # Ajouter les fichiers sources ici au fur et à mésure de leur création.
)


# --- Configuration de compilation ---

# Ajouter les répertoires d'en-têtes
include_directories(${HEADERS_DIR}) # Ton répertoire d'en-têtes (où sont les headers)

# Options de compilation standard
add_compile_options(-std=c++20 -Wall -Wextra -pedantic -lm -O0) # Options de warning strictes + options maths


# --- Création des exécutables ---

# Ajout des executables
add_executable(Main.exe ${MAIN_SRC})


# --- Création des cibles (personalisés) ---

# Cible pour nettoyer le répertoire de build
add_custom_target(my_clean
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}
    COMMENT "Cleaning all contents in the build directory (${CMAKE_BINARY_DIR})"
)
# usage (go to the build directory): cmake --build . --target my_clean && cd ../ 


# --- Ajout des test automatiques: commande ctest ---

enable_testing()
add_test(Test Main.exe)