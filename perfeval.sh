#!/bin/bash
# perfeval.sh

# on se deplace dans le repertoire qui nous interesse
cd build-release

# on compile le projet
cmake --build .

# on modifie le format de sortie pour la fonction bash time
TIMEFORMAT=%R

# ici on mesure le temps de plusieurs executions, afin d'obtenir une moyenne.
# Attention: la focntion bash time renvoie le resultat dans stderr, pour le recuperer, on fait comme ceci:
t=$( { time for i in `seq 1 10`; do ./Main.exe > log.txt; done } 2>&1 )

# on oublie pas de calculer la moyenne et on affiche le resultat
echo "Nombre d'executions: 10"
echo "Temps de toutes les executions: $t"
echo "(à diviser par 10 pour obtenir la moyenne d'une execution)"