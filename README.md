# Modèle de prédation cyclique "Pierre Papier Ciseaux"

## Résumé

Ce travail porte sur un modèle de prédation cyclique de type "pierre papier ciseaux". L'aspect étudié concerne l'extinction des espèces en compétition, qui se produisent tôt ou tard avec un nombre fini d'individus. Le temps d'extinction est vu comme un temps de première sortie d'une région de l'espace de phases et calculé en utilisant une description en équation de Fokker-Planck. 

Le rapport écrit de ce travail est 'TP/TP.pdf'

![alt text](TP/figures/illustration.png?raw=true "Une trajectoire dans l'espace des phases")

## Contenu du code source

Le dossier 'code' contient les fichiers suivants:

| Fichier | Contenu |
|---------|---------|
| **FP\_escape\_problem.h** | Calcul du temps d'échappement d'un intervalle à partir de l'équation de Fokker-Planck d'un processus 1D |
| **Integration\_D\_rho.h** | Calcul numérique du coefficient de diffusion effectif D(rho) |
| **Interpolation.h** | Objet encapsulant l'interpolation d'une fonction à une variable |
| **graph\_D\_rho.cpp** | Code source pour produire le graphe de D(rho) |
| **graph\_extc\_time.cpp** | Code source pour produire le graphe du temps d'extinction moyen en fonction de rho |
| **phase\_portrait.cpp** | Code source pour produire le portrait de phase |
| **gillespie\_trajectory.cpp** | Simulation d'une trajectoire |
| **gillespie\_extc\_time.cpp** | Calcul du temps moyen d'extinction à partir des simulations |


## Installation

Le code servant à reproduire les données numériques présentées dans le rapport se trouve dans le dossier "code". Il a été développé sous Ubuntu 18.04 avec les programmes suivants

> g++ 7.5.0 <br>
> gnuplot 5.2 <br>
> pdflatex (from texlive-latex-base) <br>

Les exécutables peuvent tous être compilés en exécutant le script 'code/compile_everything.sh'

> cd code <br>
> ./compile\_everything.sh <br>

## Reproduction des données numériques

Portraits de phases:

> cd code/portrait <br>
> ./run\_portrait.sh <br>
> ./run\_simulation.sh <br>
> ./plot.sh <br>

Graphique D(rho):

> cd code <br>
> ./graph\_D\_rho <br>
> cd D\_rho <br>

Graphique du temps d'extinction moyen (théorie + simulation) (prend environ 20min):

> cd code/comparaison\_extc <br>
> ./run\_theory.sh <br>
> ./run\_simulations.sh <br>
> ./plot.sh <br>

Graphiques du temps d'extinction moyen pour plusieurs N (prend environ 1h):

> cd code/quality\_theo\_with\_N <br>
> ./run\_several\_N.sh



