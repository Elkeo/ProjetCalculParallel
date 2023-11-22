#include "main.h"


void remplissageVariables(std::string parametersFilePath);

/* Déclaration du terme source */

double f(int i, int j, double t);

/* Déclaration des fonctions de bord */

double g(int i, int j, double t);

double h(int i, int j, double t);

double exactSolution(int i, int j);

/* Déclaration de la fonction qui sauvegarde
la solution dans un fichier au temps t */

void saveSolution(std::vector<double>& U);

void saveErrorFile(const std::vector<double>& U);

void concatenateSolutionToRootProcessAndSave(const std::vector<double>& U, const procData& procDF);

/* Déclaration de la fonction qui créé un fichier script et qui sort un gif */

void createGnuplotScriptAndShowPlot();