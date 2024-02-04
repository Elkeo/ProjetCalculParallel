#include "main.h"


void remplissageVariables(const std::string parametersFilePath, SpaceTimeDomain& dom);

/* Déclaration du terme source */

double f(const int i, const int j, const double t, const SpaceTimeDomain& dom);

/* Déclaration des fonctions de bord */

double g(const int i, const int j, const double t, const SpaceTimeDomain& dom);

double h(const int i, const int j, const double t, const SpaceTimeDomain& dom);

/* Déclaration de la fonction qui sauvegarde
la solution dans un fichier au temps t */

void saveSolution(const std::vector<double> U, const int timeIteration, const procData& proc, const SpaceTimeDomain& dom);
void saveErrorFile(const std::vector<double>& U, const int timeIteration, const procData& proc, const SpaceTimeDomain& dom);

/* Déclaration de la fonction qui créé un fichier script et qui sort un gif */

void createGnuplotScriptAndShowPlot(const SpaceTimeDomain& dom, const int timeIteration);