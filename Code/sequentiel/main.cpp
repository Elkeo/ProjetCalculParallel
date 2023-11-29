#include "fonctions.h"
#include "prodMatVect_BC.h"
#include "conjGrad.h"
#include "main.h"
#include <chrono>
#include <string>

using namespace std::chrono;
using namespace std;

/* Déclaration des fonctions locales */
void calculateSolutionStep(vector<double>& U, vector<double>& U0);

double t(0.0);
int timeIteration(0);

/*
      FONCTION MAIN
*/

int main(int argc, char const* argv[])
{
   /* Initialisation des variables */
   remplissageVariables("parameters.txt");

   /* Déclaration des variables */
   vector<double> U(Nx * Ny, 1.0), U0(Nx * Ny, 1.0);
   double elapsedTime(0.0);

   /* Sauvegarde de l'état initial */
   saveSolution(U);

   /* Calcul de la solution et sortie */
   while (t < tmax)
   {
      auto start = high_resolution_clock::now();
      calculateSolutionStep(U, U0);
      auto end = high_resolution_clock::now();
      saveSolution(U);
      saveErrorFile(U);
      elapsedTime += duration_cast<microseconds>(end - start).count() * 1e-6;
   }
   cout << "Temps de calcul : " << elapsedTime << endl;
   createGnuplotScriptAndShowPlot();

   return 0;
}

/*
   FONCTIONS ANNEXES
*/

void calculateSolutionStep(vector<double>& U, vector<double>& U0)
{
   U0 = U;
   vector<double> RHS = calculateRightHandSide(U0, t);

   U = conjugateGradient(U0, RHS);

   t += dt;
   timeIteration++;
}
