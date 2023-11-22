#include "fonctions.h"
#include "prodMatVect_BC.h"
#include "conjGrad.h"
#include "main.h"

using namespace std;


/* Déclaration des fonctions locales */
void calculateSolutionStep(vector<double>& U, vector<double>& U0, const procData& procDF);

double t(0.0);
int timeIteration(0);

/*
      FONCTION MAIN
*/

int main(int argc, char** argv)
{

   MPI_Init(&argc, &argv);

   /* Initialisation des variables */
   remplissageVariables("parameters.txt");

   procData procDF; //"Fichier" de données du processeur me
   procDF.n = Nx * Ny;

   MPI_Comm_rank(MPI_COMM_WORLD, &procDF.me);
   MPI_Comm_size(MPI_COMM_WORLD, &procDF.nproc);

   charge(procDF.me, procDF.n, procDF.nproc, procDF.iBeg, procDF.iEnd);
   procDF.locSize = procDF.iEnd - procDF.iBeg + 1;

   /* Déclaration des variables */
   vector<double> U(procDF.locSize, 1.0), U0(procDF.locSize, 1.0);
   double elapsedTime(0.0), startTime(0.0), endTime(0.0);

   /* Sauvegarde de l'état initial */
   concatenateSolutionToRootProcessAndSave(U, procDF);

   /* Calcul de la solution et sortie */
   while (t < tmax)
   {
      startTime = MPI_Wtime();
      calculateSolutionStep(U, U0, procDF);
      endTime = MPI_Wtime();
      if (procDF.me == 0)
      {
         cout << "Itération : " << timeIteration << endl;
         elapsedTime += endTime - startTime;
      }
      concatenateSolutionToRootProcessAndSave(U, procDF);
   }
   if (procDF.me == 0)
   {
      cout << "Temps de calcul : " << elapsedTime << endl;
      createGnuplotScriptAndShowPlot();
   }

   MPI_Finalize();
   return 0;
}

/*
   FONCTIONS ANNEXES
*/

void charge(int me, int n, int np, int& iBeg, int& iEnd)
{
   // Attention il faut bien écrire (n/np) qui renvoie la partie entière de la division
   int r = n % np;
   if (me < r)
   {
      iBeg = me * ((n / np) + 1);
      iEnd = (me + 1) * ((n / np) + 1) - 1;
   }
   if (me >= r)
   {
      iBeg = me * (n / np) + r;
      iEnd = iBeg + (n / np) - 1;
   }
}

void calculateSolutionStep(vector<double>& U, vector<double>& U0, const procData& procDF)
{
   U0 = U;
   vector<double> RHS = calculateRightHandSide(U0, t, procDF);

   U = conjugateGradient(U0, RHS, procDF);

   t += dt;
   timeIteration++;
}
