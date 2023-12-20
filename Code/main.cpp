#include "fonctions.h"
#include "prodMatVect_BC.h"
#include "conjGrad.h"
#include "main.h"

using namespace std::chrono;
using namespace std;

/* Déclaration des fonctions locales */
void calculateSolutionStep(valarray<double>& U, double t, procData& proc, const SpaceTimeDomain& dom, int maxIterSchwarz = 10000, double tolSchwarz = 1e-8);
void charge(int me, int n, int np, int& iBeg, int& iEnd);
void communicationStencil(const valarray<double>& U, procData& proc, const SpaceTimeDomain& dom);

/*
      FONCTION MAIN
*/

int main(int argc, char* argv[])
{
   MPI_Init(&argc, &argv);

   /* Initialisation des variables */
   procData proc;                // Variables locales au processeur
   SpaceTimeDomain dom;          // Variables globales du problème
   remplissageVariables("parameters.txt", dom);

   double t(0.0);
   int timeIteration(0);

   MPI_Comm_rank(MPI_COMM_WORLD, &proc.me);
   MPI_Comm_size(MPI_COMM_WORLD, &proc.nbProc);

   /* Remplissage des informations sur le processeur */
   charge(proc.me, dom.Ny, proc.nbProc, proc.iBeg, proc.iEnd);

   // On étend le domaine de calcul en bas si on n'est pas le premier processeur
   proc.iBeg = proc.iBeg - dom.r * (proc.me != 0);
   // On étend le domaine de calcul en haut si on n'est pas le dernier processeur
   proc.iEnd = proc.iEnd + dom.r * (proc.me != proc.nbProc - 1);
   // On calcule le nombre total de noeuds du processeur
   proc.nbElem_y = proc.iEnd - proc.iBeg + 1;

   proc.neighborsToMe[0] = (proc.me == proc.nbProc - 1) ? MPI_PROC_NULL : proc.me + 1;
   proc.neighborsToMe[1] = (proc.me == 0) ? MPI_PROC_NULL : proc.me - 1;

   proc.stencilLower.resize(3 * dom.Nx);
   proc.stencilUpper.resize(3 * dom.Nx);
   proc.stencilLower = 0.0;
   proc.stencilUpper = 0.0;

   /* Déclaration des variables */
   valarray<double> U(1.0, dom.Nx * proc.nbElem_y);
   double elapsedHumanTime(0.0), maxElapsedHumanTime(0.0);

   /* Sauvegarde de l'état initial */
   saveSolution(U, timeIteration, proc, dom);
   if (proc.me == 0)
   {
      std::string solutionFilePath = "solution" + std::to_string(dom.testCase) + "/solutionFile_" + std::to_string(timeIteration) + ".dat";
      for (size_t i = 0; i < proc.nbProc; i++)
      {
         std::string solutionFilePathProc = "solution" + std::to_string(dom.testCase) + "/solutionFile_" + std::to_string(timeIteration) + "_" + std::to_string(i) + ".dat";
         system(("cat " + solutionFilePathProc + " >> " + solutionFilePath).c_str());
         system(("rm " + solutionFilePathProc).c_str());
      }
      //system(("sort -nk1 " + solutionFilePath + " -o " + solutionFilePath).c_str());
   }

   /* Calcul de la solution et sortie */
   while (t < dom.tmax)
   {
      auto start = MPI_Wtime();
      calculateSolutionStep(U, t, proc, dom);
      auto end = MPI_Wtime();
      elapsedHumanTime += end - start;
      t += dom.dt;
      timeIteration++;
      saveSolution(U, timeIteration, proc, dom);
      saveErrorFile(U, timeIteration, proc, dom);


      /* Post-processing */
      if (proc.me == 0)
      {
         std::string solutionFilePath = "solution" + std::to_string(dom.testCase) + "/solutionFile_" + std::to_string(timeIteration) + ".dat";
         for (size_t i = 0; i < proc.nbProc; i++)
         {
            std::string solutionFilePathProc = "solution" + std::to_string(dom.testCase) + "/solutionFile_" + std::to_string(timeIteration) + "_" + std::to_string(i) + ".dat";
            system(("cat " + solutionFilePathProc + " >> " + solutionFilePath).c_str());
            system(("rm " + solutionFilePathProc).c_str());
         }
         //system(("sort -nk1 " + solutionFilePath + " -o " + solutionFilePath).c_str());
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }
   MPI_Reduce(&elapsedHumanTime, &maxElapsedHumanTime, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
   if (proc.me == 0)
   {
      cout << "Temps de calcul total : " << elapsedHumanTime << endl;
      createGnuplotScriptAndShowPlot(dom, timeIteration);
   }
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

void calculateSolutionStep(valarray<double>& U, double t, procData& proc, const SpaceTimeDomain& dom, int maxIterSchwarz, double tolSchwarz)
{
   int iterSchwarz(0);
   double error = 1 + tolSchwarz;
   /* Boucle d'itérations de Schwarz */
   while (iterSchwarz <= maxIterSchwarz && error > tolSchwarz)
   {
      communicationStencil(U, proc, dom);
      valarray<double> RHS = U;
      calculateRightHandSide(RHS, t, proc, dom);
      U = conjugateGradient(RHS, proc, dom);
      MPI_Barrier(MPI_COMM_WORLD);
      iterSchwarz++;
      /* FAIRE UN CHOIX SUR LA CONDITION D'ARRÊT
      error = ? */

   }

}

void communicationStencil(const valarray<double>& U, procData& proc, const SpaceTimeDomain& dom)
{
   // On envoie les recouvrements en bas
   MPI_Send(&U[0], 3 * dom.Nx, MPI_DOUBLE, proc.neighborsToMe[proc.Lower], proc.tag, MPI_COMM_WORLD);

   // On envoie les recouvrements en haut
   MPI_Send(&U[((proc.nbElem_y - 3)) * dom.Nx - 1], 3 * dom.Nx, MPI_DOUBLE, proc.neighborsToMe[proc.Upper], proc.tag, MPI_COMM_WORLD);

   // On reçoit les recouvrements du bas
   MPI_Recv(&proc.stencilLower[0], 3 * dom.Nx, MPI_DOUBLE, proc.neighborsToMe[proc.Lower], proc.tag, MPI_COMM_WORLD, &proc.status);

   // On reçoit les recouvrements du haut
   MPI_Recv(&proc.stencilUpper[0], 3 * dom.Nx, MPI_DOUBLE, proc.neighborsToMe[proc.Upper], proc.tag, MPI_COMM_WORLD, &proc.status);
};