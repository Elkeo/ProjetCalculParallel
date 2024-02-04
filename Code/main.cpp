#include "fonctions.h"
#include "prodMatVect_BC.h"
#include "conjGrad.h"
#include "main.h"

using namespace std::chrono;
using namespace std;

/* Déclaration des fonctions locales */
void calculateSolutionStep(vector<double>& U, double t, procData& proc, const SpaceTimeDomain& dom);
void charge(int me, int n, int np, int& iBeg, int& iEnd);

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

   charge(proc.me, dom.Ny, proc.nbProc, proc.iBeg, proc.iEnd);

   // On étend le domaine de calcul en bas si on n'est pas le premier processeur
   proc.iBeg = proc.iBeg - dom.r * (proc.me != 0);
   // On étend le domaine de calcul en haut si on n'est pas le dernier processeur
   proc.iEnd = proc.iEnd + dom.r * (proc.me != proc.nbProc - 1);
   // On calcule le nombre total de noeuds du processeur
   proc.nbElem_y = proc.iEnd - proc.iBeg + 1;
   proc.neighborsToMe = { MPI_PROC_NULL, MPI_PROC_NULL };
   proc.neighborsToMe[0] = (proc.me == proc.nbProc - 1) ? MPI_PROC_NULL : proc.me + 1;
   proc.neighborsToMe[1] = (proc.me == 0) ? MPI_PROC_NULL : proc.me - 1;
   proc.stencilOver = vector<double>(3 * dom.Nx, 0.0);
   proc.stencilUnder = vector<double>(3 * dom.Nx, 0.0);

   /* Déclaration des variables */
   vector<double> U(dom.Nx * proc.nbElem_y, 1.0);
   double elapsedHumanTime(0.0), maxElapsedHumanTime(0.0);

   /* Sauvegarde de l'état initial */
   saveSolution(U, timeIteration, proc, dom);
   if (proc.me == 0)
   {
      std::string solutionFilePath = "solution" + std::to_string(dom.testCase) + "/solutionFile_" + std::to_string(timeIteration) + ".dat";
      for (int i = 0; i < proc.nbProc; i++)
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

      if (proc.me == 0)
      {
         std::string solutionFilePath = "solution" + std::to_string(dom.testCase) + "/solutionFile_" + std::to_string(timeIteration) + ".dat";
         for (int i = 0; i < proc.nbProc; i++)
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
      cout << "Temps de calcul total : " << maxElapsedHumanTime << endl;
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

void calculateSolutionStep(vector<double>& U, double t, procData& proc, const SpaceTimeDomain& dom)
{
   vector<double> U1 = U;
   vector<double> RHS = U1;
   vector<double> actualValuesUnder(3 * dom.Nx, 0.0), actualValuesOver(3 * dom.Nx, 0.0);
   bool convergedUnder(false), convergedOver(false);
   int k(0), kSchwarz(100);
   // Début de la boucle de Schwarz
   while ((convergedUnder == false and convergedOver == false) or k <= kSchwarz) {
      RHS = U1;
      // Stockage des valeurs actuelles des recouvrements en bas et en haut
      for (int i = 0; i < 3 * dom.Nx; i++)
      {
         actualValuesUnder.at(i) = U1.at(i + (dom.r - 1) * dom.Nx);
         actualValuesOver.at(i) = U1.at(i + (proc.iEnd - 3 - (dom.r - 1)) * dom.Nx);
      }

      if (actualValuesUnder - proc.stencilUnder > dom.epsilon)
      {
         // On envoie les recouvrements en bas
         MPI_Send(&actualValuesUnder[0], 3 * dom.Nx, MPI_DOUBLE, proc.neighborsToMe[0], 0, MPI_COMM_WORLD);
         // On reçoit les recouvrements du bas
         MPI_Recv(&proc.stencilUnder[0], 3 * dom.Nx, MPI_DOUBLE, proc.neighborsToMe[0], MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         if (proc.me == 0)
         {
            proc.stencilUnder = actualValuesUnder;
         }
         convergedUnder = false;
      }
      else
      {
         convergedUnder = true;
      }

      if (actualValuesOver - proc.stencilOver > dom.epsilon)
      {
         // On envoie les recouvrements en haut
         MPI_Send(&actualValuesOver[0], 3 * dom.Nx, MPI_DOUBLE, proc.neighborsToMe[1], 1, MPI_COMM_WORLD);
         // On reçoit les recouvrements du haut
         MPI_Recv(&proc.stencilOver[0], 3 * dom.Nx, MPI_DOUBLE, proc.neighborsToMe[1], MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         if (proc.me == proc.nbProc - 1)
         {
            proc.stencilOver = actualValuesOver;
         }
         convergedOver = false;
      }
      else
      {
         convergedOver = true;
      }

      calculateRightHandSide(RHS, t, proc, dom);
      U1 = conjugateGradient(RHS, proc, dom);
      k++;
   }
   U = U1;
}
