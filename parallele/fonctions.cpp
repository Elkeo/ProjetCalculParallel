#include "fonctions.h"

int Nx(0), Ny(0), testCase(0);
double Lx(0.0), Ly(0.0), D(0.0), dt(0.0), tmax(0.0), dx(0.0), dy(0.0);

std::string parametersLine;

/* Récupération des variables dans le fichier parameters.dat */

void remplissageVariables(std::string parametersFilePath)
{
   /* Cette fonction remplit les variables
   Nx, Ny, Lx, Ly, D, dt
   à partir du fichier parameters.dat et affecte dx, dy */

   // On ouvre le fichier
   std::ifstream parametersFile(parametersFilePath);

   if (not parametersFile.is_open())
   {
      printf("Le fichier de paramètres n'a pas pu être ouvert.");
      abort();
   }

   getline(parametersFile, parametersLine);

   parametersFile >> Nx;

   parametersFile >> Ny;

   parametersFile >> Lx;

   parametersFile >> Ly;

   parametersFile >> D;

   parametersFile >> dt;

   parametersFile >> tmax;

   parametersFile >> testCase;

   dx = Lx / (Nx + 1);

   dy = Ly / (Ny + 1);


   if ((testCase != 1) and
      (testCase != 2) and
      (testCase != 3))
   {
      printf("Attention à choisir un cas test");
      abort();
   }

   parametersFile.close();
}

/* Implémentation du terme source */

double f(int i, int j, double t)
{
   double x(i * dx), y(j * dy);
   switch (testCase)
   {
   case 1:
      return 2 * (y - y * y + x - x * x);
   case 2:
      return sin(x) + cos(y);
   case 3:
      return exp(-(pow((x - Lx / 2), 2) + pow((y - Ly / 2), 2))) * cos(M_PI * t / 2);
   default:
      abort();
   }
}

/* Implémentation des fonctions de bord */

double g(int i, int j, double t)
{
   double x(i * dx), y(j * dy);
   switch (testCase)
   {
   case 1:
      return 0.0;
   case 2:
      return sin(x) + cos(y);
   case 3:
      return 0.0;
   default:
      abort();
   }
}

double h(int i, int j, double t)
{
   double x(i * dx), y(j * dy);
   switch (testCase)
   {
   case 1:
      return 0.0;
   case 2:
      return sin(x) + cos(y);
   case 3:
      return 1.0;
   default:
      abort();
   }
}

/* Implémentation de la solution exacte */

double exactSolution(int i, int j)
{
   double x(i * dx), y(j * dy);
   switch (testCase)
   {
   case 1:
      return (1.0 / D) * x * (1 - x) * y * (1 - y);
   case 2:
      return (1.0 / D) * (sin(x) + cos(y));
   default:
      return 0.0;
      break;
   }
};

/* fonction qui sauvegarde la solution au temps n dans le fichier solutionFile_n.dat */

void saveSolution(std::vector<double>& U, int timeIteration)
{

   std::string solutionFilePath = "solution" + std::to_string(testCase) + "/solutionFile_" + std::to_string(timeIteration) + ".dat";
   std::ofstream solutionFile(solutionFilePath);

   if (not solutionFile.is_open())
   {
      printf("Le fichier de solution n'a pas pu être ouvert.");
      abort();
   }

   solutionFile << " x, \t y, \t U" << std::endl;

   // bords bas et haut
   for (int i = 0; i <= Nx + 1; i++)
   {
      if ((i == 0) || (i == Nx + 1))
      {
         solutionFile << i * dx << "\t" << 0.0 << "\t" << (g(i, 0, timeIteration * dt) + h(i, 0, timeIteration * dt)) / 2 << std::endl;
      }
      else
         solutionFile << i * dx << "\t" << 0.0 << "\t" << g(i, 0, timeIteration * dt) << std::endl;
   }

   // bords gauche et droite
   for (int j = 0; j < Ny; j++)
   {
      solutionFile << 0.0 << "\t" << (j + 1) * dy << "\t" << h(0, j + 1, timeIteration * dt) << std::endl;
      for (int i = 0; i < Nx; i++)
      {
         solutionFile << (i + 1) * dx << "\t" << (j + 1) * dy << "\t" << U[i + j * Nx] << std::endl;
      }
      solutionFile << Lx << "\t" << (j + 1) * dy << "\t" << h(Nx + 1, j + 1, timeIteration * dt) << std::endl;
   }

   // intérieur
   for (int i = 0; i <= Nx + 1; i++)
   {
      if ((i == 0) || (i == Nx + 1))
      {
         solutionFile << i * dx << "\t" << (Ny + 1) * dy << "\t" << (g(i, Ny + 1, timeIteration * dt) + h(i, Ny + 1, timeIteration * dt)) / 2 << std::endl;
      }
      else
         solutionFile << i * dx << "\t" << (Ny + 1) * dy << "\t" << g(i, Ny + 1, timeIteration * dt) << std::endl;
   }

   solutionFile.close();
};

void saveErrorFile(const std::vector<double>& U)
{
   /* fonction qui sauvegarde la solution au temps n dans le fichier solutionFile_n.dat */

   std::string errorFilePath = "solution" + std::to_string(testCase) + "/solutionError.dat";
   std::ofstream errorFile(errorFilePath, std::fstream::app);

   if (not errorFile.is_open())
   {
      printf("Le fichier de solution n'a pas pu être ouvert.");
      abort();
   }

   double error_L2norm(0.0);

   for (int I = 0; I < Nx * Ny; I++)
   {
      int i = I % Nx;
      int j = I / Nx;

      error_L2norm += pow(U[I] - exactSolution(i + 1, j + 1), 2);
   }
   error_L2norm = sqrt(std::max(dx, dy) * error_L2norm);
   errorFile << timeIteration << "\t" << error_L2norm << std::endl;

   errorFile.close();
};

void concatenateSolutionToRootProcessAndSave(const std::vector<double>& U, const procData& procDF)
{
   // Collecte des tailles de vecteur sur tous les processus
   std::vector<int> vectorSizes(procDF.nproc);
   MPI_Allgather(&procDF.locSize, 1, MPI_INT, vectorSizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

   // Calcul des décalages pour la collecte des vecteurs
   std::vector<int> displacements(procDF.nproc);
   int totalSize = 0;
   for (int i = 0; i < procDF.nproc; i++) {
      displacements[i] = totalSize;
      totalSize += vectorSizes[i];
   }

   // Allocation de mémoire pour les vecteurs concaténés sur le processus racine
   std::vector<double> concatenatedVector(totalSize);

   // Collecte des vecteurs locaux sur le processus racine
   MPI_Gatherv(U.data(), procDF.locSize, MPI_DOUBLE, concatenatedVector.data(),
      vectorSizes.data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

   // Concaténation des vecteurs sur le processus racine 0
   if (procDF.me == 0) {
      std::vector<double> finalVector;
      for (int i = 0; i < procDF.nproc; i++) {
         finalVector.insert(finalVector.end(), concatenatedVector.begin() + displacements[i],
            concatenatedVector.begin() + displacements[i] + vectorSizes[i]);
      }

      saveSolution(finalVector, timeIteration);
      saveErrorFile(finalVector);
   }

};

void createGnuplotScriptAndShowPlot()
{
   /* On écrit un script Gnuplot à partir de la solution */
   std::string scriptPath = "Gnuplot_script/solution" + std::to_string(testCase) + ".gp";
   std::ofstream scriptFile(scriptPath);

   if (not scriptFile.is_open())
   {
      printf("Le fichier de script Gnuplot n'a pas pu être ouvert.");
      abort();
   }

   scriptFile << "set palette defined (-5 0 0 1, 0 1 1 1, 5 1 0 0)" << std::endl;
   scriptFile << "set terminal gif enhanced font Arial 30 animate delay 50 loop 1 optimize size 2000,2000" << std::endl;
   scriptFile << "set output \"solution" + std::to_string(testCase) + "/solution" + std::to_string(testCase) + ".gif\"" << std::endl;
   scriptFile << "set xr[0:1]" << std::endl;
   scriptFile << "set yr[0:1]" << std::endl;
   //+ std::to_string(dt) + "  "
   scriptFile << "dt=" << dt << std::endl;
   scriptFile << "do for [i=0:" + std::to_string(timeIteration) + "] {" << std::endl
      << "t=i*dt" << std::endl
      << "set title \"t = \".sprintf(\"%f\", t).\" s\"" << std::endl
      << "plot \"solution" + std::to_string(testCase) + "/solutionFile_\".i.\".dat\" u 1:2:3 with image" << std::endl
      << "}" << std::endl
      << "set output";

   scriptFile.close();
   /* On demande à Gnuplot de tracer les fichiers à partir du script précédent correspondant à chaque solution */
   std::string plotCommand = "gnuplot Gnuplot_script/solution" + std::to_string(testCase) + ".gp";
   system(plotCommand.c_str());
}
