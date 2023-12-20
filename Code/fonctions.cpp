#include "fonctions.h"

std::string parametersLine;

/* Récupération des variables dans le fichier parameters.dat */

void remplissageVariables(const std::string parametersFilePath, SpaceTimeDomain& dom)
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
   std::cout << parametersLine << ", dx, dy" << std::endl;

   parametersFile >> dom.Nx;
   printf("%i, ", dom.Nx);

   parametersFile >> dom.Ny;
   printf("%i, ", dom.Ny);

   parametersFile >> dom.Lx;
   printf("%g, ", dom.Lx);

   parametersFile >> dom.Ly;
   printf("%g, ", dom.Ly);

   parametersFile >> dom.D;
   printf("%g, ", dom.D);

   parametersFile >> dom.dt;
   printf("%g, ", dom.dt);

   parametersFile >> dom.tmax;
   printf("%g, ", dom.tmax);

   parametersFile >> dom.testCase;
   printf("%i, ", dom.testCase);

   parametersFile >> dom.r;
   printf("%i, ", dom.r);

   dom.dx = dom.Lx / (dom.Nx + 1);
   printf("%g, ", dom.dx);

   dom.dy = dom.Ly / (dom.Ny + 1);
   printf("%g", dom.dy);

   printf("\n");

   if ((dom.testCase != 1) and
      (dom.testCase != 2) and
      (dom.testCase != 3))
   {
      printf("Attention à choisir un cas test");
      abort();
   }

   parametersFile.close();
}

/* Implémentation du terme source */

double f(const int i, const int j, const double t, const SpaceTimeDomain& dom)
{
   double x(i * dom.dx), y(j * dom.dy);
   switch (dom.testCase)
   {
   case 1:
      return 2 * (y - y * y + x - x * x);
   case 2:
      return sin(x) + cos(y);
   case 3:
      return exp(-(pow((x - dom.Lx / 2), 2) + pow((y - dom.Ly / 2), 2))) * cos(M_PI * t / 2);
   default:
      abort();
   }
}

/* Implémentation des fonctions de bord */

double g(const int i, const int j, const double t, const SpaceTimeDomain& dom)
{
   double x(i * dom.dx), y(j * dom.dy);
   switch (dom.testCase)
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

double h(const int i, const int j, const double t, const SpaceTimeDomain& dom)
{
   double x(i * dom.dx), y(j * dom.dy);
   switch (dom.testCase)
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

double exactSolution(const int i, const int j, const SpaceTimeDomain& dom)
{
   double x(i * dom.dx), y(j * dom.dy);
   switch (dom.testCase)
   {
   case 1:
      return (1.0 / dom.D) * x * (1 - x) * y * (1 - y);
   case 2:
      return (1.0 / dom.D) * (sin(x) + cos(y));
   default:
      return 0.0;
      break;
   }
};

void saveSolution(const std::valarray<double> U, const int timeIteration, const procData& proc, const SpaceTimeDomain& dom)
{
   /* fonction qui sauvegarde la solution au temps n dans le fichier solutionFile_n.dat */
   double t(timeIteration * dom.dt);
   std::string solutionFilePath = "solution" + std::to_string(dom.testCase) + "/solutionFile_" + std::to_string(timeIteration) + "_" + std::to_string(proc.me) + ".dat";
   std::ofstream solutionFile(solutionFilePath);

   if (not solutionFile.is_open())
   {
      printf("Le fichier de solution n'a pas pu être ouvert.");
      abort();
   }
   if (proc.me == 0)
   {
      solutionFile << " x, \t y, \t U" << std::endl;
   }


   for (int i = 0; i <= dom.Nx + 1; i++)
   {
      if ((i == 0) || (i == dom.Nx + 1))
      {
         solutionFile << i * dom.dx << "\t" << 0.0 << "\t" << (g(i, 0, t, dom) + h(i, 0, t, dom)) / 2 << std::endl;
      }
      else
         solutionFile << i * dom.dx << "\t" << 0.0 << "\t" << g(i, 0, t, dom) << std::endl;
   }

   for (int j = 0; j < dom.Ny; j++)
   {
      solutionFile << 0.0 << "\t" << (j + 1) * dom.dy << "\t" << h(0, j + 1, t, dom) << std::endl;
      for (int i = 0; i < dom.Nx; i++)
      {
         solutionFile << (i + 1) * dom.dx << "\t" << (j + 1) * dom.dy << "\t" << U[i + j * dom.Nx] << std::endl;
      }
      solutionFile << dom.Lx << "\t" << (j + 1) * dom.dy << "\t" << h(dom.Nx + 1, j + 1, t, dom) << std::endl;
   }

   for (int i = 0; i <= dom.Nx + 1; i++)
   {
      if ((i == 0) || (i == dom.Nx + 1))
      {
         solutionFile << i * dom.dx << "\t" << (dom.Ny + 1) * dom.dy << "\t" << (g(i, dom.Ny + 1, t, dom) + h(i, dom.Ny + 1, t, dom)) / 2 << std::endl;
      }
      else
         solutionFile << i * dom.dx << "\t" << (dom.Ny + 1) * dom.dy << "\t" << g(i, dom.Ny + 1, t, dom) << std::endl;
   }

   solutionFile.close();
};

void saveErrorFile(const std::valarray<double>& U, const int timeIteration, const procData& proc, const SpaceTimeDomain& dom)
{
   /* fonction qui sauvegarde la solution au temps n dans le fichier solutionFile_n.dat */

   std::string errorFilePath = "solution" + std::to_string(dom.testCase) + "/solutionError.dat";
   std::ofstream errorFile(errorFilePath, std::fstream::app);

   if (not errorFile.is_open())
   {
      printf("Le fichier de solution n'a pas pu être ouvert.");
      abort();
   }

   double error_L2norm(0.0);

   for (int I = 0; I < dom.Nx * dom.Ny; I++)
   {
      int i = I % dom.Nx;
      int j = I / dom.Nx;

      error_L2norm += pow(U[I] - exactSolution(i + 1, j + 1, dom), 2);
   }
   error_L2norm = sqrt(std::max(dom.dx, dom.dy) * error_L2norm);
   errorFile << timeIteration << "\t" << error_L2norm << std::endl;

   errorFile.close();
};

void createGnuplotScriptAndShowPlot(const SpaceTimeDomain& dom, const int timeIteration)
{
   /* On écrit un script Gnuplot à partir de la solution */
   std::string scriptPath = "Gnuplot_script/solution" + std::to_string(dom.testCase) + ".gp";
   std::ofstream scriptFile(scriptPath);

   if (not scriptFile.is_open())
   {
      printf("Le fichier de script Gnuplot n'a pas pu être ouvert.");
      abort();
   }

   scriptFile << "set palette defined (-5 0 0 1, 0 1 1 1, 5 1 0 0)" << std::endl;
   scriptFile << "set terminal gif enhanced font Arial 30 animate delay 50 loop 1 optimize size 1300,1000" << std::endl;
   scriptFile << "set output \"solution" + std::to_string(dom.testCase) + "/solution" + std::to_string(dom.testCase) + ".gif\"" << std::endl;
   scriptFile << "set xr[0:1]" << std::endl;
   scriptFile << "set yr[0:1]" << std::endl;
   //+ std::to_string(dt) + "  "
   scriptFile << "dt=" << dom.dt << std::endl;
   scriptFile << "do for [i=0:" + std::to_string(timeIteration) + "] {" << std::endl
      << "t=i*dt" << std::endl
      << "set title \"t = \".sprintf(\"%f\", t).\" s\"" << std::endl
      << "plot \"solution" + std::to_string(dom.testCase) + "/solutionFile_\".i.\".dat\" u 1:2:3 with image" << std::endl
      << "}" << std::endl
      << "set output";

   scriptFile.close();
   /* On demande à Gnuplot de tracer les fichiers à partir du script précédent correspondant à chaque solution */
   std::string plotCommand = "gnuplot Gnuplot_script/solution" + std::to_string(dom.testCase) + ".gp";
   system(plotCommand.c_str());
}
