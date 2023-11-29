#include "fonctions.h"
#include "prodMatVect_BC.h"

using namespace std;

vector<double> widenU(const vector<double>& U, const procData& procDF);

/* Définition des fonctions */
vector<double> prodMatvect(const vector<double>& U, const procData& procDF)
{
   double alpha = 1 + D * dt * 2 * (1.0 / (dx * dx) + 1.0 / (dy * dy));
   double beta = -D * dt * (1.0 / (dx * dx));
   double gamma = -D * dt * (1.0 / (dy * dy));

   // Solution à retourner et solution élargie
   vector<double> F(procDF.locSize, 0.0), wide_Uloc(procDF.locSize + 2 * Nx);

   /*

   La matrice s'écrit :
   [alpha beta           gamma             ]
   |beta alpha beta         gamma          |
   |                                       |
   |gamma         0 alpha 0        gamma   |
   |                                       |
   |           gamma        beta alpha beta|
   [               gamma         beta alpha]

    */

    // On élargit U
   wide_Uloc = widenU(U, procDF);

   for (int I = procDF.iBeg; I <= procDF.iEnd; I++)
   {
      // Indices locaux pour se déplacer dans la matrice
      int i = I % Nx;
      int j = I / Nx;

      // Décalage d'indice pour se déplacer dans F et U
      int J = I - procDF.iBeg;
      int K = J + Nx;

      /* Pour chaque cas, on ajoute avec U récupéré d'un autre processeur ou non */

      if (j > 0)
      {
         F[J] += gamma * wide_Uloc[K - Nx];
      }

      if (i > 0)
      {
         F[J] += beta * wide_Uloc[K - 1];
      }

      F[J] += alpha * wide_Uloc[K];

      if (i < Nx - 1)
      {
         F[J] += beta * wide_Uloc[K + 1];
      }
      if (j < Ny - 1)
      {
         F[J] += gamma * wide_Uloc[K + Nx];
      }
   }

   return F;
}

/* Fonction qui élargit la solution U grâce à des communications avec d'autres processeurs */

vector<double> widenU(const vector<double>& U, const procData& procDF)
{
   vector<double> U1(Nx, 0.0), U2(Nx, 0.0), wide_Uloc(procDF.locSize + 2 * Nx, 0.0);

   if (procDF.me == 0)
   {
      // On reçoit les Nx premiers termes de me+1
      MPI_Recv(&U1[0], Nx, MPI_DOUBLE, procDF.me + 1, procDF.tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // On envoie les Nx derniers termes à me+1
      MPI_Send(&U[procDF.locSize - Nx], Nx, MPI_DOUBLE, procDF.me + 1, procDF.tag, MPI_COMM_WORLD);
   }
   else if (procDF.me == procDF.nproc - 1)
   {
      // On envoie les Nx premiers termes à me-1
      MPI_Send(&U[0], Nx, MPI_DOUBLE, procDF.me - 1, procDF.tag, MPI_COMM_WORLD);

      // On reçoit les Nx derniers termes de me-1
      MPI_Recv(&U2[0], Nx, MPI_DOUBLE, procDF.me - 1, procDF.tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   }
   else
   {
      // Si le proc est pair
      if (procDF.me % 2 == 0)
      {
         // On reçoit les Nx premiers termes de me+1
         MPI_Recv(&U1[0], Nx, MPI_DOUBLE, procDF.me + 1, procDF.tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         // On reçoit les Nx derniers termes de me-1
         MPI_Recv(&U2[0], Nx, MPI_DOUBLE, procDF.me - 1, procDF.tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

         // On envoie les Nx premiers termes à me-1
         MPI_Send(&U[0], Nx, MPI_DOUBLE, procDF.me - 1, procDF.tag, MPI_COMM_WORLD);
         // On envoie les Nx derniers termes à me+1
         MPI_Send(&U[procDF.locSize - Nx], Nx, MPI_DOUBLE, procDF.me + 1, procDF.tag, MPI_COMM_WORLD);
      }
      // Si le proc est impair
      else
      {
         // On envoie les Nx premiers termes à me-1
         MPI_Send(&U[0], Nx, MPI_DOUBLE, procDF.me - 1, procDF.tag, MPI_COMM_WORLD);
         // On envoie les Nx derniers termes à me+1
         MPI_Send(&U[procDF.locSize - Nx], Nx, MPI_DOUBLE, procDF.me + 1, procDF.tag, MPI_COMM_WORLD);

         // On reçoit les Nx premiers termes de me+1
         MPI_Recv(&U1[0], Nx, MPI_DOUBLE, procDF.me + 1, procDF.tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
         // On reçoit les Nx derniers termes de me-1
         MPI_Recv(&U2[0], Nx, MPI_DOUBLE, procDF.me - 1, procDF.tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
   }

   // On assemble le vecteur total
   for (int i = 0; i < Nx; i++)
   {
      wide_Uloc[i] = U2[i];
   }
   for (int i = 0; i < procDF.locSize; i++)
   {
      wide_Uloc[i + Nx] = U[i];
   }
   for (int i = 0; i < Nx; i++)
   {
      wide_Uloc[i + procDF.locSize + Nx] = U1[i];
   }

   return wide_Uloc;
};

/* Définition du second membre */
vector<double> calculateRightHandSide(vector<double>& U, double t, const procData& procDF)
{
   double beta = -D * dt * (1.0 / (dx * dx));
   double gamma = -D * dt * (1.0 / (dy * dy));

   for (int I = procDF.iBeg; I <= procDF.iEnd; I++)
   {
      // Indices locaux pour se déplacer dans la matrice
      int i = I % Nx;
      int j = I / Nx;

      // Décalage d'indice pour se déplacer dans F et U
      int J = I - procDF.iBeg;

      U[J] += dt * f(i + 1, j + 1, t + dt);

      if (j == 0)
      {
         U[J] -= gamma * g(i + 1, 0, t + dt);
      }
      if (j == Ny - 1)
      {
         U[J] -= gamma * g(i + 1, Ny + 1, t + dt);
      }
      if (i == 0)
      {
         U[J] -= beta * h(0, j + 1, t + dt);
      }
      if (i == Nx - 1)
      {
         U[J] -= beta * h(Nx + 1, j + 1, t + dt);
      }
   }

   return U;
};