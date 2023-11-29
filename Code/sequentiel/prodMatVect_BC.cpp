#include "fonctions.h"
#include "prodMatVect_BC.h"

using namespace std;

vector<double> prodMatvect(const vector<double>& U)
{
   double alpha = 1 + D * dt * 2 * (1.0 / (dx * dx) + 1.0 / (dy * dy));
   double beta = -D * dt * (1.0 / (dx * dx));
   double gamma = -D * dt * (1.0 / (dy * dy));

   vector<double> F(Nx * Ny, 0.0);
   int I;
   /*

   La matrice s'écrit :
   [alpha beta           gamma             ]
   |beta alpha beta         gamma          |
   |                                       |
   |gamma     beta alpha 0       gamma  |
   |                                       |
   |           gamma        beta alpha beta|
   [               gamma         beta alpha]

    */

   for (int I = 0; I <= Nx * Ny; I++)
   {
      // Indices locaux pour se déplacer dans la matrice
      int i = I % Nx;
      int j = I / Nx;

      I = i + j * Nx;
      if (j > 0)
      {
         F[I] += gamma * (U[I - Nx]);
      }

      if (i > 0)
      {
         F[I] += beta * (U[I - 1]);
      }

      F[I] += alpha * U[I];

      if (i < Nx - 1)
      {
         F[I] += beta * (U[I + 1]);
      }
      if (j < Ny - 1)
      {
         F[I] += gamma * (U[I + Nx]);
      }
   }
   return F;
}


vector<double> calculateRightHandSide(vector<double>& U, double t)
{
   int I;
   double alpha = 1 + D * dt * 2 * (1.0 / (dx * dx) + 1.0 / (dy * dy));
   double beta = -D * dt * (1.0 / (dx * dx));
   double gamma = -D * dt * (1.0 / (dy * dy));

   for (int j = 0; j < Ny; j++)
   {
      for (int i = 0; i < Nx; i++)
      {
         I = i + j * Nx;
         U[I] += dt * f(i + 1, j + 1, t + dt);

         if (j == 0)
         {
            U[I] -= gamma * g(i + 1, 0, t + dt);
         }
         if (j == Ny - 1)
         {
            U[I] -= gamma * g(i + 1, Ny + 1, t + dt);
         }
         if (i == 0)
         {
            U[I] -= beta * h(0, j + 1, t + dt);
         }
         if (i == Nx - 1)
         {
            U[I] -= beta * h(Nx + 1, j + 1, t + dt);
         }
      }
   }

   return U;
}