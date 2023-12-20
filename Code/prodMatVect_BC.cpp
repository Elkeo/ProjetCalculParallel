#include "fonctions.h"
#include "prodMatVect_BC.h"

valarray<double> prodMatvect(const valarray<double>& U, const procData proc, const SpaceTimeDomain& dom)
{
   size_t N = U.size();
   double alpha = 1 + dom.D * dom.dt * 2 * (1.0 / (dom.dx * dom.dx) + 1.0 / (dom.dy * dom.dy));
   double beta = -dom.D * dom.dt * (1.0 / (dom.dx * dom.dx));
   double gamma = -dom.D * dom.dt * (1.0 / (dom.dy * dom.dy));
   valarray<double> F(0.0, N);
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

   for (size_t I = 0; I <= N; I++)
   {
      // Indices locaux pour se déplacer dans la matrice
      int i = I % dom.Nx;
      int j = I / dom.Nx;

      if (j > 0)
      {
         F[I] += gamma * (U[I - dom.Nx]);
      }

      if (i > 0)
      {
         F[I] += beta * (U[I - 1]);
      }

      F[I] += alpha * U[I];

      if (i < dom.Nx - 1)
      {
         F[I] += beta * (U[I + 1]);
      }
      if (j < proc.nbElem_y - 1)
      {
         F[I] += gamma * (U[I + dom.Nx]);
      }
   }
   return F;
}


void calculateRightHandSide(valarray<double>& U, const double t, const procData& proc, const SpaceTimeDomain& dom)
{
   int I;
   double alpha = 1 + dom.D * dom.dt * 2 * (1.0 / (dom.dx * dom.dx) + 1.0 / (dom.dy * dom.dy));
   double beta = -dom.D * dom.dt * (1.0 / (dom.dx * dom.dx));
   double gamma = -dom.D * dom.dt * (1.0 / (dom.dy * dom.dy));

   for (int j = 0; j < proc.nbElem_y; j++)
   {
      for (int i = 0; i < dom.Nx; i++)
      {
         I = i + j * dom.Nx;
         U[I] += dom.dt * f(i + 1, j + 1, t + dom.dt, dom);

         if (j == 0)
         {
            U[I] -= gamma * g(i + 1, 0, t + dom.dt, dom);
         }
         if (j == proc.nbElem_y - 1)
         {
            U[I] -= gamma * g(i + 1, proc.nbElem_y + 1, t + dom.dt, dom);
         }
         if (i == 0)
         {
            U[I] -= beta * h(0, j + 1, t + dom.dt, dom);
         }
         if (i == dom.Nx - 1)
         {
            U[I] -= beta * h(dom.Nx + 1, j + 1, t + dom.dt, dom);
         }
      }
   }
}