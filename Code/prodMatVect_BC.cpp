#include "fonctions.h"
#include "prodMatVect_BC.h"

vector<double> prodMatvect(const vector<double>& U, const procData proc, const SpaceTimeDomain& dom)
{
   int N = U.size();
   double a = 1 + dom.D * dom.dt * 2 * (1.0 / (dom.dx * dom.dx) + 1.0 / (dom.dy * dom.dy));
   double b = -dom.D * dom.dt * (1.0 / (dom.dx * dom.dx));
   double c = -dom.D * dom.dt * (1.0 / (dom.dy * dom.dy));
   double a_mod = a + 2 * dom.D * dom.dt * (b / (a * dom.dy));
   vector<double> F(N, 0.0);
   int I(0);

   for (int j = 0; j < proc.nbElem_y; j++)
   {
      for (int i = 0; i < dom.Nx; i++)
      {
         I = i + j * dom.Nx;

         if (j > 0)
         {
            if (proc.me == 0 or proc.me == proc.nbProc - 1)
               F[I] += a_mod * U[I];
            else
               F[I] += a * U[I];
         }
         if (i > 0)
         {
            F[I] += b * U[I - 1];
         }

         // Diagonale de la matrice globale
         if (proc.me == 0 or proc.me == proc.nbProc - 1)
            F[I] += a_mod * U[I];
         else
            F[I] += a * U[I];

         if (i < dom.Nx - 1)
         {
            F[I] += b * U[I + 1];
         }
         if (j < proc.nbElem_y - 1)
         {
            /* code */
         }

      }
   }
   return F;
}


void calculateRightHandSide(vector<double>& U, const double t, const procData& proc, const SpaceTimeDomain& dom)
{
   int I;
   double b = -dom.D * dom.dt * (1.0 / (dom.dx * dom.dx));
   double c = -dom.D * dom.dt * (1.0 / (dom.dy * dom.dy));
   double lambda = 2 * dom.D * dom.beta * dom.dt / (dom.alpha * dom.dy);
   double nu = dom.D * dom.dt / (dom.dy * dom.dy);

   for (int j = 0; j < proc.nbElem_y; j++)
   {
      for (int i = 0; i < dom.Nx; i++)
      {
         I = i + j * dom.Nx;
         U[I] += dom.dt * f(i + 1, j + 1, t + dom.dt, dom);

         if (j == 0)
         {
            if (proc.me != 0)
               U[I] += nu * (proc.stencilUnder[i + 2 * dom.Nx] - proc.stencilUnder[i]) + lambda * proc.stencilUnder[i + 1 * dom.Nx];
            else
               U[I] -= c * g(i + 1, 0, t + dom.dt, dom);
         }
         if (j == proc.nbElem_y - 1)
         {
            if (proc.me != proc.nbProc - 1)
               U[I] += nu * (proc.stencilOver[i] - proc.stencilOver[i + 2 * dom.Nx]) + lambda * proc.stencilOver[i + 1 * dom.Nx];
            else
               U[I] -= c * g(i + 1, dom.Ny + 1, t + dom.dt, dom);
         }
         if (i == 0)
         {
            U[I] -= b * h(0, proc.iBeg + j + 1, t + dom.dt, dom);
         }
         if (i == dom.Nx - 1)
         {
            U[I] -= b * h(dom.Nx + 1, proc.iBeg + j + 1, t + dom.dt, dom);
         }
      }
   }
}