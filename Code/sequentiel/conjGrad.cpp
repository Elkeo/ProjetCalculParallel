#include "conjGrad.h"
#include "prodMatVect_BC.h"
#include "main.h"

using namespace std;

vector<double> sommeVect(vector<double> A, double scalar, vector<double> B)
{
   // calcule A + sB
   vector<double> result(A.size(), 0.0);
   for (int i = 0; i < A.size(); i++)
   {
      result[i] = A[i] + scalar * B[i];
   }

   return result;
}

double innerProduct(vector<double> U, vector<double> V)
{
   double product(0.0);
   for (int i = 0; i < U.size(); i++)
   {
      product += U[i] * V[i];
   }
   return product;
}

vector<double> conjugateGradient(vector<double> x, vector<double> RightHandSide, double epsilon, int kmax)
{
   /* Calcule la solution par la méthode du gradient conjugué, à epsilon près */

   vector<double> z(Nx * Ny, 0.0), d(Nx * Ny, 0.0), r(Nx * Ny, 0.0);
   double alpha(0.0), beta(0.0), gamma(0.0);
   int k(0);

   // r0 = b - Ax0
   r = sommeVect(RightHandSide, -1.0, prodMatvect(x));

   //direction initiale
   d = r;

   //norme du résidu initial
   gamma = innerProduct(r, r);
   beta = sqrt(gamma);


   while ((beta > epsilon) and (k < kmax))
   {
      z = prodMatvect(d);

      // On a supprimé un produit scalaire en trop...
      alpha = gamma / innerProduct(z, d);
      x = sommeVect(x, alpha, d);

      r = sommeVect(r, -alpha, z);
      beta = innerProduct(r, r);

      gamma = beta / gamma;
      d = sommeVect(r, gamma, d);

      // ...grâce à cette ligne
      gamma = beta;
      beta = sqrt(gamma);

      k++;
   }
   if (k > kmax)
   {
      printf("Tolérance non-atteinte : %g.", beta);
   }
   return x;
}
