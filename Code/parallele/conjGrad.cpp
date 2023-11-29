#include "conjGrad.h"
#include "main.h"
#include "prodMatVect_BC.h"

using namespace std;

vector<double> vectorSum(const vector<double>& A, double scalar, const vector<double>& B)
{
   // performs A + sB
   vector<double> result(A.size(), 0.0);
   for (int i = 0; i < A.size(); i++)
   {
      result[i] = A[i] + scalar * B[i];
   }
   return result;
}

double innerProduct(const vector<double>& U, const vector<double>& V)
{
   /* performs the inner product U•V */
   double product(0.0);
   for (int i = 0; i < U.size(); i++)
   {
      product += U[i] * V[i];
   }
   return product;
}

vector<double> conjugateGradient(const vector<double>& U0, const vector<double>& RightHandSide, const procData& procDF, double epsilon, int kmax)
{
   /* Calcule la solution par la méthode du gradient conjugué, à epsilon près */

   vector<double> z(procDF.locSize, 0.0), d(procDF.locSize, 0.0), r(procDF.locSize, 0.0), sol = U0;
   double alpha(0.0), beta(0.0), gamma(0.0), delta(0.0);
   int k(0);

   // r0 = b - AU0
   r = vectorSum(RightHandSide, -1.0, prodMatvect(U0, procDF));

   //direction initiale
   d = r;

   //norme du résidu initial
   gamma = innerProduct(r, r);
   MPI_Allreduce(MPI_IN_PLACE, &gamma, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   beta = sqrt(gamma);

   while ((beta > epsilon) and (k < kmax))
   {
      z = prodMatvect(d, procDF);

      // On a supprimé un produit scalaire en trop...
      delta = innerProduct(z, d);
      MPI_Allreduce(MPI_IN_PLACE, &delta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      alpha = gamma / delta;
      sol = vectorSum(sol, alpha, d);

      r = vectorSum(r, -alpha, z);
      beta = innerProduct(r, r);
      MPI_Allreduce(MPI_IN_PLACE, &beta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      gamma = beta / gamma;
      d = vectorSum(r, gamma, d);

      // ...grâce à cette ligne
      gamma = beta;
      beta = sqrt(gamma);

      k++;
   }
   if (k >= kmax)
   {
      printf("Tolérance non-atteinte : %g. \n", beta);
   }
   return sol;
}
