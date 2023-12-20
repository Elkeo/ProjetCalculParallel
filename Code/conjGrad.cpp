#include "conjGrad.h"

using namespace std;

valarray<double> conjugateGradient(const valarray<double>& RightHandSide, const procData& proc, const SpaceTimeDomain& dom, double epsilon, int kmax)
{
   /* Calcule la solution par la méthode du gradient conjugué, à epsilon près */

   size_t N = RightHandSide.size();
   valarray<double> x(0.0, N), z(0.0, N), d(0.0, N), r(0.0, N);
   double alpha(0.0), beta(0.0), gamma(0.0);
   int k(0);

   // r0 = b - Ax0
   r = RightHandSide;

   //direction initiale
   d = r;

   //norme du résidu initial
   gamma = (r * r).sum();
   beta = sqrt(gamma);

   while ((beta > epsilon) and (k < kmax))
   {
      z = prodMatvect(d, proc, dom);

      // On a supprimé un produit scalaire en trop...
      alpha = gamma / (z * d).sum();
      x = x + alpha * d;

      r = r - alpha * z;
      beta = (r * r).sum();

      gamma = beta / gamma;
      d = r + gamma * d;

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
