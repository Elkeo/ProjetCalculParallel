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



valarray<double> BiCGstab(const valarray<double>& RightHandSide, const procData& proc, const SpaceTimeDomain& dom, double epsilon, int kmax)
{
   /*Calcul de la solutiona avec le solveur Bi gradient conjugué stabilisé*/

   size_t N = RightHandSide.size();
   valarray<double> x(0.0, N), r(0.0, N), r0(0.0, N), p(0.0, N), v(0.0, N), s(0.0, N), t(0.0, N);
   double alpha(1.0), beta(0.0), omega(1.0), rho(1.0), rhoPrev(1.0);
   int k(0);

   // Calcule r0 = b - Ax0
   r = RightHandSide;
   r0 = r;
   // Initialise p, v, s, t
   p = r;
   v = r;
   s = r;
   t = r;

   // Norme du résidu initial
   rho = (r0 * r).sum();
   beta = sqrt(rho);

   while ((beta > epsilon) and (k < kmax))
   {
      rho = (r0 * r).sum();
      beta = (rho / rhoPrev) * (alpha / omega);
      p = r + beta * (p - omega * v);

      v = prodMatvect(p, proc, dom);

      
      alpha = rho / (r0 * v).sum();

      s = r - alpha * v;

      t = prodMatvect(s, proc, dom);

      omega = (t * s).sum() / (t * t).sum();

      // Mise à jour de la solution x, du résidu r
      x = x + alpha * p + omega * s;
      r = s - omega * t;


      rhoPrev = rho;
      beta = sqrt(rho);

      k++;
   }


   if (k >= kmax)
   {
      printf("Tolérance non-atteinte : %g.", beta);
   }

   return x;


}
