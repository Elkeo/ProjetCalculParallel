#include "conjGrad.h"

using namespace std;

vector<double> conjugateGradient(const vector<double>& RightHandSide, const procData& proc, const SpaceTimeDomain& dom, double epsilon, int kmax)
{
   /* Calcule la solution par la méthode du gradient conjugué, à epsilon près */

   int N = RightHandSide.size();
   vector<double> x(N, 0.0), z(N, 0.0), d(N, 0.0), r(N, 0.0);
   double alpha(0.0), beta(0.0), gamma(0.0);
   int k(0);

   // r0 = b - Ax0
   r = RightHandSide;

   //direction initiale
   d = r;

   //norme du résidu initial
   gamma = sum(r * r);
   beta = sqrt(gamma);

   while ((beta > epsilon) and (k < kmax))
   {
      z = prodMatvect(d, proc, dom);

      // On a supprimé un produit scalaire en trop...
      alpha = gamma / sum(z * d);
      x = x + alpha * d;

      r = r - alpha * z;
      beta = sum(r * r);

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

vector<double> BiCGstab(const vector<double>& RightHandSide, const procData& proc, const SpaceTimeDomain& dom, double epsilon, int kmax)
{
   /*Calcul de la solutiona avec le solveur Bi gradient conjugué stabilisé*/

   size_t N = RightHandSide.size();
   vector<double> x(N, 0.0), r(N, 0.0), r0(N, 0.0), p(N, 0.0), v(N, 0.0), s(N, 0.0), t(N, 0.0);
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
   rho = sum(r0 * r);
   beta = sqrt(rho);

   while ((beta > epsilon) and (k < kmax))
   {
      rho = sum(r0 * r);
      beta = (rho / rhoPrev) * (alpha / omega);
      p = r + beta * (p - omega * v);

      v = prodMatvect(p, proc, dom);


      alpha = rho / sum(r0 * v);

      s = r - alpha * v;

      t = prodMatvect(s, proc, dom);

      omega = sum(t * s) / sum(t * t);

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
