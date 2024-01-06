#include "prodMatVect_BC.h"
#include "main.h"

std::valarray<double> conjugateGradient(const std::valarray<double>& RightHandSide, const procData& proc, const SpaceTimeDomain& dom, double epsilon = 1.0e-8, int kmax = 1000);
std::valarray<double> BiCGstab(const valarray<double>& RightHandSide, const procData& proc, const SpaceTimeDomain& dom, double epsilon, int kmax);