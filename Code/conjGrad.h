#include "prodMatVect_BC.h"
#include "main.h"

std::vector<double> conjugateGradient(const std::vector<double>& RightHandSide, const procData& proc, const SpaceTimeDomain& dom, double epsilon = 1.0e-8, int kmax = 1000);
std::vector<double> BiCGstab(const vector<double>& RightHandSide, const procData& proc, const SpaceTimeDomain& dom, double epsilon, int kmax);