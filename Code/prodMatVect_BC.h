#include "main.h"

std::valarray<double> prodMatvect(const std::valarray<double>& U, const procData proc, const SpaceTimeDomain& dom);

void calculateRightHandSide(std::valarray<double>& U, const double t, procData& proc, const SpaceTimeDomain& dom);