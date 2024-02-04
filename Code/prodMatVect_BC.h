#include "main.h"

std::vector<double> prodMatvect(const std::vector<double>& U, const procData proc, const SpaceTimeDomain& dom);

void calculateRightHandSide(std::vector<double>& U, const double t, const procData& proc, const SpaceTimeDomain& dom);