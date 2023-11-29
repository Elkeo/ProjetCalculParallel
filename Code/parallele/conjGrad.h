#ifndef CONJGRAD_H
#define CONJGRAD_H

#include "main.h"

std::vector<double> conjugateGradient(const std::vector<double>& x, const std::vector<double>& RightHandSide, const procData& procDatafile, double epsilon = 1.0e-8, int kmax = 1000000000);

#endif