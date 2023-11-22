#include <vector>

std::vector<double> conjugateGradient(std::vector<double> U0, std::vector<double> RightHandSide, double epsilon = 1.0e-8, int kmax = 1000);