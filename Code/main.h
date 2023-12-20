#ifndef MAIN_H
#define MAIN_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <valarray>
#include <chrono>
#include "mpi.h"

using namespace std;

/* Déclaration d'un type contenant les informations sur un processeur */
struct procData {
   int me, nbProc, iBeg, iEnd, nbElem_y, tag, Upper = 1, Lower = 0;
   MPI_Status status;
   valarray<int> neighborsToMe = { MPI_PROC_NULL, MPI_PROC_NULL };
   valarray<int> stencilLower, stencilUpper;
};

struct SpaceTimeDomain {
   int Nx, Ny, testCase, nbProc, r;
   double Lx, Ly, D, dt, tmax, dx, dy, a, b;
};

#endif 