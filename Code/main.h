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
   int me, nbProc, iBeg, iEnd, nbElem_y, tag;
   MPI_Status status;
};

struct SpaceTimeDomain {
   int Nx, Ny, testCase, nbProc, r;
   double Lx, Ly, D, dt, tmax, dx, dy;
};

#endif 