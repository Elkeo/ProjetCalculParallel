#ifndef MAIN_H
#define MAIN_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>

#include <mpi.h>

/* Déclaration des variables globales */

extern int Nx, Ny, testCase, timeIteration;
extern double Lx, Ly, D, dt, tmax, dx, dy, t;

extern void charge(int me, int n, int np, int& iBeg, int& iEnd);

/* Déclaration d'un type contenant les informations sur un processeur */
struct procData {
   int me, n, nproc, iBeg, iEnd, locSize, tag = 10;
   MPI_Status status;
};

#endif