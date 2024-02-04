#ifndef MAIN_H
#define MAIN_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <vector>
#include <chrono>
#include "mpi.h"

using namespace std;


/* Déclaration d'un type contenant les informations sur un processeur */
struct procData {
   int me, iBeg, iEnd, nbElem_y, tag;
   int nbProc;
   MPI_Status status;
   vector<int> neighborsToMe;
   vector<double> stencilUnder, stencilOver;
};

struct SpaceTimeDomain {
   int Nx, Ny, testCase, nbProc, r;
   double Lx, Ly, D, dt, tmax, dx, dy, alpha, beta, epsilon = 1e-6;
};

/* Surcharge d'opérateur */
// Addition
template <typename T>
vector<T> operator+(const vector<T>& a, const vector<T>& b)
{
   vector<T> c(a.size());
   for (size_t i = 0; i < a.size(); i++)
   {
      c[i] = a[i] + b[i];
   }
   return c;
}

// Soustraction
template <typename T>
vector<T> operator-(const vector<T>& a, const vector<T>& b)
{
   vector<T> c(a.size());
   for (size_t i = 0; i < a.size(); i++)
   {
      c[i] = a[i] - b[i];
   }
   return c;
}

// Multiplication par un scalaire
template <typename T>
vector<T> operator*(const T& scalar, const vector<T>& a)
{
   vector<T> result(a.size());
   for (size_t i = 0; i < a.size(); i++)
   {
      result[i] = scalar * a[i];
   }
   return result;
}

template <typename T>
vector<T> operator*(const vector<T>& a, const T& scalar)
{
   return scalar * a;
}

//Division par un scalaire
template <typename T>
vector<T> operator/(const vector<T>& a, const T& scalar)
{
   vector<T> result(a.size());
   for (size_t i = 0; i < a.size(); i++)
   {
      result[i] = a[i] / scalar;
   }
   return result;
}

// Produit de deux vecteurs
template <typename T>
vector<T> operator*(const vector<T>& a, const vector<T>& b)
{
   vector<T> result(a.size());
   for (size_t i = 0; i < a.size(); i++)
   {
      result[i] = a[i] * b[i];
   }
   return result;
}

// Fonction pour sommer les éléments d'un vecteur
template <typename T>
T sum(const vector<T>& vec)
{
   T result = 0;
   for (const auto& element : vec)
   {
      result += element;
   }
   return result;
}

// Fonction pour comparer un vecteur à un scalaire
template <typename T>
bool operator>(const vector<T>& a, const T& b)
{
   bool result = true;
   for (size_t i = 0; i < a.size(); i++)
   {
      result *= (a[i] > b);
   }
   return result;
}



#endif 