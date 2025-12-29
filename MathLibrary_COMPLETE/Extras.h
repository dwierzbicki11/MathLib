#ifndef EXTRAS_H
#define EXTRAS_H
#include "Common.h"
DLLEXPORT int SolveQuadraticInequality(double a, double b, double c, double* x1, double* x2);
DLLEXPORT double CalculateMean(int n, double d[]);
DLLEXPORT double CalculateStandardDeviation(int n, double d[]);
DLLEXPORT int IsPrime(int n);
#endif
