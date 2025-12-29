#ifndef ANALYSIS_H
#define ANALYSIS_H
#include "Common.h"
DLLEXPORT double PolynomialDerivativeAtPoint(int degree, double coefficients[], double x);
DLLEXPORT int OptimizeQuadraticFunction(double a, double b, double c, double* xOpt, double* yOpt);
DLLEXPORT void GetTangentLine(int degree, double coefficients[], double x0, double* a, double* b);
// Uniwersalna pochodna (Studia)
DLLEXPORT double CalculateNumericalDerivative(int type, double x, double params[]);
#endif
