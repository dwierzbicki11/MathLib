#ifndef NUMERICAL_H
#define NUMERICAL_H
#include "Common.h"
DLLEXPORT double Integrate_Trapezoidal(int n, double y[], double h);
DLLEXPORT int LinearRegression(int n, double x[], double y[], double* a, double* b, double* r_squared);
DLLEXPORT int MatrixInverse3x3(double A[], double Result[]);
#endif
