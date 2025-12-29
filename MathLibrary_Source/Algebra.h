#ifndef ALGEBRA_H
#define ALGEBRA_H
#include "Common.h"
DLLEXPORT int SolveLinearEquation(double a, double b, double* x);
DLLEXPORT int SolveLinearSystem2x2(double a1, double b1, double c1, double a2, double b2, double c2, double* x, double* y);
DLLEXPORT int SolveQuadraticEquation(double a, double b, double c, double* r1, double* r2);
DLLEXPORT void GetQuadraticVertex(double a, double b, double c, double* p, double* q);
DLLEXPORT double VieteSum(double a, double b);
DLLEXPORT double VieteProduct(double a, double c);
DLLEXPORT double ArithmeticSequenceTerm(double a1, double r, int n);
DLLEXPORT double ArithmeticSequenceSum(double a1, double an, int n);
DLLEXPORT double GeometricSequenceTerm(double a1, double q, int n);
DLLEXPORT double GeometricSequenceSum(double a1, double q, int n);
DLLEXPORT double InfiniteGeometricSum(double a1, double q);
DLLEXPORT double EvaluatePolynomialHorner(int degree, double coefficients[], double x);
DLLEXPORT double LogarithmArbitraryBase(double base, double value);
DLLEXPORT double CompoundInterest(double capital, double rate, int years);
DLLEXPORT void MatrixAdd3x3(double A[], double B[], double R[]);
DLLEXPORT void MatrixMultiply3x3(double A[], double B[], double R[]);
DLLEXPORT double MatrixDeterminant3x3(double A[]);
#endif
