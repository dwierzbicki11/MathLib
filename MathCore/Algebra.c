#include "Algebra.h"
#include <math.h>

DLLEXPORT int SolveLinearEquation(double a, double b, double* x) {
    if (a == 0) return (b == 0) ? 2 : 0;
    *x = -b / a;
    return 1;
}

DLLEXPORT int SolveLinearSystem2x2(double a1, double b1, double c1, double a2, double b2, double c2, double* x, double* y) {
    double W = (a1 * b2) - (b1 * a2);
    double Wx = (c1 * b2) - (b1 * c2);
    double Wy = (a1 * c2) - (c1 * a2);
    if (W == 0) return (Wx == 0 && Wy == 0) ? 2 : 0;
    *x = Wx / W; *y = Wy / W;
    return 1;
}

DLLEXPORT int SolveQuadraticEquation(double a, double b, double c, double* root1, double* root2) {
    if (a == 0.0) return -1;
    double delta = (b * b) - (4.0 * a * c);
    if (delta < 0) return 0;
    else if (delta == 0) { *root1 = -b / (2.0 * a); return 1; }
    else {
        double sqrtDelta = sqrt(delta);
        *root1 = (-b - sqrtDelta) / (2.0 * a);
        *root2 = (-b + sqrtDelta) / (2.0 * a);
        return 2;
    }
}

DLLEXPORT void GetQuadraticVertex(double a, double b, double c, double* p, double* q) {
    if (a == 0) return;
    *p = -b / (2.0 * a);
    double delta = (b * b) - (4.0 * a * c);
    *q = -delta / (4.0 * a);
}

DLLEXPORT double VieteSum(double a, double b) { return (a == 0) ? 0 : -b / a; }
DLLEXPORT double VieteProduct(double a, double c) { return (a == 0) ? 0 : c / a; }

DLLEXPORT double ArithmeticSequenceTerm(double a1, double r, int n) { return a1 + (n - 1) * r; }
DLLEXPORT double ArithmeticSequenceSum(double a1, double an, int n) { return ((a1 + an) / 2.0) * n; }
DLLEXPORT double GeometricSequenceTerm(double a1, double q, int n) { return a1 * pow(q, n - 1); }
DLLEXPORT double GeometricSequenceSum(double a1, double q, int n) {
    if (q == 1.0) return a1 * n;
    return a1 * ((1.0 - pow(q, n)) / (1.0 - q));
}
DLLEXPORT double InfiniteGeometricSum(double a1, double q) {
    if (fabs(q) >= 1.0) return 0;
    return a1 / (1.0 - q);
}

DLLEXPORT double EvaluatePolynomialHorner(int degree, double coefficients[], double x) {
    double result = coefficients[0];
    for (int i = 1; i <= degree; i++) result = result * x + coefficients[i];
    return result;
}

DLLEXPORT double LogarithmArbitraryBase(double base, double value) {
    if (base <= 0 || base == 1 || value <= 0) return -1.0;
    return log(value) / log(base);
}

DLLEXPORT double CompoundInterest(double capital, double rate, int years) {
    return capital * pow(1.0 + (rate / 100.0), years);
}