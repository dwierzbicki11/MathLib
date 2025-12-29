#include "Geometry.h"
#include <math.h>

DLLEXPORT double GetTriangleAreaHeron(double a, double b, double c) {
    if (a + b <= c || a + c <= b || b + c <= a) return -1.0;
    double p = (a + b + c) / 2.0;
    return sqrt(p * (p - a) * (p - b) * (p - c));
}

DLLEXPORT double GetCircleArea(double r) { return (r < 0) ? -1 : M_PI * r * r; }
DLLEXPORT double GetCircleCircumference(double r) { return (r < 0) ? -1 : 2 * M_PI * r; }

DLLEXPORT double GetDistance2D(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

DLLEXPORT void GetMidPoint(double x1, double y1, double x2, double y2, double* mx, double* my) {
    *mx = (x1 + x2) / 2.0; *my = (y1 + y2) / 2.0;
}

DLLEXPORT int ArePointsCollinear(double x1, double y1, double x2, double y2, double x3, double y3) {
    double det = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
    return (fabs(det) < 1e-9) ? 1 : 0;
}

DLLEXPORT double DistancePointToLine(double A, double B, double C, double x0, double y0) {
    double den = sqrt(A * A + B * B);
    return (den == 0) ? -1 : fabs(A * x0 + B * y0 + C) / den;
}

DLLEXPORT int AreLinesParallel(double a1, double a2) { return (fabs(a1 - a2) < 1e-9); }
DLLEXPORT int AreLinesPerpendicular(double a1, double a2) { return (fabs(a1 * a2 + 1.0) < 1e-9); }

DLLEXPORT double GetLineSlope(double x1, double y1, double x2, double y2) {
    return (x1 == x2) ? 0 : (y2 - y1) / (x2 - x1);
}