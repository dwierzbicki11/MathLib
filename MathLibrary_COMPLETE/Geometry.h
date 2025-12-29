#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "Common.h"
DLLEXPORT double GetTriangleAreaHeron(double a, double b, double c);
DLLEXPORT double GetCircleArea(double r);
DLLEXPORT double GetCircleCircumference(double r);
DLLEXPORT double GetDistance2D(double x1, double y1, double x2, double y2);
DLLEXPORT void GetMidPoint(double x1, double y1, double x2, double y2, double* mx, double* my);
DLLEXPORT int ArePointsCollinear(double x1, double y1, double x2, double y2, double x3, double y3);
DLLEXPORT double DistancePointToLine(double A, double B, double C, double x0, double y0);
DLLEXPORT int AreLinesParallel(double a1, double a2);
DLLEXPORT int AreLinesPerpendicular(double a1, double a2);
DLLEXPORT double GetLineSlope(double x1, double y1, double x2, double y2);
#endif
