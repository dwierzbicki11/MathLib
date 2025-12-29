#ifndef TRIGONOMETRY_H
#define TRIGONOMETRY_H
#include "Common.h"
DLLEXPORT double LawOfCosines_GetSide(double a, double b, double gammaDeg);
DLLEXPORT double LawOfSines_GetSide(double side, double ang1, double ang2);
DLLEXPORT double GetCircumradius(double a, double alphaDeg);
DLLEXPORT double Cotangent(double deg);
#endif
