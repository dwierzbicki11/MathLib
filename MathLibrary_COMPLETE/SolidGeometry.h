#ifndef SOLID_GEOMETRY_H
#define SOLID_GEOMETRY_H
#include "Common.h"
DLLEXPORT double GetCubeVolume(double a); DLLEXPORT double GetCubeArea(double a);
DLLEXPORT double GetSphereVolume(double r); DLLEXPORT double GetSphereArea(double r);
DLLEXPORT double GetCylinderVolume(double r, double h); DLLEXPORT double GetCylinderArea(double r, double h);
DLLEXPORT double GetConeVolume(double r, double h); DLLEXPORT double GetConeArea(double r, double h);
#endif
