#include "SolidGeometry.h"
#include <math.h>
DLLEXPORT double GetCubeVolume(double a){return pow(a,3);}
DLLEXPORT double GetCubeArea(double a){return 6*pow(a,2);}
DLLEXPORT double GetSphereVolume(double r){return (4.0/3.0)*M_PI*pow(r,3);}
DLLEXPORT double GetSphereArea(double r){return 4*M_PI*pow(r,2);}
DLLEXPORT double GetCylinderVolume(double r, double h){return M_PI*r*r*h;}
DLLEXPORT double GetCylinderArea(double r, double h){return 2*M_PI*r*(r+h);}
DLLEXPORT double GetConeVolume(double r, double h){return (1.0/3.0)*M_PI*r*r*h;}
DLLEXPORT double GetConeArea(double r, double h){return M_PI*r*(r+sqrt(r*r+h*h));}
