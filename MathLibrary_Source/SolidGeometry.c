#include "SolidGeometry.h"
#include <math.h>

// Graniastoslupy
DLLEXPORT double GetCubeVolume(double a){return pow(a,3);}
DLLEXPORT double GetCubeArea(double a){return 6*pow(a,2);}
DLLEXPORT double GetCuboidVolume(double a, double b, double c){return a*b*c;}
DLLEXPORT double GetCuboidArea(double a, double b, double c){return 2*(a*b + a*c + b*c);}
DLLEXPORT double GetCylinderVolume(double r, double h){return M_PI*r*r*h;}
DLLEXPORT double GetCylinderArea(double r, double h){return 2*M_PI*r*(r+h);}
DLLEXPORT double GetRegularPrismVolume(double baseArea, double h){return baseArea*h;}

// Ostroslupy
DLLEXPORT double GetConeVolume(double r, double h){return (M_PI*r*r*h)/3.0;}
DLLEXPORT double GetConeArea(double r, double h){double l=sqrt(r*r+h*h); return M_PI*r*(r+l);}
DLLEXPORT double GetPyramidVolume(double baseArea, double h){return (baseArea*h)/3.0;}

// Wielościany Foremne (Wzory z tablic)
DLLEXPORT double GetTetrahedronVolume(double a){return (pow(a,3)*sqrt(2))/12.0;}
DLLEXPORT double GetTetrahedronArea(double a){return pow(a,2)*sqrt(3);}
DLLEXPORT double GetOctahedronVolume(double a){return (pow(a,3)*sqrt(2))/3.0;}

// Bryly Sciete
DLLEXPORT double GetFrustumConeVolume(double R, double r, double h){
    // V = (PI * h / 3) * (R^2 + R*r + r^2)
    return (M_PI * h / 3.0) * (R*R + R*r + r*r);
}
DLLEXPORT double GetFrustumConeArea(double R, double r, double h){
    double l = sqrt(pow(R-r, 2) + h*h); // tworzaca stozka scietego
    return M_PI * (R*R + r*r + (R+r)*l);
}
DLLEXPORT double GetFrustumPyramidVolume(double A1, double A2, double h){
    // V = (h/3) * (A1 + A2 + sqrt(A1*A2))
    return (h/3.0) * (A1 + A2 + sqrt(A1*A2));
}

// Kula i Elementy
DLLEXPORT double GetSphereVolume(double r){return (4.0/3.0)*M_PI*pow(r,3);}
DLLEXPORT double GetSphereArea(double r){return 4*M_PI*r*r;}
DLLEXPORT double GetSphericalCapVolume(double r, double h){
    // V = (PI * h^2 / 3) * (3r - h)
    return (M_PI * h * h / 3.0) * (3*r - h);
}
DLLEXPORT double GetSphericalCapArea(double r, double h){
    // Pole samej czaszy (bez podstawy): 2 * PI * r * h
    return 2 * M_PI * r * h;
}
DLLEXPORT double GetEllipsoidVolume(double a, double b, double c){
    return (4.0/3.0) * M_PI * a * b * c;
}
