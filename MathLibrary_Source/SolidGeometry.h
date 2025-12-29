#ifndef SOLID_GEOMETRY_H
#define SOLID_GEOMETRY_H
#include "Common.h"

// --- GRANIASTOSLUPY I WALCE ---
DLLEXPORT double GetCubeVolume(double a);
DLLEXPORT double GetCubeArea(double a);
DLLEXPORT double GetCuboidVolume(double a, double b, double c); // Prostopadloscian
DLLEXPORT double GetCuboidArea(double a, double b, double c);
DLLEXPORT double GetCylinderVolume(double r, double h);
DLLEXPORT double GetCylinderArea(double r, double h);
DLLEXPORT double GetRegularPrismVolume(double baseArea, double h); // Graniastoslup ogolny

// --- OSTROSLUPY I STOZKI ---
DLLEXPORT double GetConeVolume(double r, double h);
DLLEXPORT double GetConeArea(double r, double h); // Pole calkowite (z pitagorasa)
DLLEXPORT double GetPyramidVolume(double baseArea, double h); // Ostroslup dowolny
DLLEXPORT double GetTetrahedronVolume(double a); // Czworoscian foremny
DLLEXPORT double GetTetrahedronArea(double a);
DLLEXPORT double GetOctahedronVolume(double a); // Osmioscian foremny

// --- BRYLY SCIETE (Tablice + Studia) ---
DLLEXPORT double GetFrustumConeVolume(double r1, double r2, double h); // Stozek sciety
DLLEXPORT double GetFrustumConeArea(double r1, double r2, double h);   // Pole boczne + podstawy
DLLEXPORT double GetFrustumPyramidVolume(double areaBase1, double areaBase2, double h); // Ostroslup sciety

// --- KULA I SFERA ---
DLLEXPORT double GetSphereVolume(double r);
DLLEXPORT double GetSphereArea(double r);
DLLEXPORT double GetSphericalCapVolume(double r, double h); // Czasza kuli (odcinek kuli)
DLLEXPORT double GetSphericalCapArea(double r, double h);   // Czasza kuli (pole sfery)

// --- INNE ---
DLLEXPORT double GetEllipsoidVolume(double a, double b, double c); // Studia
#endif
