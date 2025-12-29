#ifndef VECTORS_H
#define VECTORS_H
#include "Common.h"
DLLEXPORT void Vector2D_FromPoints(double xA, double yA, double xB, double yB, double* vx, double* vy);
DLLEXPORT double Vector2D_Length(double vx, double vy);
DLLEXPORT void Vector2D_Add(double ux, double uy, double vx, double vy, double* resX, double* resY);
DLLEXPORT double Vector2D_DotProduct(double ux, double uy, double vx, double vy);
DLLEXPORT double Vector2D_Angle(double ux, double uy, double vx, double vy);
DLLEXPORT double Vector3D_Length(double vx, double vy, double vz);
DLLEXPORT double Vector3D_DotProduct(double ux, double uy, double uz, double vx, double vy, double vz);
DLLEXPORT void Vector3D_CrossProduct(double ux, double uy, double uz, double vx, double vy, double vz, double* rx, double* ry, double* rz);
#endif
