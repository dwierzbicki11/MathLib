#include "Vectors.h"
#include <math.h>
DLLEXPORT void Vector2D_FromPoints(double xA, double yA, double xB, double yB, double* vx, double* vy) { *vx = xB - xA; *vy = yB - yA; }
DLLEXPORT double Vector2D_Length(double vx, double vy) { return sqrt(vx*vx + vy*vy); }
DLLEXPORT void Vector2D_Add(double ux, double uy, double vx, double vy, double* rx, double* ry) { *rx = ux + vx; *ry = uy + vy; }
DLLEXPORT double Vector2D_DotProduct(double ux, double uy, double vx, double vy) { return ux*vx + uy*vy; }
DLLEXPORT double Vector2D_Angle(double ux, double uy, double vx, double vy) {
    double L1=Vector2D_Length(ux,uy), L2=Vector2D_Length(vx,vy); if(L1==0||L2==0)return 0;
    double c=Vector2D_DotProduct(ux,uy,vx,vy)/(L1*L2); if(c>1)c=1; if(c<-1)c=-1; return acos(c)*(180.0/M_PI);
}
DLLEXPORT double Vector3D_Length(double x, double y, double z){return sqrt(x*x+y*y+z*z);}
DLLEXPORT double Vector3D_DotProduct(double ux, double uy, double uz, double vx, double vy, double vz){return ux*vx+uy*vy+uz*vz;}
DLLEXPORT void Vector3D_CrossProduct(double ux, double uy, double uz, double vx, double vy, double vz, double* rx, double* ry, double* rz){
    *rx=uy*vz-uz*vy; *ry=uz*vx-ux*vz; *rz=ux*vy-uy*vx;
}
