#ifndef COMMON_H
#define COMMON_H
#ifdef _WIN32
    #define DLLEXPORT __declspec(dllexport)
#else
    #define DLLEXPORT
#endif
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#endif
