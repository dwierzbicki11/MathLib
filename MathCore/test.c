#include <stdlib.h>
#define DLLEXPORT __declspec(dllexport)

DLLEXPORT int Dodaj(int a, int b) {
    return a + b;
}