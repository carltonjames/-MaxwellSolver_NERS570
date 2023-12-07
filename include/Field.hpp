#ifndef FIELD_H
#define FIELD_H
#include <array>

struct Field {
    std::array<float, 3> E; // Electric field components (Ex, Ey, Ez)
    std::array<float, 3> B; // Magnetic flux density components (Bx, By, Bz)
    std::array<float, 3> H; // Magnetic field strength components (Hx, Hy, Hz)
    std::array<float, 3> D; // Displacement field components (Dx, Dy, Dz)
    Field() : E({ 0.0f, 0.0f, 0.0f }),
        B({ 0.0f, 0.0f, 0.0f }),
        H({ 0.0f, 0.0f, 0.0f }),
        D({ 0.0f, 0.0f, 0.0f }) {}
};
#endif // FIELD_H