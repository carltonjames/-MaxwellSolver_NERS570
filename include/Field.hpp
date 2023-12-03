#include <array>

struct Field {
    std::array<float, 3> E; // Electric field components (Ex, Ey, Ez)
    std::array<float, 3> B; // Magnetic field components (Bx, By, Bz)

    Field() : E({ 0.0f, 0.0f, 0.0f }), B({ 0.0f, 0.0f, 0.0f }) {}
};