#ifndef FIELD_H
#define FIELD_H
#include <array>
#include <stdexcept>

struct Field {
    std::array<float, 3> E; // Electric field components (Ex, Ey, Ez)
    std::array<float, 3> B; // Magnetic flux density components (Bx, By, Bz)
    std::array<float, 3> H; // Magnetic field strength components (Hx, Hy, Hz)
    std::array<float, 3> D; // Displacement field components (Dx, Dy, Dz)

    Field() : E({ 0.0f, 0.0f, 0.0f }),
        B({ 0.0f, 0.0f, 0.0f }),
        H({ 0.0f, 0.0f, 0.0f }),
        D({ 0.0f, 0.0f, 0.0f }) {}

    // Operator overloads
    Field operator+(const Field& other) const {
        Field result;
        for (int i = 0; i < 3; ++i) {
            result.E[i] = this->E[i] + other.E[i];
            result.B[i] = this->B[i] + other.B[i];
            result.H[i] = this->H[i] + other.H[i];
            result.D[i] = this->D[i] + other.D[i];
        }
        return result;
    }

    Field operator-(const Field& other) const {
        Field result;
        for (int i = 0; i < 3; ++i) {
            result.E[i] = this->E[i] - other.E[i];
            result.B[i] = this->B[i] - other.B[i];
            result.H[i] = this->H[i] - other.H[i];
            result.D[i] = this->D[i] - other.D[i];
        }
        return result;
    }

    Field operator*(float scalar) const {
        Field result;
        for (int i = 0; i < 3; ++i) {
            result.E[i] = this->E[i] * scalar;
            result.B[i] = this->B[i] * scalar;
            result.H[i] = this->H[i] * scalar;
            result.D[i] = this->D[i] * scalar;
        }
        return result;
    }

    Field operator/(float scalar) const {
        if (scalar == 0.0f) {
            throw std::runtime_error("Division by zero in Field struct");
        }
        Field result;
        for (int i = 0; i < 3; ++i) {
            result.E[i] = this->E[i] / scalar;
            result.B[i] = this->B[i] / scalar;
            result.H[i] = this->H[i] / scalar;
            result.D[i] = this->D[i] / scalar;
        }
        return result;
    }
};
#endif // FIELD_H
