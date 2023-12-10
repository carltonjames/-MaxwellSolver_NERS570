#include "Dipole.hpp"
#include "Source.hpp"
#include "globals.hpp"
#include <cmath>
#include <complex>
using namespace std;
Dipole::Dipole() : Source() {
    charge = e * 10000000.0f; // Absurd charge to get interesting results in less time.
    location = { 0.0f, 0.0f, 0.0f };
    distribution = false;
    dynamic = false;
    boolFree = false;
    mass = m_e * 0.00000001f;
    name = "point";
    special = 'x';
    neutral = true;
}

Field Dipole::getField(int T) override {
    Field field{};
    double frequency = 5.0e12.0;
    complex<double> I(0, -1);
    if (special == 'x') {
        //field.D[0] = expf(-0.5f * powf(static_cast<float>((t0 - T) / spread), 2.0f)); // Single gaussian pulse
        field.D[0] = static_cast <float>(real(exp(I * frequency)); // oscillating dipole
    } else if (special == 'y') {
       // field.D[1] = expf(-0.5f * powf(static_cast<float>((t0 - T) / spread), 2.0f)); // Single gaussian pulse
        field.D[1] = static_cast <float>(real(exp(I * frequency)); // oscillating dipole
    } else if (special == 'z') {
        //field.D[2] = expf(-0.5f * powf(static_cast<float>((t0 - T) / spread), 2.0f)); // Single gaussian pulse
        field.D[2] = static_cast <float>(real(exp(I * frequency)); // oscillating dipole
    }
    return field;
}