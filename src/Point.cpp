#include "Source.hpp"
#include "Point.hpp"
#include "globals.hpp"
#include <cmath>

Point::Point() : Source() {
    charge = e * 10000000.0f; // Absurd charge to get interesting results in less time.
    location = { 0.0f, 0.0f, 0.0f };
    distribution = false;
    dynamic = true;
    boolFree = true;
    mass = m_e * 0.00000001f;
    name = "point";
}