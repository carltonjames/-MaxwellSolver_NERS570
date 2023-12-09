#include "Source.hpp"
#include "Point.hpp"
extern double delta;
extern const float timeStep; // sec
extern const float c; // cm/sec
extern const float e; // statcoloumbs
extern const float m_e; // MeV-cm

Point::Point() : Source() {
    charge = e;
    location = { 0.0f, 0.0f, 0.0f };
    distribution = false;
    dynamic = true;
    boolFree = true;
    mass = 1.0f;
    name = "point";
}