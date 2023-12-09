#include "Source.hpp"
#include "Point.hpp"
#include "globals.hpp"

Point::Point() : Source() {
    charge = e;
    location = { 0.0f, 0.0f, 0.0f };
    distribution = false;
    dynamic = true;
    boolFree = true;
    mass = 1.0f;
    name = "point";
}