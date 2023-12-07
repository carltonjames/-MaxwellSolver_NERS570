#ifndef POINT_H
#define POINT_H

#include "Source.h"

extern double delta;
extern float timeStep;
extern float c;
extern float e;
extern float m_e;

class Point : public Source {
public:
    Point() {
        charge = 1.0f;
        location = { 50.0f, 50.0f, 50.0f };
        velocity = { 0.0f, 0.0f, 0.0f };
        distribution = false;
        dynamic = true;
        free = true;
        mass = 1.0f;
    }

    Point(std::array<float, 3> loc, q, v, m) {
        location = loc;
        charge = q;
        velocity = v;
        mass = m;
    }
};

#endif // POINT_H
