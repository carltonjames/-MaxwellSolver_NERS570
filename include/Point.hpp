#ifndef POINT_H
#define POINT_H

#include "Source.h"

class Point : public Source {
public:
    Point() {
        charge = 1.0f;
        location = { 50, 50, 50 };
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
