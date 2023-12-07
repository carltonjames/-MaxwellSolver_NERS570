#include "Source.h"
#include <array>
#include <vector>
#include <iostream>

class Point : public Source {
public:
    Point() : Source() {
        charge = e;
        location = { 0.0f, 0.0f, 0.0f };
        distribution = false;
        dynamic = true;
        isFree = true;
        mass = 1.0f;
        name = "point";
    }

    float getCharge(std::array<float, 3> loc) const override {
        if (loc[0] == location[0] && loc[1] == location[1] && loc[2] == location[2]) {
            return charge;
        }
        else {
            return 0.0f;
        }
    }

    std::array<float, 3> getVelocity(std::array<float, 3> loc) const override {
        if (loc[0] == location[0] && loc[1] == location[1] && loc[2] == location[2]) {
            return velocity;
        }
        else {
            return { 0.0f, 0.0f, 0.0f };
        }
    }

};
