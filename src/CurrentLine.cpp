#include "Source.h"
#include <vector>
#include <iostream>

extern double delta;
extern float timeStep;
extern float c;
extern float e;
extern float m_e;

class CurrentLine : public Source {
public:
    CurrentLine() {
        charge = 1.0f;
        location = { 0.0f, 0.0f, 0.0f };
        distribution = true;
        dynamic = false;
        isFree = false;
        mass = 1.0f;
        name = "CurrentLine";
    }

    currentLine() {
        cout << "Hello from the CurrentLine default constructor!" << endl;
    };


    float getCharge(int i, int j, int k) const override {
        std::array<int, 3> charge;
        if (distribution) {
            if (i == 50 && j == 50) {
                return charge;
            }
            else {
                return 0.0f;
            }
        }
        else {
            cout << "Not a distribution!" << endl;
        }
    }

    std::array<float, 3> velocity(std::array<int, 3> loc) const override {
        std::array<float, 3> vel = { 0.0f, 0.0f, 0.0f };
        return vel;
    }
};
