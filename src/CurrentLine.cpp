#include "Source.h"
#include <vector>
#include <iostream>

class CurrentLine : public Source {
public:
    CurrentLine() {
        charge = 1.0f;
        location = { 0, 0, 0 };
        distributionFlag = true;
        dynamic = true;
        free = true;
        mass = 1.0f;
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

    std::array<float, 3> velocity(int i, int j, int k) const override {
        std::array<float, 3> vel = { 0.0f, 0.0f, 0.0f };
        return vel;
    }
};
