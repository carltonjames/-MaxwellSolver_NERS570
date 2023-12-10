#include "CurrentLine.hpp"
CurrentLine::CurrentLine() : Source() {
        charge = 1.0f;
        location = { 0.0f, 0.0f, 0.0f };
        distribution = true;
        dynamic = false;
        boolFree = false;
        mass = 1.0f;
        name = "CurrentLine";
}

