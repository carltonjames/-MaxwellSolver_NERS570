#ifndef FIELD_H
#define FIELD_H

#include "VectorField.hpp"

class Field {
public:
    Field();
    // Methods to manage fields

private:
    VectorField E_Field;
    VectorField B_Field;
};

#endif // FIELD_H
