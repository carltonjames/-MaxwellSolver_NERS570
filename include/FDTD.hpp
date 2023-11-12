#ifndef FDTD_H
#define FDTD_H

#include "Mesh.hpp"
#include "Field.hpp"
#include "Source.hpp"

class FDTD {
public:
    FDTD();
    void runSimulation();
    // Other necessary methods

private:
    Mesh* mesh;
    Field* field;
    Source* source;
    // Other necessary attributes
};

#endif // FDTD_H
