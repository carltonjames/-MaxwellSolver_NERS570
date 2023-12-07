#ifndef FDTD_H
#define FDTD_H

#include "Mesh.hpp"
#include "Field.hpp"
#include "Source.hpp"

class FDTD {
public:
    FDTD();
    void runSimulation();

private:
    Mesh* mesh;
    Field* field;
    Source* source;
};

#endif // FDTD_H
