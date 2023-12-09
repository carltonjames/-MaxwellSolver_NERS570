#ifndef FDTD_H
#define FDTD_H

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include <memory>
#include <array>
#include "Mesh.hpp"
#include "Field.hpp"
#include "Source.hpp"
#include "Point.hpp"
#include "CurrentLine.hpp"

class FDTD {
public:
    int xDim;
    int yDim;
    int zDim;

    FDTD(int dimensions, int divisions);
    FDTD(std::array<int, 3> dimensions, int divisions);
    ~FDTD();

    void setSources(std::vector<std::unique_ptr<Source>>&& sourceList);
    int runSimulation();

private:
    Mesh mesh;
    float timeStep;
    std::vector<std::unique_ptr<Source>> sources;

    double fieldEnergy();
    void progressFreeParticles();
    void outputSourceData(int T);
    std::array<float, 3> crossProductFloat(std::array<float, 3> v1, std::array<float, 3> v2);
    std::array<int, 3> crossProductInt(std::array<int, 3> v1, std::array<int, 3> v2);
    std::array<float, 3> subtract(std::array<float, 3> v1, std::array<float, 3> v2);
    float norm(std::array<float, 3> v);
    std::array<float, 3> calcForce(Source& source1, Source& source2);
    Field calcField(Source& source, std::array<float, 3>& loc);
};

#endif // FDTD_H
