#ifndef MESH_H
#define MESH_H

#include <cmath>
#include <array>
#include <vector>
#include "Source.hpp"
#include "Field.hpp"

using namespace std;
class Mesh {
private:
    int n;
    std::array<float, 3> dimensions;
    double scale;
    std::vector<std::unique_ptr<Source>> sources; // list of sources within simulation
    vector<vector<vector<Field>>> map; // Defines EM field at all points in space

    void zeroMap();

public:
    Field getField(std::array<float, 3> loc) const;
    vector<vector<vector<Field>>>& getMap();
    void setSources((std::vector<std::unique_ptr<Source>>& sourceList);
    void addSource(std::unique_ptr<Source> source);
    void removeSource(size_t sourceID);
    float* physicalLoc(std::array<int, 3> loc) const;
    int* discreteLoc(std::array<float, 3> loc) const;

    Mesh(std::array<float, 3> dimensions, int divs) : dimensions(dimensions) {}

    ~Mesh(){}
};


#endif // MESH_H
