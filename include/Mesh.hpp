#ifndef MESH_H
#define MESH_H

#include <vector>
#include "Source.hpp"
#include "Mesh.hpp"

class Mesh {
public:
    Mesh(float dimensions, int divs);
    float* physicalLoc(int i, int j, int k) const;
    int* discreteLoc(float i, float j, float k) const;
    void addSource(const Source& source);
    void removeSource(int sourceID);
    float* calculateField(int i, int j, int k); // calculate the field at (i, j, k) and return it
    void calculateField(); // Calculate the field everywhere within the mesh
    std::vector<Source> sources;
    
    void Mesh::addSource(const Source& source);

    void Mesh::removeSource(int sourceID);

    void Mesh::calculateField();

    void Mesh::calculateField(int i, int j, int k);

private:
    int xDim, yDim, zDim, n;
    std::vector<std::vector<std::vector<std::array<float, 6>>>> fList;
};

#endif // MESH_H
