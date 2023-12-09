#ifndef MESH_H
#define MESH_H

#include <array>
#include <vector>
#include <memory>
#include "Source.hpp"
#include "Field.hpp"

class Mesh {
private:
    int n, X, Y, Z;
    double scale;
    std::vector<std::unique_ptr<Source>> sources;
    std::vector<std::vector<std::vector<Field>>> map;

    void zeroMap();

public:
    Mesh(std::array<int, 3> dimensions, int divs);
    ~Mesh();

    std::vector<std::vector<std::vector<Field>>>& getMap();
    Field getField(std::array<int, 3> loc) const;
    void setSources(std::vector<std::unique_ptr<Source>>& sourceList);
    void addSource(std::unique_ptr<Source> source);
    void removeSource(size_t sourceID);
    std::array<float, 3> physicalLoc(std::array<int, 3> loc) const;
    std::array<int, 3> discreteLoc(std::array<float, 3> loc) const;
};

#endif // MESH_H
