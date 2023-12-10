#include "Mesh.hpp"
#include <cmath>
#include "globals.hpp"

using namespace std;

Mesh::Mesh(std::array<int, 3> dimensions, int divs) {
    scale = 1;
    n = divs;
    X = dimensions[0];
    Y = dimensions[1];
    Z = dimensions[2];

    map.resize(X);
    for (int i = 0; i < X; i++) {
        map[i].resize(Y);
        for (int j = 0; j < Y; j++) {
            map[i][j].resize(Z);
        }
    }
    zeroMap();
}

Mesh::~Mesh() = default;

// Sets all the fields within this object to zero
void Mesh::zeroMap() {
    for (int i = 0; i < X; i++) {
        for (int j = 0; j < Y; j++) {
            for (int k = 0; k < Z; k++) {
                map[i][j][k].E.fill(0);
                map[i][j][k].B.fill(0);
                map[i][j][k].D.fill(0);
                map[i][j][k].H.fill(0);
            }
        }
    }
}

// Returns a reference to the vector map of the fields
// within this Mesh object.
vector<vector<vector<Field>>>& Mesh::getMap() {
    return map;
}

// Returns the Field object at the specified location.
// Will return zero if location is outside the mesh.
Field Mesh::getField(std::array<int, 3> loc) const {
    // Checks if loc is a valid location within this Mesh
    if (loc[0] < X && loc[1] < Y && loc[2] < Z &&
        loc[0] >= 0 && loc[1] >= 1 && loc[2] >= 0) {
        return map[loc[0]][loc[1]][loc[2]];
    }
    else {
        return Field();
    }
}

// Converts a given discrete location (int array) to a 
// physical location (float array). *Note imprecise calculation!
// High mesh precision necessary for accurate conversions!
std::array<float, 3> Mesh::physicalLoc(std::array<int, 3> loc) const {
    return { static_cast<float>(loc[0] * scale),
             static_cast<float>(loc[1] * scale),
             static_cast<float>(loc[2] * scale) };
}

// Converts a given physical location (float array) to a physical
// location (int array).
std::array<int, 3> Mesh::discreteLoc(std::array<float, 3> loc) const {
    return { static_cast<int>(loc[0] / scale),
             static_cast<int>(loc[1] / scale),
             static_cast<int>(loc[2] / scale) };
}
