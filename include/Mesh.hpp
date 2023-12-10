#ifndef MESH_H
#define MESH_H

#include <array>
#include <vector>
#include <memory>
#include "Source.hpp"
#include "Field.hpp"

// The "Mesh" object, serves as a map for each other object within a
// simulation to a discrete location in space. Information of the 'Yee'
// cells of the simulation are stored in the 'map' variable.
class Mesh {
private:
    int n, X, Y, Z; // Discrete dimensions of mesh
    double scale; // Conversion factor between discrete and physical locations of the mesh. units cell/cm
    std::vector<std::vector<std::vector<Field>>> map; // Defines the EM field at every point of the mesh

    // Sets all the fields within this object to zero
    void zeroMap();

public:
    Mesh(std::array<int, 3> dimensions, int divs);
    ~Mesh();

    // Returns a reference to the vector map of the fields
    // within this Mesh object.
    std::vector<std::vector<std::vector<Field>>>& getMap();

    // Returns the Field object at the specified location.
    // Will return zero if location is outside the mesh.
    Field getField(std::array<int, 3> loc) const;

    // Converts a given discrete location (int array) to a 
    // physical location (float array). *Note imprecise calculation!
    // High mesh precision necessary for accurate conversions!
    std::array<float, 3> physicalLoc(std::array<int, 3> loc) const;

    // Converts a given physical location (float array) to a physical
    // location (int array).
    std::array<int, 3> discreteLoc(std::array<float, 3> loc) const;
};

#endif // MESH_H
