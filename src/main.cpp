#include <iostream>
#include <memory>
#include <array>
#include "Mesh.hpp"
#include "Field.hpp"
#include "Source.hpp"
#include "Point.hpp"
#include "FDTD.hpp"
#include "CurrentLine.hpp"
#include "globals.hpp"

using namespace std;
int main(int argc, char** argv) {
    std::array<int, 3> dimensions = { 100, 100, 100 };
    int divisions = 100;

    // Create unique_ptr instances of Source (or derived classes like Point)
    auto source1 = std::make_unique<Point>();
    auto source2 = std::make_unique<Point>();

    // Set properties for source1 and source2
    source1->setLocation(std::array<float, 3>{10.0f, 10.0f, 10.0f});
    source1->setVelocity(std::array<float, 3>{1.5f, 1.5f, 1.5f});
    source2->setLocation(std::array<float, 3>{-10.0f, -10.0f, -10.0f});
    source2->setVelocity(std::array<float, 3>{-1.5f, -1.5f, -1.5f});

    // Create a vector to hold the sources
    std::vector<std::unique_ptr<Source>> sources;

    // Move the unique_ptr instances into the vector
    sources.push_back(std::move(source1));
    sources.push_back(std::move(source2));

    // Create an instance of FDTD and set up the simulation
    FDTD simulation(dimensions, divisions);

    // Move the sources into the simulation
    // Note: This will leave the 'sources' vector empty
    simulation.setSources(std::move(sources));

    // Run the simulation
    simulation.runSimulation();

    return 0;
}