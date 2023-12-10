#include <iostream>
#include <memory>
#include <array>
#include <chrono>
#include "Mesh.hpp"
#include "Field.hpp"
#include "Source.hpp"
#include "Point.hpp"
#include "FDTD.hpp"
#include "CurrentLine.hpp"
#include "globals.hpp"

using namespace std;
using namespace std::chrono;  // Use the std::chrono namespace

int main(int argc, char** argv) {
    std::array<int, 3> dimensions = { 100, 100, 100 };
    int divisions = 100;

    // Create unique_ptr instances of Source (or derived classes like Point)
    auto source1 = std::make_unique<Point>();
    auto source2 = std::make_unique<Point>();

    // Set properties for source1 and source2
    source1->setLocation(std::array<float, 3>{10.0f, 10.0f, 10.0f});
    source1->setVelocity(std::array<float, 3>{1.5f, 0, 0});
    source2->setLocation(std::array<float, 3>{-10.0f, -10.0f, -10.0f});
    source2->setVelocity(std::array<float, 3>{-1.5f, 0, 0});
    source2->setCharge(-1.0f * source2->getCharge());
    // Create a vector to hold the sources
    std::vector<std::unique_ptr<Source>> sources;

    // Move the unique_ptr instances into the vector
    sources.push_back(std::move(source1));
    sources.push_back(std::move(source2));

    // Create an instance of FDTD and set up the simulation
    FDTD simulation(dimensions, divisions);

    std::cout << "FDTD Created" << endl;

    // Move the sources into the simulation
    simulation.setSources(std::move(sources));

    // Start timing
    auto start = high_resolution_clock::now();

    // Run the simulation
    simulation.runSimulation();

    // Stop timing
    auto stop = high_resolution_clock::now();

    // Calculate duration
    auto duration = duration_cast<microseconds>(stop - start);

    std::cout << "runSimulation() completed in " << duration.count() << " microseconds." << endl;

    return 0;
}
