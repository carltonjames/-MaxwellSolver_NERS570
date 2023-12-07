#include <iostream>
#include <vector>
#include "Mesh.hpp"
#include "Field.hpp"
#include "Source.hpp"

class Mesh;
class Source;

// TODO: Add templates to implemented sources to allow looping over the 'sources' vector
// TODO: replace functions with coordinate inputs to take in array<int, 3> instead
//using namespace std;
class FDTD {
public:
    FDTD();
    FDTD(int xDim, int yDim, int zDim, int n, float timeStep);
    void runSimulation();
    void outputResults();

    FDTD() {
        cout << "This is the FDTD constructor!" << endl;
    }
    ~FDTD() {
        cout << "This is the FDTD destructor!" << endl;
    }


private:
    Mesh mesh;
    float timeStep;
    void progressElements(); // Move free sources according to fields
    std::vector<Source> sources;
    void runSimulation() {
        for (float currentTime = 0; /* condition to end simulation */; currentTime += timeStep) {
            mesh.updateFields();

            // Output data to file for the current time step
            outputFieldsToFile("output_" + std::to_string(currentTime) + ".csv");
        }
    }

    void progressElements(){
        for (Source s : sources) {
            if (s.isFree()) {
                std::array<int, 3> loc = s.getLocation();
                float mass = s.getMass();
                std::array<int, 3> v = s.getVelocity();
            }
        }
    }
};

void FDTD::outputFieldsToFile(const std::string& filename) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    // Output the E field
    for (int i = 0; i < mesh.xDim; ++i) {
        for (int j = 0; j < mesh.yDim; ++j) {
            for (int k = 0; k < mesh.zDim; ++k) {
                file << mesh.field.E[i][j][k] << ",";
            }
            file << "\n";
        }
    }

    file.close();
}

int main() {
    int xDim = 100; // Example values
    int yDim = 100;
    int zDim = 100;
    int n = 10;
    float timeStep = 1.0f; // Example timestep

    FDTD simulation(xDim, yDim, zDim, n, timeStep);
    simulation.runSimulation();

    return 0;
}
