#include <iostream>
#include <vector>
#include <cmath>
#include "Mesh.hpp"
#include "Field.hpp"
#include "Source.hpp"

double delta = 0;
const float timeStep = 0.001f;

class FDTD {
public:
    std::vector getSources() const;
    void setSources();
    int xDim;
    int yDim;
    int zDim;

    FDTD();
    FDTD(int dimensions);
    

    FDTD() {
        cout << "This is the FDTD constructor!" << endl;
    }

    FDTD(int dimensions) {
        delta = static_cast <double> (1.0f / dimensions);
        xDim = dimensions; yDim = dimensions; zDim = Dimensions;
    }

    ~FDTD() {
        cout << "This is the FDTD destructor!" << endl;
        sources.clear(); // Destroys all objects within sources
        sources.shrink_to_fit(); // Shrink the capacity of sources to zero
        delete mesh;
    }

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

    std::vector getSources() const {
        return sources;
    }
    void setSources() {
        
    }

private:
    Mesh mesh;
    float timeStep;
    std::vector<Source> sources;
    void runSimulation();
    void outputResults();
    double fieldEnergy();

    void FDTD::runSimulation() {
        for (float currentTime = 0; currentTime < 1; currentTime += timeStep) {
            mesh.updateFields();

            // Output data to file for the current time step
            // outputFieldsToFile("output_" + std::to_string(currentTime) + ".csv");
        }
    }

    // Calculates the energy stored inthe electric field for the entire mesh.
    // Returns the calculated energy
    double FDTD::fieldEnergy() {
        double totalEnergy = 0;
        Field temp;

        for (int i = 0; i < X; i++) {
            for (int j = 0; j < Y; j++) {
                for (int k = 0; k < Z; k++) {
                    temp = mesh.getField();
                    double E_magnitude = sqrt(pow(temp.E, 2));
                    double B_magnitude = sqrt(pow(temp.B, 2));
                    totalEnergy += E_magnitude * B_magnitude * pow(delta, 3) * timeStep;
                }
            }
        }

        return totalEnergy;
    }

    void outputResults() {

    }
};


int main() {
    int dimensions = 1000;

    FDTD simulation(dimensions);
    simulation.runSimulation();

    return 0;
}