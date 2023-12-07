#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <memory>
#include <array>
#include "Mesh.hpp"
#include "Field.hpp"
#include "Source.hpp"

using namespace std;
double delta = 0;
const float timeStep = 0.001f;
const float c = 2.9979f * pow(10, 10); // cm/sec
const float e = 4.8029f * pow(10, -10); // statcoloumbs
const float m_e = 9.1083f * pow(10, -13); // MeV-cm

class FDTD;
int main(int argc, char** argv) {
    int dimensions = 1000;

    FDTD simulation(dimensions);
    simulation.runSimulation();

    return 0;
}

class FDTD {
public:
    std::vector<std::unique_ptr<Source>> getSources() const;
    void setSources();
    int xDim;
    int yDim;
    int zDim;

    FDTD() {
        cout << "This is the FDTD constructor!" << endl;
    }

    FDTD(int dimensions) {
        delta = static_cast <double> (1.0f / dimensions);
        xDim = dimensions; yDim = dimensions; zDim = dimensions;
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
    std::vector<std::unique_ptr<Source>> sources;
    std::array<float, 3> crossProduct(std::array<float, 3> v1, std::array<float, 3> v2);
    std::array<int, 3> crossProduct(std::array<int, 3> v1, std::array<int, 3> v2);
    
    // Calculates the energy stored inthe electric field for the entire mesh.
    // Returns the calculated energy
    double fieldEnergy() {
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

    // double counts sources, fix later.
    void progressFreeParticles() {
        for (size_t i = 0; i < sources.size(); ++i) {
            for (size_t j = 0; j < sources.size(); ++j) {
                if (i != j && sources[i]->isfree()) {
                    float mass, charge;
                    std::array<float, 3> location, velocity, force, force_sources;
                    std::array<int, 3> discreteloc;
                    Field fields;

                    mass = sources[i]->getMass();
                    charge = sources[i]->getCharge();
                    location = sources[i]->getLocation();
                    discreteloc = mesh.discreteLoc(location);
                    velocity = sources[i]->getVelocity();
                    fields = mesh[discreteloc[0]][discreteloc[1]][discreteloc[2]]);
                    force_sources = calcForce(sources[i], sources[j]);
                    force = charge * (fields.E + (1.0f / c) * crossProduct(velocity, field.B));
                    location = location + velocity * timeStep + 0.5f * (force / mass) * pow(timeStep, 2);
                    velocity = velocity + (force / mass) * timeStep;

                    sources[i]->setLocation(location);
                    sources[i]->setVelocity(velocity);
                }
            }
        }
    }

    void outputSourceData(int T) {
        ofstream outputFile("source_output.csv");
        if (T == 0) {
            outputFile << "T,i,j,k,id,source\n";
        }
        std::array<int, 3> loc;
        for (int i = 0; i < sources.size(); i++) {
            loc = sources[i].getLocation();
            outputFile << loc[0] << "," 
                << loc[1] << "," << loc[2] 
                << i << "," << sources[i].name;
        }
        outputFile.close();
    }

    std::array<float, 3> crossProduct(std::array<float, 3> v1, std::array<float, 3> v2) {
        std::array<float, 3> cross = v1[0] * v2[1] - v1[0] * v2[2]
            + v1[1] * v2[2] - v1[1] * v2[0]
            + v1[2] * v2[0] - v1[2] * v2[1];
        return cross;
    }

    std::array<int, 3> crossProduct(std::array<int, 3> v1, std::array<int, 3> v2) {
        std::array<int, 3> cross = { v1[1] * v2[2] - v1[2] * v2[1] ,
                                     v1[2] * v2[0] - v1[0] * v2[2],
                                     v1[0] * v2[1] - v1[1] * v2[0] }
        return cross;
    }

    std::array<float, 3> subtract(std::array<float, 3> v1, std::array<float, 3> v2) {
        return { v2[0] - v1[0], v2[1] - v1[1] , v2[2] - v1[2] };
    }

    float norm(std::array<float, 3> v) {
        return sqrtf(powf(v[0], 2) + powf(v[1], 2) + powf(v[2], 2));
    }

    std::array<float, 3> calcForce(Source source1, Source source2) {
        Field F;
        std::array<float, 3> r, rhat, velocity, force;
        r = subtract(source1.getLocation(), source2.getLocation());
        rhat = r / norm(r);
        velocity = subtract(source1.getVelocity(), source2.getVelocity());
        F.E = rhat * source2.getCharge() / powf(norm(r), 2);
        F.B = (1 / c) * crossProduct(source2.getVelocity(), rhat) / powf(norm(r), 2);

        force = source1.getCharge() * (F.E + (1 / c) * crossProduct(velocity, F.B));
        return force;
    }

/* Implemented as described in Electromagntic Simulation Using The FDTD Method by Dennis M. Sullivan */
    int runSimulation() {
        FILE* fp;
        int ic = xDim / 2; // center of x
        int jc = yDim / 2; // center of y
        int kc = zDim / 2; // center of z
        float ddx = 0.01f; // Cell Size
        float dt = ddx / 6.0E8f; // Time Step
        float epsz = 8.8E-12f;
        float pi = 3.14159f;
        float tolerance = 1.0e-8f;
        char dipole_orientation = 'y';

        vector<vector<vector<Field>>>* field = mesh.getMap();

        for (int i = 0; i < xDim; ++i) {
            for (int j = 0; j < yDim; ++j) {
                for (int k = 0; k < zDim; ++k) {
                    gax[i][j][k] = 1.0f;
                    gay[i][j][k] = 1.0f;
                    gaz[i][j][k] = 1.0f;
                }
            }
        }

        // Initialize the gaz array
        for (int k = 11; k < 30; k++) {
            gaz[ic][jc][k] = 0.0f;
        }
        gaz[ic][jc][kc] = 0.0f;

        float t0 = 0.0f;
        float spread = 6.0f;
        int T = 0;
        int nsteps = 1;

        ofstream outputFile("fdtd_output.csv");
        outputFile << "T,i,j,k,Ex,Ey,Ez,Hx,Hy,Hz\n";

        while (nsteps > 0) {
            cout << "nsteps --> ";
            cin >> nsteps;
            cout << nsteps << endl;

            for (int n = 1; n <= nsteps; n++) {
                T++;

                // Calculate the Dx, Dy, Dz fields
                for (int k = 1; k < zDim; k++) {
                    for (int j = 1; j < yDim; j++) {
                        for (int i = 1; i < xDim; i++) {
                            map[i][j][k].D = { map[i][j][k].D[0] + 0.5f * (map[i][j][k].H[2] - map[i][j - 1][k].H[2] - map[i][j][k].H[1] + map[i][j][k - 1].H[1]),
                                               map[i][j][k].D[1] + 0.5f * (map[i][j][k].H[0] - map[i][j][k - 1].H[0] - map[i][j][k].H[2] + map[i - 1][j][k].H[2]),
                                               map[i][j][k].D[2] + 0.5f * (map[i][j][k].H[1] - map[i - 1][j][k].H[1] - map[i][j][k].H[0] + map[i][j - 1][k].H[0]) };
                        }
                    }
                }

                // Calculate the E field from D field
                for (int k = 1; k < zDim - 1; k++) {
                    for (int j = 1; j < yDim - 1; j++) {
                        for (int i = 1; i < xDim - 1; i++) {
                            map[i][j][k].E = { gax[i][j][k] * map[i][j][k].D[0],
                                               gay[i][j][k] * map[i][j][k].D[1],
                                               gaz[i][j][k] * map[i][j][k].D[2] };
                        }
                    }
                }

                // Calculate the Hx, Hy, Hz fields
                for (int k = 1; k < zDim - 1; k++) {
                    for (int j = 1; j < yDim - 1; j++) {
                        for (int i = 1; i < xDim; i++) {
                            map[i][j][k].H = { map[i][j][k].H[0] + 0.5f * (map[i][j][k + 1].E[1] - map[i][j][k].E[1] - map[i][j + 1][k].E[2] + map[i][j][k].E[2]),
                                               map[i][j][k].H[1] + 0.5f * (map[i + 1][j][k].E[2] - map[i][j][k].E[2] - map[i][j][k + 1].E[0] + map[i][j][k].E[0]),
                                               map[i][j][k].H[2] + 0.5f * (map[i][j + 1][k].E[0] - map[i][j][k].E[0] - map[i + 1][j][k].E[1] + map[i][j][k].E[1]) };
                        }
                    }
                }


                // Source
                float pulse = exp(-0.5f * pow((t0 - T) / spread, 2.0f));
                if (dipole_orientation == 'x') {
                    map[ic][jc][kc].E[0] = pulse;
                    map[ic][jc][kc].D[0] = pulse;
                }
                else if (dipole_orientation == 'y') {
                    map[ic][jc][kc].E[1] = pulse;
                    map[ic][jc][kc].D[1] = pulse;
                }
                else if (dipole_orientation == 'z') {
                    map[ic][jc][kc].E[2] = pulse;
                    map[ic][jc][kc].D[2] = pulse;
                }

                // Absorbing Boundary Conditions for E fields
                for (int j = 0; j < yDim; j++) {
                    for (int i = 0; i < xDim; i++) {
                        map[i][j][0].E[0] = ex[i][j][1].E[0];
                        map[i][j][zDim - 1].E[0] = ex[i][j][zDim - 2].E[0];
                        map[i][0][j].E[1] = map[i][1][j].E[1];
                        map[i][yDim - 1][j].E[1] = map[i][yDim - 2][j].E[1];
                        map[0][i][j].E[2] = map[1][i][j].E[2];
                        map[xDim - 1][i][j].E[2] = map[xDim - 2][i][j].E[2];

                        map[i][j][0].H[0] = ex[i][j][1].H[0];
                        map[i][j][zDim - 1].H[0] = ex[i][j][zDim - 2].H[0];
                        map[i][0][j].H[1] = map[i][1][j].H[1];
                        map[i][yDim - 1][j].H[1] = map[i][yDim - 2][j].H[1];
                        map[0][i][j].H[2] = map[1][i][j].H[2];
                        map[xDim - 1][i][j].H[2] = map[xDim - 2][i][j].H[2];
                    }
                }

                progressFreeParticles();

                // Output data to file
                for (int i = 0; i < xDim; i++) {
                    for (int j = 0; j < yDim; j++) {
                        for (int k = 0; k < zDim; k++) {
                            if (abs(map[i][j][k].E[0] + map[i][j][k].E[1] + map[i][j][k].E[2]
                                + map[i][j][k].H[0] + map[i][j][k].H[1] + map[i][j][k].H[2]) < tolerance) {
                                outputFile << T << "," << i << "," << j << "," << k << ","
                                    << map[i][j][k].E[0] << "," << map[i][j][k].E[1] << "," << map[i][j][k].E[2] << ","
                                    << map[i][j][k].H[0] << "," << map[i][j][k].H[1] << "," << map[i][j][k].H[2] << "\n";
                            }
                        }
                    }
                }

            }
        }

        outputFile.close();
        return 0;
    }

};

