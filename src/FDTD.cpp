#include "FDTD.hpp"

using namespace std;

FDTD::FDTD(int dimensions, int divisions) : mesh({ dimensions, dimensions, dimensions }, divisions) {
    delta = static_cast<double>(1.0f / dimensions);
    xDim = dimensions; yDim = dimensions; zDim = dimensions;
}

FDTD::FDTD(std::array<int, 3> dimensions, int divisions) : mesh(dimensions, divisions) {
    delta = static_cast<double>(1.0f / dimensions[0]);
    xDim = dimensions[0]; yDim = dimensions[1]; zDim = dimensions[2];
}

FDTD::~FDTD() {
    std::cout << "This is the FDTD destructor!" << std::endl;
}

void FDTD::setSources(std::vector<std::unique_ptr<Source>>&& sourceList) {
    sources = std::move(sourceList);
}

//std::vector<std::unique_ptr<Source>> getSources() const {
//    return sources;
//}


   
// Calculates the energy stored in the electromagnetic fields across the entire mesh.
// Returns the calculated energy
double FDTD::fieldEnergy() {
    double totalEnergy = 0;
    Field temp;

    for (int i = 0; i < xDim; i++) {
        for (int j = 0; j < yDim; j++) {
            for (int k = 0; k < zDim; k++) {
                temp = mesh.getField(array<int, 3>{ i, j, k });
                double E_magnitude = sqrt(pow(temp.E[0], 2) + pow(temp.E[1], 2) + pow(temp.E[2], 2));
                double B_magnitude = sqrt(pow(temp.B[0], 2) + pow(temp.B[1], 2) + pow(temp.B[2], 2));
                totalEnergy += E_magnitude * B_magnitude * pow(delta, 3) * timeStep;
            }
        }
    }

    return totalEnergy;
}

// nonphysical: this function assumes instantaneous interactions between sources.
// For this reason, this particular implementation should be avoided.
void FDTD::old_progressFreeParticles() {
    std::cout << "Starting progressFreeParticles" << std::endl;
    for (auto& source : sources) {
        std::cout << "Processing source" << std::endl;
        if (source->isFree()) {
            float mass = source->getMass();
            float charge = source->getCharge();
            std::array<float, 3> location = source->getLocation();
            std::array<int, 3> discreteloc = mesh.discreteLoc(location);
            std::array<float, 3> velocity = source->getVelocity();
            Field fields = mesh.getField(discreteloc);

            // Calculate the net force on the particle
            std::array<float, 3> netForce = { 0.0f, 0.0f, 0.0f };
            for (auto& otherSource : sources) {
                if (otherSource != source) {
                    std::array<float, 3> forceFromOther = calcForce(*source, *otherSource);
                    std::cout << "Force from other source: " << forceFromOther[0] << ", " << forceFromOther[1] << ", " << forceFromOther[2] << std::endl;
                    netForce[0] += forceFromOther[0];
                    netForce[1] += forceFromOther[1];
                    netForce[2] += forceFromOther[2];
                }
            }

            // Add Lorentz force from the electromagnetic field
            std::array<float, 3> lorentzForce;
            std::array<float, 3> crossProd = crossProductFloat(velocity, fields.B);
            for (int i = 0; i < 3; ++i) {
                lorentzForce[i] = charge * (fields.E[i] + (1.0f / c) * crossProd[i]);
            }
            netForce[0] += lorentzForce[0];
            netForce[1] += lorentzForce[1];
            netForce[2] += lorentzForce[2];

            // Update position and velocity
            std::array<float, 3> newLocation = {
                location[0] + velocity[0] * timeStep + 0.5f * (netForce[0] / mass) * powf(timeStep, 2),
                location[1] + velocity[1] * timeStep + 0.5f * (netForce[1] / mass) * powf(timeStep, 2),
                location[2] + velocity[2] * timeStep + 0.5f * (netForce[2] / mass) * powf(timeStep, 2)
            };

            std::array<float, 3> newVelocity = {
                velocity[0] + (netForce[0] / mass) * timeStep,
                velocity[1] + (netForce[1] / mass) * timeStep,
                velocity[2] + (netForce[2] / mass) * timeStep
            };

            // Check for overflow
            for (int i = 0; i < 3; ++i) {
                if (!std::isfinite(newLocation[i]) || !std::isfinite(newVelocity[i])) {
                    std::cerr << "Overflow detected in progressFreeParticles" << std::endl;
                    std::cout << "Updated source velocity: " << newVelocity[0] << ", " << newVelocity[1] << ", " << newVelocity[2] << std::endl;
                    std::cout << "Updated source location: " << newLocation[0] << ", " << newLocation[1] << ", " << newLocation[2] << std::endl;
                    return; // Exit the function to prevent further calculations
                }
            }

            source->setLocation(newLocation);
            source->setVelocity(newVelocity);

            std::cout << "Updated source location: " << newLocation[0] << ", " << newLocation[1] << ", " << newLocation[2] << std::endl;
            std::cout << "Updated source velocity: " << newVelocity[0] << ", " << newVelocity[1] << ", " << newVelocity[2] << std::endl;
        }
    }

    std::cout << "Finished progressFreeParticles" << std::endl;
}

// Advances particles forward according to fields propagating in the Mesh
void FDTD::progressFreeParticles() {
    std::cout << "Starting progressFreeParticles" << std::endl;
    for (auto& source : sources) {
        std::cout << "Processing source" << std::endl;
        if (source->isFree()) {
            float mass = source->getMass();
            float charge = source->getCharge();
            std::array<float, 3> location = source->getLocation();
            std::array<int, 3> discreteloc = mesh.discreteLoc(location);
            std::array<float, 3> velocity = source->getVelocity();
            Field fields = mesh.getField(discreteloc);

            // Calculate the net force on the particle
            std::array<float, 3> netForce = { 0.0f, 0.0f, 0.0f };

            // Add Lorentz force from the electromagnetic field
            std::array<float, 3> lorentzForce;
            std::array<float, 3> crossProd = crossProductFloat(velocity, fields.B);
            for (int i = 0; i < 3; ++i) {
                lorentzForce[i] = charge * (fields.E[i] + (1.0f / c) * crossProd[i]);
            }
            netForce[0] += lorentzForce[0];
            netForce[1] += lorentzForce[1];
            netForce[2] += lorentzForce[2];

            // Update position and velocity
            std::array<float, 3> newLocation = {
                location[0] + velocity[0] * timeStep + 0.5f * (netForce[0] / mass) * powf(timeStep, 2),
                location[1] + velocity[1] * timeStep + 0.5f * (netForce[1] / mass) * powf(timeStep, 2),
                location[2] + velocity[2] * timeStep + 0.5f * (netForce[2] / mass) * powf(timeStep, 2)
            };

            std::array<float, 3> newVelocity = {
                velocity[0] + (netForce[0] / mass) * timeStep,
                velocity[1] + (netForce[1] / mass) * timeStep,
                velocity[2] + (netForce[2] / mass) * timeStep
            };

            // Check for overflow
            for (int i = 0; i < 3; ++i) {
                if (!std::isfinite(newLocation[i]) || !std::isfinite(newVelocity[i])) {
                    std::cerr << "Overflow detected in progressFreeParticles" << std::endl;
                    std::cout << "Updated source velocity: " << newVelocity[0] << ", " << newVelocity[1] << ", " << newVelocity[2] << std::endl;
                    std::cout << "Updated source location: " << newLocation[0] << ", " << newLocation[1] << ", " << newLocation[2] << std::endl;
                    return; // Exit the function to prevent further calculations
                }
            }

            source->setLocation(newLocation);
            source->setVelocity(newVelocity);

            std::cout << "Updated source location: " <<
                newLocation[0] << ", " << newLocation[1] <<
                ", " << newLocation[2] << std::endl;

            std::cout << "Updated source velocity: " <<
                newVelocity[0] << ", " << newVelocity[1] <<
                ", " << newVelocity[2] << std::endl;
        }
    }

    std::cout << "Finished progressFreeParticles" << std::endl;
}

// Update position of dynamic source based on its velocity
void FDTD::updateDynamicSource() {
    std::cout << "Starting updateDynamicSource" << std::endl;
    for (auto& source : sources) {
        std::cout << "Processing source" << std::endl;
        if (source->isDynamic() && !(source->isFree())) {
            std::array<float, 3> location = source->getLocation();
            std::array<float, 3> velocity = source->getVelocity();

            std::cout << "Dynamic source location before update: " << location[0] <<
                ", " << location[1] << ", " << location[2] << std::endl;

            // Update position according to velocity
            location[0] += velocity[0] * timeStep;
            location[1] += velocity[1] * timeStep;
            location[2] += velocity[2] * timeStep;

            std::cout << "Dynamic source location after update: " << location[0] <<
                ", " << location[1] << ", " << location[2] << std::endl;

            source->setLocation(location);
        }
    }
}

// Outputs the data for each source file at the current time of the simulation
void FDTD::outputSourceData() {
    ios_base::openmode mode = ios_base::out;
    if (T > 1) {
        mode |= ios_base::app;
    }

    ofstream outputFile("source_output.csv", mode);

    if (T == 1) {
        outputFile << "T,i,j,k,id,source\n";
    }

    std::array<float, 3> loc;
    for (int i = 0; i < sources.size(); i++) {
        loc = sources[i]->getLocation();
        outputFile << T << ","
            << std::fixed << std::setprecision(20) // Sets precision of floats in output
            << loc[0] << ","
            << loc[1] << ","
            << loc[2] << ","
            << i << ","
            << sources[i]->getName() << endl;
    }

    outputFile.close();
}

// Takes the cross product of two vectors (array of floats) from eachother and returns the result
std::array<float, 3> FDTD::crossProductFloat(std::array<float, 3> v1, std::array<float, 3> v2) {
    std::array<float, 3> cross = { v1[1] * v2[2] - v1[2] * v2[1],
                                    v1[2] * v2[0] - v1[0] * v2[2],
                                    v1[0] * v2[1] - v1[1] * v2[0] };
    return cross;
}

// Takes the cross product of two vectors (array of ints) from eachother and returns the result
std::array<int, 3> FDTD::crossProductInt(std::array<int, 3> v1, std::array<int, 3> v2) {
    std::array<int, 3> cross = { v1[1] * v2[2] - v1[2] * v2[1],
                                    v1[2] * v2[0] - v1[0] * v2[2],
                                    v1[0] * v2[1] - v1[1] * v2[0] };
    return cross;
}

// Subtracts two vectors (arrays) from eachother and returns the result
std::array<float, 3> FDTD::subtract(std::array<float, 3> v1, std::array<float, 3> v2) {
    return { v2[0] - v1[0], v2[1] - v1[1] , v2[2] - v1[2] };
}

// Calculates the norm of an input vector (array) and returns it.
float FDTD::norm(std::array<float, 3> v) {
    return sqrtf(powf(v[0], 2) + powf(v[1], 2) + powf(v[2], 2));
}

// Calculates the force between two sources. Takes references to two
// Source objects as inputs
std::array<float, 3> FDTD::calcForce(Source& source1, Source& source2) {
    Field F;
    std::array<float, 3> r, rhat, velocity, force;
    r = subtract(source1.getLocation(), source2.getLocation());

    float distance = norm(r);
    if (distance == 0.0f) {
        std::cerr << "Error: Division by zero in calcForce" << std::endl;
        return { 0.0f, 0.0f, 0.0f };
    }

    for (int i = 0; i < 3; ++i) {
        rhat[i] = r[i] / distance;
    }

    velocity = subtract(source1.getVelocity(), source2.getVelocity());

    float chargeFactor = source2.getCharge() / (distance * distance);
    for (int i = 0; i < 3; ++i) {
        F.E[i] = rhat[i] * chargeFactor;
    }

    std::array<float, 3> crossProd = crossProductFloat(source2.getVelocity(), rhat);
    float factor = (1 / c) / (distance * distance);
    for (int i = 0; i < 3; ++i) {
        F.B[i] = crossProd[i] * factor;
    }

    crossProd = crossProductFloat(velocity, F.B);
    float charge = source1.getCharge();
    for (int i = 0; i < 3; ++i) {
        force[i] = charge * (F.E[i] + (1 / c) * crossProd[i]);
    }

    std::cout << "Force calculated: " << force[0] << ", " << force[1] << ", " << force[2] << std::endl;
    return force;
}

// Calculates the field generated by a given source at a
// particular location. Take a reference to a source and a
// reference to location (array of floats) as input
Field FDTD::calcField(Source& source, array<float, 3>& loc) {
    Field F{};

    std::array<float, 3> r, rhat, velocity, force;
    r = subtract(loc, source.getLocation());

    for (int i = 0; i < 3; ++i) {
        rhat[i] = r[i] / norm(r);
    }

    velocity = source.getVelocity();

    float chargeFactor = source.getCharge() / powf(norm(r), 2);
    for (int i = 0; i < 3; ++i) {
        F.E[i] = rhat[i] * chargeFactor;
    }

    std::array<float, 3> crossProd = crossProductFloat(velocity, rhat);
    float factor = (1 / c) / powf(norm(r), 2);
    for (int i = 0; i < 3; ++i) {
        F.B[i] = crossProd[i] * factor;
    }

    return F;
}

// Discretizes electromagnetic field generated by sources.
Field FDTD::discreteCalcField(Source& source) {
    Field F{};
    float cellVolume = xDim * yDim * zDim;
    float charge = source.getCharge();
    array<float, 3> velocity = source.getVelocity();
    float regularizationFactor = sqrtf(powf(xDim / 2, 2) + powf(xDim / 2, 2) + powf(xDim / 2, 2));

    if (source.isNeutral()) {
        return F;
    }

    float chargeFactor = source.getCharge() / powf(regularizationFactor, 2);
    for (int i = 0; i < 3; ++i) {
        F.E[i] = regularizationFactor * chargeFactor;
    }

    std::array<float, 3> crossProd =
        crossProductFloat(velocity,
            array<float, 3>{regularizationFactor, regularizationFactor, regularizationFactor});

    float factor = (1 / c) / powf(norm(
        array<float, 3>{regularizationFactor, regularizationFactor, regularizationFactor}), 2);
    for (int i = 0; i < 3; ++i) {
        F.B[i] = crossProd[i] * factor;
    }

    return F;
}

/* Implemented as described in Electromagntic Simulation Using The FDTD Method by Dennis M. Sullivan */
int FDTD::runSimulation() {
    std::cout << "FDTD::runSimulation() called" << endl;
    int ic = xDim / 2; // center of x
    int jc = yDim / 2; // center of y
    int kc = zDim / 2; // center of z
    float ddx = 0.01f; // Cell Size
    float dt = ddx / 6.0E4f; // Time Step
    float epsz = 8.8E-12f;
    float pi = 3.14159f;
    float tolerance = 1.0e-8f;
    char dipole_orientation = 'y';
    std::cout << "1. variables initialized in runSimulation()" << endl;
    array<float, 3> center = { static_cast<float>(ic),
                                static_cast<float>(jc),
                                static_cast<float>(kc) };

    // Setting gauge factors
    std::vector<std::vector<std::vector<float>>> gax(xDim, std::vector<std::vector<float>>(yDim, std::vector<float>(zDim, 1.0f)));
    std::vector<std::vector<std::vector<float>>> gay(xDim, std::vector<std::vector<float>>(yDim, std::vector<float>(zDim, 1.0f)));
    std::vector<std::vector<std::vector<float>>> gaz(xDim, std::vector<std::vector<float>>(yDim, std::vector<float>(zDim, 1.0f)));

    std::cout << "2. variables initialized in runSimulation()" << endl;
    vector<vector<vector<Field>>>& map = mesh.getMap();
    timeStep = dt;

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

    float t0 = 10.0f;
    float spread = 6.0f * 100.0f;
    int T = 0;
    nsteps = 100;

    ofstream outputFile("fdtd_output.csv");
    outputFile << "T,i,j,k,Ex,Ey,Ez,Hx,Hy,Hz\n";

    std::cout << "fdtd_output.csv declared, written" << endl;

    cout << "nsteps --> " << nsteps << endl;

    // Primary FDTD loop
    for (int n = 1; n <= nsteps; n++) {
        std::cout << "Top of for-loop" << endl;
        T++;

        cout << "Time --> " << T << endl;

        std::cout << "Calculate the Dx, Dy, Dz fields" << endl;
        // Calculate the Dx, Dy, Dz fields
        for (int k = 1; k < zDim; k++) {
            for (int j = 1; j < yDim; j++) {
                for (int i = 1; i < xDim; i++) {
                    map[i][j][k].D = { map[i][j][k].D[0] + 0.5f * (map[i][j][k].H[2] - map[i][j - 1][k].H[2] - map[i][j][k].H[1] + map[i][j][k - 1].H[1]),
                                        map[i][j][k].D[1] + 0.5f * (map[i][j][k].H[0] - map[i][j][k - 1].H[0] - map[i][j][k].H[2] + map[i - 1][j][k].H[2]),
                                        map[i][j][k].D[2] + 0.5f * (map[i][j][k].H[1] - map[i - 1][j][k].H[1] - map[i][j][k].H[0] + map[i][j - 1][k].H[0]) };
                    map[i][j][k].E = { gax[i][j][k] * map[i][j][k].D[0],
                                        gay[i][j][k] * map[i][j][k].D[1],
                                        gaz[i][j][k] * map[i][j][k].D[2] };
                }
            }
        }

        std::cout << "Calculate the Hx, Hy, Hz fields" << endl;
        // Calculate the Hx, Hy, Hz fields
        for (int k = 1; k < zDim - 1; ++k) {
            for (int j = 1; j < yDim - 1; ++j) {
                for (int i = 1; i < xDim - 1; ++i) {
                    map[i][j][k].H = { map[i][j][k].H[0] + 0.5f * (map[i][j][k + 1].E[1] - map[i][j][k].E[1] - map[i][j + 1][k].E[2] + map[i][j][k].E[2]),
                                        map[i][j][k].H[1] + 0.5f * (map[i + 1][j][k].E[2] - map[i][j][k].E[2] - map[i][j][k + 1].E[0] + map[i][j][k].E[0]),
                                        map[i][j][k].H[2] + 0.5f * (map[i][j + 1][k].E[0] - map[i][j][k].E[0] - map[i + 1][j][k].E[1] + map[i][j][k].E[1]) };
                    map[i][j][k].B = { gax[i][j][k] * map[i][j][k].H[0],
                                        gay[i][j][k] * map[i][j][k].H[1],
                                        gaz[i][j][k] * map[i][j][k].H[2] };
                }
            }
        }

        std::cout << "Setting pulse" << endl;
        // Source setting dipole source
        // TODO: Replace with dipole source object        
        // float pulse = exp(-0.5f * powf(static_cast<float>((t0 - T) / spread), 2.0f));
        Field field{};
        for (auto& source : sources) {
            std::cout << "iterating over source " << source->getName() <<endl;
            array<int, 3> loc;
            loc = mesh.discreteLoc(source->getLocation());
            field = discreteCalcField(*source);
            if (source->getName().compare("dipole") == 0) {
                map[loc[0]][loc[1]][loc[2]] = source->getField(T);
            } else if (loc[0] < xDim && loc[1] < yDim && loc[2] < zDim &&
                loc[0] >= 0 && loc[1] >= 1 && loc[2] >= 0) {
                map[loc[0]][loc[1]][loc[2]].E[0] += gax[loc[0]][loc[1]][loc[2]] * field.D[0];
                map[loc[0]][loc[1]][loc[2]].E[1] += gay[loc[0]][loc[1]][loc[2]] * field.D[1];
                map[loc[0]][loc[1]][loc[2]].E[2] += gaz[loc[0]][loc[1]][loc[2]] * field.D[2];
                map[loc[0]][loc[1]][loc[2]].B[0] += field.B[0];
                map[loc[0]][loc[1]][loc[2]].B[1] += field.B[1];
                map[loc[0]][loc[1]][loc[2]].B[2] += field.B[2];

                map[loc[0]][loc[1]][loc[2]].D[0] += field.E[0];
                map[loc[0]][loc[1]][loc[2]].D[1] += field.E[1];
                map[loc[0]][loc[1]][loc[2]].D[2] += field.E[2];
                map[loc[0]][loc[1]][loc[2]].H[0] += field.B[0];
                map[loc[0]][loc[1]][loc[2]].H[1] += field.B[1];
                map[loc[0]][loc[1]][loc[2]].H[2] += field.B[2];
            }
        }

        std::cout << "Absorption" << endl;
        // Absorbing Boundary Conditions for E fields
        // Perfectly-matched Layer
        for (int j = 0; j < yDim; j++) {
            for (int i = 0; i < xDim; i++) {
                map[i][j][0].E[0] = map[i][j][1].E[0];
                map[i][j][zDim - 1].E[0] = map[i][j][zDim - 2].E[0];
                map[i][0][j].E[1] = map[i][1][j].E[1];
                map[i][yDim - 1][j].E[1] = map[i][yDim - 2][j].E[1];
                map[0][i][j].E[2] = map[1][i][j].E[2];
                map[xDim - 1][i][j].E[2] = map[xDim - 2][i][j].E[2];

                map[i][j][0].H[0] = map[i][j][1].H[0];
                map[i][j][zDim - 1].H[0] = map[i][j][zDim - 2].H[0];
                map[i][0][j].H[1] = map[i][1][j].H[1];
                map[i][yDim - 1][j].H[1] = map[i][yDim - 2][j].H[1];
                map[0][i][j].H[2] = map[1][i][j].H[2];
                map[xDim - 1][i][j].H[2] = map[xDim - 2][i][j].H[2];
            }
        }

        outputSourceData(); // Output source data
        progressFreeParticles(); // Advances free particles according to subjected forces
        updateDynamicSource(); // Advances dynamic, non-free, particles

        std::cout << "Outputting data" << endl;
        // Output data to file
        for (int i = 0; i < xDim; i++) {
            for (int j = 0; j < yDim; j++) {
                for (int k = 0; k < zDim; k++) {
                    if (abs(map[i][j][k].E[0]) + abs(map[i][j][k].E[1]) + abs(map[i][j][k].E[2])
                        + abs(map[i][j][k].H[0]) + abs(map[i][j][k].H[1]) + abs(map[i][j][k].H[2]) > tolerance) {
                        outputFile << T << "," << i << "," << j << "," << k << ","
                            << map[i][j][k].E[0] << "," << map[i][j][k].E[1] << "," << map[i][j][k].E[2] << ","
                            << map[i][j][k].H[0] << "," << map[i][j][k].H[1] << "," << map[i][j][k].H[2] << "\n";
                    }
                }
            }
        }
        std::cout << "Bottom of for-loop" << endl;
    }

    outputFile.close();

    return 0;
}