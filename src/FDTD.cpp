#include "FDTD.hpp"
#include "globals.hpp"

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

/* Implemented as described in Electromagntic Simulation Using The FDTD Method by Dennis M. Sullivan */
int FDTD::runSimulation() {
    int ic = xDim / 2; // center of x
    int jc = yDim / 2; // center of y
    int kc = zDim / 2; // center of z
    float ddx = 0.01f; // Cell Size
    float dt = ddx / 6.0E8f; // Time Step
    float epsz = 8.8E-12f;
    float pi = 3.14159f;
    float tolerance = 1.0e-8f;
    char dipole_orientation = 'y';

    array<float, 3> center = { static_cast<float>(ic),
                                static_cast<float>(jc),
                                static_cast<float>(kc) };

    vector<vector<vector<Field>>>& field = mesh.getMap();

    float gax[xDim][yDim][zDim];
    float gay[xDim][yDim][zDim];
    float gaz[xDim][yDim][zDim];

    vector<vector<vector<Field>>> map = mesh.getMap();

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
    float spread = 6.0f * 100.0f;
    int T = 0;
    int nsteps = 100;

    ofstream outputFile("fdtd_output.csv");
    outputFile << "T,i,j,k,Ex,Ey,Ez,Hx,Hy,Hz\n";

    while (nsteps > 0) {
        cout << "nsteps --> " << nsteps;
        //cin >> nsteps;
        //cout << nsteps << endl;

        for (int n = 1; n <= nsteps; n++) {
            T++;

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

            // Calculate the Hx, Hy, Hz fields
            for (int k = 1; k < zDim - 1; k++) {
                for (int j = 1; j < yDim - 1; j++) {
                    for (int i = 1; i < xDim; i++) {
                        map[i][j][k].H = { map[i][j][k].H[0] + 0.5f * (map[i][j][k + 1].E[1] - map[i][j][k].E[1] - map[i][j + 1][k].E[2] + map[i][j][k].E[2]),
                                            map[i][j][k].H[1] + 0.5f * (map[i + 1][j][k].E[2] - map[i][j][k].E[2] - map[i][j][k + 1].E[0] + map[i][j][k].E[0]),
                                            map[i][j][k].H[2] + 0.5f * (map[i][j + 1][k].E[0] - map[i][j][k].E[0] - map[i + 1][j][k].E[1] + map[i][j][k].E[1]) };
                        map[i][j][k].B = { map[i][j][k].H[0],
                                            map[i][j][k].H[1],
                                            map[i][j][k].H[2] };
                    }
                }
            }


            // Source
            float pulse = exp(-0.5f * pow((t0 - T) / spread, 2.0f));
            Field field{};
            for (int i = 0; i < sources.size(); i++) {
                Field temp{};
                array<int, 3> loc = { ic, jc, kc };
                std::array<float, 3> physicalLoc = mesh.physicalLoc(loc);
                temp = calcField(*sources[i], physicalLoc);
                field.E[0] += temp.E[0]; field.E[1] += temp.E[1]; field.E[2] += temp.E[2];
                field.B[0] += temp.B[0]; field.B[1] += temp.B[1]; field.B[2] += temp.B[2];
            }
            if (dipole_orientation == 'x') {
                map[ic][jc][kc].E[0] = pulse + field.E[0];
                map[ic][jc][kc].D[0] = pulse + field.E[0];
            }
            else if (dipole_orientation == 'y') {
                map[ic][jc][kc].E[1] = pulse + field.E[1];
                map[ic][jc][kc].D[1] = pulse + field.E[1];
            }
            else if (dipole_orientation == 'z') {
                map[ic][jc][kc].E[2] = pulse + field.E[2];
                map[ic][jc][kc].D[2] = pulse + field.E[2];
            }

            // Absorbing Boundary Conditions for E fields
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
        nsteps = 0;
    }

    outputFile.close();

    return 0;
}

   
// Calculates the energy stored inthe electric field for the entire mesh.
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

// source fields should also propagate
void FDTD::progressFreeParticles() {
    for (auto& source : sources) {
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
            location[0] += velocity[0] * timeStep + 0.5f * (netForce[0] / mass) * pow(timeStep, 2);
            location[1] += velocity[1] * timeStep + 0.5f * (netForce[1] / mass) * pow(timeStep, 2);
            location[2] += velocity[2] * timeStep + 0.5f * (netForce[2] / mass) * pow(timeStep, 2);

            velocity[0] += (netForce[0] / mass) * timeStep;
            velocity[1] += (netForce[1] / mass) * timeStep;
            velocity[2] += (netForce[2] / mass) * timeStep;

            source->setLocation(location);
            source->setVelocity(velocity);
        }
        else if (source->isDynamic()) {
            std::array<float, 3> location = source->getLocation();
            std::array<float, 3> velocity = source->getVelocity();

            // Update position according to velocity
            location[0] += velocity[0] * timeStep;
            location[1] += velocity[1] * timeStep;
            location[2] += velocity[2] * timeStep;

            source->setLocation(location);
        }
    }
}

void FDTD::outputSourceData(int T) {
    ofstream outputFile("source_output.csv");
    if (T == 0) {
        outputFile << "T,i,j,k,id,source\n";
    }
    std::array<float, 3> loc;
    for (int i = 0; i < sources.size(); i++) {
        loc = sources[i]->getLocation();
        outputFile << loc[0] << "," 
            << loc[1] << "," << loc[2] << ","
            << i << "," << sources[i]->getName();
    }
    outputFile.close();
}

std::array<float, 3> FDTD::crossProductFloat(std::array<float, 3> v1, std::array<float, 3> v2) {
    std::array<float, 3> cross = { v1[1] * v2[2] - v1[2] * v2[1],
                                    v1[2] * v2[0] - v1[0] * v2[2],
                                    v1[0] * v2[1] - v1[1] * v2[0] };
    return cross;
}

std::array<int, 3> FDTD::crossProductInt(std::array<int, 3> v1, std::array<int, 3> v2) {
    std::array<int, 3> cross = { v1[1] * v2[2] - v1[2] * v2[1],
                                    v1[2] * v2[0] - v1[0] * v2[2],
                                    v1[0] * v2[1] - v1[1] * v2[0] };
    return cross;
}

std::array<float, 3> FDTD::subtract(std::array<float, 3> v1, std::array<float, 3> v2) {
    return { v2[0] - v1[0], v2[1] - v1[1] , v2[2] - v1[2] };
}

float FDTD::norm(std::array<float, 3> v) {
    return sqrtf(powf(v[0], 2) + powf(v[1], 2) + powf(v[2], 2));
}

std::array<float, 3> FDTD::calcForce(Source& source1, Source& source2) {
    Field F;
    std::array<float, 3> r, rhat, velocity, force;
    r = subtract(source1.getLocation(), source2.getLocation());

    for (int i = 0; i < 3; ++i) {
        rhat[i] = r[i] / norm(r);
    }

    velocity = subtract(source1.getVelocity(), source2.getVelocity());

    float chargeFactor = source2.getCharge() / powf(norm(r), 2);
    for (int i = 0; i < 3; ++i) {
        F.E[i] = rhat[i] * chargeFactor;
    }

    std::array<float, 3> crossProd = crossProductFloat(source2.getVelocity(), rhat);
    float factor = (1 / c) / powf(norm(r), 2);
    for (int i = 0; i < 3; ++i) {
        F.B[i] = crossProd[i] * factor;
    }

    crossProd = crossProductFloat(velocity, F.B);
    float charge = source1.getCharge();
    for (int i = 0; i < 3; ++i) {
        force[i] = charge * (F.E[i] + (1 / c) * crossProd[i]);
    }

    return force;
}

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
