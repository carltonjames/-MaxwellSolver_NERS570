#include "Mesh.hpp"
#include <cmath>

extern double delta;
extern float timeStep;
extern float c;
extern float e;
extern float m_e;

using namespace std;

Mesh::Mesh(std::array<int, 3> dimensions, int divs) {
    scale = static_cast<double>(dimensions[0]) / divs;
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

vector<vector<vector<Field>>>& Mesh::getMap() {
    return map;
}

Field Mesh::getField(std::array<int, 3> loc) const {
    return map[loc[0]][loc[1]][loc[2]];
}

void Mesh::setSources(std::vector<std::unique_ptr<Source>>& sourceList) {
    sources = std::move(sourceList);
}

void Mesh::addSource(std::unique_ptr<Source> source) {
    sources.push_back(std::move(source));
}

void Mesh::removeSource(size_t sourceID) {
    if (sourceID < sources.size()) {
        sources.erase(sources.begin() + sourceID);
    }
}

std::array<float, 3> Mesh::physicalLoc(std::array<int, 3> loc) const {
    return { static_cast<float>(loc[0] * scale),
             static_cast<float>(loc[1] * scale),
             static_cast<float>(loc[2] * scale) };
}

std::array<int, 3> Mesh::discreteLoc(std::array<float, 3> loc) const {
    return { static_cast<int>(loc[0] / scale),
             static_cast<int>(loc[1] / scale),
             static_cast<int>(loc[2] / scale) };
}
