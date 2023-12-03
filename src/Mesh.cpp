#ifndef MESH_H
#define MESH_H

#include <cmath>
#include <array>
#include <vector>
#include "Source.hpp"
#include "Field.hpp"

extern float delta;
extern float timeStep;

/* Structural changes:
* - removed flist (deprecated by Field struct)
* 
*/
using namespace std;
class Mesh {
private:
    int X, Y, Z, n;
    double scale;
    vector<Source> sources; // list of sources within simulation
    vector<vector<vector<Field>>> map; // Defines EM field at all points in space

    void Mesh::zeroMap() {
        for (int i; i < X; i++) {
            for (int j; j < Y; j++) {
                for (int k; k < Z; k++) {
                    for (int l; l < 3; l++) {
                        map[i][j][k].E[l] = 0;
                        map[i][j][k].B[l] = 0;
                    }
                }
            }
        }
    }

    void Mesh::calcE(Field &E, Source s, std::array<int, 3> at) {
        std::array<int, 3> r = at - s.getLocation();
        int rsquared = static_cast <int>(pow(r[0], 2) + pow(r[1], 2) + pow(r[2], 2));

        if (rsquared == 0) {
            return;
        }

        std::array<float, 3> contribution;
        std::array<float, 3> rhat = r / sqrt(static_cast <double>rsquared);
        float q = s.getCharge();
        contribution = (q / rsquared) * rhat;
        E += contribution;
    }

    void Mesh::calcB(Field &B, Source s, std::array<int, 3> at) {
        std::array<int, 3> r = at - s.getLocation();
        int rsquared = static_cast <int>(pow(r[0], 2) + pow(r[1], 2) + pow(r[2], 2));

        if (rsquared == 0) {
            return;
        }

        std::array<float, 3> contribution;
        std::array<int, 3> v = s.getVelocity();
        std::array<float, 3> rhat = r / sqrt(static_cast <double>(rsquared));
        std::array<int, 3> cross = v[0] * rhat[1] - v[0] * rhat[2]
            + v[1] * rhat[2] - v[1] * rhat[0]
            + v[2] * rhat[0] - v[2] * rhat[1];
        float q = s.getCharge();
        contribution = (q / rsquared) * cross;
        B += contribution;
    }

public:
    Mesh(float dimensions, int n)
    void calculateField();
    Field calculateField(int i, int j, int k);
    Field getField(int i, int j, int k) const;
    void propagateFields();
    void setSources(std::vector &sourceList);
    void addSource(Source &source);
    void removeSource(int sourceID);
    float* physicalLoc(int i, int j, int k) const;
    int* discreteLoc(float i, float j, float k) const;

    /* Mesh Constructor.*/
    Mesh(float dimensions, int divs) {
        scale = dimensions / divs;
        n = divs;
        X = dimensions;
        Y = dimensions;
        Z = dimensions;        

        // initialize the field map
        for (int i; i < X; i++) {
            for (int j; j < Y; j++) {
                for (int k; k < Z; k++) {
                    map[i][j][k] = new Field();
                }
            }
        }
    }

    ~Mesh() {
        for (int i; i < X; i++) {
            for (int j; j < Y; j++) {
                for (int k; k < Z; k++) {
                        delete [] map[i][j][k].E;
                        delete [] map[i][j][k].B;
                }
            }
        }
    }

    // Calculates the field at every point in the Mesh instantaneously.
    void Mesh::calculateField() {
        zeroMap();
        std::array<int, 3> at;

        for (Source s : sources) {
            if (!s.isDistribution) {
                for (int i = 0; i < X; i++) {
                    for (int j = 0; j < Y; j++) {
                        for (int k = 0; k < Z; k++) {
                            at = { i, j, k };
                            calcE(map[i][j][k].E, s, at);
                            calcB(map[i][j][k].B, s, at);
                        }
                    }
                }
            }
        }
    }

    Field Mesh::calculateField(int i, int j, int k) {
        std::array<int, 3> at = { i, j, k };
        map[i][j][k].E[0] = 0; map[i][j][k].E[1] = 0; map[i][j][k].E[2] = 0;
        map[i][j][k].B[0] = 0; map[i][j][k].B[1] = 0; map[i][j][k].B[2] = 0;
        for (Source s : sources) {
                calcE(map[i][j][k].E, s, at);
                calcB(map[i][j][k].B, s, at);
        }

        return map[i][j][k];
    }

    Field Mesh::getField(int i, int j, int k) const {
        return map[i][j][k];
    }

    // Uses Yee's method to advance fields forward in time.
    void Mesh::propagateFields() {
        const float deltaT = timeStep; // time discretization

        for (int i = 1; i < xDim - 1; ++i) {
            for (int j = 1; j < yDim - 1; ++j) {
                for (int k = 1; k < zDim - 1; ++k) {
                    field.E[i][j][k].x += (deltaT / delta) * (field.B[i][j][k].z - field.B[i][j - 1][k].z) - (deltaT / delta) * (field.B[i][j][k].y - field.B[i][j][k - 1].y);
                    field.E[i][j][k].y += (deltaT / delta) * (field.B[i][j][k].x - field.B[i][j][k - 1].x) - (deltaT / delta) * (field.B[i][j][k].z - field.B[i - 1][j][k].z);
                    field.E[i][j][k].z += (deltaT / delta) * (field.B[i][j][k].y - field.B[i - 1][j][k].y) - (deltaT / delta) * (field.B[i][j][k].x - field.B[i][j - 1][k].x);
                }
            }
        }

        for (int i = 0; i < xDim - 1; ++i) {
            for (int j = 0; j < yDim - 1; ++j) {
                for (int k = 0; k < zDim - 1; ++k) {
                    field.B[i][j][k].x += (deltaT / delta) * (field.E[i][j + 1][k].z - field.E[i][j][k].z) - (deltaT / delta) * (field.E[i][j][k + 1].y - field.E[i][j][k].y);
                    field.B[i][j][k].y += (deltaT / delta) * (field.E[i][j][k + 1].x - field.E[i][j][k].x) - (deltaT / delta) * (field.E[i + 1][j][k].z - field.E[i][j][k].z);
                    field.B[i][j][k].z += (deltaT / delta) * (field.E[i + 1][j][k].y - field.E[i][j][k].y) - (deltaT / delta) * (field.E[i][j + 1][k].x - field.E[i][j][k].x);
                }
            }
        }
    }

    void Mesh::setSources(std::vector &sourceList) {
        sources = sourceList;
    }

    void Mesh::addSource(Source &source) {
        sources.push_back(source);
    }

    void Mesh::removeSource(int sourceID) {
        delete sources[sourceID];
        sources.pop_back[sourceID];
    }

    float* Mesh::physicalLoc(int i, int j, int k) const {
        float* location [3] = { i * scale, j * scale, k * scale };
        return location;
    }

    int* Mesh::discreteLoc(float i, float j, float k) const {
        int x = static_cast <int> i;
        int y = static_cast <int> j;
        int z = static_cast <int> k;
        int* location [3] = { x / scale, y / scale, z / scale };
        return location;
    }
};

#endif // MESH_H
