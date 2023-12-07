#include <cmath>
#include <array>
#include <vector>
#include "Source.hpp"
#include "Field.hpp"

extern double delta;
extern float timeStep;
extern float c;
extern float e;
extern float m_e;

using namespace std;
class Mesh {
private:
    int n;
    std::array<float, 3> dimensions;
    double scale;
    std::vector<std::unique_ptr<Source>> sources; // list of sources within simulation
    vector<vector<vector<Field>>> map; // Defines EM field at all points in space

    void zeroMap() {
        
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

public:
    /* Mesh Constructor.*/
    Mesh(std::array<float, 3> dimensions, int divs) : dimensions(dimensions) {
        scale = dimensions[0] / divs; // Assuming uniform scaling for simplicity
        n = divs;
        X = static_cast<int>(dimensions[0]);
        Y = static_cast<int>(dimensions[1]);
        Z = static_cast<int>(dimensions[2]);

        map.resize(X);
        for (int i = 0; i < X; i++) {
            map[i].resize(Y);
            for (int j = 0; j < Y; j++) {
                map[i][j].resize(Z);
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

    vector<vector<vector<Field>>>& getMap() {
        return map;
    }

    Field getField(std::array<float, 3> loc) const {
        return map[loc[0]][loc[1]][loc[2]];
    }

    void setSources(std::vector<std::unique_ptr<Source>>& sourceList) {
        sources = std::move(sourceList);
    }

    void addSource(std::unique_ptr<Source> source) {
        sources.push_back(std::move(source));
    }

    void removeSource(size_t sourceID) {
        if (sourceID < sources.size()) {
            sources.erase(sources.begin() + sourceID);
        }
    }

    std::array<float, 3> physicalLoc(std::array<int, 3> loc) const {
        return { loc[0] * scale, loc[1] * scale, loc[2] * scale };
    }

    std::array<int, 3> discreteLoc(std::array<float, 3> loc) const {
        return { static_cast<int>(loc[0] / scale), static_cast<int>(loc[1] / scale), static_cast<int>(loc[2] / scale) };
    }
};
