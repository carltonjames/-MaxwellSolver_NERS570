/* Fd3d_4.1.cpp. 3D FDTD. Dipole in free space. */
/* Implemented as described in Electromagntic Simulation Using The FDTD Method by Dennis M. Sullivan */
#include <cmath>
#include <iostream>
#include <fstream>
#include <array>

using namespace std;

const int IE = 40;
const int JE = 40;
const int KE = 40;

template<typename T>
using Array3D = array<array<array<T, KE>, JE>, IE>;

int main() {
    Array3D<float> gax{};
    Array3D<float> gay{};
    Array3D<float> gaz{};
    Array3D<float> dx{};
    Array3D<float> dy{};
    Array3D<float> dz{};
    Array3D<float> ex{};
    Array3D<float> ey{};
    Array3D<float> ez{};
    Array3D<float> hx{};
    Array3D<float> hy{};
    Array3D<float> hz{};

    FILE* fp;
    int ic = IE / 2;
    int jc = JE / 2;
    int kc = KE / 2;
    float ddx = 0.01f; // Cell Size
    float dt = ddx / 6.0E8f; // Time Step
    float epsz = 8.8E-12f;
    float pi = 3.14159f;
    float tolerance = 1.0e-8f;
    char dipole_orientation = 'y';

    for (int i = 0; i < IE; ++i) {
        for (int j = 0; j < JE; ++j) {
            for (int k = 0; k < KE; ++k) {
                gax[i][j][k] = 1.0f;
                gay[i][j][k] = 1.0f;
                gaz[i][j][k] = 1.0f;
                dx[i][j][k] = 0.0f;
                dy[i][j][k] = 0.0f;
                dz[i][j][k] = 0.0f;
                ex[i][j][k] = 0.0f;
                ey[i][j][k] = 0.0f;
                ez[i][j][k] = 0.0f;
                hx[i][j][k] = 0.0f;
                hy[i][j][k] = 0.0f;
                hz[i][j][k] = 0.0f;
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
            for (int k = 1; k < KE; k++) {
                for (int j = 1; j < JE; j++) {
                    for (int i = 1; i < IE; i++) {
                        dx[i][j][k] = dx[i][j][k] + 0.5f * (hz[i][j][k] - hz[i][j - 1][k] - hy[i][j][k] + hy[i][j][k - 1]);
                        dy[i][j][k] = dy[i][j][k] + 0.5f * (hx[i][j][k] - hx[i][j][k - 1] - hz[i][j][k] + hz[i - 1][j][k]);
                        dz[i][j][k] = dz[i][j][k] + 0.5f * (hy[i][j][k] - hy[i - 1][j][k] - hx[i][j][k] + hx[i][j - 1][k]);
                    }
                }
            }

            // Calculate the E field from D field
            for (int k = 1; k < KE - 1; k++) {
                for (int j = 1; j < JE - 1; j++) {
                    for (int i = 1; i < IE - 1; i++) {
                        ex[i][j][k] = gax[i][j][k] * dx[i][j][k];
                        ey[i][j][k] = gay[i][j][k] * dy[i][j][k];
                        ez[i][j][k] = gaz[i][j][k] * dz[i][j][k];
                    }
                }
            }

            // Calculate the Hx, Hy, Hz fields
            for (int k = 1; k < KE - 1; k++) {
                for (int j = 1; j < JE - 1; j++) {
                    for (int i = 1; i < IE; i++) {
                        hx[i][j][k] = hx[i][j][k] + 0.5f * (ey[i][j][k + 1] - ey[i][j][k] - ez[i][j + 1][k] + ez[i][j][k]);
                        hy[i][j][k] = hy[i][j][k] + 0.5f * (ez[i + 1][j][k] - ez[i][j][k] - ex[i][j][k + 1] + ex[i][j][k]);
                        hz[i][j][k] = hz[i][j][k] + 0.5f * (ex[i][j + 1][k] - ex[i][j][k] - ey[i + 1][j][k] + ey[i][j][k]);
                    }
                }
            }


            // Source
            float pulse = exp(-0.5f * pow((t0 - T) / spread, 2.0f));
            if (dipole_orientation == 'x') {
                ex[ic][jc][kc] = pulse;
                dx[ic][jc][kc] = pulse;
            }
            else if (dipole_orientation == 'y') {
                ey[ic][jc][kc] = pulse;
                dy[ic][jc][kc] = pulse;
            }
            else if (dipole_orientation == 'z') {
                ez[ic][jc][kc] = pulse;
                dz[ic][jc][kc] = pulse;
            }


            // Absorbing Boundary Conditions for E fields
            for (int j = 0; j < JE; j++) {
                for (int i = 0; i < IE; i++) {
                    ex[i][j][0] = ex[i][j][1];
                    ex[i][j][KE - 1] = ex[i][j][KE - 2];
                    ey[i][0][j] = ey[i][1][j];
                    ey[i][JE - 1][j] = ey[i][JE - 2][j];
                    ez[0][i][j] = ez[1][i][j];
                    ez[IE - 1][i][j] = ez[IE - 2][i][j];
                }
            }

            // Absorbing Boundary Conditions for H fields
            for (int j = 0; j < JE; j++) {
                for (int i = 0; i < IE; i++) {
                    hx[i][j][0] = hx[i][j][1];
                    hx[i][j][KE - 1] = hx[i][j][KE - 2];
                    hy[i][0][j] = hy[i][1][j];
                    hy[i][JE - 1][j] = hy[i][JE - 2][j];
                    hz[0][i][j] = hz[1][i][j];
                    hz[IE - 1][i][j] = hz[IE - 2][i][j];
                }
            }

            // Output data to file
            for (int i = 0; i < IE; i++) {
                for (int j = 0; j < JE; j++) {
                    for (int k = 0; k < KE; k++) {
                        if (abs(ex[i][j][k] + ey[i][j][k] + ez[i][j][k]
                            + hx[i][j][k] + hy[i][j][k] + hz[i][j][k]) < tolerance) {
                            outputFile << T << "," << i << "," << j << "," << k << ","
                                << ex[i][j][k] << "," << ey[i][j][k] << "," << ez[i][j][k] << ","
                                << hx[i][j][k] << "," << hy[i][j][k] << "," << hz[i][j][k] << "\n";
                        }
                    }
                }
            }

        }
    }

    outputFile.close();
    return 0;
}
