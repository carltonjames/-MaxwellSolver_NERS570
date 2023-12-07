/* Fd3d_4.1.c. 3D FDTD. Dipole in free space. */
/* Implemented as described in Electromagntic Simulation Using The FDTD Method by Dennis M. Sullivan */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define IE 40
#define JE 40
#define KE 40

int main() {
	float gax[IE][JE][KE], gay[IE][JE][KE], gaz[IE][JE][KE];
	float dx[IE][JE][KE], dy[IE][JE][KE], dz[IE][JE][KE];
	float ex[IE][JE][KE], ey[IE][JE][KE], ez[IE][JE][KE];
	float hx[IE][JE][KE], hy[IE][JE][KE], hz[IE][JE][KE];
	int l, n, i, j, k, ic, jc, kc, nsteps, npml;
	float ddx, dt, T, epsz, pi, epsilon, sigma, eaf;
	float xn, xxn, xnum, xd, curl_e;
	float t0, spread, pulse;
	char dipole_orientation;
	FILE* fp, * fopen();

	ic = IE / 2;
	jc = JE / 2;
	kc = KE / 2;
	ddx = 0.01; /* Cell Size*/
	dt = ddx / 6.0E8; /* Time Steps */
	epsz = 8.8E-12;
	pi = 3.14159;
	dipole_orientation = 'y';

	/* Initialize the arrays */
	for ( k=0; k < KE; k++ ) {
		for ( j=0; j < JE; j++ ) {
			for ( i=0; i < IE; i++ ) {
				ex[i][j][k] = 0.0;
				ey[i][j][k] = 0.0;
				ez[i][j][k] = 0.0;
				dx[i][j][k] = 0.0;
				dy[i][j][k] = 0.0;
				dz[i][j][k] = 0.0;
				hx[i][j][k] = 0.0;
				hy[i][j][k] = 0.0;
				hz[i][j][k] = 0.0;
				gax[i][j][k] = 1.0;
				gay[i][j][k] = 1.0;
				gaz[i][j][k] = 1.0;
			}
		} 
	}

	for (k = 11; k < 30; k++) {
		gaz[ic][jc][k] = 0.0;
	}
	gaz[ic][jc][kc] = 0.0;

	t0 = 0.0;
	spread = 6.0;
	T = 0;
	nsteps = 1;

	while (nsteps > 0) {
		printf("nsteps --> ");
		scanf("%d", &nsteps);
		printf("%d \n", nsteps);
		for (n = 1; n <= nsteps; n++) {
			T = T + 1;
			/*Start of the Main FDTD loop*/
			/* Calculate the Dx field */
			for (k = 1; k < KE; k++) {
				for (j = 1; j < JE; j++) {
					for (i = 1; i < IE; i++) {
						dx[i][j][k] = dx[i][j][k]
							+ .5 * (hz[i][j][k] - hz[i][j - 1][k]
								- hy[i][j][k] + hy[i][j][k - 1]);
					}
				}
			}

			/* Calculate the Dy field */
			for (k = 1; k < KE; k++) {
				for (j = 1; j < JE; j++) {
					for (i = 1; i < IE; i++) {
						dy[i][j][k] = dy[i][j][k]
							+ .5 * (hx[i][j][k] - hx[i][j][k - 1]
								- hz[i][j][k] + hz[i - 1][j][k]);
					}
				}
			}

			/* Calculate the Dz field */
			for (k = 1; k < KE; k++) {
				for (j = 1; j < JE; j++) {
					for (i = 1; i < IE; i++) {
						dz[i][j][k] = dz[i][j][k]
							+ .5 * (hy[i][j][k] - hy[i - 1][j][k]
								- hx[i][j][k] + hx[i][j - 1][k]);
					}
				}
			}

			/* Source */
			pulse = exp(-0.5 * pow((t0 - T) / (spread), 2.0));
			if (dipole_orientation == 'x') {
				ex[ic][jc][kc] = pulse;
			}
			else if (dipole_orientation == 'y') {
				ey[ic][jc][kc] = pulse;
			}
			else if (dipole_orientation == 'z') {
				ez[ic][jc][kc] = pulse;
			}

			/* Calculate the E from D field */
			for (k = 1; k < KE - 1; k++) {
				for (j = 1; j < JE - 1; j++) {
					for (i = 1; i < IE - 1; i++) {
						ex[i][j][k] = gax[i][j][k] * dx[i][j][k];
						ey[i][j][k] = gay[i][j][k] * dy[i][j][k];
						ez[i][j][k] = gaz[i][j][k] * dz[i][j][k];
					}
				}
			}

			/* Calculate the Hx field*/
			for (k = 1; k < KE - 1; k++) {
				for (j = 1; j < JE - 1; j++) {
					for (i = 1; i < IE; i++) {
						hx[i][j][k] = hx[i][j][k]
							+ .5 * (ey[i][j][k + 1] - ey[i][j][k]
								- ez[i][j + 1][k] + ez[i][j][k]);
					}
				}
			}

			/* Calculate the Hy field*/
			for (k = 1; k < KE - 1; k++) {
				for (j = 1; j < JE; j++) {
					for (i = 1; i < IE - 1; i++) {
						hy[i][j][k] = hx[i][j][k]
							+ .5 * (ez[i + 1][j][k] - ez[i][j][k]
								- ex[i][j][k + 1] + ex[i][j][k]);
					}
				}
			}

			/* Calculate the Hz field*/
			for (k = 1; k < KE; k++) {
				for (j = 1; j < JE - 1; j++) {
					for (i = 1; i < IE - 1; i++) {
						hz[i][j][k] = hx[i][j][k]
							+ .5 * (ex[i][j + 1][k] - ex[i][j][k]
								- ey[i + 1][j][k] + ey[i][j][k]);
					}
				}
			}

			/* End of the main FDTD loop */

			/*printf("HZ \n");
			for (k = 1; k < KE; k++) {
				printf("%2d ", k);
				for (i = 1; i < IE - 1; i++) {
					printf("%6.3f", hz[i][jc][k]);
				}
				printf(" \n");
			}*/

			/* Write the E field out to a file "Hz" */
			fp = fopen("HZ", "a");
			float num = 0;
			for (j = 0; j < JE; j++) {
				for (i = 0; i < IE; i++) {
					num = hz[i][j][k];
					if (num < 0.000000001) {
						fprintf(fp, "%1.7f, ", num);
					}
					else {
						fprintf(fp, "         , ");
					}
				}
				fprintf(fp, " \n");
			}
			fprintf(fp, "\n");
			fclose(fp);
		}
	}
}