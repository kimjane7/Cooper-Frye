//#include "hypersurface.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>

const double pi = 4.0*atan(1.0);

int main(int argc, char *argv[]){

	int nmax = 1000;
	double R, t, phi;
	double R0 = 10.0, v_R = 1.0;
	double tmax = R0/v_R, dt = tmax/nmax, dR = v_R*dt;

	// place holder for volume element
	std::vector<double> dOmega; 
	dOmega.resize(4);

	// list of volume elements for each time step
	std::vector<std::vector<std::vector<double>>> Omega;
	Omega.resize(nmax);

	// open file
	FILE *fptr;
	fptr = fopen("volume_elements.dat", "w");
	fprintf(fptr,"# t\tR\tx\ty\tdOmega[0]\tdOmega[1]\tdOmega[2]\tdOmega[3]\n");

	// make list of volume elements around circle
	for(int n = 0; n < nmax; n++) {

		t = (n+0.5)*dt;
		R = Rmax-v_R*t;
		v_R = -1.0;

		for(x = -Rmax+0.5*dR; x < Rmax; x+=dR) {

			for(y = -Rmax+0.5*dR; y < Rmax; y+=dR) {

				r = sqrt(x*x+y*y);
				phi = atan(y/x);
				tolerance = abs(0.5*dR/cos(phi));

				if(abs(r-R)<tolerance) {

					dOmega[0] = dphi*dR* //delta-phi Delta-r tau (think of tau as length)
					dOmega[1] = -dR*dR*dt*cos(phi)/v_R; //Delta Omega_x=cos(phi) *Delta-Omega_r = cos(phi) Omega_0/v
					dOmega[2] = -dR*dR*dt*sin(phi)/v_R;
					dOmega[3] = 0.0;

					Omega[n].push_back(dOmega);
					fprintf(fptr,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
						t, R, x, y, dOmega[0], dOmega[1], dOmega[2], dOmega[3]);

				}

			}
		}

	}

	fclose(fptr);

	return 0;
}

// DeltaOmega_r = (delT/delr)/(delT/delt)

//pick temp, mass, outward velocity, 
// u_r = gamma*v = umax R(tau)/R_0

// DONT FORGET ABOUT 2PI'S

// calculate number of particles per phi (here it won't depend on phi
// M, T, v, Rmax, umax
// for each Delta Omega, DeltaN = e^(-E/T) d^3p (p . Delta Omega) Theta(p dot Delta Omega/ (R*(2pihbar)^3))
// integral over 2D p, t, phi