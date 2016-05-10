//	\author Abhijit Bhattacharyya (abhijit.bhattacharyya@cern.ch)
#include <iostream>
#include <cstdio>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <ctime>
#include <climits>
#include <string>
#include <limits>
#include <algorithm>
#include <random>
#include <chrono>
#include <vector>
#include <Vc/Vc>
#include <Vc/IO>
#include <Vc/vector.h>
// #include <Vc/limits>
// #include <Vc/common/macros.h>

#define PI 2.0*asin(1.0)	// Pi
#define sqr(x) (x*x)
#define sqrt(x) (pow(x,0.5))
#define Max(x,y) ((x>y)?x:y)
#define Min(x,y) ((x<y)?x:y)

using namespace std;
using Vc::float_v;
using Vc::double_v;
using Vc::int_v;

struct particle {
	int p_type;		// particle type
	long p_id;		// particle id number
	float_v p_pos;		// particle position  3 vectors     [  x   y   z NULL]
	float_v p_vel;		// particle velocity  3 vectors     [ Vx  Vy  Vz NULL]
	float_v p_avvel;	// particle average velocity        [AVx AVy AVz NULL]
	float_v p_accln;	// particle acceleration 3 vectors  [ Ax  Ay  Az NULL]
	int countiac;		// number of neighboring particles
	double hsml;		// smoothing length h
	double p_density;	// particle density
    double p_dRhodt;    // Density change rate of each particle
	double p_mass;		// particle mass;
    double eta;         // coeff of viscosity
    float_v stress_1;   // stress elements with 1i  [rr     rphi      rz   NULL]
    float_v stress_2;   // stress elements with 2i  [phir   phiphi  phiz   NULL]
    float_v stress_3;   // stress elements with 3i  [zr     zphi      zz   NULL]
    float_v thermod;    // thermodynamical vars     [TdSdT  dedT       p   NULL]
};


//	define file pointers
ifstream iFile1;
ofstream oFile1;

//	define file names
string pathData		= "./";
string inFile1		= pathData+"SeismoIN.dat";
string inFile2		= pathData+"SeismoUR.dat";

//	define some global params
double Rdist		= 60.0;		// 60 km
double rMaxGeom		= Rdist;	// 60 km
double rMinGeom		= 1.0e-10;	// 1.0e-7 km
double phiMin		= 0.0;		// 0 degree
double phiMax		= 359.99;	// 360 degree
double zMaxGeom		= 40.0;		// 40 km inside elastic half space
double zMinGeom		= -20.0;	// -20 km in air
double maxIter		= 10e+6;	// number of iterations
double airVolume	= PI * sqr(rMaxGeom - rMinGeom) * abs(zMinGeom);	// in meter
double soilVolume	= PI * sqr(rMaxGeom - rMinGeom) * zMaxGeom;		// in meter
int approx;				// particle approximation for SPH algorithm
int nnps;				// nearest neighbor approx 1=direct 2=grid linked 3=tree algorithm
int sle;				// sm_len evol. 0:unchanged 1:h=fac*(m/rho)^(1/dim) 2:dhdt=(-1/dim)*(h/rho)*(droh/dt) 3:h=h0*(rho_0/rho)**(1/dim)
int skf;				// smoothing kernel 1=cubic splineby W4 2=gauss kernel 3=Quintic kernel

int summ_den, av_vel, config_input, virt_part, vp_input, visc, ex_force, visc_art, heat_art, self_grav, nor_den, symm_type;
int dim;					//  dimension 
long nVirt, totNum;					// total number of particles

