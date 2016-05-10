// \author Abhijit Bhattacharyya (abhijit.bhattacharyya@cern.ch)
#include "vegaSPH.hh"

//	Get info about Virtual Particles
void virtPart() {
	return;
}

//	Search on the basis of Linked List
void linkedList() {
	return;
}

//	Search on the basis of Tree search Algorithm
void treeSearch() {
	return;
}


//	Compute smoothing kernel W_ij for SPH
void sphKernel(double r, float_v &dist, double hsml, double &kernel, float_v &tdwdx) {
	double q, w, fac, tmp1, tmp2;
	int i;

	tdwdx = float_v::Zero();
	q = (r/hsml);
    cout << "q=" << q << endl;

	if (skf == 1) {
		fac = (dim == 1) ? 1.0/hsml : ((dim == 2) ? 15.0/(7.0*PI*sqr(hsml)) : 3.0/(2.0*PI*(sqr(hsml)*hsml)));
		if ((q >= 0.0) && (q <= 1.0)) {
			kernel = fac*((2.0/3.0)-sqr(q)+0.5*(sqr(q)*q));
			for (i = 0; i < float_v::Size; i++) tdwdx[i] = fac * (-2.0 +(3.0/2.0)*q)/sqr(hsml) * dist[i];
		} else if ((q > 1.0) && (q <= 2.0)) {
			kernel = fac * (1.0/6.0) * pow((2.0 - q), 3);
			for (i = 0; i < float_v::Size; i++) tdwdx[i] = -fac * (1.0/6.0) * 3.0 * sqr(2.0 - q)/hsml * dist[i]/r; 
		} else {
			kernel = 0.0;
			for (i = 0; i < float_v::Size; i++) tdwdx[i] = 0.0;
		}
	} else if (skf == 2) {
		fac = 1.0 / (pow(hsml, dim) * pow(PI,(0.5*dim)));
	       if ((q >= 0.0) && (q <= 3.0)) {
		       kernel = fac * exp(-q*q);
		       for (i = 0; i < float_v::Size; i++) tdwdx[i] = kernel * (-2.0 * dist[i]/sqr(hsml));
	       } else {
		       kernel = 0.0;
		       for(i = 0; i < float_v::Size; i++) tdwdx[i] = 0.0;
	       }
	} else if (skf ==3) {
		fac = (dim==1) ? 1.0/(120.0*hsml) : ((dim==2) ? 7.0/(478.0*PI*sqr(hsml)) : 1.0/(120.0*PI*sqr(hsml)*hsml));
		if ((q >= 0.0) && ( q <= 1.0)) {
			kernel = fac * (pow((3.0 - q), 5) - 6.0 * pow((2.0 - q), 5) + 15.0 * pow((1.0 - q), 5));
			for (i = 0; i < float_v::Size; i++) tdwdx[i] = fac * ((-120.0 + 120.0 * q - 50.0 * q * q)/sqr(hsml) * dist[i]);
		} else if ((q > 1.0) && (q <= 2.0)) {
			kernel = fac * (pow((3.0 - q), 5) - 6.0 * pow((2.0 - q), 5));
			for (i = 0; i < float_v::Size; i++) tdwdx[i] = (-5.0 * pow((3.0 - q), 4) + 30.0 * pow((2.0 - q), 4)/hsml) * dist[i] / r;
		} else if ((q > 2.0) && ( q <= 3.0)) {
			kernel = fac * pow((3.0 - q), 5);
			for (i = 0; i < float_v::Size; i++) tdwdx[i] = fac * (-5.0 * pow((3.0 - q), 4))/hsml * dist[i] / r;
		} else {
			kernel = 0.0;
			for (i = 0; i < float_v::Size; i++) tdwdx[i] = 0.0;
		}
	}
	return;
}


//	Direct Search to find nearest neighbor
void direct_find(vector<particle> &sphPar, long maxInterac, long &niac, vector<long> &pair_i, vector<long> &pair_j, vector<double> &kernel, vector<float_v> &dwdx) {
	int i, j, scale_k;
	float_v dist, tdwdx;
	double driac, mhsml, r, minr=1.0e+10, maxr=-minr;
	
	driac = 0.0;	niac = 0; 	dist = float_v::Zero();
	pair_i.resize(maxInterac); pair_j.resize(maxInterac);  dwdx.resize(maxInterac);

	scale_k = (skf==1) ? 2 : 3;
	for (i = 0; i < totNum - 1; i++) {
		for (j = i+1; j < totNum; j++) {
			dist = (sphPar[i].p_pos	- sphPar[j].p_pos);
			for (int k = 0; k < float_v::Size; k++)  driac += dist[k] * dist[k];
			mhsml = 0.5 * (sphPar[i].hsml + sphPar[j].hsml);
			r = sqrt(driac); maxr=Max(maxr,r); minr = Min(minr,r);
			if ((r < (scale_k * mhsml)) && (niac < maxInterac)){
				niac++;
				pair_i[niac] = i; pair_j[niac] = j;
				sphPar[i].countiac++;  sphPar[j].countiac++;
				sphKernel(r, dist, mhsml, kernel[niac], tdwdx);
				for (int l=0; l < float_v::Size; l++) {
					dwdx[niac][l] = tdwdx[l];  
                    cout << "niac = " << niac << " dwdx(niac)(l) = " << dwdx[niac][l] << endl;
				}
			}
		}
	}
    // cout << " Max r = " << maxr << " Min r = " << minr << endl;
	return;
}

//  Compute density with SPH continuity approach
void conDensity(vector<particle> &sphPar, long niac, vector<long> &pair_i, vector<long> &pair_j, vector<double> &kernel, vector<float_v> &dwdx) {
    int i, j, k, l;
    long total = totNum + nVirt;
    float_v dvx;
    double vcc;

    for (i = 0; i < total; i++) sphPar[i].p_dRhodt = 0.0;
    for (k = 0; k < niac; k++) {
        i = pair_i[k];      j = pair_j[k];
        for (l = 0; l < float_v::Size; l++) dvx[l] = sphPar[i].p_vel[l] - sphPar[j].p_vel[l];
        for (l = 0; l < float_v::Size; l++) vcc += dvx[l] * dwdx[niac][l];
        sphPar[i].p_dRhodt += sphPar[j].p_mass * vcc;
        sphPar[j].p_mass += sphPar[i].p_mass * vcc;
    }
    return;
}


//	Compute density with SPH summation algorithm
void sumDensity(vector<particle> &sphPar, long niac, vector<long> &pair_i, vector<long> &pair_j, vector<double> &kernel) {
	int i, j, k;
	long total = totNum+nVirt;
	double selfDens, r, mhsml;
	float_v dist = float_v::Zero();
	vector<double> wi;

	wi.resize(total);

	//	self density of each particle Wii i.e. kernel for distance 0 and take contribution of particle itself
	r = 0.0;

	//	Compute integration of kernel over the space
	for (i=0; i < total; i++) {
		mhsml = sphPar[i].hsml;
		sphKernel(r, dist, mhsml, selfDens, dist);
		wi[i] = selfDens*sphPar[i].p_mass/sphPar[i].p_density;
	}
	for (i=0; i < niac; i++) {
		j = pair_i[i];    k = pair_j[i];
		wi[j] += sphPar[k].p_mass/(sphPar[k].p_density*kernel[i]); 
		wi[k] += sphPar[j].p_mass/(sphPar[j].p_density*kernel[i]);
    }


	//	Compute density integration over the entire space
    for (i = 0; i < total; i++) {
        mhsml = sphPar[i].hsml;
        sphKernel(r, dist, mhsml, selfDens, dist);
        sphPar[i].p_density = selfDens * sphPar[i].p_mass;
    }
	for (i=0; i < niac; i++) {
		j = pair_i[i];    k = pair_j[i];
		wi[j] += sphPar[k].p_mass * kernel[i]; 
		wi[k] += sphPar[j].p_mass * kernel[i];
    }

    //  Compute normalized density
    if (nor_den) {
        for (i = 0; i < total; i++) sphPar[i].p_density /= kernel[i];
    }
	return;
}


//  Compute dynamic viscosity (mu) in N.s/m^2 e.g. water ~ 1.0e-3 N.s/m^2 ==  kinematic viscosity ~ 1.0e-6 m^2/s
void viscosity(vector<particle> &sphPar) {
    long i;

    for (i = 0; i < totNum; i++) {
        sphPar[i].eta = (sphPar[i].p_type == 0) ? 0.0 : ((sphPar[i].p_type ==1) ? 1.0e-3 : 1.0e+10);
    }
    return;
}

//  Compute internal force on the RHS of the Navier-Stokes i.e. pressure gradient, stress tensor etc.
void int_force(vector<particle> &sphPar, long niac, vector<long> &pair_i, vector<long> &pair_j, vector<float_v> &dwdx) {
    long i, j, k, l;
    float_v dvx;
    vector<float_v> htensor;

    htensor.resize(3);
    for (i = 0; i < 3; i++) htensor[i] = float_v::Zero();
    for (i = 0; i < totNum; i++) {
        sphPar[i].stress_1 = float_v::Zero();       sphPar[i].stress_2 = float_v::Zero();       sphPar[i].stress_3 = float_v::Zero();
        sphPar[i].thermod = float_v::Zero();        sphPar[i].p_accln = float_v::Zero();
    }

    if (visc) {
        for (k = 0; k < niac; k++) {
            i = pair_i[k];  j = pair_j[k];
            for (l = 0; l < float_v::Size; l++) dvx[l] = sphPar[j].p_vel[l] - sphPar[i].p_vel[l];
            if (dim ==1) {
                cout << "htensor "<< endl;
            }
        }
    }
    return;
}


//	Do Single step Computation
void stepSingle(vector<particle> &sphPar, long MaxInterac) {
	long nVirt = 0, niac = 0;
	vector<long> pair_i, pair_j;
	vector<double> kernel;
    vector<float_v> dwdx;

    kernel.resize(MaxInterac);

	if (virt_part) virtPart();				// if virtual particle is considered by user then determine their status

	//	interaction
	if (nnps==1) direct_find(sphPar, MaxInterac, niac, pair_i, pair_j, kernel, dwdx);		// nearest neighbor direct search
	if (nnps==2) linkedList();				// nearest neighbor grid Linked List
	if (nnps==3) treeSearch();				// nearest neighbor Tree Algorithm

    //  Density
	if (summ_den) {
        sumDensity(sphPar, niac, pair_i, pair_j, kernel);
    } else {
        conDensity(sphPar, niac, pair_i, pair_j, kernel, dwdx); 
    }

    //  Viscosity
    if (visc) viscosity(sphPar);               // for soil viscosity is very very high

    //  Internal force
    int_force(sphPar, niac, pair_i, pair_j, dwdx);
    
	return;
}


//	Calculate time Integration
void time_integ(long MaxInterac, double MaxTime, double dt, vector<particle> &sphPar) {
	long i, niter;

	niter = (long)(MaxTime/dt+0.5);
	for (i = 0; i < niter; i++) {
		stepSingle(sphPar, MaxInterac);
	}
	return;
}


//	Calculate artificial Heat
void art_heat(){
	return;
}

//	Read User Data
void ReadDat(void) {
	int nchr = 100;
	char line[nchr], *thestrptr;
	iFile1.open(inFile1.c_str(), ifstream::in);	// read the input file

	if(iFile1) {
		iFile1.getline(line, nchr);
		totNum				= (long)strtod(line, &thestrptr);

		iFile1.getline(line, nchr);
		summ_den			= strtod(line, &thestrptr);

		iFile1.getline(line, nchr);
		av_vel				= strtod(line, &thestrptr);

		iFile1.getline(line, nchr);
		config_input		= strtod(line, &thestrptr);

		iFile1.getline(line, nchr);
		virt_part			= strtod(line, &thestrptr);

		iFile1.getline(line, nchr);
		vp_input			= strtod(line, &thestrptr);

		iFile1.getline(line, nchr);
		visc				= strtod(line, &thestrptr);

		iFile1.getline(line, nchr);
		ex_force			= strtod(line, &thestrptr);

		iFile1.getline(line, nchr);
		visc_art			= strtod(line, &thestrptr);

		iFile1.getline(line, nchr);
		heat_art			= strtod(line, &thestrptr);

		iFile1.getline(line, nchr);
		self_grav			= strtod(line, &thestrptr);

		iFile1.getline(line, nchr);
		nor_den				= strtod(line, &thestrptr);

		iFile1.getline(line, nchr);
		symm_type			= strtod(line, &thestrptr);
	}
	iFile1.close();

	//	Display initial simulation algorithm conditions for the SPH
	cout << "Total Number of particles = " << totNum << endl;
	cout << "    Other options selected : " << endl;
	cout << "         Summation density = " << ((summ_den) ? " True <----" : " False ") << endl <<
		    "         Average Velocity = "  << ((av_vel) ? " True <----" : " False ")   << endl <<
			"         Virtual Particles = " << ((virt_part) ? " True <----" : " False ") << endl <<
			"         External Force = "    << ((ex_force) ? " True <----" : " False ") << endl <<
			"         Artificial Viscosity = " << ((visc_art) ? " True <----" : " False ") << endl <<
			"         Artificial heat = " << ((heat_art) ? " True <----" : " False ") << endl <<
			"         Self gravity = " << ((self_grav) ? " True <----" : " False ") << endl <<
			"         Normalized Density = " << ((nor_den) ? " True <----" : " False ") << endl <<
			"         Symmetry = " << ((!symm_type) ? " No Symmetry " : ((symm_type==1) ? " Axial Symmetry <----" : " Central Symmetry <----")) <<
			endl;
	return;
}


//	Initialize Particle parameters defined in the vector of structure :: sphPar
void defPar(vector<particle> &sphPar) {
	unsigned seed;
	long i;
	float temprnd;
	double volTake;

	dim = 3;					// 2D axially symm case cylindrical
	seed = chrono::system_clock::now().time_since_epoch().count(); default_random_engine generator (seed); srand(seed);

	for (i = 0; i < totNum; i++) {
		particle temp;				// define temporary struct for each particle
	    	temp.p_id   	= i;			// normal id i.e. roll no.
		temp.p_pos	= float_v::Zero();
		temp.p_vel	= float_v::Zero();	// initialize to Zero()
		temp.p_avvel	= float_v::Zero();
		temp.p_accln	= float_v::Zero();
		temprnd = rMinGeom + (static_cast<double>(rand())/static_cast<double>(RAND_MAX+1.0) * (rMaxGeom - rMinGeom)); temp.p_pos[0] = temprnd;
		temprnd = phiMin + (static_cast<double>(rand())/static_cast<double>(RAND_MAX+1.0) * (phiMax - phiMin)); temp.p_pos[1] = temprnd;
		temprnd = zMinGeom + (static_cast<double>(rand())/static_cast<double>(RAND_MAX+1.0) * (zMaxGeom - zMinGeom)); temp.p_pos[2] = temprnd;
		temprnd = temp.p_pos[2];
		temp.p_type 	= (temprnd<=0.0) ? 0 : 2;		// gas = 0, liquid = 1, solid = 2 :: ver.1.0
		temp.countiac	= 0;			// nearest neighbor initalize to zero
		temp.hsml 	= 2.0*(zMaxGeom-zMinGeom);			// smoothing length h
		temp.p_density	= (temprnd<=0.0) ? 1.225 : 2300.0;		// air:1.225 kg/m^3 soil:2300 kg/m^3
		volTake		= (temprnd<=0.0) ? airVolume : soilVolume;
		temp.p_mass 	= (temp.p_density*volTake)/(double)totNum; 	// in kg
		sphPar.push_back(temp);			// Pack data into the vector of particles
	}
	return;
}

//	Shutdown the Code
int shutDown(vector<particle> &sphPar) {
	sphPar.clear();
	// vector<particle>().swap(sphPar);
	return 0;
}


//	MAIN SPH  Simulation code for the simulation of 
//			detection of signal in air and under ground for 
//			wave propagation after initiation of wave
int main() {
	double dt, MaxTime;
       long maxInterac;

	//	Declare variables
	vector<particle> sphPar;	// declare particles
//	sphPar.resize(totNum);		// Initialize all particles
	dt = 0.5;			// Fix the time step dt
	approx = 1;  nnps = 1;  sle = 1;   skf = 1;  dim = 3;

	//	Initialize Particle parameters defined in the vector of particles
	ReadDat();		// Read user input and display SPH algorithm conditions
	defPar(sphPar);		// Initialize particle properties of each SPH particle

	// testing any particle property
	cout << " TESTING some parameters..............." << endl;
	cout << " particle id = " << sphPar[1].p_id << "  Position = " << sphPar[1].p_pos << " Velocity = " << sphPar[1].p_vel << endl;
	cout << endl;

	// 	Run Processes
	MaxTime = 60.0;    // 60 minutes rMaxGeom/330.0;
	MaxTime *= 60.0;
	MaxTime = 5.0;
	maxInterac = 100*totNum;
	cout << "Max time to simulate = " << MaxTime << " secs. " << endl << endl;
	time_integ(maxInterac, MaxTime, dt, sphPar);		// Start time integration


	//	Shut Down the Simulation
	if (!shutDown(sphPar)) {return EXIT_SUCCESS;} else {return EXIT_FAILURE;}
}

