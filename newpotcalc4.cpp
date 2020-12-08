#include<iostream>
#include<fstream>
#include<string>
#include<functional>
#include<cmath>
#include<random>
#include<thread>
#include<mutex>
#include<algorithm>
#include"./GalPot-master/src/GalPot.h"


#include <math.h>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_spline.h>

#define DIMENSION 3
#define THREADNO 6


#define PI 3.14159265358979323846264338328
#define parsec 3.08567759756e+16
#define VC 239.1
#define RSUN 8.29

#define NBASE   500
#define NPERBASE  10
#define TIMESTEPFAC 0.005        //0.0002   //the time-step factor
#define MINSTEP 0.2  //the minimum step
#define MAXSTEP 200.0



using namespace std;
using std::setprecision;

// parameters of normalization array
#define CVALN 40
#define GVALN 27
#define ZVALN 23                //50
#define CVL 0.05
#define GVL 0.8
#define ZVL 9.0


#define VPHIPLACE 0.215
#define VPHIPLACESIG 0.045
#define VRPLACEFAC 1.6
#define VZPLACEFAC 2.0
#define RPLACEL 4.0
#define ZPLACEL 0.6
#define RMIN 6.79
#define RMAX 9.79
#define CONTIME 10000.0 

#define VCG 0.235 //in 1000km/s
#define RCR 6.0 //in kpc

// real and imag for fft
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])


#define FROMFILE1 "Halo1.Tpot"
#define FROMFILE2 "Disc1.Tpot"




//global variable problem, instantiating Potential
std::ifstream fromfile("./GalPot-master/pot/PJM17_best.Tpot");
GalaxyPotential PhiGalPot(fromfile);
/*
std::ifstream fromGalPot("./barPJM16_best.Tpot");
GalaxyPotential PhiGalPotBar(fromGalPot);

std::ifstream fromBarsmall("./barsmalla2212Normal.Tpot");
GalaxyPotential PhiBarSmall(fromBarsmall);

std::ifstream fromBarlarge("./barsmalla2212Normal.Tpot");
GalaxyPotential PhiBarLarge(fromBarlarge);

std::ifstream fromBarsmallgrow("./barsmalla2212Normal.Tpot");
GalaxyPotential PhiBarSmallGrow(fromBarsmallgrow);

std::ifstream fromBarlargegrow("./barsmalla2212Normal.Tpot");
GalaxyPotential PhiBarLargeGrow(fromBarlargegrow);
*/


std::mutex mutex_var;

std::random_device random_var;
std::mt19937 randomGen(random_var());


const double delta = 10e-7;
const double divDelta = 1.0 / delta;


double timefunc(double time) {
    return time/CONTIME;//1.0 by Simon ;
}

//function for returning the azimuthal angle of the bar
double barphi(double time){
    double Omegabar = VCG/RCR; //in Myr //angular velocity of bar
    return Omegabar*time; //gives the azimuthal angle of bar at the specific time: phi_b
}

//function returning the omega_bar
double barOm(double time){
    const double tinytime = 0.01; //10 thousand years //remember 1 unit is 1 Myr
    //note to me: in the thesis, write down actual timelengths of typical dynamical motions of the milky way
    
    //gives the angular velocity omega in terms of tinytime (and angular displacement)
    return (barphi(time+tinytime) - barphi(time))/tinytime;
    }

//so that amplitude grows with time
double baramp(double A, double time){
//    double growtime = 500.0; //the growth time in units of Myr, can set however i want
    //to slowly grow the amplitude, mimicking a growing bar
    //choose any type of growing function i want
    return 0.0; //by Simon ;      //*atan(time*time/(growtime*growtime)); //*2.0/PI;
}


void fPot(double coord[3], double grad[], GalaxyPotential &pH) {
	double R = std::sqrt(coord[0] * coord[0] + coord[1] * coord[1]);
	double dR; double dz;
	pH(R, coord[2], dR, dz);
		grad[0] = -dR * coord[0] / R;
		grad[1] = -dR * coord[1] / R;
		grad[2] = -dz;
}


void fPot(double coord[3], double grad[], GalaxyPotential &pH, GalaxyPotential &pD, double time) {
	double R = std::sqrt(coord[0] * coord[0] + coord[1] * coord[1]);
	double phi = std::atan2(coord[1], coord[0]);
	double dR; double dz;
    double f2 = timefunc(time);
	pH(R, coord[2], dR, dz);
    double dR2; double dz2;
    pD(R, coord[2], dR2, dz2); 
    grad[0] = -(dR + f2*dR2)*coord[0]/R;
	grad[1] = -(dR + f2*dR2)*coord[1]/R;
	grad[2] = -(dz + f2*dz2);
    //double pot_unper = (pH(R,coord[2]) + f2*pD(R,coord[2])); 
    
/*    const double A = 0.017; const double b =  0.28; const double vc = VCG; const double m = 2.0; const double Rcr = RCR;    //gives the azimuthal angle of bar //in rads
    double phib = barphi(time);
    //gives the angle difference between bar and coordinate frame //in rads
    double deltaphi = phi - phib;  
    double Phim = - (baramp(A,time)*(vc*vc/m)) * (R/Rcr)*(R/Rcr) * pow(((b + 1.0)/(b + R/Rcr)),5.0);
    double potper = Phim*cos(deltaphi*m);
    double gradfR = -potper * ((2.0/R) - (5.0)*((b + R/Rcr)/(b + 1.0))*(1.0/Rcr)*(b+1.0)/((b + R/Rcr)*(b+R/Rcr))); //         (1.0/(b+(R/R_cr))));
    double gradfphi = Phim * (m/R) * sin(deltaphi*m);
    double gradfx = (gradfR *coord[0]/R) - gradfphi * coord[1]/R; // x = R*cos(phi) - phi*sin(phi)
    double gradfy = (gradfR *coord[1]/R) + gradfphi * coord[0]/R; // y = R*sin(phi) + phi*cos(phi)
    grad[0] += gradfx;
    grad[1] += gradfy;
*/
//    grad[2] += gradf_z;
}

double gPot(double coord[3], GalaxyPotential &pH1, GalaxyPotential &pH2, double time) {
    double R = std::sqrt(coord[0] * coord[0] + coord[1] * coord[1]);
    double phi = std::atan2(coord[1], coord[0]); 
    double f2 = timefunc(time);
    const double A = 0.017; const double b =  0.28; const double vc = VCG; const double m = 2.0; const double Rcr = RCR; //in kpc    
    double phib = barphi(time);
    double deltaphi = phi - phib;
    double Phim = - (baramp(A,time)*(vc*vc/m)) * (R/Rcr)*(R/Rcr) * pow(((b + 1.0)/(b + R/Rcr)),5.0);
    double potunper = (pH1(R,coord[2]) + f2*pH2(R,coord[2]));
    double potper = Phim* cos(deltaphi*m);
    return (potper + potunper);
}


class particle {
public:
    int stepnumber = 0;
    double weight = 1.0;
    double timestep = 1.0;
    double htimestep = 0.5;
    double time = 0.0;
    double timefine = 0.0;
    double x[DIMENSION];
    double p[DIMENSION];
    double f[DIMENSION];
    void drift();
    double getR();
    double getz();
    double getr();
    double gettime();
    void settime();
    void timefineadd();
    double getE(GalaxyPotential &pH);
    void kick(GalaxyPotential &pH);    
    void step(GalaxyPotential &pH);
    void zstep(GalaxyPotential &pH, GalaxyPotential &pD);
    void zkick(GalaxyPotential &pH, GalaxyPotential &pD);
    void zdrift();
    void zleapfrog(GalaxyPotential &pH, GalaxyPotential &pD);
    void getJz(GalaxyPotential &pH, GalaxyPotential &pH2, double JrJz[2]);
    void leapfrog(GalaxyPotential &pH);
    void print(GalaxyPotential &pH);
    void print(GalaxyPotential &pH, ostringstream &ostr);
//    void print(GalaxyPotential &pH1, GalaxyPotential &pH2, ostringstream &ostr);
    void getsimpleperiods(GalaxyPotential &pH,double periods[2], double maxtime);
    void integrate(double integtime, GalaxyPotential &pH);
    void kick(GalaxyPotential &pH, GalaxyPotential &pD);    
    void step(GalaxyPotential &pH, GalaxyPotential &pD);
    void leapfrog(GalaxyPotential &pH, GalaxyPotential &pD);
    void print(GalaxyPotential &pH, GalaxyPotential &pD);
    void print(GalaxyPotential &pH, GalaxyPotential &pD, ostringstream &ostr);
//  void Nprint(GalaxyPotential &pH, GalaxyPotential &pD, ostringstream &ostr);
//    void getsimplestartperiods(GalaxyPotential &pH, GalaxyPotential &pD, periods[2], double maxtime);
    void integrate(double integtime, GalaxyPotential &pH, GalaxyPotential &pD);    
    void phishift(std::uniform_real_distribution<double> phasedist);
    void phishift(double delphi);
    void phiset(double R, double phi, double vR, double vphi);
    void phitrans(double &R, double &phi, double &vR, double &vphi);
    void clone(particle &pp);
    void fourier(GalaxyPotential &pH1, GalaxyPotential &pH2, ostringstream &ostr1, ostringstream &ostr2, ostringstream &ostr3, ostringstream &ostr22);
    void gaussianvelocity(GalaxyPotential &pH1, GalaxyPotential &pH2, ostringstream &ostr,double steps);
private:
    ;
};

double particle::getz() {
    return x[2];
}

void particle::getJz(GalaxyPotential &pH1, GalaxyPotential &pH2, double JRJz[2]) {
    JRJz[0] = 0.0; JRJz[1] = 0.0;
    double z0 = x[2];
    double z = x[2];
    double zold = z0 - p[2];
    double zold2 = zold - p[2];
    double delz = z- zold;
    double timestepsave = timestep;
    for (int i = 0; i < 2; i++) {
      zleapfrog(pH1, pH2);
        zold2 = zold;
        z = x[2];
        delz = z - zold;
        JRJz[1] += delz*p[2];
        zold = z;
//        cout << "1 " << " " << htimestep << " " << zold << " " << zold2 << " " << x[2] << " " << p[2] << "  " << JRJz[1] << "\n";
    }        
    while (((zold - z0)*(zold2 - z0)) >= 0.0) {
        zleapfrog(pH1, pH2);
        zold2 = zold;
        z = x[2];
        delz = z - zold;
        JRJz[1] += delz*p[2];
        zold = z;
//        cout << "1 " << " " << htimestep << " " << zold << " " << zold2 << " " << x[2] << " " << p[2] << "  " << JRJz[1] << "\n";
    }
    for (int i = 0; i < 3; i++) {
       zleapfrog(pH1, pH2);
        zold2 = zold;
        z = x[2];
        delz = z - zold;
        JRJz[1] += delz*p[2];
        zold = z;
//        cout << "2 " << " " << x[2] << " " << p[2] << "  " << JRJz[1] << "\n";
    }        
    while (((zold - z0)*(zold2 - z0)) >= 0.0) {
        zleapfrog(pH1, pH2);
        zold2 = zold;
        z = x[2];
        delz = z - zold;
        JRJz[1] += delz*p[2];
        zold = z;
//        cout << "3 " <<  " " << x[2] << " " << p[2] << "  " << JRJz[1] << "\n";
    }
    timestep = timestepsave;
}

void particle::zleapfrog(GalaxyPotential &pH, GalaxyPotential &pD) {
	double grad[3];
	for (int i = 0; i < 3; ++i) {
		grad[i] = 0.0;
	}
	fPot(x, grad, pH, pD, gettime());
	double acc = std::sqrt(grad[2] * grad[2]);
    double df = (f[2] - grad[2])*(f[2] - grad[2]);
    df = 10.0*sqrt(df)/acc;
    df = 1.0/(df + 1.0e-14);
//    cout << 0.0001 << " " << df << " " << sqrt(q2)/acc << " " << 1.0/acc << "\n";
    if (1.0/acc < df) {df = 1.0/acc;} 
//        else {
//            cout << "dfrule " << df << " " << 1.0/acc << "\n";
//        }
    if (df > MINSTEP) {
//        cout << "WARNING " << df << " " << x[0] << " " << x[1] << " " << x[2] << " " << 1.0/acc << "\n";
        df = MINSTEP;
    }
    double tnew = TIMESTEPFAC*df;
    if ((tnew > timestep*2.0) || (tnew < htimestep)){
        if (stepnumber > 2) {
            if (tnew < htimestep) {
//                cout << "WARNING tnew " << tnew << " " << timestep << " " << stepnumber << "\n";
            }
        }
        if (tnew > timestep) {
            tnew = timestep*1.2;
        }
    }
	timestep = tnew;    
    htimestep = 0.5*timestep;
    for (int j = 0; j < 4; j++) {
        zstep(pH, pD);
    }
}

void particle::timefineadd() {
    if (timefine > 0.00001*time) {
        time += timefine; //adding timefine to time
        timefine = 0.0;
    }
}

double particle::getR() {
    return sqrt(x[0]*x[0] + x[1]*x[1]);
}

double particle::getr() {
    return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}


double particle::gettime() {
    return (time + timefine);
}

void particle::settime() {
    time = 0.0;
    timefine = 0.0;
    stepnumber = 0;
}

void particle::kick(GalaxyPotential &pH) {
    fPot(x, f, pH);
    for (int i = 0; i < DIMENSION; i++) {
        p[i] += timestep*f[i];
    }
}    

void particle::kick(GalaxyPotential &pH, GalaxyPotential &pD) {
    fPot(x, f, pH, pD, gettime());
    for (int i = 0; i < DIMENSION; i++) {
        p[i] += timestep*f[i];
    }
}    

void particle::zkick(GalaxyPotential &pH, GalaxyPotential &pD) {
    fPot(x, f, pH, pD, time);
        p[2] += timestep*f[2];
}    


void particle::zdrift() {
        x[2] += htimestep*p[2];
}


void particle::drift() {
    for (int i = 0; i < DIMENSION; i++) {
        x[i] += htimestep*p[i];
    }
}

void particle::step(GalaxyPotential &pH) {
    drift();
    timefine += htimestep;
    kick(pH);
    drift();
    timefine += htimestep;
}

void particle::zstep(GalaxyPotential &pH, GalaxyPotential &pD) {
    zdrift();
    timefine += htimestep;
    zkick(pH, pD);
    zdrift();
    timefine += htimestep;
}


void particle::step(GalaxyPotential &pH, GalaxyPotential &pD) {
    drift();
    timefine += htimestep;
    kick(pH, pD);
    drift();
    timefine += htimestep;
}

void particle::leapfrog(GalaxyPotential &pH) {
    stepnumber++;
	double grad[3];
	for (int i = 0; i < 3; ++i) {
		grad[i] = 0.0;
	}
    
        //difference here is the fPot taking 5 arguments
    //taking in the coordinates, force, pH, pD, and gettime
	fPot(x, grad, pH);
        //the total acceleration of all 3 space-dimension
	double acc = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);
    double sacc = sqrt(acc);
    
    double q2 = 0.0;//.............................................................
    double df = 0.0;
    for (int j = 0; j < DIMENSION; j++) { //for each coordinates
        q2 += x[j]*x[j]; //position^2
        df += (f[j] - grad[j])*(f[j] - grad[j]);//force difference ^2
    }//.............................................................................
    df = sqrt(df/acc); //some sort of normalisation
    double tssug = (1.0/(10.0*df + 1.0e-14 + sacc + (1.0/MAXSTEP))) + MINSTEP;

    double tnew = TIMESTEPFAC*tssug;
    if ((tnew > timestep*2.0) || (tnew < htimestep)){ // if ... OR ...
        if (stepnumber > 2) {
            if (tnew < htimestep) {
                double saccold = sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);
                cout << "WARNING tnew " << tnew << " " << timestep << " " << tssug << " " << 1.0/sacc << " " << 1.0/saccold << " " << 1.0/(10.0*df) << " " << stepnumber << "\n";
            }
        }
        if (tnew > timestep) {
            tnew = timestep*1.2;
        }
    }
	timestep = tnew;    
    htimestep = 0.5*timestep;
    //cout << gettime() << " "  << timestep << " " << df << " " << acc << " " << sqrt(x[0]*x[0] + x[1]*x[1]) << "\n";
//......................................................................................................
    //implementing the step function
    for (int j = 0; j < 8; j++) {
        step(pH);  //here calling the step function with 2 potentials
    }
    if (timefine > 0.00001*time) {
        time += timefine;
        timefine = 0.0;
    }
}

void particle::leapfrog(GalaxyPotential &pH, GalaxyPotential &pD) {
    stepnumber++;
	double grad[3];
	for (int i = 0; i < 3; ++i) {
		grad[i] = 0.0;
	}
    
        //difference here is the fPot taking 5 arguments
    //taking in the coordinates, force, pH, pD, and gettime
	fPot(x, grad, pH, pD, gettime());
        //the total acceleration of all 3 space-dimension
	double acc = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);
    double sacc = sqrt(acc);
    
    double q2 = 0.0;//.............................................................
    double df = 0.0;
    for (int j = 0; j < DIMENSION; j++) { //for each coordinates
        q2 += x[j]*x[j]; //position^2
        df += (f[j] - grad[j])*(f[j] - grad[j]);//force difference ^2
    }//.............................................................................
    df = sqrt(df/acc); //some sort of normalisation
    double tssug = (1.0/(10.0*df + 1.0e-14 + sacc + (1.0/MAXSTEP))) + MINSTEP;

    double tnew = TIMESTEPFAC*tssug;
    if ((tnew > timestep*2.0) || (tnew < htimestep)){ // if ... OR ...
        if (stepnumber > 2) {
            if (tnew < htimestep) {
                double saccold = sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);
                cout << "WARNING tnew " << tnew << " " << timestep << " " << tssug << " " << 1.0/sacc << " " << 1.0/saccold << " " << 1.0/(10.0*df) << " " << stepnumber << "\n";
            }
        }
        if (tnew > timestep) {
            tnew = timestep*1.2;
        }
    }
	timestep = tnew;    
    htimestep = 0.5*timestep;
    //cout << gettime() << " "  << timestep << " " << df << " " << acc << " " << sqrt(x[0]*x[0] + x[1]*x[1]) << "\n";
//......................................................................................................
    //implementing the step function
    for (int j = 0; j < 8; j++) {
        step(pH, pD);  //here calling the step function with 2 potentials
    }
    timefineadd();

}//=======================================================================================================


double particle::getE(GalaxyPotential &pH) {
    double Ekin = 0.0;
    for (int i = 0; i < 3; i++) {
        Ekin += 0.5*p[i]*p[i];
    }
    return (Ekin + pH(sqrt(x[0]*x[0] + x[1]*x[1]), x[2]));
}


void particle::print(GalaxyPotential &pH, ostringstream &ostr) {
    double Ekin = 0.0;
    for (int i = 0; i < 3; i++) {
        ostr << x[i] << " " << p[i] << " " << f[i] << " ";
        Ekin += 0.5*p[i]*p[i];
    }
    double R,phi,vR,vphi;
    phitrans(R,phi,vR,vphi);
    ostr << R << " " << phi << " " << vR << " " << vphi << " " << stepnumber << " ";
    ostr << timestep << " " << htimestep << " " << Ekin << " " << pH(R, x[2]) << " " << pH(R, 0.0) << " " << time+timefine;
    ostr << "\n";

    
}

void particle::print(GalaxyPotential &pH1, GalaxyPotential &pH2, ostringstream &ostr) {
    double Ekin = 0.0;
    for (int i = 0; i < 3; i++) {
        ostr << x[i] << " " << p[i] << " " << f[i] << " ";
        Ekin += 0.5*p[i]*p[i];
    }
    double R,phi,vR,vphi;
    phitrans(R,phi,vR,vphi);
    double coord[3]; 
    coord[0] = x[0]; coord[1] = x[1]; coord[2] = 0.0;
    ostr << R << " " << phi << " " << vR << " " << vphi << " " << stepnumber << " ";
    ostr << timestep << " " << htimestep << " " << Ekin << " " << gPot(x, pH1, pH2, gettime()) << " " << gPot(coord, pH1, pH2, gettime()) << " " << gettime();
    ostr << "\n";
}

void particle::fourier(GalaxyPotential &pH1, GalaxyPotential &pH2, ostringstream &ostr1, ostringstream &ostr2, ostringstream &ostr3, ostringstream &ostr22){    
    //Instructions for the gnuplot: 
    //(printorbit1) 1:5=time againsts r; 1:6=time againsts z
    //(printfourier) 1:3=frequencyr; 2:4=frequencyz
    #define M 8192
    //generate the data of an orbital trajectory
    double t[M];
    double rinitial[M];
    double zinitial[M];
    double periods[2];
    
    double maxtime = 10000.0;
    double magnitude=std::sqrt(M);
    x[0] = 4.0; x[1] = 0.0; x[2] = 0.0; 
    p[0] = 0.15; p[1] =0.22; p[2] = 0.02;
    t[0]=0.0;
    rinitial[0]=x[0];
    zinitial[0]=0.0;
    integrate(maxtime, pH1, pH2);
    for (int i = 1; i < M; i++) {
        for (int ii = 0; ii < 1; ii++) {
            leapfrog(pH1, pH2);
        }
        t[i] = gettime()-10000.0;
        ostr1<<gettime()-10000.0<<" "<<x[0]<<" "<< x[1]<<" "<<x[2]<<" ";
        rinitial[i] = std::sqrt(x[0] * x[0] + x[1] * x[1]);
        zinitial[i] = x[2];
        ostr1<<rinitial[i]<<" "<<zinitial[i]<<" ";
        ostr1 << "\n";
    }
    
    //interpolate between points to have fixed time intervals
    int i=0;
    double totaltime = gettime()-10000.0;
    double delta=(totaltime)/M;
    double ti, ri, zi;
    double rfinal[2*M];
    double zfinal[2*M];
    
    gsl_interp_accel *accr
      = gsl_interp_accel_alloc ();
    gsl_spline *spliner
      = gsl_spline_alloc (gsl_interp_cspline, M);
    gsl_interp_accel *accz
      = gsl_interp_accel_alloc ();
    gsl_spline *splinez
      = gsl_spline_alloc (gsl_interp_cspline, M);

    gsl_spline_init (spliner, t, rinitial, M);
    gsl_spline_init (splinez, t, zinitial, M);
    
    for (ti = t[0]; ti < t[M-1]; ti += delta)
      {
        ri = gsl_spline_eval (spliner, ti, accr);
        zi = gsl_spline_eval (splinez, ti, accz);
      }
    //find rmax and zmax
    double rmax=0;
    double zmax=0;
    double nr=0;//number of nr for Rmax to occur
    double nz=0;//number of nz for zmax to occur
    for (int i = 0; i < M; i++) {
        ti=1.0*delta*i;
        REAL(rfinal,i) = gsl_spline_eval (spliner, ti, accr); IMAG(rfinal,i) = 0.0;
        REAL(zfinal,i) = gsl_spline_eval (splinez, ti, accz); IMAG(zfinal,i) = 0.0;
        ostr2<<REAL(rfinal,i)<<" "<<REAL(zfinal,i)<<" ";
        ostr2 << "\n";
        if (REAL(rfinal,i) > rmax) {
                rmax = REAL(rfinal,i);
                nr = i;
        } 
        
        if (fabs(REAL(zfinal,i)) > zmax) {
                zmax = fabs(REAL(zfinal,i));
                nz = i;
        }
        if (i==(M-1)){
            ostr22 << nr+1 << " "<< nz+1 << " "<< rmax << " " << zmax << " ";
            ostr22 << "\n";
        }
        
      }
    //free the storage
    gsl_spline_free (spliner);
    gsl_interp_accel_free (accr);
    gsl_spline_free (splinez);
    gsl_interp_accel_free (accz);
  
    
    //Fourier analysis
    gsl_fft_complex_radix2_forward (rfinal, delta, M);
    gsl_fft_complex_radix2_forward (zfinal, delta, M);
    
    // store fr, fz, magnitudes 
    double frequencyr[M];
    double frequencyz[M];
    double sqrtmagr[M];
    double sqrtmagz[M];

    //Fourier shift (centre)
    int j;
    for (i = 0; i < M; i++)
    {
        if (i<M/2){
            j = i+M/2;
            sqrtmagr[j]=sqrt(REAL(rfinal,j)*REAL(rfinal,j)+IMAG(rfinal,j)*IMAG(rfinal,j));
            sqrtmagz[j]=sqrt(REAL(zfinal,j)*REAL(zfinal,j)+IMAG(zfinal,j)*IMAG(zfinal,j));
            frequencyr[j]=-1.0+j*1.0/(M);
            frequencyz[j]=-1.0+j*1.0/(M);
        }
        else{
            j= i-M/2;
            sqrtmagr[j]=sqrt(REAL(rfinal,j)*REAL(rfinal,j)+IMAG(rfinal,j)*IMAG(rfinal,j));
            sqrtmagz[j]=sqrt(REAL(zfinal,j)*REAL(zfinal,j)+IMAG(zfinal,j)*IMAG(zfinal,j));
            frequencyr[j]=j*1.0/(M);
            frequencyz[j]=j*1.0/(M);
        }
        ostr3<<frequencyr[j]<<" "<<frequencyz[j]<<" "<<sqrtmagr[j]/M<<" "<<sqrtmagz[j]/M<<" ";
        ostr3 << "\n";
    }

}

void particle::gaussianvelocity(GalaxyPotential &pH1, GalaxyPotential &pH2, ostringstream &ostr,double steps){
    //Instructions for the gnuplot: (printorbitNa) 4:5=heated; 6:7=unperturbed; 6:8=levitation
    int a;
    int b;
    a= rand() % 300+1;
    b= rand() % 300+1;
    // For the adiabaticly heated particle
    double settime(); 
    x[0] = 4.0; x[1] = 0.0; x[2] = 0.0; 
    p[0] = 0.15; p[1] =0.22; p[2] = 0.0002*a;
    ostr<< p[0]<< " "<< p[1]<< " "<< p[2]<< " ";
    double periods[2];
    double maxtime = 10000.0;
    integrate(maxtime, pH1, pH2);
    double RmaxA=0;
    double zmaxA=0;
    double RA = getR();
    double zA = getz();
    double nrA=0;//number of nr for Rmax to occur
    double nzA=0;//number of nz for zmax to occur
    for (int i = 0; i < steps; i++) {
        for (int ii = 0; ii < 1; ii++) {
            leapfrog(pH1, pH2);
        }
        //print(pH1, pH2, ostr);
        double RA,phiA,vRA,vphiA;
        phitrans(RA,phiA,vRA,vphiA);
        //extract z step by step
        double zA;
        zA = fabs(x[2]);
        
        if (RA > RmaxA) {
                RmaxA = RA;
                nrA = i;
        } 
        
        if (zA > zmaxA) {
                zmaxA = zA;
                nzA = i;
        }
        if (i==(steps-1)){
            ostr << nrA+1 << " "<< nzA+1 << " "<< RmaxA << " " << zmaxA << " ";
        }
    }
    

    double settime();  
    // For the contracted particle
    x[0] = 4.0; x[1] = 0.0; x[2] = 0.0; 
    p[0] = 0.15; p[1] =0.22; p[2] = 0.0002*a;
    integrate(maxtime, pH1);
    double RmaxB=0;
    double zmaxB=0;
    double RB = getR();
    double zB = getz();
    for (int i = 0; i < steps; i++) {
        for (int ii = 0; ii < 1; ii++) {
            leapfrog(pH1);
        }
        double RB,phiB,vRB,vphiB;
        phitrans(RB,phiB,vRB,vphiB);
        //extract z step by step
        double zB;
        zB = fabs(x[2]);
        
        if (RB > RmaxB) {
                RmaxB = RB;
        }        
        if (zB > zmaxB) {
                zmaxB = zB;
        }
        if (i==(steps-1)){
            ostr << RmaxB << " " << zmaxB << " "<< zmaxA-zmaxB << " " << "\n";
        }
    }

}

void particle::print(GalaxyPotential &pH) {
    double Ekin = 0.0;
    for (int i = 0; i < 3; i++) {
        cout << x[i] << " " << p[i] << " " << f[i] << " ";
        Ekin += 0.5*p[i]*p[i];
    }
    cout << timestep << " " << Ekin << " " << pH(sqrt(x[0]*x[0] + x[1]*x[1]), x[2]) << " " << time+timefine;
    cout << "\n";
}


void particle::print(GalaxyPotential &pH, GalaxyPotential &pD) {
    double Ekin = 0.0;
    for (int i = 0; i < 3; i++) {
        cout << x[i] << " " << p[i] << " " << f[i] << " ";
        Ekin += 0.5*p[i]*p[i];
    }
    cout << timestep << " " << Ekin << " " << pH(sqrt(x[0]*x[0] + x[1]*x[1]), x[2]) + pD(sqrt(x[0]*x[0] + x[1]*x[1]), x[2]) << " " << time+timefine;
    cout << "\n";
}


void particle::phitrans(double &R, double &phi, double &vR, double &vphi) {
    phi = 0.0;
    if (x[1] != 0.0 || x[0] != 0.0) {
        phi = atan2(x[1],x[0]);
        R = sqrt(x[0]*x[0] + x[1]*x[1]);
//    double vs = p[0]*p[0] + p[1]*p[1];
        vR = (x[0]*p[0] + x[1]*p[1])/R;
        vphi = (-x[1]*p[0] + x[0]*p[1])/R;
    }
    else {
       R = 0.0; phi = 0.0; vR = 0.0; vphi = 0.0;
       cout << "WARNING do not use origin in phitrans() \n";
    }
}


void particle::phiset(double R,double phi,double vR, double vphi) {
    x[0] = R*cos(phi);
    x[1] = R*sin(phi);
    p[0] = vR*cos(phi) - vphi*sin(phi);
    p[1] = vR*sin(phi) + vphi*cos(phi);
}

void particle::phishift(std::uniform_real_distribution<double> phasedist) {
    double delphi = phasedist(randomGen);
    double R; double phi; double vR; double vphi;
    phitrans(R,phi,vR,vphi);
    phi += delphi;
    phiset(R,phi,vR,vphi);
}

void particle::phishift(double delphi) {
    double R; double phi; double vR; double vphi;
    phitrans(R,phi,vR,vphi);
    phi += delphi;
    phiset(R,phi,vR,vphi);
}


void particle::getsimpleperiods(GalaxyPotential &pH, double periods[2], double maxtime) {
    double Rmax = getR();
    double t0 = gettime();
    double zmax = 0.0;
    double zmin = 0.0;
    double Rmin = Rmax;
    double tRmin = 0.0;
    double tRmax = 0.0;
    double tzmin = 0.0;
    double tzmax = 0.0;
    double R = getR();
    double z = getz();
    double Rsave[5];
    double zsave[5];
    for (int i = 0; i < 5; i++) {
        Rsave[i] = 0.0; zsave[i] = 0.0;
    }
    int Rrising = 0;
    int zrising = 0;
    int Rrisingold = 0;
    int zrisingold = 0;
    int icount = 0;
    int Rcount = 0;
    int zcount = 0;
//    cout << "time - t0 " << (gettime()-t0) << " " << time <<  " " << maxtime << "\n";
    while (((gettime() - t0) < maxtime) && ((Rcount < 2) || (zcount < 2))) {
        icount++;
        leapfrog(pH);
        R = getR();
        z = getz();
        Rrisingold = Rrising;
        zrisingold = zrising;
        if (Rsave[4] > R) {Rrising = -1;} else {Rrising = 1;}
        if (zsave[4] > z) {zrising = -1;} else {zrising = 1;}
        for (int ii = 4; ii >= 1; ii--) {
            Rsave[ii] = Rsave[ii-1];
            zsave[ii] = zsave[ii-1];
        }
        Rsave[0] = R;
        zsave[0] = z;
//        cout << "time - t0 " << icount << " " << (gettime()-t0) << " " << time <<  " " << maxtime << "\n";

        if (icount > 6) {        
            if (Rrising*Rrisingold < 0) {Rcount++;}
            if (zrising*zrisingold < 0) {zcount++;}
            if (R < Rmin) {
                Rmin = R;
                tRmin = time;
            }
            if (R > Rmax) {
                Rmax = R;
                tRmax = time;
            }
            if (z > zmax) {
                zmax = z;
                tzmax = time;
            }
            if (z < zmin) {
                zmin = z;
                tzmin = time;
            }
            cout << "time " << time << " " <<  R << " " << z << " " << Rmin << " " << Rmax << " " << zmin << " " << zmax << " " << Rcount << " " << zcount << "\n";
        }
    }
    
    periods[0] = 2.0*fabs(tRmax - tRmin);
    periods[1] = 2.0*fabs(tzmax - tzmin);
//  cout << "periods " << periods[0] << " " << periods[1] << "\n";
}


void particle::clone(particle &pp) {
    for (int i = 0; i < DIMENSION; i++) {
        x[i] = pp.x[i];
        p[i] = pp.p[i];
        f[i] = pp.f[i];
    }
    weight = pp.weight;
    time = pp.time;
    timestep = pp.timestep; htimestep = pp.htimestep;
}

void particle::integrate(double integtime, GalaxyPotential &pH) {
    stepnumber = 0;
    double ttarget = gettime() + integtime;
    while (gettime() < ttarget) {
        leapfrog(pH);
    }
}

void particle::integrate(double integtime, GalaxyPotential &pH, GalaxyPotential &pD) {
    stepnumber = 0;
    double ttarget = gettime() + integtime;
    while (gettime() < ttarget) {
        leapfrog(pH, pD);
    }
}
    
    



void compint(int i, std::uniform_int_distribution<int> intdist, int *savearray) {
	for (int j = 0; j < i; ++j) {
		if (savearray[j] == savearray[i]) {
			savearray[i] = intdist(randomGen);
			compint(i, intdist, savearray);
		}
	}
}

void diffinteg(int no, std::uniform_int_distribution<int> intdist, int *savearray) {
	for (int i = 0; i < no; ++i) {
		savearray[i] = intdist(randomGen);
		compint(i, intdist, savearray);
	}
}


//--------------------Phiweight---------------------
double intfunc(double c, double Rg, double Ez, double x, double kalpha, double scalelength) { //actual function to be integrated
	double res = 0.0;
	double expo = c * (2.0*log(Rg / x) + 1.0 - (Rg*Rg) / (x*x) - (2.0*Ez / (VC*VC))*(exp(-((x - Rg) / scalelength)*kalpha) - 1.0));
	if ((expo < 350.0) && (expo > -350.0)) {
		res += exp(expo);
	}
	return res;
}

double integ(double c, double Rg, double Ez, double kalpha, double scalelength) { //brute force integrator to obtain normalization
																				  //cout << "enter integ \n";
	double sum = 0.0;
	double sumb = 0.0;
	double dn = 0.00001;
	double x = 0.000001;
	double valold = intfunc(c, Rg, Ez, x, kalpha, scalelength);
	double valnew = intfunc(c, Rg, Ez, x + dn, kalpha, scalelength);
	double valm = intfunc(c, Rg, Ez, x + 0.5*dn, kalpha, scalelength);
	while (x < 700000000.0) {
		valold = intfunc(c, Rg, Ez, x, kalpha, scalelength);
		valnew = intfunc(c, Rg, Ez, x + dn, kalpha, scalelength);
		valm = intfunc(c, Rg, Ez, x + 0.5*dn, kalpha, scalelength);
		//    if ((valold*valold > 1.0e-80) && (((valold + valnew)*0.5*dn) > (1.0e-20 + (1.0e-8)*sum))) {
		while ((((valold - valnew)*(valold - valnew) / (valold*valold)) > 0.002) && (dn > 1.0e-9) && ((valold + valnew)*0.5*dn) > (1.0e-11*sum)) {
			dn *= 0.4;
			valnew = intfunc(c, Rg, Ez, x + dn, kalpha, scalelength);
			valm = intfunc(c, Rg, Ez, x + 0.5*dn, kalpha, scalelength);
		}
		//    }
		sumb += (0.25*(valold + valnew + valm + valm))*dn;
		if (sumb > 0.000001*sum) {
			sum += sumb;
			sumb = 0.0;
		}
		x += dn;

		//cout<< "integ " << x << "  " << dn << "  " << ((valold - valnew)*(valold - valnew)/(valold*valold)) << "  " <<  0.5*(valold+valnew) << "  " << sum << "  " << sumb << "\n";

		if ((valnew < valold) || (dn < 0.05*x)) {
			dn *= 1.8;
		}
	}
	sum += sumb;
	// cout << "leave integ \n";
	return sum;
}

double getintval(double ***values, double c, double Rg, double Ez) { //read from the normalization array
																	 //cout << "enter getintval \n";

	int i0 = ((int)(sqrt((c - 0.5) / CVL)));
	//cout << "i0 " << i0 << "\n";
	int g0 = ((int)(Rg / GVL));
	int zz0 = ((int)((sqrt(Ez)) / ZVL));
	//int iu = i0 + 1;
	if (i0 < 1)
	{
		i0 = 1;
	}
	if (i0 >(CVALN - 2))
	{
		i0 = CVALN - 2;
	}
	double i0d = ((double)(i0));

	double iud = i0d + 1.0;
	if (g0 < 1)
	{
		g0 = 1;
	}
	if (g0 >(GVALN - 2))
	{
		g0 = GVALN - 2;
	}
	double g0d = ((double)(g0))*GVL;
	if (zz0 < 0)
	{
		zz0 = 0;
	}
	if (zz0 > ZVALN - 2)
	{
		zz0 = ZVALN - 2;
	}
	double zz0d = ((double)(zz0));
	double il = iud - i0d;
	double gl = GVL;
	double zl = ZVL;
	double ifac = (sqrt((c - 0.5) / CVL) - i0d) / il;
	double gfac = (Rg - g0d) / gl;
	double zzfac = ((sqrt(Ez)) - zz0d * ZVL) / zl;
	//
	//cout << "getinterval " << i0 << "  " << g0 << "  " << zz0 << "  " << ifac << "  " << gfac << "  " << zzfac << "\n";
	//
	if (ifac < -0.0) {
		ifac = -0.0;
	}
	if (gfac < -0.0) {
		gfac = -0.0;
	}
	if (zzfac < -0.0) {
		zzfac = -0.0;
	}
	if (ifac > 1.0) {
		ifac = 1.0;
	}
	if (gfac > 1.0) {
		gfac = 1.0;
	}
	if (zzfac > 1.0) {
		zzfac = 1.0;
	}
	double back = 0.0;
	double backi00 = values[i0][g0][zz0] + (values[i0 + 1][g0][zz0] - values[i0][g0][zz0])*ifac;
	double backi01 = values[i0][g0][zz0 + 1] + (values[i0 + 1][g0][zz0 + 1] - values[i0][g0][zz0 + 1])*ifac;
	double backi10 = values[i0][g0 + 1][zz0] + (values[i0 + 1][g0 + 1][zz0] - values[i0][g0 + 1][zz0])*ifac;
	double backi11 = values[i0][g0 + 1][zz0 + 1] + (values[i0 + 1][g0 + 1][zz0 + 1] - values[i0][g0 + 1][zz0 + 1])*ifac;
	double backg0 = backi00 + (backi10 - backi00)*gfac;
	double backg1 = backi01 + (backi11 - backi01)*gfac;
	back = backg0 + (backg1 - backg0)*zzfac;
	//
	//cout << "leave getintval \n";
	//
	return back;
}

//gives c, approximates it if necessary  
double getnormfact(double c)
{
	if (c < 70.0) {
		return exp(-c)*pow(c, (c - 0.5)) / (tgamma(c - 0.5));
	}
	else {
		return sqrt(c - 0.913) / (sqrt(2.0*3.1415));
	}
}

void getint(double ***values, double alpha, double scalelength) {
	cout << "enter getint \n";
	for (int i = 1; i < CVALN; i++)
	{
		for (int g = 1; g < GVALN; g++)
		{
			double Rg = ((double)(g))*GVL;
			for (int zz = 0; zz < ZVALN; zz++)
			{
				double Ez = ((double)(zz))*ZVL;
				Ez *= Ez;

				double c = ((double)(i))*((double)(i))*CVL + 0.5;

				double kalpha = 2.0 / (2.0 + alpha);
				double intres = integ(c, Rg, Ez, kalpha, scalelength);  //integrate
				double fac = getnormfact(c)*2.0 / Rg;
				values[i][g][zz] = 1.0 / intres;						//set normalisation constants
				values[i][g][zz] = (1.0 / intres) / fac;
			}
		}
	}
	cout << "leave getint \n";
}


void getintpar(double ***values, double alpha, double scalelength, int i) {
	//cout << "enter getint \n";
	//	for (int i = 1; i < CVALN; i++)
	{
		for (int g = 1; g < GVALN; g++)
		{
			double Rg = ((double)(g))*GVL;
			for (int zz = 0; zz < ZVALN; zz++)
			{
				double Ez = ((double)(zz))*ZVL;
				Ez *= Ez;

				double c = ((double)(i))*((double)(i))*CVL + 0.5;

				double kalpha = 2.0 / (2.0 + alpha);
				double intres = integ(c, Rg, Ez, kalpha, scalelength);  //integrate
				double fac = getnormfact(c)*2.0 / Rg;
				values[i][g][zz] = 1.0 / intres;						//set normalisation constants
				values[i][g][zz] = (1.0 / intres) / fac;
			}
		}
	}
	//cout << "leave getint \n";
}


double func(double ***inttab, double *param, double v) {
	if (param[8] < -1.0) {
		param[8] = -1.0;
	}

	if (param[8] > 2.0) {
		param[8] = 2.0;
	}  //limits alpha

	double a = param[0];  //normalization
	double s0 = param[1];  //horizontal vel. dispersion at local galactocentric radius
	double h0 = param[2];  //scale height
	double vc = param[3];  //circular velocity
	double l = param[4];   //horizontal vel. dispersion scale length in km/s
//	double m = param[5];   //vertical vel. dispersion scale length
	double nx = param[6];   //disc scale length
	double hh = param[7];   //altitude of measurement
	double alpha = param[8];  //potential adiabatic index alpha
							  ///double Rg = param[9];		///added myself

	double n = exp(-v / nx);
	double f = 0.0;
	if (v > 0.0)		///seems strange condition but well
	{
		double s = s0 * exp(-(v - vc) / l);
		//cout << "s " << s << "\n";
		double c = vc * vc / (2.0*s*s);
		//cout << "c " << c << "\n";
		double spi = PI;
		double sMsol = 1.989E+30;
		double sG = 6.67e-11;
		//        double a = 1.7051;
		double spars = 3.0856776E+16;			///1 pc = x meters
		double m0 = 30.0;    ///no clue what the following stuff means
		double mg = 12.0;
		//double hg = 10.0;
		double h1 = 300.0;			///scale lenght inner disc
		double h2 = 1000.0;			///scale length outer disc
		double m3 = 70.0;
		double h3 = 4000.0;

		double mactb = mg + m0 * (1.0 - exp(-h0 / h1)) + m0 * 0.34*(1.0 - exp(-h0 / h2)) + m3 * (1.0 - exp(-h0 / h3));
		double mbfact = 0.000001*2.0*spi*sG*sMsol / spars;
		double Epotz = mg * hh + (m0 + m0 * 0.34 + m3)*hh + (h1*m0*(exp(-hh / h1) - 1.0)) + (m0*0.34*h2*(exp(-hh / h2) - 1.0)) + (m3*h3*(exp(-hh / h3) - 1.0));
		Epotz *= mbfact;
		double sz0 = sqrt(0.5*h0*(mactb + mactb)*mbfact);

		if (sz0 < 10.0) {  // limits the local vertical dispersion as derived from the local scale height
			sz0 = 10.0;
		}
		if (sz0 > 50.0) {
			sz0 = 50.0;
		}

		double sz = sz0 * exp(-(v - vc) / (2.0*nx));

		double kalpha = 2.0 / (2.0 + alpha);
		double Ez = (0.5*sz*sz) + Epotz * exp(-((v - vc)*kalpha / nx));
		double b = exp((v - vc) / ((2.0 + alpha)*nx))*exp((hh / h0)*(1.0 - exp((v - vc) / ((2.0 + alpha)*nx))));  // weight difference by origin dependent scaleheight
		double fac = 1.0;
		double Rg = v * RSUN / vc;		///sense??
		fac = getintval(inttab, c, Rg, Ez);  //difference to normalization without adiabatic potential
											 //cout << "fac in func " << fac << "\n";
		double facold = 0.0;  //normalization of the normal effective potential

		if (c < 70.0) {
			facold = exp(-c)*pow(c, (c - 0.5)) / (tgamma(c - 0.5));
		}

		else {
			facold = sqrt(c - 0.913) / (sqrt(2.0*3.1415));
		}
		//actual function
		//  cout << a << "  " << n << "  " << b << "  " << fac << "  " << facold << "  " << nx << "  "<< l << "  " << kalpha << "  " << sz << "  " << sz0 << "  "  << Epotz << "  " <<  Ez << "\n";

		f = a * n*b*fac*facold*exp(-(c*(v*v / (vc*vc) - log(v*v / (vc*vc)) - 1.0 + (2.0*Ez / (vc*vc))*(exp(-(vc - v)*kalpha / nx) - 1.0))));
		if (f < 0.0) {
			f = 0.0;
		}
	}
	return f;
}

double *funcal(double ***inttab, double *param, double vmin, double vmax) { //calculates moments, integrates between vmin and vmax
	double vvv = vmin;
	//double vc = param[3];	///added
	double vc = VC;
	double l = param[4];
	double m = param[5];
	double nx = param[6];
	double alpha = param[8];
	double ww = 0.0;
	double suu = 0.0;
	double sww = 0.0;
	while (vvv < vmax) {
		double weight = func(inttab, param, vvv);
		//double sig = param[1] * exp(-(vvv - vc) / l);
		double sig = param[1] * exp(-(vvv - vc) / l);
		//double sigW = exp(-(vvv - vc) / m)*sqrt(exp(((vvv - vc) / nx)*(2.0 / (2.0 + alpha))));
		double sigW = param[10] * exp(-(vvv - vc) / m)*sqrt(exp(((vvv - vc) / nx)*(2.0 / (2.0 + alpha)))); ///added param[10]

		ww += weight;
		suu += sig * sig*weight;
		sww += sigW * sigW*weight;
		vvv += 0.25;
	}
	double *back = new double[2];
	back[0] = sqrt(suu / ww); // U velocity dispersion r direction
	back[1] = sqrt(sww / ww); // W velocity dispersion/dispersion of local population
							  //cout << back[0] << "  " << back[1] << "\n";
	return back;
}



void weighting(int nset, int nper, particle *pp, std::uniform_real_distribution<double> phasedist, std::normal_distribution<double> pphidistribution, double sigphi, std::normal_distribution<double> prdistribution, std::normal_distribution<double> pzdistribution, std::uniform_real_distribution<double> rInidist, std::uniform_real_distribution<double> zInidist, double ***intval, double *param, double scalelength, double rInird, double zInizd, int threadno) {
//    std::ifstream fromfile("./GalPot-master/pot/PJM17_best.Tpot");
    std::ifstream fromfile1(FROMFILE1);
    GalaxyPotential pH(fromfile1);

    double vphiplace = VPHIPLACE;
	double paramthread[11];
	for (int n = 0; n < 11; n++) {
		paramthread[n] = param[n];
	}

	double sigphig = 1.0 / (2.0 * sigphi*sigphi);

	for (int i = 0; i < nset; i++)
	{
		double xi = rInidist(randomGen);
		double ri = -rInird * log(exp(-RMIN / rInird) - xi / rInird);
		xi = zInidist(randomGen);
		double zstart = xi > 0.0 ? (-zInizd * log(1.0 - xi / zInizd)) : (zInizd * log(1.0 + xi / zInizd));

		pp[i].x[0] = ri;
		pp[i].x[1] = phasedist(randomGen);
		pp[i].x[2] = zstart;

		paramthread[7] = std::fabs(zstart) * 1000.0;

		double zweight = exp(-fabs(zstart / 0.3)) / exp(-fabs(zstart / zInizd));
		double Rweight = exp(-ri / 2.9) / exp(-ri / rInird);

		double pphi = -1.0;
		while (pphi < 0.01 || pphi > 0.45) {
			pphi = pphidistribution(randomGen);
		}
		pp[i].p[1] = pphi;

		double vphigaussweight = /*(1.0 / (sigmavPhi * std::sqrt(2.0 * PI)) * */ std::exp(-(pp[i].p[1] - vphiplace)*(pp[i].p[1] - vphiplace) *sigphig);

		double vphiweight = func(intval, paramthread, (pp[i].p[1])*1000.0) / vphigaussweight;

		double *back = funcal(intval, paramthread, (pp[i].p[1])*1000.0 - 5.0, (pp[i].p[1])*1000.0 + 5.0);

		//std::normal_distribution<double> prdistribution(0.0, sigmaVr);
		//std::normal_distribution<double> pzdistribution(0.0, sigmaVz);


		double pr = -8000.0; //*back[0];
		while ((pr <-0.007*back[0] || pr > 0.007*back[0]) && (fabs(pr) > 0.45)) {
			pr = prdistribution(randomGen);
		}
		pp[i].p[0] = pr;

		double pz = -8000.0;   //*back[1];
		while ((pz <-0.007*back[1] || pz > 0.007*back[1]) && (fabs(pz) > 0.4)) {
			pz = pzdistribution(randomGen);
		}
		pp[i].p[2] = 0.2*pz;
        pp[i].x[2] *= 0.3;
		//cout << pr << " " << pz << "\n";

		double r3 = sqrt(ri*ri + zstart * zstart);
		double prtilt = (pp[i].p[0] * (ri / r3) - pp[i].p[2] * (zstart / r3))* 1000.0;
		double pztilt = (pp[i].p[2] * (ri / r3) + pp[i].p[0] * (zstart / r3))* 1000.0;
		//cout << PSvector[1][i][0] * 1000.0 << " " << prtilt << " " << PSvector[1][i][2] * 1000.0 << " " << pztilt << "\n";

		//double prweightGauss = (1.0 / (2.0 * back[0] * std::sqrt(2.0 * PI))  * std::exp(-(pztilt*pztilt) / (2.0 *4.0 *back[0] * back[0])));
		//double prweightReal = (1.0 / (back[0] * std::sqrt(2.0 * PI))  * std::exp(-(pztilt*pztilt) / (2.0 *back[0] * back[0])));
		//double prweight = prweightReal / prweightGauss;

		//double pzweightGauss = (1.0 / (2.0 * back[1] * std::sqrt(2.0 * PI))  * std::exp(-(prtilt*prtilt) / (2.0 *4.0 *back[1] * back[1])));
		//double pzweightReal = (1.0 / (back[1] * std::sqrt(2.0 * PI))  * std::exp(-(prtilt*prtilt) / (2.0 *back[1] * back[1])));
		//double pzweight = pzweightReal / pzweightGauss;

		double sigVrd = 1.0 / (paramthread[1] * paramthread[1] * VRPLACEFAC*VRPLACEFAC);
		double sigVzd = 1.0 / (paramthread[10] * paramthread[10] * VZPLACEFAC*VZPLACEFAC);
		double back0d = 1.0 / (back[0] * back[0]);
		double back1d = 1.0 / (back[1] * back[1]);
		//cout << back[0] << " " << back[1]<< "\n";

		double prweight = (1.0 / back[0]) /* /sigmaVr */ * exp(-((prtilt*prtilt) / 2.0) *(back0d - sigVrd)); //for pr/pz distributed with 
		double pzweight = (1.0 / back[1]) /* /sigmaVz */ * exp(-((pztilt*pztilt) / 2.0)*(back1d - sigVzd)); //for pr/pz distributed with 
																											//double prweight = back[0] /* /sigmaVr */ * exp(-((prtilt*prtilt) / (2.0)) *(back0d- sigVrd)); //for pr/pz distributed with 
																											//double pzweight = back[1] /* /sigmaVz */ * exp(-((pztilt*pztilt) / (2.0) )*(back1d - sigVzd)); //for pr/pz distributed with 


																											//std::cout << "vphiweight " << vphiweight << " pzweight " << pzweight << " prweight " << prweight << "\n";

		pp[i].weight = vphiweight * pzweight * prweight * Rweight * zweight;
        double RR = pp[i].x[0];
        double phiphi = pp[i].x[1];
//        double pphi = pp[i].p[1];
        double xx = RR*cos(phiphi);
        double yy = RR*sin(phiphi);
        double px = pr*cos(phiphi) - pphi*sin(phiphi);
        double py = pr*sin(phiphi) + pphi*cos(phiphi);
        pp[i].x[0] = xx; pp[i].x[1] = yy;
        pp[i].p[0] = px; pp[i].p[1] = py;
        pp[i].time = 0.0; pp[i].timefine = 0.0; pp[i].htimestep = 0.05; pp[i].timestep = 0.1;
        double E = pp[i].getE(pH);
        if (E > 0) {
            i--;
        }
		delete[] back;
	}
	
	
	for (int i = 0; i < nset; i++) {
        double periods[2];
        double maxrun = 4000.0;
        double minrun = 1000.0;
        pp[i].getsimpleperiods(pH, periods, maxrun);
        double testperiod = periods[0];
        double masterperiod = periods[0];
        double target = 2.0;
        while (testperiod < maxrun) {
            double xtest = 1.0 - cos(2.0*PI*testperiod/periods[1]);
            if (xtest < target) {
                target = xtest;
                masterperiod = testperiod;
            }
            testperiod += periods[0];
        }
        while (masterperiod < minrun) {
            masterperiod *= 2.0;
        }
        cout << "periods " << periods[0] << " " << periods[1] << " " << target << " " << masterperiod << " " << testperiod << "\n";
        double integtime = masterperiod/((double)(nper));
        pp[i + nset].clone(pp[i]);
        pp[i + nset].settime();
        for (int ii = 1; ii < nper-1; ii++) {
            pp[i].integrate(integtime, pH);
            pp[i + (nset*(ii + 1))].clone(pp[i]);
            pp[i + (nset*(ii + 1))].settime();
            pp[i+(nset*(ii+1))].phishift(phasedist);
        }
        pp[i].integrate(integtime, pH);
        pp[i].settime();
    }
}

void integset(particle *ppa, int nsub, double integtime, int threadno) {
   std::ifstream fromfile1(FROMFILE1); std::ifstream fromfile2(FROMFILE2);
    GalaxyPotential pH1(fromfile1);
    GalaxyPotential pH2(fromfile2);    
    cout << "integset potentials set " << threadno << "\n";
    for (int i = 0; i < nsub; i++) {
        cout << "integset " << threadno << " " << i << "\n";
        ppa[i].integrate(integtime, pH1, pH2);
    }
}

void DoDisk() {
    std::ifstream fromfile1(FROMFILE1); std::ifstream fromfile2(FROMFILE2);
    GalaxyPotential pH1(fromfile1);
    GalaxyPotential pH2(fromfile2);
        double sigmaVrIn = 0.04;
	double sigmaVzIn = 0.03;
//    double frequarray[500];
//	frequarray[0] = 0.0;
    
	double alpha = 1.0;
	double rd = 2.9;

	//normalization array for weights
    
	double ***intval = new double**[CVALN];
	for (int n = 0; n < CVALN; n++) {
		intval[n] = new double*[GVALN];
		for (int m = 0; m < GVALN; m++) {
			intval[n][m] = new double[ZVALN];
			for (int p = 0; p < ZVALN; p++) {
				intval[n][m][p] = 0.0;
			}
		}
	}

	thread *persei = new thread[THREADNO];
	int perbes[THREADNO];
	for (int i = 0; i < THREADNO; i++) {
//		persei[i];
		perbes[i] = 0;
	}
	int countper = 0;
	for (int i = 1; i < CVALN; i++) {
		if (perbes[countper] == 1) {
			if (persei[countper].joinable()) {
				persei[countper].join();
				perbes[countper] = 0;
			}
		}
		persei[countper] = thread(getintpar, intval, alpha, rd, i);
		perbes[countper] = 1;
		countper++;
		if (countper >= THREADNO) {
			countper = 0;
		}
	}
	for (int i = 0; i < THREADNO; i++) {
		if (perbes[i] > 0) {
			if (persei[i].joinable()) {
				persei[i].join();
				perbes[i] = 0;
			}
		}
	}
	delete[] persei;
    
    
 	cout << "done caluclating normarray\n";
	double sigmaVr = sigmaVrIn;
	double sigmaVz = sigmaVzIn;

	double param[11];
	for (int n = 0; n < 11; n++) {
		param[n] = 0.0;
	}

	/*
	double a = param[0];  //normalization
	double s0 = param[1];  //horizontal vel. dispersion at local galactocentric radius
	double h0 = param[2];  //scale height
	double vc = param[3];  //circular velocity
	double l = param[4];   //horizontal vel. dispersion scale length in km/s
	double m = param[5];   //vertical vel. dispersion scale length
	double nx = param[6];   //disc scale length
	double hh = param[7];   //altitude of measurement
	double alpha = param[8];  //potential adiabatic index alpha
	*/

	double vcs = VC;
	double Rsun = RSUN;
	param[0] = 0.1;
	param[1] = sigmaVr * 1000.0;		//35.0;		
	param[2] = 360.0;						/// changed to scale height for thin disc
	param[3] = vcs;
	param[4] = 7.5 * vcs / Rsun;			///7.5*lz ~ 200
	param[5] = 6.0 * vcs / Rsun;			///~150
	double scalelength = 2.9;   //disc scalelength
	param[6] = scalelength * vcs / Rsun;  /// get rd in velocity units
	param[7] = 0.0; //defined in loop	
	param[8] = 1.0;							///alpha

	param[10] = sigmaVz * 1000.0;
	double sigphi = VPHIPLACESIG;

	double rInird = RPLACEL;
	double zInizd = ZPLACEL;

	std::uniform_real_distribution<double> rInidist(0.0, rInird* (exp(-RMIN / rInird) - exp(-RMAX / rInird)));
	std::uniform_real_distribution<double> zInidist(-zInizd * (1.0 - exp(-2.0 / zInizd)), zInizd * (1.0 - exp(-2.0 / zInizd)));

	std::uniform_real_distribution<double>  phasedist(0.0, 2.0 * PI);
	std::normal_distribution<double> pphidistribution(VPHIPLACE, sigphi);
	std::normal_distribution<double> pzdistribution(0.0, /*2.0* */ param[10] * VZPLACEFAC / 1000.0);
	std::normal_distribution<double> prdistribution(0.0, /*2.0* */ param[1] * VRPLACEFAC / 1000.0);

    particle **ppa = new particle*[THREADNO];
    particle **ppa0 = new particle*[THREADNO];
    int nbb = NBASE;
    int nperb = NPERBASE;
    std::thread *weightthreads = new std::thread[THREADNO];
    for (int i = 0; i < THREADNO; i++) {
        ppa0[i] = new particle[NBASE*NPERBASE];
        ppa[i] = new particle[NBASE*NPERBASE];
    }
    
    cout << "launch the threads \n";
    	for (int u = 0; u < THREADNO; u++)
	{
		weightthreads[u] = std::thread(weighting, nbb, nperb, ppa[u], phasedist, pphidistribution, sigphi, prdistribution, pzdistribution, rInidist, zInidist, intval, param, scalelength, rInird, zInizd, u);
	}
//	weighting(THREADNO*particlePerThread, particleNumber, PSvector, phasedist, pphidistribution, sigphi, prdistribution, pzdistribution, rInidist, zInidist, intval, param, scalelength, rInird, zInizd, THREADNO);

    cout << "join weighting \n";
	for (int u = 0; u < THREADNO; u++)
	{
		weightthreads[u].join();
        cout << "thread " << u << " joined \n";
	}
	
    cout << "delete weighting threads \n";
	delete[] weightthreads;

	for (int n = 0; n < CVALN; n++) {
		for (int m = 0; m < GVALN; m++) {
			delete[] intval[n][m];
			}
		delete[] intval[n];
		}
	delete[] intval;

	//------------------------------end intital conditions-----------------------------
	cout << "done weighing basisorbits\n";

   int nsub = NBASE*NPERBASE;
   for (int i = 0; i < THREADNO; i++) {
       for (int jj = 0; jj < nsub; jj++) {
           ppa0[i][jj].clone(ppa[i][jj]);
       }
   }
   cout << "cloned particles \n";
   std::thread *perseint = new std::thread[THREADNO];   
   double integtime = 10000.0;
   
   for (int i = 0; i < THREADNO; i++) {
        perseint[i] = std::thread(integset, ppa[i], nsub, integtime, i);
   }
   for (int i = 0; i < THREADNO; i++) {
       perseint[i].join();
       cout << "thread " <<  i << " joined \n";
   }
   delete [] perseint;
   
   ostringstream ostr;
   for (int i = 0; i < THREADNO; i++) {
       for (int jj = 0; jj < nsub; jj++) {
        ppa[i][jj].print(pH1, pH2, ostr); 
        ppa0[i][jj].print(pH1, ostr);
       }
   }
   
   FILE *allout = fopen("allparticles4", "w");
   string stri = ostr.str();
   char *ch = &stri[0];
   fputs(ch, allout);
   fflush(allout);
   fclose(allout);
   
   for (int i = 0; i < THREADNO; i++) {
       delete [] ppa[i];
       delete [] ppa0[i];
   }
   delete [] ppa;
   delete [] ppa0;

}

void Dotest1() {
    std::ifstream fromfile("Halo1.Tpot");
    GalaxyPotential pH(fromfile);
    std::ifstream fromfile2("Disc1.Tpot");
    GalaxyPotential pD(fromfile2);
    particle pp;
    ostringstream ostr1;
    ostringstream ostr2;
    ostringstream ostr22;
    ostringstream ostr3;
    //fourier analysis
    pp.fourier(pH,pD,ostr1,ostr2,ostr3,ostr22);
    FILE *fout1 = fopen("printorbit1", "w");
    string stri1 = ostr1.str();
    char *ch1 = &stri1[0];
    fputs(ch1, fout1);
    fclose(fout1);
    FILE *fout2 = fopen("interpolation", "w");
    string stri2 = ostr2.str();
    char *ch2 = &stri2[0];
    fputs(ch2, fout2);
    fclose(fout2);
    FILE *fout22 = fopen("max", "w");
    string stri22 = ostr22.str();
    char *ch22 = &stri22[0];
    fputs(ch22, fout22);
    fclose(fout22);
    FILE *fout3 = fopen("printfourier", "w");
    string stri3 = ostr3.str();
    char *ch3 = &stri3[0];
    fputs(ch3, fout3);
    fclose(fout3);
   
} 

void Dotest2() {
    std::ifstream fromfile("Halo1.Tpot");
    GalaxyPotential pH(fromfile);
    std::ifstream fromfile2("Disc1.Tpot");
    GalaxyPotential pD(fromfile2);
    ostringstream ostrNa;
    // making particles
    for (int i = 0; i < 500; i++) {
        particle pp;
        pp.gaussianvelocity(pH,pD,ostrNa,8192);
    }
    FILE *foutNa = fopen("printorbitNa", "w");
    string striNa = ostrNa.str();
    char *chNa = &striNa[0];
    fputs(chNa, foutNa);
    fclose(foutNa);
}

void Dotest3() {
    //fourier transform of a SHO
    #define T 10 //periods
    #define TIME 64
    ostringstream ostrf1;
    double positionx[TIME];
    double omg=2.0*PI/T;
    double sqrtmag[TIME];
    double freq[TIME];
    
    int i; double data[2*TIME];
    for (i = 0; i < TIME; i++)
    {
       REAL(data,i) = cos(1.0*i*omg); IMAG(data,i) = 0.0;
    }

    gsl_fft_complex_radix2_forward (data, 1, TIME);
    for (i = 0; i < TIME; i++)
    {
      printf ("%d %e %e\n", i,
              REAL(data,i)/sqrt(TIME),
              IMAG(data,i)/sqrt(TIME));
    }
/*
    for (i = TIME/2; i < TIME; i++){
        positionx[i]=sin(1.0*i*omg);
        sqrtmag[i]=sqrt(REAL(data,i)*REAL(data,i)+IMAG(data,i)*IMAG(data,i));
        freq[i]=-1.0/(T)+i*1.0/(T*TIME);
        //ostrf1<<freq[i]<<" "<<sqrtmag[i]<<" ";
        //ostrf1 << "\n";
    }
*/
    int j;
    for (i = 0; i < TIME; i++)
    {
        if (i<TIME/2){
            j = i+TIME/2;
            positionx[j]=cos(1.0*i*omg);
            sqrtmag[j]=sqrt(REAL(data,j)*REAL(data,j)+IMAG(data,j)*IMAG(data,j));
            freq[j]=-1.0+j*1.0/(TIME);
        }
        else{
            j= i-TIME/2;
            positionx[j]=cos(1.0*i*omg);
            sqrtmag[j]=sqrt(REAL(data,j)*REAL(data,j)+IMAG(data,j)*IMAG(data,j));
            freq[j]=j*1.0/(TIME);
        }
        
        //ostrf1<<i*1.0/(T*TIME)<<" "<<positionx[i]<<" "<<REAL(data,i)/sqrt(TIME)<<" "<<IMAG(data,i)/sqrt(TIME)<<" ";
        //ostrf1<<positionx[j]<<" ";
        ostrf1<<freq[j]<<" "<<sqrtmag[j]/TIME<<" ";
        ostrf1 << "\n";
    }
    
    FILE *foutf1 = fopen("printf1", "w");
    string strif1 = ostrf1.str();
    char *chf1 = &strif1[0];
    fputs(chf1, foutf1);
    fclose(foutf1);
}


int main() {
    Dotest1();

	return 0;
}

/*
#define S 100
int
main (void)
{
  int i;
  double xi, yi, x[S], y[S];
  ostringstream ostrint1;
  ostringstream ostrint2;
  printf ("#m=0,S=17\n");
  
  for (i = 0; i < S; i++)
    {
      x[i] = i + 0.5 * sin (i);
      y[i] = sin(1.0*i)+cos(1.0*i);
      printf ("%g %g\n", x[i], y[i]);
      ostrint1<<x[i]<<" "<<y[i]<<" ";
      ostrint1<< "\n";
    }

  printf ("#m=1,S=0\n");

  
    gsl_interp_accel *acc
      = gsl_interp_accel_alloc ();
    gsl_spline *spline
      = gsl_spline_alloc (gsl_interp_cspline, S);

    gsl_spline_init (spline, x, y, S);
    //double delta=gettime()/M;
    for (xi = x[0]; xi < x[99]; xi += 0.1345632)
      {
        yi = gsl_spline_eval (spline, xi, acc);
        ostrint2<<xi<<" "<<yi<<" ";
        ostrint2<< "\n";
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  
  
    FILE *foutint1 = fopen("printint1", "w");
    string striint1 = ostrint1.str();
    char *chint1 = &striint1[0];
    fputs(chint1, foutint1);
    fclose(foutint1);
    FILE *foutint2 = fopen("printint2", "w");
    string striint2 = ostrint2.str();
    char *chint2 = &striint2[0];
    fputs(chint2,foutint2);
  return 0;
}
*/
