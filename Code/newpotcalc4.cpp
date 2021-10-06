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
#define TIMESTEPFAC 0.008       //0.008   //the time-step factor
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
#define VCG 0.235 //in 1000km/s
#define RCR 6.0 //in kpc

// real and imag for fft
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

#define FROMFILE1 "Halo1.Tpot"
#define FROMFILE2 "Disc1.Tpot"

//Fix the dis potential
#define CONTIME 100000.0 //50000.0
#define FIXTIME 10000.0
#define TIMEXPRESS (FIXTIME/CONTIME) 


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
    return 0.1+time/CONTIME;//Growing potential 
    //return 0.1;//TIMEXPRESS;//Fixed potential, time in units of Myr, when TIMEXPRESS=1.0, we achieve the full MW2017 potential
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
    double growtime = 500.0; //the growth time in units of Myr, can set however i want
    //to slowly grow the amplitude, mimicking a growing bar
    //choose any type of growing function i want
    return 0.0;//A*2.0/PI; //by Simon ;      //*atan(time*time/(growtime*growtime)); //*2.0/PI;
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
    
    const double A = 0.017; const double b =  0.28; const double vc = VCG; const double m = 2.0; const double Rcr = RCR;    //gives the azimuthal angle of bar //in rads
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
    //grad[2] += gradf_z; // make the bar 3-D
}

void fPotGMC(double coord[3], double grad[], GalaxyPotential &pH, GalaxyPotential &pD, double mass, double time) {
	double R = std::sqrt(coord[0] * coord[0] + coord[1] * coord[1]);
	double phi = std::atan2(coord[1], coord[0]);
	double dR; double dz;
    double f2 = timefunc(time);
	pH(R, coord[2], dR, dz);
    double dR2; double dz2;
    pD(R, coord[2], dR2, dz2); 
    grad[0] = -((dR + f2*dR2)*coord[0]/R)*mass;
	grad[1] = -((dR + f2*dR2)*coord[1]/R)*mass;
	grad[2] = -(dz + f2*dz2)*mass;
    //double pot_unper = (pH(R,coord[2]) + f2*pD(R,coord[2])); 
    
}

void fPot2(double coord[3], double grad[], GalaxyPotential &pH, GalaxyPotential &pD2) {
    
    double R = std::sqrt(coord[0] * coord[0] + coord[1] * coord[1]);
	double phi = std::atan2(coord[1], coord[0]);
	double dR; double dz;
	pH(R, coord[2], dR, dz);
    double dR2; double dz2;
    pD2(R, coord[2], dR2, dz2); 
    grad[0] = -(dR+dR2 )*coord[0]/R;
	grad[1] = -(dR+dR2 )*coord[1]/R;
	grad[2] = -(dz+dz2 );
}

void fPotL(double coord[3], double grad[], GalaxyPotential &pH, GalaxyPotential &pD, double time) {
	double R = std::sqrt(coord[0] * coord[0] + coord[1] * coord[1]);
	double phi = std::atan2(coord[1], coord[0]);
	double dR; double dz;
	pH(R, coord[2], dR, dz);
    double dR2; double dz2;
    pD(R, coord[2], dR2, dz2); 
    double f2 = 0.1+time/CONTIME; // 10 per cent initial mass
    grad[0] = -(dR + f2*dR2)*coord[0]/R;
	grad[1] = -(dR + f2*dR2)*coord[1]/R;
	grad[2] = -(dz + f2*dz2);
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

double gPotconst(double coord[3], double R, GalaxyPotential &pH1, GalaxyPotential &pH2) {
    double f2 = 0.1;// 0.1 for drag; 1.0 for epicycle;
    double potunper = (pH1(R,coord[2]) + f2*pH2(R,coord[2]));
    return potunper;
}

double gPottime(double coord[3], double R, GalaxyPotential &pH1, GalaxyPotential &pH2, double time) {
    double f2 = 0.1 + time/CONTIME;//full potential
    double potunper = (pH1(R,coord[2]) + f2*pH2(R,coord[2]));
    return potunper;
}

double effPot(double coord[3], double R, double lz, GalaxyPotential &pH1, GalaxyPotential &pH2) {
    double f2 = 1.0;//full potential
    double potunper = (pH1(R,coord[2]) + f2*pH2(R,coord[2]));
    double poteff = potunper + lz*lz/(2*R*R); 
    return poteff;
}

double effPottime(double coord[3], double R, double lz, GalaxyPotential &pH1, GalaxyPotential &pH2, double time) {
    double f2 = 0.1 + time/CONTIME;//full potential
    double potunper = (pH1(R,coord[2]) + f2*pH2(R,coord[2]));
    double poteff = potunper + lz*lz/(2*R*R); 
    return poteff;
}

class clouds {
public:
    double weight = 10000000.0;
    double x[DIMENSION];
    double p[DIMENSION];
    double f[DIMENSION];
private:
    ;
};


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
    void kickGMC(GalaxyPotential &pH, GalaxyPotential &pD, double mass);    
    void stepGMC(GalaxyPotential &pH, GalaxyPotential &pD, double mass);
    void kick2(GalaxyPotential &pH, GalaxyPotential &pD2);    
    void step2(GalaxyPotential &pH, GalaxyPotential &pD2);
    void kickL(GalaxyPotential &pH, GalaxyPotential &pD, double time);    
    void stepL(GalaxyPotential &pH, GalaxyPotential &pD, double time);
    void leapfrog(GalaxyPotential &pH, GalaxyPotential &pD);
    void leapfrogGMC(GalaxyPotential &pH, GalaxyPotential &pD, double mass);
    void leapfrog2(GalaxyPotential &pH, GalaxyPotential &pD2);
    void leapfrogL(GalaxyPotential &pH, GalaxyPotential &pD, double time);
    void print(GalaxyPotential &pH, GalaxyPotential &pD);
    void print(GalaxyPotential &pH, GalaxyPotential &pD, ostringstream &ostr);
//  void Nprint(GalaxyPotential &pH, GalaxyPotential &pD, ostringstream &ostr);
//    void getsimplestartperiods(GalaxyPotential &pH, GalaxyPotential &pD, periods[2], double maxtime);
    void integrate(double integtime, GalaxyPotential &pH, GalaxyPotential &pD);   
    void integrateGMC(double integtime, GalaxyPotential &pH, GalaxyPotential &pD, double mass);  
    void integrate2(double integtime, GalaxyPotential &pH, GalaxyPotential &pD2);
    void integrateL(double integtime, GalaxyPotential &pH, GalaxyPotential &pD, double time);  
    void phishift(std::uniform_real_distribution<double> phasedist);
    void phishift(double delphi);
    void phiset(double R, double phi, double vR, double vphi);
    void phitrans(double &R, double &phi, double &vR, double &vphi);
    void clone(particle &pp);
    void epicycle(GalaxyPotential &pH1, GalaxyPotential &pH2, ostringstream &ostr1, ostringstream &ostr21, ostringstream &ostr22, ostringstream &ostr31, ostringstream &ostr32);
    void epicycletime(GalaxyPotential &pH1, GalaxyPotential &pH2, ostringstream &ostr1, ostringstream &ostr21, ostringstream &ostr22, ostringstream &ostr31, ostringstream &ostr32);
    void fourier(GalaxyPotential &pH1, GalaxyPotential &pH2, ostringstream &ostr1, ostringstream &ostr11, ostringstream &ostr12, ostringstream &ostr2, ostringstream &ostr21, ostringstream &ostr3);
    void sos(GalaxyPotential &pH1, GalaxyPotential &pH2, ostringstream &ostr1, ostringstream &ostr11, ostringstream &ostr12, ostringstream &ostr2);
    void dragging(GalaxyPotential &pH1, GalaxyPotential &pH2, GalaxyPotential &pH3, ostringstream &ostr1,ostringstream &ostr11,ostringstream &ostr2,ostringstream &ostr21,ostringstream &ostr3);
    void simplelevitation(GalaxyPotential &pH1, GalaxyPotential &pH2, GalaxyPotential &pH3, ostringstream &ostr);
    void GMClevitation(GalaxyPotential &pH1, GalaxyPotential &pH2, GalaxyPotential &pH3, ostringstream &ostr1, ostringstream &ostr11);
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

void particle::kickGMC(GalaxyPotential &pH, GalaxyPotential &pD, double mass) {
    fPotGMC(x, f, pH, pD, mass, gettime());
    for (int i = 0; i < DIMENSION; i++) {
        p[i] += timestep*f[i];
    }
}    

void particle::kick2(GalaxyPotential &pH, GalaxyPotential &pD2) {
    fPot2(x, f, pH, pD2);
    for (int i = 0; i < DIMENSION; i++) {
        p[i] += timestep*f[i];
    }
}    

void particle::kickL(GalaxyPotential &pH, GalaxyPotential &pD, double time) {
    fPotL(x, f, pH, pD, time);
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

void particle::stepGMC(GalaxyPotential &pH, GalaxyPotential &pD, double mass) {
    drift();
    timefine += htimestep;
    kickGMC(pH, pD, mass);
    drift();
    timefine += htimestep;
}

void particle::step2(GalaxyPotential &pH, GalaxyPotential &pD2) {
    drift();
    timefine += htimestep;
    kick2(pH, pD2);
    drift();
    timefine += htimestep;
}

void particle::stepL(GalaxyPotential &pH, GalaxyPotential &pD, double time) {
    drift();
    timefine += htimestep;
    kickL(pH, pD, time);
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

void particle::leapfrogGMC(GalaxyPotential &pH, GalaxyPotential &pD,double mass) {
    stepnumber++;
	double grad[3];
	for (int i = 0; i < 3; ++i) {
		grad[i] = 0.0;
	}
    
        //difference here is the fPot taking 5 arguments
    //taking in the coordinates, force, pH, pD, and gettime
	fPotGMC(x, grad, pH, pD, mass,gettime());
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

void particle::leapfrog2(GalaxyPotential &pH, GalaxyPotential &pD2) {
    stepnumber++;
	double grad[3];
	for (int i = 0; i < 3; ++i) {
		grad[i] = 0.0;
	}
    
        //difference here is the fPot taking 5 arguments
    //taking in the coordinates, force, pH, pD, and gettime
	fPot2(x, grad, pH, pD2);
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
        step2(pH, pD2);  //here calling the step function with 2 potentials
    }
    timefineadd();

}//=======================================================================================================

void particle::leapfrogL(GalaxyPotential &pH, GalaxyPotential &pD, double time) {
    stepnumber++;
	double grad[3];
	for (int i = 0; i < 3; ++i) {
		grad[i] = 0.0;
	}
    
        //difference here is the fPot taking 5 arguments
    //taking in the coordinates, force, pH, pD, and gettime
	fPotL(x, grad, pH, pD, time);
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
        stepL(pH, pD, time);  //here calling the step function with 2 potentials
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

void particle::epicycle(GalaxyPotential &pH1, GalaxyPotential &pH2, ostringstream &ostr1, ostringstream &ostr21, ostringstream &ostr22, ostringstream &ostr31, ostringstream &ostr32){ 
    #define NR 1500 // number of R
    #define NZ 200// number of Z
    #define NRG 12 // number of Rg

   //double Rg = 8.0;
    double delR = 0.001;
    double delz = 0.0001;
    double delRg = 0.1;
/*    // for circular velocity
    double pot [NRG][NR];
    double vc [NRG];
    double dpotdR [NRG][NR];
    // for epicycle frequencies
    double Rg [NRG];
    double Lz [NRG];
    double R [NRG][NR];
    double effpot [NRG][NR];
    double dphidR [NRG][NR];
    double kap [NRG][NR];
    
    // not working
    double z [NRG][NZ];
    double effpotz [NRG][NZ];
    double dphidz [NRG][NZ];
    double nu [NRG][NZ];
*/
    double **pot = new double*[NRG];   //the stars mean pointers, the **lvplane is a pointer to pointers which gets initialised as two-hundred pointers
    double *vc = new double[NRG];
    double **dpotdR = new double*[NRG];
    double *Rg = new double[NRG];
    double *Lz = new double[NRG];
    double **R = new double*[NRG];
    double **effpot = new double*[NRG];
    double **dphidR = new double*[NRG];
    double **kap = new double*[NRG];
    double **z = new double*[NRG];
    double **potz = new double*[NRG];
    double **dphidz = new double*[NRG];
    double **nu = new double*[NRG];
    for (int i = 0; i < NRG; i++) { //this is a loop using the index i running from 0 to 199
        pot[i] = new double[NR];
        dpotdR[i] = new double[NR];
        R[i] = new double[NR];
        effpot[i] = new double[NR];
        dphidR[i] = new double[NR];
        kap[i] = new double[NR];
        z[i] = new double[NZ];
        potz[i] = new double[NZ];
        dphidz[i] = new double[NZ];
        nu[i] = new double[NZ];
    }
    // find circular velocities vc(Rg)
    int i; int j;
    for (i = 0; i < NRG; i++){
        Rg[i]=i*delRg;
        for (j = 0; j < NR; j++){
            settime();
            x[0] = 0.0 + delR* j; x[1] = 0.0; x[2] = 0.0;
            double r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
            pot[i][j] = gPotconst(x, r, pH1, pH2);
        }
    }
    for (i = 0; i < NRG; i++){
        for (j = 0; j < NR-1; j++){
            dpotdR[i][j] = (pot[i][j+1] - pot[i][j])/delR;

        }
    }
    for (i = 0; i < NRG; i++){
        int nRg = Rg[i]/delR;
        vc[i]=std::sqrt(Rg[i]*dpotdR[i][nRg]);
        ostr1 << Rg[i] << " "<< vc[i] << " ";
        ostr1<<"\n";
    }

    /////////////////////////////////////////////
    // radial epicycle frequency
    for (i = 0; i < NRG; i++){
        Lz[i]= Rg[i] * vc[i];//(8.29kpc*(233.14/977.77)kpc/Myr)
        for (j = 0; j < NR; j++){
            settime();
            x[0] = 0.0 + delR* j; x[1] = 0.0; x[2] = 0.0;
            double r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
            double lz = Lz [i];
            R[i][j] = r;
            effpot[i][j] = effPot(x, r, lz, pH1, pH2);
            ostr21<<R[i][j]<<" "<<effpot[i][j]<<" ";
            ostr21<<"\n";
        }
    }


    //differentiation to find kap
    for (i = 0; i < NRG; i++){
        for (j = 0; j < NR-1; j++){
            dphidR[i][j] = (effpot[i][j+1] - effpot[i][j])/delR;
            //ostr2 << dphidR [i][j] << " "; 
            //ostr2<<"\n";
        }
    }
    for (i = 0; i < NRG; i++){
        for (j = 0; j < NR-1; j++){
            kap[i][j] = std::sqrt((dphidR[i][j+1] - dphidR[i][j])/delR);
            //ostr2 << kap[i][j] << " "; 
            //ostr2<<"\n";
        }
    }
    for (i = 0; i < NRG; i++){
        int nRg = Rg[i]/delR;
        ostr22 << Rg[i] << " "<< kap[i][nRg] << " ";
        ostr22<<"\n";
    }
    
    /////////////////////////////////////////////
    // vertical epicycle frequency
    //double coord[3];
    for (i = 0; i < NRG; i++){
        for (j = 0; j < NZ; j++){
            settime();
            x[0] = Rg[i]; x[1] = 0.0; x[2] = -0.1+delz * ((double)(j));
            double r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
            z[i][j] = x[2];
            //coord[0] = Rg[i]; coord[1] = 0.0; coord[2] =0.0;
            potz[i][j] = gPotconst(x, r, pH1, pH2);//-gPotconst(coord, r, pH1, pH2);
            ostr31<<z[i][j]<<" "<<potz[i][j]<<" ";
            ostr31<<"\n";
        }
    }

    //differentiation to find nu
    //#define NR 1500 // number of R #define NZ 1000 // number of Z #define NRG 120
    for (i = 0; i < NRG; i++){
        for (j = 0; j < NZ-1; j++){
            dphidz[i][j] = (potz[i][j+1] - potz[i][j])/delz;
            //ostr32 << dphidz [i][j] << " "; 
            //ostr32<<"\n";
        }
    }

    
    for (i = 0; i < NRG; i++){
        for (j = 0; j < NZ-1; j++){
            nu[i][j] =std::sqrt(fabs(dphidz[i][j+1] - dphidz[i][j])/delz);

            //ostr32 << nu[i][j] << " "; 
            //ostr32<<"\n";
        }
    }

    
    int nz = NZ*delz/2; // find nz when z=0
    
    for (i = 0; i < NRG; i++){
        int nRg = Rg[i]/delR;
        ostr32 << Rg[i] << " "<< kap[i][nRg]<< " "<<nu[i][nz] << " ";
        ostr32<<"\n";
    }


    for (i = 0; i < NRG; i++){
       
        delete [] pot[i];      //free all memory at all pointers in lvplane
        delete [] dpotdR[i];
        delete [] R[i];
        delete [] effpot[i];
        delete [] dphidR[i];
        delete [] kap[i];
        delete [] z[i];
        delete [] potz[i];
        delete [] dphidz [i];
        delete [] nu[i];
            
        }

    delete [] pot;      //free all memory at all pointers in lvplane
    delete [] vc;
    delete [] dpotdR;
    delete [] Rg;
    delete [] Lz;
    delete [] R;
    delete [] effpot;
    delete [] dphidR;
    delete [] kap;
    delete [] z;
    delete [] potz;
    delete [] dphidz ;
    delete [] nu;

}


void particle::epicycletime(GalaxyPotential &pH1, GalaxyPotential &pH2, ostringstream &ostr1, ostringstream &ostr21, ostringstream &ostr22, ostringstream &ostr31, ostringstream &ostr32){ 
    #define A 1000
    
/*
    
    double T [NRG][A];
    // for circular velocity
    double pot [NRG][A][NR];
    double vc [NRG][A];
    double dpotdR [NRG][A][NR];
    // for epicycle frequencies

    double Lz [NRG][A];
    double R [NRG][A][NR];
    double effpot [NRG][A][NR];
    double dphidR [NRG][A][NR];
    double kap [NRG][A][NR];

    double z [NRG][A][NZ];
    double effpotz [NRG][A][NZ];
    double dphidz [NRG][A][NZ];
    double nu [NRG][A][NZ];
*/    
    //#define NR 1500 // number of R #define NZ 200 // number of Z
    //double Rg = 8.0;
    double delRg = 1.0;
    double delR = 0.01;
    double delz = 0.001;
    double deltime = 2.0;
    
    double *Rg = new double[NRG];
    double ***pot = new double**[NRG];   //the stars mean pointers, the **lvplane is a pointer to pointers which gets initialised as two-hundred pointers
    double **vc = new double*[NRG];
    double ***dpotdR = new double**[NRG];
    double **T = new double*[NRG];
    double **Lz = new double*[NRG];
    double ***R = new double**[NRG];
    double ***effpot = new double**[NRG];
    double ***dphidR = new double**[NRG];
    double ***kap = new double**[NRG];
    double ***z = new double**[NRG];
    double ***effpotz = new double**[NRG];
    double ***dphidz = new double**[NRG];
    double ***nu = new double**[NRG];
    int i; int ii; int j;
    for (int i = 0; i < NRG; i++) {
        pot[i] = new double*[A];
        vc[i] = new double[A];
        dpotdR[i] = new double*[A];
        T[i] = new double[A];
        Lz[i] = new double[A];
        R[i] = new double*[A];
        effpot[i] = new double*[A];
        dphidR[i] = new double*[A];
        kap[i] = new double*[A];
        z[i] = new double*[A];
        effpotz[i] = new double*[A];
        dphidz[i] = new double*[A];
        nu[i] = new double*[A];
        for (int ii = 0; ii < A; ii++) { 
            pot[i][ii] = new double[NR];
            dpotdR[i][ii] = new double[NR];
            R[i][ii] = new double[NR];
            effpot[i][ii] = new double[NR];
            dphidR[i][ii] = new double[NR];
            kap[i][ii] = new double[NR];
            z[i][ii] = new double[NZ];
            effpotz[i][ii] = new double[NZ];
            dphidz[i][ii] = new double[NZ];
            nu[i][ii] = new double[NZ];
        }
    }
    
    // find circular velocities vc(Rg)
    
    for (i = 0; i < NRG; i++){
        Rg[i]=i*delRg;
        for (ii = 0; ii < A; ii++){
            int time = deltime*ii;
            T[i][ii] = time;
            for (j = 0; j < NR; j++){
                settime();
                x[0] = 0.0 + delR* j; x[1] = 0.0; x[2] = 0.0;
                double r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
                pot[i][ii][j] = gPottime(x, r, pH1, pH2, time);
                //ostr1 << pot[i][ii][j] << " ";
                //ostr1<<"\n";
            }
        }
    }

     
    for (i = 0; i < NRG; i++){
        for (ii = 0; ii < A; ii++){
            for (j = 0; j < NR-1; j++){
                dpotdR[i][ii][j] = (pot[i][ii][j+1] - pot[i][ii][j])/delR;

            }
        }
    }
    for (i = 0; i < NRG; i++){
        Rg[i]=i*delRg;
        for (ii = 0; ii < A; ii++){
            int nRg = Rg[i]/delR;
            vc[i][ii]=std::sqrt(Rg[i]*dpotdR[i][ii][nRg]);
            ostr1 << Rg[i] << " "<< T[i][ii] << " "<< vc[i][ii] << " ";
            ostr1<<"\n";
        }
    }
   
     /////////////////////////////////////////////
    // radial epicycle frequency
    for (i = 0; i < NRG; i++){
        Rg[i]=i*delRg;
        for (ii = 0; ii < A; ii++){
            Lz[i][ii]= Rg[i] * vc[i][ii];//(8.29kpc*(233.14/977.77)kpc/Myr)
            int time = deltime*ii;
            for (j = 0; j < NR; j++){
                settime();
                x[0] = 0.0 + delR* j; x[1] = 0.0; x[2] = 0.0;
                double r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
                double lz = Lz [i][ii];
                R[i][ii][j] = r;
                effpot[i][ii][j] = effPottime(x, r, lz, pH1, pH2, time);
                //ostr21<<R[i][ii][j]<<" "<<effpot[i][ii][j]<<" ";
                //ostr21<<"\n";
            }
        }

    }
    //differentiation to find kap
    for (i = 0; i < NRG; i++){
        for (ii = 0; ii < A; ii++){
            for (j = 0; j < NR-1; j++){
                dphidR[i][ii][j] = (effpot[i][ii][j+1] - effpot[i][ii][j])/delR;
                //ostr2 << dphidR [i][ii][j] << " "; 
                //ostr2<<"\n";
            }
        }
    }
    for (i = 0; i < NRG; i++){
        for (ii = 0; ii < A; ii++){
            for (j = 0; j < NR-1; j++){
                kap[i][ii][j] = std::sqrt((dphidR[i][ii][j+1] - dphidR[i][ii][j])/delR);
                //ostr2 << kap[i][ii][j] << " "; 
                //ostr2<<"\n";
            }
        }
    }
    for (i = 0; i < NRG; i++){
        Rg[i]=i*delRg;
        for (ii = 0; ii < A; ii++){
            int time = deltime*ii;
            T[i][ii] = time;
            int nRg = Rg[i]/delR;
            ostr22 << T[i][ii] << " "<< kap[i][ii][nRg] << " ";
            ostr22<<"\n";
        }
    }
    

    /////////////////////////////////////////////
    // vertical epicycle frequency
    for (i = 0; i < NRG; i++){
        Rg[i]=i*delRg;
        for (ii = 0; ii < A; ii++){
            Lz[i][ii]= Rg[i] * vc[i][ii];//(8.29kpc*(233.14/977.77)kpc/Myr)
            int time = deltime*ii;
            for (j = 0; j < NZ; j++){
                settime();
                x[0] = Rg[i]; x[1] = 0.0; x[2] = -0.1+delz * j;
                z[i][ii][j]= x[2];
                double lz = Lz [i][ii];
                double r = x[0];
                effpotz[i][ii][j] = effPottime(x, r, lz, pH1, pH2, time);
                //ostr21<<z[i][ii][j]<<" "<<effpotz[i][ii][j]<<" ";
                //ostr21<<"\n";
            }
        }

    }
/*
    for (i = 0; i < 1; i++){
        Rg[i]=8.0;
        for (ii = 0; ii < A; ii++){
            Lz[i][ii]= Rg[i] * vc[i][ii];//(8.29kpc*(233.14/977.77)kpc/Myr)
            int time = deltime*ii;
            T[i][ii] = time;
            for (j = 0; j < NZ; j++){
                settime();
                x[0] = Rg[i]; x[1] = 0.0; x[2] = -0.1+delz * j;
                z[i][ii][j]= x[2];
                double lz = Lz [i][ii];
                double r = x[0];
                effpotz[i][ii][j] = effPottime(x, r, lz, pH1, pH2, time);
                ostr31<<T[i][ii]<<" "<<z[i][ii][j]<<" "<<effpotz[i][ii][j]<<" ";
                ostr31<<"\n";
            }
        }

    }

*/
    //differentiation to find nu
    //#define NR 1500 // number of R #define NZ 1000 // number of Z #define NRG 120
    for (i = 0; i < NRG; i++){
        for (ii = 0; ii < A; ii++){
            for (j = 0; j < NZ-1; j++){
                dphidz[i][ii][j] = (effpotz[i][ii][j+1] - effpotz[i][ii][j])/delz;
                //ostr2 << dphidz[i][ii][j] << " "; 
                //ostr2<<"\n";
            }
        }
    }
    for (i = 0; i < NRG; i++){
        for (ii = 0; ii < A; ii++){
            for (j = 0; j < NZ-1; j++){
                nu[i][ii][j] = std::sqrt((dphidz[i][ii][j+1] - dphidz[i][ii][j])/delz);
                //ostr2 << nu[i][ii][j] << " "; 
                //ostr2<<"\n";
            }
        }
    }

    
    int nz = NZ*delz/2; // find nz when z=0
    for (i = 0; i < NRG; i++){
        Rg[i]=i*delRg;
        for (ii = 0; ii < A; ii++){
            int nRg = Rg[i]/delR;
            ostr31 << Rg[i] << " "<< T[i][ii] << " "<< kap[i][ii][nRg]<< " "<<nu[i][ii][nz] << " ";
            ostr31<<"\n";
        }
    }
    
    for (i = 0; i < NRG; i++){
        Rg[i]=8.0;
        for (ii = 0; ii < A; ii++){
            int nRg = Rg[i]/delR;
            ostr32 << Rg[i] << " "<< T[i][ii] << " "<< kap[8][ii][nRg]<< " "<<nu[8][ii][nz] << " ";
            ostr32<<"\n";
        }
    }


    for (i = 0; i < NRG; i++){
        for (ii = 0; ii < A; ii++){
            delete [] pot[i][ii];      //free all memory at all pointers in lvplane
            delete [] dpotdR[i][ii];
            delete [] R[i][ii];
            delete [] effpot[i][ii];
            delete [] dphidR[i][ii];
            delete [] kap[i][ii];
            delete [] z[i][ii];
            delete [] effpotz[i][ii];
            delete [] dphidz[i][ii];
            delete [] nu[i][ii];
            
        }
        delete [] pot[i];      //free all memory at all pointers in lvplane
        delete [] vc[i]; 
        delete [] dpotdR[i];
        delete [] T[i];
        delete [] Lz[i];
        delete [] R[i];
        delete [] effpot[i];
        delete [] dphidR[i];
        delete [] kap[i];
        delete [] z[i];
        delete [] effpotz[i];
        delete [] dphidz[i];
        delete [] nu[i];
        
    }
   
    delete [] pot;      //free all memory at all pointers in lvplane
    delete [] vc;
    delete [] dpotdR;
    delete [] T;
    delete [] Rg;
    delete [] Lz;
    delete [] R;
    delete [] effpot;
    delete [] dphidR;
    delete [] kap;
    delete [] z;
    delete [] effpotz;
    delete [] dphidz ;
    delete [] nu;

}
void particle::fourier(GalaxyPotential &pH1, GalaxyPotential &pH2, ostringstream &ostr1, ostringstream &ostr11, ostringstream &ostr12,ostringstream &ostr2, ostringstream &ostr21, ostringstream &ostr3){    
    //Instructions for the gnuplot: 
    //(printorbit1) 1:5=time againsts r; 1:6=time againsts z
    //(resonance) plot "printfourier" u 1:2 title "r" with lines, "printfourier" u 1:($3*5.0) title "z" with lines
    //(Phase space) plot "printorbit1" u 7:9 title "" with lines 
    //(surface of sections) plot "printorbit1" u 5:7:9 title "" with lines palette, "SOS" u 5:10:9 title "" with points palette
    #define M 16384
    double settime();
    //generate the data of an orbital trajectory
    double t[M];
    double rinitial[M];
    double zinitial[M];
    double periods[2];
    double maxtime = 10000.0;
    double magnitude=std::sqrt(M);
    
    //Potential
    double coord[3]; 
    coord[0] = 8.0; coord[1] = 0.0; coord[2] = 0.0;//The velocities of an orbit required to find the potential
    double potential = gPot(coord, pH1, pH2, gettime());
    //Energy
    double px = 0.1; double py = 0.22; double pz = 0.0; // The velocities of an orbit required to find E
    double E = 0.5*(px*px+py*py+pz*pz) + potential; //Energy is constant. E=0.5*(vx^2+vy^2+vz^2)+potential
        
    double number = 1.0;//100.0;
    double del = 0.000125;// delta is the difference between each orbit in p[0]
    //near-circular orbit in 10% of McMillian17 disc, x is in kpc, p is in Mm/s
    //fix the x position and vary the radial velocity 
    x[0] = 8.0; x[1] = 0.0; x[2] = 0.0; 
    
    p[0] = 0.0; p[1] = 0.177861; p[2] = 0.01; //resonant
/*
    double potential1=gPot(x, pH1, pH2, gettime());
    
    p[1] = py; 
    p[2] = sqrt((E-potential1)*2.0 -p[1]*p[1]-0.0)-del*145.0;//-del*145.0 (1:1 MOR) -del*585.0 (0.664921 2:3 MOR)
    //-del*585.0 (0.664921 2:3 MOR), -del*578.0 (0.668421), -del*566.0 (0.679144)
    p[0] = sqrt((E-potential1)*2.0 -p[1]*p[1] -p[2]*p[2]); //sqrt((E-potential)*2.0 -p[1]*p[1])
*/

   
    // intital conditions for the surface of sections
    double oldtime=0.0;
    double oldx=x[0]; 
    double oldy=x[1];
    double oldz=x[2];
    double oldvx = p[0];
    double oldvy = p[1];
    double oldvz = p[2];
    double oldR,oldphi,oldvR,oldvphi;
    phitrans(oldR,oldphi,oldvR,oldvphi);
    double oldlibratephi = 4.0;//atan2(1000.0*oldvR,(oldR-8.0));
    
    t[0]=0.0;
    rinitial[0]=x[0];
    zinitial[0]=0.0;

    double Ri,phii,vRi,vphii;
    phitrans(Ri,phii,vRi,vphii);
    ostr1<<0.0<<" "<<x[0]<<" "<< 0.0<<" "<<0.0<<" ";
    ostr1<<Ri<<" "<<x[2]<<" "<<vRi<<" "<<p[0]<<" "<<p[2]<<" "<<vphii<<" ";
    ostr1 << "\n";
    
    integrate(maxtime, pH1, pH2);
    for (int i = 1; i < M; i++) {
        for (int ii = 0; ii < 1; ii++) {
            leapfrog(pH1, pH2);
        }
        
        double R,phi,vR,vphi;
        phitrans(R,phi,vR,vphi);
        double libratephi = atan2(1000.0*vR,(R-8.0)); // angle of libration
        double vx = p[0];
        double vy = p[1];
        double vz = p[2];
        t[i] = gettime()-10000.0;
        rinitial[i] = R;
        zinitial[i] = x[2];
        
        // ostri1 is the printorbit1
        ostr1<<gettime()-10000.0<<" "<<x[0]<<" "<< x[1]<<" "<<x[2]<<" ";
        ostr1<<rinitial[i]<<" "<<zinitial[i]<<" "<<vR<<" "<<vx<<" "<<vz<<" "<<vphi<<" ";
        ostr1 << "\n";
    
        // detect the change of sign
        if (oldz*x[2] < 0) { 
            if (x[2] >= 0){ // Select z with positive sign
            
            //ostr11 is the SOS, surface of sections
            ostr11<<gettime()-10000.0<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<R<<" "<<phi<<" "; // positions
            ostr11<<vx<<" "<<vy<<" "<<vz<<" "<<vR<<" "<<vphi<<" "<<libratephi<<" "; // velocities
            ostr11 << "\n";
            }
                if (p[2] >= 0){
                        //ostr12 is the SOS of positive/rising vz
                        ostr12<<gettime()-10000.0<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<R<<" "<<phi<<" "; // positions
                        ostr12<<vx<<" "<<vy<<" "<<vz<<" "<<vR<<" "<<vphi<<" "<<libratephi<<" "; // velocities
                        ostr12<<i+1.0<<" "<<E<<" "<<potential<<" "<< "\n";
                    }    
            else {
                ostr11<<oldtime<<" "<<oldx<<" "<<oldy<<" "<<oldz<<" "<<oldR<<" "<<oldphi<<" "; // positions
                ostr11<<oldvx<<" "<<oldvy<<" "<<oldvz<<" "<<oldvR<<" "<<oldvphi<<" "<<oldlibratephi<<" "; // velocities
                ostr11 << "\n";
            }           
                if (p[2] >= 0){
                        //ostr12 is the SOS of positive/rising vz
                        ostr12<<gettime()-10000.0<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<R<<" "<<phi<<" "; // positions
                        ostr12<<vx<<" "<<vy<<" "<<vz<<" "<<vR<<" "<<vphi<<" "<<oldlibratephi<<" "; // velocities
                        ostr12<<i+1.0<<" "<<E<<" "<<potential<<" "<< "\n";
                    }    
        } 
        oldtime=gettime()-10000.0;oldx=x[0]; oldy=x[1]; oldz=x[2]; oldR=R; oldphi=phi;
        oldvx=vx; oldvy=vy; oldvz=vz; oldvR=vR; oldvphi=vphi;
    }
    
    //interpolate between points to have fixed time intervals
    int i=0;
    double totaltime = gettime()-10000.0;
    double delta=(totaltime)/M;
    double ti, ri, zi;
    double rfinal[2*M];
    double zfinal[2*M];
    double thetafinal[2*M];//theta=arctan(z/r)
    
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
    
    // test for interpolations for r and z
    for (ti = t[0]; ti < t[M-1]; ti += delta)
      {
        ri = gsl_spline_eval (spliner, ti, accr);
        zi = gsl_spline_eval (splinez, ti, accz);
      }
    
    
    //ostr2 gives the interpolation
    for (int i = 0; i < M; i++) {
        ti=1.0*delta*i;
        REAL(rfinal,i) = gsl_spline_eval (spliner, ti, accr); IMAG(rfinal,i) = 0.0;
        REAL(zfinal,i) = gsl_spline_eval (splinez, ti, accz); IMAG(zfinal,i) = 0.0;
        REAL(thetafinal,i) = atan2(REAL(zfinal,i),REAL(rfinal,i)); IMAG(thetafinal,i) = 0.0;
        //ostr2<<rfinal<<" "<<zfinal<<" ";
        ostr2<<i*delta<<" "<<REAL(rfinal,i)<<" "<<REAL(zfinal,i)<<" "<<REAL(thetafinal,i)<<" ";
        ostr2 << "\n";
      }
    //free the storage
    gsl_spline_free (spliner);
    gsl_interp_accel_free (accr);
    gsl_spline_free (splinez);
    gsl_interp_accel_free (accz);
    
    
    //Fourier analysis
    //test1: fourier transform the fist half, then fourier transform the late half.
    //Fourier analysis
    gsl_fft_complex_radix2_forward (rfinal, delta, M);
    gsl_fft_complex_radix2_forward (zfinal, delta, M);
    gsl_fft_complex_radix2_forward (thetafinal, delta, M);
    
    // store fr, fz, magnitudes 
    double frequency[M];
    double sqrtmagr[M];
    double sqrtmagz[M];
    double sqrtmagtheta[M];
    
    //find maximum magnitude
    double magr[M];
    double magz[M];
    double magtheta[M];
    double magrmax=0.0;
    double magzmax=0.0;
    double magthetamax=0.0;
    double nr=0.0;//the location when magrmax occurs
    double nz=0.0;//the location when magzmax occurs
    double ntheta=0.0;//the location when magthetamax occurs

    //Find the maximum magnitudes and corresponding frequencies to understand the MOR in the first half
    for (i = 1; i < M/2; i++){
        magr[i]=sqrt(REAL(rfinal,i)*REAL(rfinal,i)+IMAG(rfinal,i)*IMAG(rfinal,i));
        magz[i]=sqrt(REAL(zfinal,i)*REAL(zfinal,i)+IMAG(zfinal,i)*IMAG(zfinal,i));
        magtheta[i]=sqrt(REAL(thetafinal,i)*REAL(thetafinal,i)+IMAG(thetafinal,i)*IMAG(thetafinal,i));
    
        if (magr[i] > magrmax) {
            magrmax = magr[i];
            nr = i;
        } 
        
        if (fabs(magz[i]) > magzmax) {
            magzmax = fabs(magz[i]);
            nz = i;
        }
        if (fabs(magtheta[i]) > magthetamax) {
            magthetamax = fabs(magtheta[i]);
            ntheta = i;
        }
        if (i==(M/2-1)){
            double fr = nr*1.0/(M);
            double fz = nz*1.0/(M);
            double ftheta = ntheta*1.0/(M);
            ostr21 << fr << " "<< ftheta<< " ";
            ostr21<< magrmax/M << " "<< magthetamax/M << " "<< (nr/nz)<< " ";
            //<< magrmax/M << " " << magzmax/M << " "<< magthetamax/M << " ";
            ostr21 << "\n";
        }
    }
        
    //Fourier shift (centre)
    int j;
    for (i = 0; i < M; i++)
    {
        // shifting the second half to the negative side
        if (i<M/2){
            j = i+M/2;
            sqrtmagr[j]=sqrt(REAL(rfinal,j)*REAL(rfinal,j)+IMAG(rfinal,j)*IMAG(rfinal,j));
            sqrtmagz[j]=sqrt(REAL(zfinal,j)*REAL(zfinal,j)+IMAG(zfinal,j)*IMAG(zfinal,j));
            sqrtmagtheta[j]=sqrt(REAL(thetafinal,j)*REAL(thetafinal,j)+IMAG(thetafinal,j)*IMAG(thetafinal,j));
            frequency[j]=-1.0+j*1.0/(M);
        }
        else{
            j= i-M/2;
            sqrtmagr[j]=sqrt(REAL(rfinal,j)*REAL(rfinal,j)+IMAG(rfinal,j)*IMAG(rfinal,j));
            sqrtmagz[j]=sqrt(REAL(zfinal,j)*REAL(zfinal,j)+IMAG(zfinal,j)*IMAG(zfinal,j));
            sqrtmagtheta[j]=sqrt(REAL(thetafinal,j)*REAL(thetafinal,j)+IMAG(thetafinal,j)*IMAG(thetafinal,j));
            frequency[j]=j*1.0/(M);
        }
        ostr3<<frequency[j]<<" "<<sqrtmagr[j]/M<<" "<<sqrtmagz[j]/M<<" "<<sqrtmagtheta[j]/M<<" ";
        ostr3 << "\n";
    }

}

void particle::sos(GalaxyPotential &pH1, GalaxyPotential &pH2, ostringstream &ostr1,ostringstream &ostr11,ostringstream &ostr12, ostringstream &ostr2){    
    //(surface of sections) plot "printorbit1" u 5:7:9 title "" with lines palette, "SOSN" u 5:10:9 title "" with points palette

    double t[M];
    double rinitial[M];
    double zinitial[M];    
    double rfinal[2*M];
    double zfinal[2*M];
    double thetafinal[2*M];//theta=arctan(z/r)
    
    //Potential
    double coord[3]; 
    coord[0] = 8.0; coord[1] = 0.0; coord[2] = 0.0;//The velocities of an orbit required to find the potential
    double potential = gPot(coord, pH1, pH2, gettime());
    //Energy
    double px = 0.0; double py = 0.177861; double pz = 0.001; // The velocities of an orbit required to find E
    double E = 0.5*(px*px+py*py+pz*pz) + potential; //Energy is constant. E=0.5*(vx^2+vy^2+vz^2)+potential
        
    int i;int j;
    double number = 100.0;//100.0;
    double del = 0.001;// delta is the difference between each orbit in p[0]

    for (j = 0; j < number; j++){
        
        settime();
        double periods[2];
        double maxtime = 10000.0;
        t[0]=0.0;
        
        //constant potential, from maximum pz to maximum pR  
        x[0] = 8.0; x[1] = 0.0; x[2] = 0.0;
        double potential1=gPot(x, pH1, pH2, gettime());
        p[0] = 0.0; 
        p[2] = sqrt((E-potential1)*2.0 -py*py)+del*j;//;//from maximum pz
        p[1] = sqrt((E-potential1)*2.0 -p[0]*p[0] -p[2]*p[2]); //sqrt((E-potential)*2.0 -p[1]*p[1])
                
/*       
        //constant potential, from maximum pz to maximum pR  
        x[0] = 8.0; x[1] = 0.0; x[2] = 0.0;
        double potential1=gPot(x, pH1, pH2, gettime());
        p[1] = py; 
        p[2] = sqrt((E-potential1)*2.0 -p[1]*p[1]-0.0)-del*j;//;//from maximum pz
        p[0] = sqrt((E-potential1)*2.0 -p[1]*p[1] -p[2]*p[2]); //sqrt((E-potential)*2.0 -p[1]*p[1])
*/        
        rinitial[0]=x[0];
        zinitial[0]=x[2];
    
        //////////////////////the surface of sections//////////////////////
        // intital conditions
        double oldtime=0.0;
        double oldx=x[0]; 
        double oldy=x[1];
        double oldz=x[2];
        double oldvx = p[0];
        double oldvy = p[1];
        double oldvz = p[2];
        double oldR,oldphi,oldvR,oldvphi;
        phitrans(oldR,oldphi,oldvR,oldvphi);
        ostr11<<E<<" "<<potential1<<" "<<p[0]<<" "<<p[1]<<" "<<p[2]<<" "<< "\n";
        
        integrate(maxtime, pH1, pH2);
        for (int i = 1; i < M; i++) {
            for (int ii = 0; ii < 1; ii++) {
                leapfrog(pH1, pH2);
            }
            
            double R,phi,vR,vphi;
            phitrans(R,phi,vR,vphi);
            double vx = p[0];
            double vy = p[1];
            double vz = p[2];
            double Lz = vphi; // Lz is the angular momentum p_phi or v_phi
            t[i] = gettime()-10000.0;
            rinitial[i]=R;
            zinitial[i]=x[2];

            // detect the change of sign
            if (oldz*x[2] < 0) { 
                if (x[2] >= 0){ // Select z with positive sign
                
                    //ostr1 is the SOS, surface of sections
                    ostr1<<gettime()-10000.0<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<R<<" "<<phi<<" "; // positions
                    ostr1<<vx<<" "<<vy<<" "<<vz<<" "<<vR<<" "<<vphi<<" "; // velocities
                    ostr1<<j+1.0<<" "<<E<<" "<<potential<<" "<< "\n";
                }
                    if (p[2] >= 0){
                        //ostr12 is the SOS of positive/rising vz
                        ostr12<<gettime()-10000.0<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<R<<" "<<phi<<" "; // positions
                        ostr12<<vx<<" "<<vy<<" "<<vz<<" "<<vR<<" "<<vphi<<" "; // velocities
                        ostr12<<j+1.0<<" "<<E<<" "<<potential<<" "<< "\n";
                    }    
                        
                else {
                    ostr1<<oldtime<<" "<<oldx<<" "<<oldy<<" "<<oldz<<" "<<oldR<<" "<<oldphi<<" "; // positions
                    ostr1<<oldvx<<" "<<oldvy<<" "<<oldvz<<" "<<oldvR<<" "<<oldvphi<<" "; // velocities
                    ostr1<<j+1.0<<" "<<E<<" "<<potential<<" "<< "\n";
                } 
                    if (p[2] >= 0){
                        //ostr12 is the SOS of positive/rising vz
                        ostr12<<gettime()-10000.0<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<R<<" "<<phi<<" "; // positions
                        ostr12<<vx<<" "<<vy<<" "<<vz<<" "<<vR<<" "<<vphi<<" "; // velocities
                        ostr12<<j+1.0<<" "<<E<<" "<<potential<<" "<< "\n";
                    }
            } 
            oldtime=gettime()-10000.0;oldx=x[0]; oldy=x[1]; oldz=x[2]; oldR=R; oldphi=phi;
            oldvx=vx; oldvy=vy; oldvz=vz; oldvR=vR; oldvphi=vphi;
        }
        
        /////////////////////Fourier analysis and find the MOR/////////////////////
       
        //Interpolation
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
        
        double totaltime = gettime()-10000.0;
        double delta=(totaltime)/M;
        for (int i = 0; i < M; i++) {
            double ti=1.0*delta*i;
            REAL(rfinal,i) = gsl_spline_eval (spliner, ti, accr); IMAG(rfinal,i) = 0.0;
            REAL(zfinal,i) = gsl_spline_eval (splinez, ti, accz); IMAG(zfinal,i) = 0.0;
            REAL(thetafinal,i) = atan2(REAL(zfinal,i),REAL(rfinal,i)); IMAG(thetafinal,i) = 0.0;
        }
        
        //free the storage
        gsl_spline_free (spliner);
        gsl_interp_accel_free (accr);
        gsl_spline_free (splinez);
        gsl_interp_accel_free (accz);
        
        
        //Fourier analysis
        //test1: fourier transform the fist half, then fourier transform the late half.
        //Fourier analysis
        gsl_fft_complex_radix2_forward (rfinal, delta, M);
        gsl_fft_complex_radix2_forward (zfinal, delta, M);
        gsl_fft_complex_radix2_forward (thetafinal, delta, M);
        
        // store fr, fz, magnitudes 
        double frequency[M];
        double sqrtmagr[M];
        double sqrtmagz[M];
        double sqrtmagtheta[M];
    
        //find maximum magnitude
        double magr[M];
        double magz[M];
        double magtheta[M];
        double magrmax=0.0;
        double magzmax=0.0;
        double magthetamax=0.0;
        double nr=0.0;//the location when magrmax occurs
        double nz=0.0;//the location when magzmax occurs
        double ntheta=0.0;//the location when magthetamax occurs

        //Find the maximum magnitudes and corresponding frequencies to understand the MOR in the first half
        for (i = 1; i < M/2; i++){
            magr[i]=sqrt(REAL(rfinal,i)*REAL(rfinal,i)+IMAG(rfinal,i)*IMAG(rfinal,i));
            magz[i]=sqrt(REAL(zfinal,i)*REAL(zfinal,i)+IMAG(zfinal,i)*IMAG(zfinal,i));
            magtheta[i]=sqrt(REAL(thetafinal,i)*REAL(thetafinal,i)+IMAG(thetafinal,i)*IMAG(thetafinal,i));
        
            if (magr[i] > magrmax) {
                magrmax = magr[i];
                nr = i;
            } 
            
            if (fabs(magz[i]) > magzmax) {
                magzmax = fabs(magz[i]);
                nz = i;
            }
            if (fabs(magtheta[i]) > magthetamax) {
                magthetamax = fabs(magtheta[i]);
                ntheta = i;
            }
            if (i==(M/2-1)){
                double fr = nr*1.0/(M);
                double fz = nz*1.0/(M);
                ostr2 << fr << " "<< fz<< " ";
                ostr2<< magrmax/M << " "<< magthetamax/M << " "<< (nr/nz)<< " ";
                //<< magrmax/M << " " << magzmax/M << " "<< magthetamax/M << " ";
                ostr2 << "\n";
            }
        }
            
      
    }
}
void particle::dragging(GalaxyPotential &pH1, GalaxyPotential &pH2, GalaxyPotential &pH3,ostringstream &ostr1,ostringstream &ostr11,ostringstream &ostr2,ostringstream &ostr21,ostringstream &ostr3){    
    
    int M1 = 16384;//16384;//2048;//8192;
    int M2 = 16384;
    int nR = 2000; // number of R
    
    ///////////////////////many orbits (Rg and vc)//////////////////
    double delR = 0.01;
    double delRg = 0.5;
    int nRG = 12.0/delRg + 1.0; // number of Rg
    int a =4.0/delRg;//
    
    double **pot = new double*[nRG];   //the stars mean pointers, the **lvplane is a pointer to pointers which gets initialised as two-hundred pointers
    //double *vc = new double[nRG];
    double **dpotdR = new double*[nRG];
    //double *Rg = new double[nRG];
    double Rg[nRG];
    double vc[nRG];
    
    int i; int j;
    for ( i = 0; i < nRG; i++) { //this is a loop using the index i running from 0 to 199
        pot[i] = new double[nR];
        dpotdR[i] = new double[nR];
    }
    // find circular velocities vc(Rg)
    for (i = 0; i < nRG; i++){
        double di1 = (double)(i);
        Rg[i]=di1*delRg;
        for (j = 0; j < nR; j++){
            settime();
            x[0] = 0.0 + delR* j; x[1] = 0.0; x[2] = 0.0;
            double r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
            pot[i][j] = gPotconst(x, r, pH1, pH2);
        }
    }
    for (i = 0; i < nRG; i++){
        for (j = 0; j < nR-1; j++){
            dpotdR[i][j] = (pot[i][j+1] - pot[i][j])/delR;

        }
    }
    for (i = 0; i < nRG; i++){
        int nRg = Rg[i]/delR;
        vc[i]=std::sqrt(Rg[i]*dpotdR[i][nRg]);
        ostr1 << Rg[i] << " "<< vc[i] << " ";
        ostr1<<"\n";
    }
    
    // this loop is not working, but I can study any orbits and its phase space/action against resonance
    double t1[M1];
    double rinitial1[M1];
    double xinitial1[M1];
    double yinitial1[M1];
    double zinitial1[M1];
    double vxinitial1[M1];
    double vyinitial1[M1];
    double vzinitial1[M1];
    double rfinal1[2*M1];
    double xfinal1[2*M1];
    double yfinal1[2*M1];
    double zfinal1[2*M1];
    double vxfinal1[2*M1];
    double vyfinal1[2*M1];
    double vzfinal1[2*M1];
    double thetafinal1[2*M1];//theta=arctan(z/r)
    
    double t[M2];
    double rinitial[M2];
    double zinitial[M2];    
    double rfinal[2*M2];
    double zfinal[2*M2];
    double thetafinal[2*M2];//theta=arctan(z/r)  
    // intital conditions
    double oldtime=0.0;
    double oldx=x[0]; 
    double oldy=x[1];
    double oldz=x[2];
    double oldvx = p[0];
    double oldvy = p[1];
    double oldvz = p[2];
    double oldR,oldphi,oldvR,oldvphi;
    phitrans(oldR,oldphi,oldvR,oldvphi);
    
    
    //number of particles
    //for (int k = a; k < nRG; k++){
    for (int k = 1; k < 20; k++){
            
            
            settime();
            //grow disc potential from 10%
            //growing disc potential
            double periods[2];
            double maxtime = 10000.0;
            t1[0]= 0.0;
            /*
            x[0] = Rg[k];//8.0;
            x[1] = 0.0; x[2] = 0.0;
            double potential1 = gPot(x, pH1, pH2, gettime());
            p[0] = 0.0;
            p[1] = vc[k]; //0.177861;
            p[2] = 0.001;
            double Rg1 = Rg[k];
            */
            double dk = (double)(k);
            x[0] = 8.0; x[1] = 0.0; x[2] = 0.0; 
            p[0] = 0.1; p[1] =0.2;p[2] = 0.005*dk;
            double Rg1 = 8.0;
            
            rinitial1[0]=x[0];
            xinitial1[0]=x[0];
            yinitial1[0]=x[1];
            zinitial1[0]=x[2];
            vxinitial1[0]=p[0];
            vyinitial1[0]=p[1];
            vzinitial1[0]=p[2];    
            integrate(maxtime, pH1, pH2);
            
            for (int i = 1; i < M1; i++) {
                for (int ii = 0; ii < 1; ii++) {
                    leapfrog(pH1, pH2);
                }
                
                double R,phi,vR,vphi;
                phitrans(R,phi,vR,vphi);
                t1[i] = gettime()-10000.0;
                rinitial1[i]=R;
                xinitial1[i]=x[0];
                yinitial1[i]=x[1];
                zinitial1[i]=x[2];
                vxinitial1[i]=p[0];
                vyinitial1[i]=p[1];
                vzinitial1[i]=p[2];
                
                
            }        

            /////////////////////Interpolation/////////////////////
            
            gsl_interp_accel *accr1
            = gsl_interp_accel_alloc ();
            gsl_spline *spliner1
            = gsl_spline_alloc (gsl_interp_cspline, M1);
                gsl_interp_accel *accx1
            = gsl_interp_accel_alloc ();
            gsl_spline *splinex1
            = gsl_spline_alloc (gsl_interp_cspline, M1);
                gsl_interp_accel *accy1
            = gsl_interp_accel_alloc ();
            gsl_spline *spliney1
            = gsl_spline_alloc (gsl_interp_cspline, M1);
            gsl_interp_accel *accz1
            = gsl_interp_accel_alloc ();
            gsl_spline *splinez1
            = gsl_spline_alloc (gsl_interp_cspline, M1);
                gsl_interp_accel *accvx1
            = gsl_interp_accel_alloc ();
            gsl_spline *splinevx1
            = gsl_spline_alloc (gsl_interp_cspline, M1);
                gsl_interp_accel *accvy1
            = gsl_interp_accel_alloc ();
            gsl_spline *splinevy1
            = gsl_spline_alloc (gsl_interp_cspline, M1);
                gsl_interp_accel *accvz1
            = gsl_interp_accel_alloc ();
            gsl_spline *splinevz1
            = gsl_spline_alloc (gsl_interp_cspline, M1);
            

            gsl_spline_init (spliner1, t1, rinitial1, M1);
            gsl_spline_init (splinex1, t1, xinitial1, M1);
            gsl_spline_init (spliney1, t1, yinitial1, M1);
            gsl_spline_init (splinez1, t1, zinitial1, M1);
            gsl_spline_init (splinevx1, t1, vxinitial1, M1);
            gsl_spline_init (splinevy1, t1, vyinitial1, M1);
            gsl_spline_init (splinevz1, t1, vzinitial1, M1);
            
            double totaltime1 = gettime()-10000.0;
            double delta1=(totaltime1)/M1;
            for (int i = 0; i < M1; i++) {
                double di = (double)(i);
                double ti1 = 1.0*delta1*di;
                REAL(rfinal1,i) = gsl_spline_eval (spliner1, ti1, accr1); IMAG(rfinal1,i) = 0.0;
                REAL(xfinal1,i) = gsl_spline_eval (splinex1, ti1, accx1); IMAG(xfinal1,i) = 0.0;
                REAL(yfinal1,i) = gsl_spline_eval (spliney1, ti1, accy1); IMAG(yfinal1,i) = 0.0;
                REAL(zfinal1,i) = gsl_spline_eval (splinez1, ti1, accz1); IMAG(zfinal1,i) = 0.0;
                REAL(vxfinal1,i) = gsl_spline_eval (splinevx1, ti1, accvx1); IMAG(vxfinal1,i) = 0.0;
                REAL(vyfinal1,i) = gsl_spline_eval (splinevy1, ti1, accvy1); IMAG(vyfinal1,i) = 0.0;
                REAL(vzfinal1,i) = gsl_spline_eval (splinevz1, ti1, accvz1); IMAG(vzfinal1,i) = 0.0;
                REAL(thetafinal1,i) = atan2(REAL(zfinal1,i),REAL(rfinal1,i)); IMAG(thetafinal1,i) = 0.0;
                x[0] = REAL(xfinal1,i); x[1] = REAL(yfinal1,i); x[2] = REAL(zfinal1,i) = REAL(zfinal1,i);
                double potgrow1 = gPot(x, pH1, pH2, ti1);
                ostr11<<ti1<<" "<<potgrow1<<" ";
                ostr11<<REAL(xfinal1,i)<<" "<<REAL(yfinal1,i)<<" "<<REAL(zfinal1,i)<<" ";
                ostr11<<sqrt(REAL(xfinal1,i)*REAL(xfinal1,i) + REAL(yfinal1,i)*REAL(yfinal1,i))<<" "<<REAL(vzfinal1,i)<<" ";
                ostr11<< "\n";
            }

            //free the storage
            gsl_spline_free (spliner1);
            gsl_interp_accel_free (accr1);
            gsl_spline_free (splinex1);
            gsl_interp_accel_free (accx1);
            gsl_spline_free (spliney1);
            gsl_interp_accel_free (accy1);    
            gsl_spline_free (splinez1);
            gsl_interp_accel_free (accz1);
                gsl_spline_free (splinevx1);
            gsl_interp_accel_free (accvx1);
                gsl_spline_free (splinevy1);
            gsl_interp_accel_free (accvy1);
                gsl_spline_free (splinevz1);
            gsl_interp_accel_free (accvz1);

            /////////////////////////////potential at specific time//////////////////////////////////  
    
            for (j = 0; j < 10; j++){ //328
                settime();
                double periods[2];
                double maxtime = 10000.0;
                double dj = (double)(j);
                double gap = 16380.0/(10.0*delta1);
                int k = dj*16380.0/10.0;
                double resolution = k*delta1; //find specific time, so that the POT is constant
                x[0] = REAL(xfinal1,k); x[1] = REAL(yfinal1,k); x[2] = REAL(zfinal1,k);
                p[0] = REAL(vxfinal1,k); p[1] = REAL(vyfinal1,k); p[2] = REAL(vzfinal1,k);
                double R = sqrt(x[0]*x[0] + x[1]*x[1]);
                rinitial[0]=sqrt(x[0]*x[0] + x[1]*x[1]);
                zinitial[0]=x[2];
                double r = sqrt(x[0]*x[0] + x[1]*x[1]);
                double potgrow = gPottime(x,r, pH1, pH2, resolution);
                //ostr3<<resolution<<" "<<potgrow<<" ";
                //ostr3<<k<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" ";
                //ostr3<< "\n";
                ostr3<< Rg1<< " "<< resolution<< " ";
                

                integrateL(maxtime, pH1, pH2,resolution);
                for (int i = 1; i < M2; i++) {
                    for (int ii = 0; ii < 1; ii++) {
                        leapfrogL(pH1, pH2, resolution);
                    }
                    
                    t[i] = gettime()-10000.0;
                    rinitial[i]=sqrt(x[0]*x[0] + x[1]*x[1]);
                    zinitial[i]=x[2];
                    // test: fine
                    //ostr3<<rinitial[i]<<" "<<zinitial[i]<<" ";
                    //ostr3<< "\n";
                }


                /////////////////////Fourier analysis and find the MOR/////////////////////
            
                //Interpolation
                gsl_interp_accel *accr
                = gsl_interp_accel_alloc ();
                gsl_spline *spliner
                = gsl_spline_alloc (gsl_interp_cspline, M2);
                gsl_interp_accel *accz
                = gsl_interp_accel_alloc ();
                gsl_spline *splinez
                = gsl_spline_alloc (gsl_interp_cspline, M2);
                

                gsl_spline_init (spliner, t, rinitial, M2);
                gsl_spline_init (splinez, t, zinitial, M2);
                
                double totaltime = gettime()-10000.0;
                double delta=(totaltime)/M2;
                for (int i = 0; i < M2; i++) {
                    double di2 = (double)(i);
                    double ti=1.0*delta*di2;
                    REAL(rfinal,i) = gsl_spline_eval (spliner, ti, accr); IMAG(rfinal,i) = 0.0;
                    REAL(zfinal,i) = gsl_spline_eval (splinez, ti, accz); IMAG(zfinal,i) = 0.0;
                    REAL(thetafinal,i) = atan2(REAL(zfinal,i),REAL(rfinal,i)); IMAG(thetafinal,i) = 0.0;
                    // test: fine
                    //ostr3<<REAL(rfinal,i)<<" "<<REAL(zfinal,i)<<" ";
                    //ostr3<< "\n";
                }
                
                //free the storage
                gsl_spline_free (spliner);
                gsl_interp_accel_free (accr);
                gsl_spline_free (splinez);
                gsl_interp_accel_free (accz);
                
                
                //Fourier analysis
                //test1: fourier transform the fist half, then fourier transform the late half.
                //Fourier analysis
                gsl_fft_complex_radix2_forward (rfinal, delta, M2);
                gsl_fft_complex_radix2_forward (zfinal, delta, M2);
                gsl_fft_complex_radix2_forward (thetafinal, delta, M2);
                
                // test: fine
                // store fr, fz, magnitudes 
                double frequency[M2];
                double sqrtmagr[M2];
                double sqrtmagz[M2];
                double sqrtmagtheta[M2];

                //find maximum magnitude
                int halfM2 = 8192;
                double magr[halfM2];
                double magz[halfM2];
                double magtheta[halfM2];
                double magrmax=0.0;
                double magzmax=0.0;
                double magthetamax=0.0;
                double nr=0.0;//the location when magrmax occurs
                double nz=0.0;//the location when magzmax occurs
                double ntheta=0.0;//the location when magthetamax occurs
                
                //Find the maximum magnitudes and corresponding frequencies to understand the MOR in the first half
            
                for (int i = 1; i < M2/2; i++){
                    magr[i]=sqrt(REAL(rfinal,i)*REAL(rfinal,i)+IMAG(rfinal,i)*IMAG(rfinal,i));
                    magz[i]=sqrt(REAL(zfinal,i)*REAL(zfinal,i)+IMAG(zfinal,i)*IMAG(zfinal,i));
                    magtheta[i]=sqrt(REAL(thetafinal,i)*REAL(thetafinal,i)+IMAG(thetafinal,i)*IMAG(thetafinal,i));
                    //ostr3 << magr[i]<< " "<< magz[i]<< " ";
                    //ostr3 << "\n";
            
            
                    double di3 = (double)(i);
                    if (magr[i] > magrmax) {
                        magrmax = magr[i];
                        nr = di3;
                    } 
                    
                    if (fabs(magz[i]) > magzmax) {
                        magzmax = fabs(magz[i]);
                        nz = di3;
                    }
                    if (fabs(magtheta[i]) > magthetamax) {
                        magthetamax = fabs(magtheta[i]);
                        ntheta = di3;
                    }
                    if (i==(M2/2-1)){
                        double fr = nr*1.0/(M2);
                        double fz = nz*1.0/(M2);
                        ostr3 << (nr/nz)<< " ";
                        ostr3 << "\n";
                    }
                }

               
                // phase space and SOS
                for (int i = 1; i < M2; i++) {
                    for (int ii = 0; ii < 1; ii++) {
                        leapfrogL(pH1, pH2, resolution);
                    }
                    
                    t[i] = gettime();
                    ////////////////SOS///////////////////////////////////
                    double R,phi,vR,vphi;
                    phitrans(R,phi,vR,vphi);
                    double vx = p[0];
                    double vy = p[1];
                    double vz = p[2];
                    double Lz = vphi; // Lz is the angular momentum p_phi or v_phi 
                    
                    //z:vz
                    ostr21<<resolution<<" "<<x[2]<<" "<<vz<<" "<<(nr/nz)<<" "<<Rg1<<"\n";
                    // detect the change of sign
                    if (oldz*x[2] < 0) { 
                        if (x[2] >= 0){ // Select z with positive sign

                            if (p[2] >= 0){
                                //ostr12 is the SOS of positive/rising vz
                                
                                ostr2<<Rg1<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<-(R-Rg1)<<" "<<phi<<" "; // positions
                                ostr2<<vx<<" "<<vy<<" "<<vz<<" "<<vR<<" "<<vphi<<" "; // velocities
                                ostr2<<gettime()<<" "<<(nr/nz)<< " "<<resolution<< " ";
                                ostr2<<R<< " "<<dk<< "\n";
                                
                            }    
                        }        
                        else {

                            if (p[2] >= 0){
                                //ostr12 is the SOS of positive/rising vz
                                
                                ostr2<<Rg1<<" "<<x[0]<<" "<<x[1]<<" "<<x[2]<<" "<<-(R-Rg1)<<" "<<phi<<" "; // positions
                                ostr2<<vx<<" "<<vy<<" "<<vz<<" "<<vR<<" "<<vphi<<" "; // velocities
                                ostr2<<gettime()<<" "<<(nr/nz)<< " "<<resolution<< " ";
                                ostr2<<R<< " "<<dk<< "\n";
                                
                            }
                            
                        }
                    } 
                    oldtime=gettime()-10000.0;oldx=x[0]; oldy=x[1]; oldz=x[2]; oldR=R; oldphi=phi;
                    oldvx=vx; oldvy=vy; oldvz=vz; oldvR=vR; oldvphi=vphi;
                }
   
            }
    }

    for (i = 0; i < nRG; i++){
       
        delete [] pot[i];      //free all memory at all pointers in lvplane
        delete [] dpotdR[i];
        }

    delete [] pot;      //free all memory at all pointers in lvplane
    //delete [] vc;
    delete [] dpotdR;
    //delete [] Rg;

}

void particle::simplelevitation(GalaxyPotential &pH1, GalaxyPotential &pH2, GalaxyPotential &pH3, ostringstream &ostr){
    //Instructions for the gnuplot: (printorbitNa) 4:5=heated; 6:7=unperturbed; 4:8=levitation
    double steps = 16384.0;
    int j;
    // For the adiabaticly heated particle
    for (int j = 0; j < 200; j++) {
        settime(); 
        x[0] = 8.0; x[1] = 0.0; x[2] = 0.0; 
        p[0] = 0.1; p[1] =0.2;p[2] = 0.0005*j;

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
                ostr << RmaxA << " " << zmaxA << " "; //<< nrA+1 << " "<< nzA+1 << " "
            }
        }
        

        settime();  
        // For the contracted particle
        x[0] = 8.0; x[1] = 0.0; x[2] = 0.0; 
        p[0] = 0.1; p[1] =0.2;p[2] = 0.0005*j;

        integrate2(maxtime, pH1, pH3);
        double RmaxB=0;
        double zmaxB=0;
        double RB = getR();
        double zB = getz();
        for (int i = 0; i < steps; i++) {
            for (int ii = 0; ii < 1; ii++) {
                leapfrog2(pH1, pH3);
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

}

void particle::GMClevitation(GalaxyPotential &pH1, GalaxyPotential &pH2, GalaxyPotential &pH3, ostringstream &ostr1, ostringstream &ostr11){
    //Instructions for the gnuplot: (printorbitNa) 4:5=heated; 6:7=unperturbed; 4:8=levitation
    
    //////////////GMC vc///////////////////////
    int M1 = 16384;//16384;//2048;//8192;
    int M2 = 16384;
    int nR = 2000; // number of R
    
    ///////////////////////many orbits (Rg and vc)//////////////////
    double delR = 0.01;
    double delRg = 0.5;
    int nGMC = 6000.0;
    
    
    //double **pot = new double*[nRG];   //the stars mean pointers, the **lvplane is a pointer to pointers which gets initialised as two-hundred pointers
    //double *vc = new double[nRG];
    //double **dpotdR = new double*[nRG];
    double N[nGMC];// mass function
    //double vc[nGMC];
    
    double mbottom = 10000.0;
    double mupper = 3000000.0;
    int numb = mupper/mbottom;
    double Ntotal = 0.0;
    double total= 1074.06;
    double mGMC[nGMC];
    int i; 
    for ( i = 1; i < numb; i++) { 
        double mass = mbottom * i;
        N[i] = pow((mass/mupper),-0.8);
        Ntotal += N[i];
        ostr1 << i << " "<< log10(mass) << " "<< log10(N[i]) << " "<< mass << " "<< Ntotal << " ";
        mGMC[i] = int(N[i]*nGMC/total); //number of GMCs in each mass region  
        ostr1 << mGMC[i]<< " ";
        ostr1<<"\n";
    }
/*    
    for ( i = 0; i < nRG; i++) { 
        pot[i] = new double[nR];
        dpotdR[i] = new double[nR];
    }
    // find circular velocities vc(Rg)
    for (i = 0; i < nRG; i++){
        double di1 = (double)(i);
        Rg[i]=di1*delRg;
        for (j = 0; j < nR; j++){
            settime();
            x[0] = 0.0 + delR* j; x[1] = 0.0; x[2] = 0.0;
            double r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
            pot[i][j] = gPotconst(x, r, pH1, pH2);
        }
    }
    for (i = 0; i < nRG; i++){
        for (j = 0; j < nR-1; j++){
            dpotdR[i][j] = (pot[i][j+1] - pot[i][j])/delR;

        }
    }
    for (i = 0; i < nRG; i++){
        int nRg = Rg[i]/delR;
        vc[i]=std::sqrt(Rg[i]*dpotdR[i][nRg]);
        ostr1 << Rg[i] << " "<< vc[i] << " ";
        ostr1<<"\n";
    }
*/
/*
    double t1[M1];
    double rinitial1[M1];
    double xinitial1[M1];
    double yinitial1[M1];
    double zinitial1[M1];
    double vxinitial1[M1];
    double vyinitial1[M1];
    double vzinitial1[M1];
    double rfinal1[2*M1];
    double xfinal1[2*M1];
    double yfinal1[2*M1];
    double zfinal1[2*M1];
    double vxfinal1[2*M1];
    double vyfinal1[2*M1];
    double vzfinal1[2*M1];
    double thetafinal1[2*M1];//theta=arctan(z/r)
    
    double t[M2];
    double rinitial[M2];
    double zinitial[M2];    
    double rfinal[2*M2];
    double zfinal[2*M2];
    double thetafinal[2*M2];//theta=arctan(z/r)  
    // intital conditions
    double oldtime=0.0;
    double oldx=x[0]; 
    double oldy=x[1];
    double oldz=x[2];
    double oldvx = p[0];
    double oldvy = p[1];
    double oldvz = p[2];
    double oldR,oldphi,oldvR,oldvphi;
    phitrans(oldR,oldphi,oldvR,oldvphi);
    
    for (int k = 0; k < nGMC; k++){
            
            
            settime();
            //grow disc potential from 10%
            //growing disc potential
            double periods[2];
            double maxtime = 10000.0;
            t1[0]= 0.0;
            
            double dk = (double)(k);
            x[0] = 8.0; x[1] = 0.0; x[2] = 0.0; 
            p[0] = 0.0; p[1] =0.177861*mass[k];
            p[2] = 0.001;
            double Rg1 = 8.0;
            rinitial1[0]=x[0];
            xinitial1[0]=x[0];
            yinitial1[0]=x[1];
            zinitial1[0]=x[2];
            vxinitial1[0]=p[0];
            vyinitial1[0]=p[1];
            vzinitial1[0]=p[2];    
            integrateGMC(maxtime, pH1, pH2, mass[k]);
            
            for (int i = 1; i < M1; i++) {
                for (int ii = 0; ii < 1; ii++) {
                    leapfrogGMC(pH1, pH2, mass[k]);
                }
                
                double R,phi,vR,vphi;
                phitrans(R,phi,vR,vphi);
                t1[i] = gettime()-10000.0;
                rinitial1[i]=R;
                xinitial1[i]=x[0];
                yinitial1[i]=x[1];
                zinitial1[i]=x[2];
                vxinitial1[i]=p[0];
                vyinitial1[i]=p[1];
                vzinitial1[i]=p[2];
                
                
            }        

            /////////////////////Interpolation/////////////////////
            
            gsl_interp_accel *accr1
            = gsl_interp_accel_alloc ();
            gsl_spline *spliner1
            = gsl_spline_alloc (gsl_interp_cspline, M1);
                gsl_interp_accel *accx1
            = gsl_interp_accel_alloc ();
            gsl_spline *splinex1
            = gsl_spline_alloc (gsl_interp_cspline, M1);
                gsl_interp_accel *accy1
            = gsl_interp_accel_alloc ();
            gsl_spline *spliney1
            = gsl_spline_alloc (gsl_interp_cspline, M1);
            gsl_interp_accel *accz1
            = gsl_interp_accel_alloc ();
            gsl_spline *splinez1
            = gsl_spline_alloc (gsl_interp_cspline, M1);
                gsl_interp_accel *accvx1
            = gsl_interp_accel_alloc ();
            gsl_spline *splinevx1
            = gsl_spline_alloc (gsl_interp_cspline, M1);
                gsl_interp_accel *accvy1
            = gsl_interp_accel_alloc ();
            gsl_spline *splinevy1
            = gsl_spline_alloc (gsl_interp_cspline, M1);
                gsl_interp_accel *accvz1
            = gsl_interp_accel_alloc ();
            gsl_spline *splinevz1
            = gsl_spline_alloc (gsl_interp_cspline, M1);
            

            gsl_spline_init (spliner1, t1, rinitial1, M1);
            gsl_spline_init (splinex1, t1, xinitial1, M1);
            gsl_spline_init (spliney1, t1, yinitial1, M1);
            gsl_spline_init (splinez1, t1, zinitial1, M1);
            gsl_spline_init (splinevx1, t1, vxinitial1, M1);
            gsl_spline_init (splinevy1, t1, vyinitial1, M1);
            gsl_spline_init (splinevz1, t1, vzinitial1, M1);
            
            double totaltime1 = gettime()-10000.0;
            double delta1=(totaltime1)/M1;
            for (int i = 0; i < M1; i++) {
                double di = (double)(i);
                double ti1 = 1.0*delta1*di;
                REAL(rfinal1,i) = gsl_spline_eval (spliner1, ti1, accr1); IMAG(rfinal1,i) = 0.0;
                REAL(xfinal1,i) = gsl_spline_eval (splinex1, ti1, accx1); IMAG(xfinal1,i) = 0.0;
                REAL(yfinal1,i) = gsl_spline_eval (spliney1, ti1, accy1); IMAG(yfinal1,i) = 0.0;
                REAL(zfinal1,i) = gsl_spline_eval (splinez1, ti1, accz1); IMAG(zfinal1,i) = 0.0;
                REAL(vxfinal1,i) = gsl_spline_eval (splinevx1, ti1, accvx1); IMAG(vxfinal1,i) = 0.0;
                REAL(vyfinal1,i) = gsl_spline_eval (splinevy1, ti1, accvy1); IMAG(vyfinal1,i) = 0.0;
                REAL(vzfinal1,i) = gsl_spline_eval (splinevz1, ti1, accvz1); IMAG(vzfinal1,i) = 0.0;
                REAL(thetafinal1,i) = atan2(REAL(zfinal1,i),REAL(rfinal1,i)); IMAG(thetafinal1,i) = 0.0;
                x[0] = REAL(xfinal1,i); x[1] = REAL(yfinal1,i); x[2] = REAL(zfinal1,i) = REAL(zfinal1,i);
                double potgrow1 = gPot(x, pH1, pH2, ti1);
                ostr11<<ti1<<" "<<potgrow1<<" ";
                ostr11<<REAL(xfinal1,i)<<" "<<REAL(yfinal1,i)<<" "<<REAL(zfinal1,i)<<" ";
                ostr11<<sqrt(REAL(xfinal1,i)*REAL(xfinal1,i) + REAL(yfinal1,i)*REAL(yfinal1,i))<<" "<<REAL(vzfinal1,i)<<" ";
                ostr11<< "\n";
            }

            //free the storage
            gsl_spline_free (spliner1);
            gsl_interp_accel_free (accr1);
            gsl_spline_free (splinex1);
            gsl_interp_accel_free (accx1);
            gsl_spline_free (spliney1);
            gsl_interp_accel_free (accy1);    
            gsl_spline_free (splinez1);
            gsl_interp_accel_free (accz1);
                gsl_spline_free (splinevx1);
            gsl_interp_accel_free (accvx1);
                gsl_spline_free (splinevy1);
            gsl_interp_accel_free (accvy1);
                gsl_spline_free (splinevz1);
            gsl_interp_accel_free (accvz1);
    }
    
*/

/*
//////////////////////////////////////////////////////////
    double steps = 16384.0;
    int j;
    // For the adiabaticly heated particle
    for (int j = 0; j < 200; j++) {
        settime(); 
        x[0] = 8.0; x[1] = 0.0; x[2] = 0.0; 
        p[0] = 0.1; p[1] =0.2;p[2] = 0.0005*j;

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
                ostr << RmaxA << " " << zmaxA << " "; //<< nrA+1 << " "<< nzA+1 << " "
            }
        }
        

        settime();  
        // For the contracted particle
        x[0] = 8.0; x[1] = 0.0; x[2] = 0.0; 
        p[0] = 0.1; p[1] =0.2;p[2] = 0.0005*j;

        integrate2(maxtime, pH1, pH3);
        double RmaxB=0;
        double zmaxB=0;
        double RB = getR();
        double zB = getz();
        for (int i = 0; i < steps; i++) {
            for (int ii = 0; ii < 1; ii++) {
                leapfrog2(pH1, pH3);
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
*/
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

void particle::integrateGMC(double integtime, GalaxyPotential &pH, GalaxyPotential &pD, double mass) {
    stepnumber = 0;
    double ttarget = gettime() + integtime;
    while (gettime() < ttarget) {
        leapfrogGMC(pH, pD, mass);
    }
}
void particle::integrate2(double integtime, GalaxyPotential &pH, GalaxyPotential &pD2) {
    stepnumber = 0;
    double ttarget = gettime() + integtime;
    while (gettime() < ttarget) {
        leapfrog2(pH, pD2);
    }
}    

void particle::integrateL(double integtime, GalaxyPotential &pH, GalaxyPotential &pD, double time) {
    stepnumber = 0;
    double ttarget = gettime() + integtime;
    while (gettime() < ttarget) {
        leapfrogL(pH, pD, time);
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




void Dotest1() {
    
    std::ifstream fromfile("Halo1.Tpot");
    GalaxyPotential pH(fromfile);
    std::ifstream fromfile1("Disc1.Tpot");
    GalaxyPotential pD(fromfile1);
    std::ifstream fromfile2("Disc2.Tpot");//MW2017
    GalaxyPotential pD2(fromfile2);
    
    particle pp;
    ostringstream ostr1;
    ostringstream ostr11;
    ostringstream ostr12;
    ostringstream ostr2;
    ostringstream ostr21;
    ostringstream ostr3;
    //fix the potential
    //fourier analysis
    pp.fourier(pH,pD,ostr1,ostr11,ostr12,ostr2,ostr21,ostr3);
    //pD is 10% of pD2
    FILE *fout1 = fopen("1printorbit", "w");
    string stri1 = ostr1.str();
    char *ch1 = &stri1[0];
    fputs(ch1, fout1);
    fclose(fout1);
    FILE *fout11 = fopen("1SOS", "w");
    string stri11 = ostr11.str();
    char *ch11 = &stri11[0];
    fputs(ch11, fout11);
    fclose(fout11);
    FILE *fout12 = fopen("1SOS+vz", "w");
    string stri12 = ostr12.str();
    char *ch12 = &stri12[0];
    fputs(ch12, fout12);
    fclose(fout12);
    FILE *fout2 = fopen("1interpolation", "w");
    string stri2 = ostr2.str();
    char *ch2 = &stri2[0];
    fputs(ch2, fout2);
    fclose(fout2);
    FILE *fout21 = fopen("1max", "w");
    string stri21 = ostr21.str();
    char *ch21 = &stri21[0];
    fputs(ch21, fout21);
    fclose(fout21);
    FILE *fout3 = fopen("1printfourier", "w");
    string stri3 = ostr3.str();
    char *ch3 = &stri3[0];
    fputs(ch3, fout3);
    fclose(fout3);
   
} 

void Dotest2() {
    std::ifstream fromfile("Halo1.Tpot");
    GalaxyPotential pH(fromfile);
    std::ifstream fromfile1("Disc1.Tpot");
    GalaxyPotential pD(fromfile1);
    std::ifstream fromfile2("Disc2.Tpot");//MW2017
    GalaxyPotential pD2(fromfile2);
    particle pp;
    ostringstream ostr1;
    ostringstream ostr11;
    ostringstream ostr12;
    ostringstream ostr2;
    //sos of 10 orbits
    pp.sos(pH,pD,ostr1,ostr11, ostr12,ostr2);
    FILE *fout1 = fopen("2SOSN", "w");
    string stri1 = ostr1.str();
    char *ch1 = &stri1[0];
    fputs(ch1, fout1);
    fclose(fout1);
    FILE *fout11 = fopen("2test", "w");
    string stri11 = ostr11.str();
    char *ch11 = &stri11[0];
    fputs(ch11, fout11);
    fclose(fout11);
    FILE *fout12 = fopen("2SOSN+vz", "w");
    string stri12 = ostr12.str();
    char *ch12 = &stri12[0];
    fputs(ch12, fout12);
    fclose(fout12);
    FILE *fout2 = fopen("2FFTN", "w");
    string stri2 = ostr2.str();
    char *ch2 = &stri2[0];
    fputs(ch2, fout2);
    fclose(fout2);
}
void Dotest3() {
    std::ifstream fromfile("Halo1.Tpot");
    GalaxyPotential pH(fromfile);
    std::ifstream fromfile1("Disc1.Tpot");
    GalaxyPotential pD(fromfile1);
    std::ifstream fromfile2("Disc2.Tpot");
    GalaxyPotential pD2(fromfile2);
    ostringstream ostrNa;
    particle pp;
    // making particles
    pp.simplelevitation(pH,pD,pD2,ostrNa);
    FILE *foutNa = fopen("3printorbitNa", "w");
    string striNa = ostrNa.str();
    char *chNa = &striNa[0];
    fputs(chNa, foutNa);
    fclose(foutNa);
}

void Dotest4() {
    //epicycle frequency and effective potential
    std::ifstream fromfile("Halo1.Tpot");
    GalaxyPotential pH(fromfile);
    std::ifstream fromfile1("Disc1.Tpot");
    GalaxyPotential pD(fromfile1);
    std::ifstream fromfile2("Disc2.Tpot");//MCM2017
    GalaxyPotential pD2(fromfile2);
    particle pp;
    ostringstream ostr1;
    ostringstream ostr21;
    ostringstream ostr22;
    ostringstream ostr31;
    ostringstream ostr32;
    pp.epicycle(pH,pD,ostr1,ostr21,ostr22,ostr31,ostr32);
    FILE *fout1 = fopen("4vc", "w");
    string stri1 = ostr1.str();
    char *ch1 = &stri1[0];
    fputs(ch1, fout1);
    fclose(fout1);
    FILE *fout21 = fopen("4effpotR", "w");
    string stri21 = ostr21.str();
    char *ch21 = &stri21[0];
    fputs(ch21, fout21);
    fclose(fout21);
    FILE *fout22 = fopen("4epicycleR", "w");
    string stri22 = ostr22.str();
    char *ch22 = &stri22[0];
    fputs(ch22, fout22);
    fclose(fout22);
    FILE *fout31 = fopen("4effpotz", "w");
    string stri31 = ostr31.str();
    char *ch31 = &stri31[0];
    fputs(ch31, fout31);
    fclose(fout31);
    FILE *fout32 = fopen("4epicyclez", "w");
    string stri32 = ostr32.str();
    char *ch32 = &stri32[0];
    fputs(ch32, fout32);
    fclose(fout32);
}

void Dotest5() {
    std::ifstream fromfile("Halo1.Tpot");
    GalaxyPotential pH(fromfile);
    std::ifstream fromfile1("Disc1.Tpot");
    GalaxyPotential pD(fromfile1);
    std::ifstream fromfile2("Disc2.Tpot");//MCM2017
    GalaxyPotential pD2(fromfile2);
    particle pp;
    ostringstream ostr1;
    ostringstream ostr11;
    ostringstream ostr2;
    ostringstream ostr21;
     ostringstream ostr3;
    //sos of 10 orbits
    pp.dragging(pH,pD,pD,ostr1,ostr11,ostr2,ostr21,ostr3);
    FILE *fout1 = fopen("5dragvc", "w");
    string stri1 = ostr1.str();
    char *ch1 = &stri1[0];
    fputs(ch1, fout1);
    fclose(fout1);
    FILE *fout11 = fopen("5dragprint", "w");
    string stri11 = ostr11.str();
    char *ch11 = &stri11[0];
    fputs(ch11, fout11);
    fclose(fout11);
    FILE *fout2 = fopen("5dragSOS", "w");
    string stri2 = ostr2.str();
    char *ch2 = &stri2[0];
    fputs(ch2, fout2);
    fclose(fout2);
    FILE *fout21 = fopen("5dragvz", "w");
    string stri21 = ostr21.str();
    char *ch21 = &stri21[0];
    fputs(ch21, fout21);
    fclose(fout21);
    FILE *fout3 = fopen("5dragfourier", "w");
    string stri3 = ostr3.str();
    char *ch3 = &stri3[0];
    fputs(ch3, fout3);
    fclose(fout3);
}

void Dotest6() {
    //epicycle frequency and effective potential with respect to time
    std::ifstream fromfile("Halo1.Tpot");
    GalaxyPotential pH(fromfile);
    std::ifstream fromfile1("Disc1.Tpot");
    GalaxyPotential pD(fromfile1);
    std::ifstream fromfile2("Disc2.Tpot");//MCM2017
    GalaxyPotential pD2(fromfile2);
    particle pp;
    ostringstream ostr1;
    ostringstream ostr21;
    ostringstream ostr22;
    ostringstream ostr31;
    ostringstream ostr32;
    pp.epicycletime(pH,pD,ostr1,ostr21,ostr22,ostr31,ostr32);
    FILE *fout1 = fopen("6Tvc", "w");
    string stri1 = ostr1.str();
    char *ch1 = &stri1[0];
    fputs(ch1, fout1);
    fclose(fout1);
    FILE *fout21 = fopen("6TeffpotR", "w");
    string stri21 = ostr21.str();
    char *ch21 = &stri21[0];
    fputs(ch21, fout21);
    fclose(fout21);
    FILE *fout22 = fopen("6TepicycleR", "w");
    string stri22 = ostr22.str();
    char *ch22 = &stri22[0];
    fputs(ch22, fout22);
    fclose(fout22);
    
    // useful
    FILE *fout31 = fopen("6Tratio", "w");
    string stri31 = ostr31.str();
    char *ch31 = &stri31[0];
    fputs(ch31, fout31);
    fclose(fout31);
    //plot "Tratio8" u ($2/1000.0):($3/$4) title "" with points pt 7 ps 0.5
    FILE *fout32 = fopen("6Tratio8", "w");
    string stri32 = ostr32.str();
    char *ch32 = &stri32[0];
    fputs(ch32, fout32);
    fclose(fout32);
}

void Dotest7() {
    std::ifstream fromfile("Halo1.Tpot");
    GalaxyPotential pH(fromfile);
    std::ifstream fromfile1("Disc1.Tpot");
    GalaxyPotential pD(fromfile1);
    std::ifstream fromfile2("Disc2.Tpot");
    GalaxyPotential pD2(fromfile2);
    ostringstream ostrNa;
    ostringstream ostr11;
    particle pp;
    clouds cc;
    // making particles
    pp.GMClevitation(pH,pD,pD2,ostrNa,ostr11);
    FILE *foutNa = fopen("7GMCmass", "w");
    string striNa = ostrNa.str();
    char *chNa = &striNa[0];
    fputs(chNa, foutNa);
    fclose(foutNa);
    FILE *fout11 = fopen("7GMCorbits", "w");
    string stri11 = ostr11.str();
    char *ch11 = &stri11[0];
    fputs(ch11, fout11);
    fclose(fout11);
}

void Dotest10() {
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
    Dotest7();

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
