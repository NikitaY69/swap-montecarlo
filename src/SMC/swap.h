#ifndef SWAP_H
#define SWAP_H

// Libraries
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <sstream> 
#include <iomanip>
#include <vector> 
#include <experimental/filesystem>
#include <boost/program_options.hpp>

namespace fs = std::experimental::filesystem;
namespace po = boost::program_options;

// Global variables
// User-defined parameters
extern int N; //Number of particles
extern int tau; //Correlation max-duration
extern int cycles; //Number of correlation cycles
extern int steps; //Monte Carlo sweeps
extern double T; //Temperature in units of 1/k_B
extern int tw; //Waiting time to start correlation calculations
extern double p_swap; //Swap-attempt probability
extern const int nr; // Number of radius calculations for the correlation lengths
extern const int ns; // Number of sigma calculations for the energy scan

// Model parameters
extern double Size; //44.721359550000003; //Size of the system
const double sigmaMax = 1.613048; //Maximum diameter of particles
const double rSkin = 1.5; //Radius of neighbours included in NL (e.g. 1.8)
const double rC = 1.25 * sigmaMax; //Cutoff radius for calculating potential
const double rNL = pow(rC+rSkin,2); //NL radius squared
const double deltaMax = 0.12; //Max particle displacement
const double deltaSMax = 0.2; //Max diameter difference for swap
const double RUpdate = pow(rSkin,2)/4; //When R2Max exceeds this, update NL
const double x_max = 1.3; // maximal value of r/s for the real neighbours

const double c0 = -28/pow(1.25,12);
const double c2 = 48/pow(1.25,14);
const double c4 = -21/pow(1.25,16);
const double pi = 3.14159265358979323846;

// Arrays
extern double *X, *Y, *S, *Sref, *X0, *Y0;
extern double *Xfull, *Yfull, *Xref, *Yref;
extern std::vector < std::vector <double>> Xtw, Ytw;
// X0 initial position at last neighbour list update
// Xfull real positions (not taking into account periodic boundaries)
// Xref positition at t=0
// Xtw position at last aging update
extern double dXCM, dYCM;
extern std::vector < std::string > allObs;

//  Neighbour Lists
extern std::vector < std::vector<int> > NL, NN;
extern std::vector < std::vector < std::vector <int>>> NN_tw, RL;
// nn_0 nearest neighbours at t=0
// nn_tw nearest neighbours at last aging update

//  Function prototypes
double bcs(double a, double b), Pshift(double a);
void UpdateAge(int cycle), UpdateNL(), UpdateNN(), UpdateRL();
double PairPotential(double x1, double y1, double s1, double x2, double y2, double s2),
       V(double xj, double yj, double rj, int j);
double VTotal(), CBLoc(int cycle, int j), CB(int cycle), MSD(), FS(int cycle),
       DispCorrLoc(int j), DispCorr(), C_sigma(), whichObs(std::string obs, int cycl);
std::vector <double> MicroDispCorrLoc(int j), MicroDispCorr(), SigmaScan(int j);
void TryDisp(int j), TrySwap(int j, int k), MC(std::string out, int n_log, int n_lin);

//  Random number between 0 and 1
#define ranf() \
    ((double)rand()/(1.0+RAND_MAX)) //check random numbers
    
#endif