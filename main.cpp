//  Project
//  Created by Kimlam Nguyen on 07/02/2022.
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <sstream> 
#include <iomanip>
#include <list>
#include <vector> 
// using namespace std;

//  Parameters
const int steps = 10000; //Monte Carlo steps
int stepCounter = 0; //Step counter
const int N = 2000; //Number of particles
const double T = 1; //Temperature in units of 1/k_B
const double Size = 44.721359550000003; //Size of the system

// const int nn = 55; //Maximum number of nearest neighbours
const double sigmaMax = 1.613048; //Maximum diameter of particles
const double rSkin = 1.25 * sigmaMax; //Radius of neighbours included in NL (e.g. 1.8)
const double rC = 1.25 * sigmaMax; //Cutoff radius for calculating potential
const double rNL = pow(rC+rSkin,2); //NL radius squared
const double deltaMax = 0.12; //Max particle displacement
const double RUpdate = pow(rSkin,2)/4; //When R2Max exceeds this, update NL
const double x_max = 1.3; // maximal value of r/s for the real neighbours

//  For log plots
const int dataPoints = 100;
double exponents = log10(steps)/dataPoints;
double samplePoints[dataPoints];

//  Constants
const double c0 = -28/pow(1.25,12);
const double c2 = 48/pow(1.25,14);
const double c4 = -21/pow(1.25,16);
const double pi = 3.14159265358979323846;

//  Arrays
double X[N],Y[N],S[N],X0[N],Y0[N],S0[N];
double Xfull[N],Yfull[N],Xref[N],Yref[N],Xtw[N],Ytw[N];
// Xref=X0 ?
// seems like Xfull is non-periodic while X is
// 

//  Neighbour List
double NL[N][N] = {0};
int numNeighbours[N];
std::vector< std::vector<int> > all_nn0;

//  Write to text file in same folder
std::ofstream log_obs, log_cfg;
std::string motherdir = "/home/allaglo/benchmarks/";

//  Function prototypes
double bcs(double a, double b);
int Find(double arr[], int len, double seek);
void UpdateList();
double PairPotential(double x1, double y1, double s1, double x2, double y2, double s2);
double V(double xj, double yj, double rj, int j);
std::vector<int> nearest_neighbours(int j, double x);
double VTotal(), CBLoc(int j), CB(), MSD(), FS(int tw, int tau, double theta);
void TryDisp(int j), TrySwap(int j, int k), MC(std::string out);

//  Random number between 0 and 1
#define ranf() \
    ((double)rand()/(1.0+RAND_MAX)) //check random numbers

//---------------------------------------------------------
//  main.cpp
int main(int argc, const char * argv[]) {
    
    srand(time(NULL)*1.0); //Random number generator
    std::string outdir = motherdir + argv[1];
    // Get sample points for log scale
    int index = 0;
    for (int x = 0; x <= dataPoints; x++){
        double value = floor(pow(10,exponents*x));
        if(Find(samplePoints, dataPoints, value) == -1){
            samplePoints[index] = value;
            index++;}
        // this if condition is actually relevent because of the floor function
    }
    
    // Read data file
    std::string line;
    std::ifstream myfile (motherdir + "config_init.txt");
    if (myfile.is_open()){
        int i = 0; // particle index
        std::vector<std::vector<double>> cfg; // array of configurations
        while (std::getline(myfile, line)){
            double value;
            std::stringstream ss(line);

            cfg.push_back(std::vector<double>());
            while (ss >> value){
                cfg[i].push_back(value);
            }
            S[i] = cfg[i][0]; X[i] = cfg[i][1]; Y[i] = cfg[i][2];
            X0[i] = X[i]; Y0[i] = Y[i]; Xfull[i] = X[i]; Yfull[i] = Y[i];
            Xref[i] = X[i]; Yref[i] = Y[i]; Xtw[i] = X[i]; Ytw[i] = Y[i];
            i++;}
        myfile.close();

    } else {
        std::cout << motherdir + "config_init.txt" << std::endl;
        return 0;
    }

    // Building list of first neighbours
    for (int i=0; i<N; i++){
        all_nn0.push_back(nearest_neighbours(i, x_max));
    }

    // for (int i=0; i<N; i++){
    //     std::cout << "##########################################" << std::endl;
    //     std::cout << "Particle " << i << std::endl;
    //     for (int n: nearest_neighbours(i, 1.3)){
    //         std::cout << n << " " << std::endl;
    //     }
    //     std::cout << std::endl;
    // }

    // Do simulation with timer
    double t0 = time(NULL); // Timer
    MC(outdir); std::cout << "Time taken: " << (time(NULL) - t0) << "s" << std::endl; 
    // Do MC simulation
    log_obs.close();
    //cout.precision(17);
    std::cout << "Done" << std::endl;
    return 0;
}
//---------------------------------------------------------

//  Calculates difference of a and b while applying periodic boundary conditions
double bcs(double a, double b) {return Size/2 - std::abs(std::abs(a-b)-Size/2);}

//  Finds index of element in array
int Find(double arr[], int len, double seek){
    for (int i = 0; i < len; ++i)
        if (arr[i] == seek) return i;
    return -1;
}

//  Creates the neighbour list for the set of particles, with another list which contains the number of nearest neighbours of each of the N particles.
//  Returns a (N x nn) matrix containing the labels of the neighbours.
void UpdateList(){
    double rij2Row[N], sortedRow[N];
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            double xij = bcs(X[j], X[i]); double yij = bcs(Y[j], Y[i]);
            rij2Row[j] = (xij*xij)+(yij*yij);
            sortedRow[j] = rij2Row[j];
        }
        std::sort(sortedRow, sortedRow+N); // neighbors sorted with respect to the distance
        
        // // finding number of effective neighbors
        // auto NL_check = [&](int k){
        //     return sortedRow[k] > rNL;
        // };
        // numNeighbours[i] = find_if(0, N, sortedRow[j]);
        
        for (int j = 0; sortedRow[j] < rNL; j++){
            numNeighbours[i] = j; // numNeighbors[i] last saved value is the last j
            NL[i][j] = Find(rij2Row, N, sortedRow[j+1]);} 
            // index of the jth neighbor of the ith particle 
        //cout << numNeighbours[i] << endl;
    }
}


//  Calculates the potential of a pair of particles
double PairPotential(double x1, double y1, double s1, double x2, double y2, double s2){
    double sigmaij = (s1+s2)*(1-0.2*abs(s1-s2))/2;
    double sigma2 = sigmaij*sigmaij;
    double rc2 = 1.25 * 1.25 * sigma2;
    double xij = bcs(x1, x2); double yij = bcs(y1, y2);
    double rij2 = (xij*xij)+(yij*yij);

    if (rij2 > rc2) return 0;
    else {
        double a2 = rij2/sigma2; double a4 = a2*a2;
        return (1/(a4*a4*a4))+c0+(c2*a2)+(c4*a4);
    }
}

//  Calculates potential of particle j
double V(double xj, double yj, double rj, int j){
    double total = 0;
    for (int i=0; i < numNeighbours[j]; i++){
        int k = NL[j][i]; // index of the ith neighbor of the jth particle
        total += PairPotential(xj, yj, rj, X[k], Y[k], S[k]);
    }
    return total;
}

//  Calculates total system energy
double VTotal(){
    double vTot = 0;
    for (int j = 0; j < N; j++)
        vTot += V(X[j], Y[j], S[j], j);
    return vTot;
}
// don't we need to divide by 2 ?

//  Calculates avg. mean square displacements
double MSD(){
    double sum = 0, deltaX, deltaY;
        for (int i = 0; i < N; i++){
            deltaX = Xfull[i]-Xref[i];
            deltaY = Yfull[i]-Yref[i];
            sum += deltaX*deltaX + deltaY*deltaY;
    }
    return sum/N;
}

//  Calculates the self scattering function
double FS(int tw, int tau, double theta){
    double dotProduct;
    double q = 2*pi/sigmaMax;
    double sum = 0, deltaX, deltaY;
    if ((stepCounter >= tw) and (stepCounter < tw + tau)){
        for (int i = 0; i < N; i++){
            if(stepCounter == tw){
                Xtw[i] = Xfull[i];
                Ytw[i] = Yfull[i];
            }
            deltaX = Xfull[i]-Xtw[i];
            deltaY = Yfull[i]-Ytw[i];
            dotProduct = q*((cos(theta*pi/180)*deltaX)+(sin(theta*pi/180)*deltaY));
            sum += cos(dotProduct);
        }
    }
    return sum/N;
}

// Computes the nearest neighbours of particle j at a given radius
std::vector<int> nearest_neighbours(int j, double x){
    std::list<int> nn {};
    for (int i=0; i<N; i++){
        double sigmaij = (S[i]+S[j])*(1-0.2*abs(S[i]-S[j]))/2;
        double xij = bcs(X[i], X[j]); double yij = bcs(Y[i], Y[j]);
        double rij = sqrt((xij*xij)+(yij*yij));
        if (rij < x*sigmaij && i != j){
            nn.push_back(i);
        }
    } return std::vector<int> { std::begin(nn), std::end(nn) } ;
}

// Computes the bond-breaking correlation function (local)
double CBLoc(int j){
    std::vector<int> intersect;
    std::vector<int> nn0 = all_nn0[j]; // neighbors at t=0
    std::vector<int> nn = nearest_neighbours(j, x_max);
    std::set_intersection(nn0.begin(), nn0.end(), nn.begin(), nn.end(),
                     std::back_inserter(intersect));

    double frac = intersect.size()/nn0.size();
    return frac;
}

// Computes the bond-breaking correlation function (averaged)
double CB(){
    double tot = 0;
    for (int j=0; j<N; j++){
        tot += CBLoc(j);
    } return tot/N;
}

//  Tries displacing one particle j by vector dr = (dx, dy)
void TryDisp(int j){
    double dx = (ranf()-0.5)*deltaMax;
    double dy = (ranf()-0.5)*deltaMax;
    double deltaE = V(X[j] + dx, Y[j] + dy, S[j], j) - V(X[j], Y[j], S[j], j);
    // why is the modulus function not in deltaE ?
    if (deltaE < 0){
        X[j] = fmod((X[j]+dx),Size); //Check modulus function
        Y[j] = fmod((Y[j]+dy),Size);
        Xfull[j] = Xfull[j]+dx;
        Yfull[j] = Yfull[j]+dy;
    }
    else if (exp(-deltaE/T) > ranf()){
        X[j] = fmod((X[j]+dx),Size);
        Y[j] = fmod((Y[j]+dy),Size);
        Xfull[j] = Xfull[j]+dx;
        Yfull[j] = Yfull[j]+dy;
    }
}

//  Tries swapping the pair of particles j, k
void TrySwap(int j, int k){
    double deltaE = V(X[j],Y[j],S[k],j)+V(X[k],Y[k],S[j],k)-V(X[j],Y[j],S[j],j)-V(X[k],Y[k],S[k],k);
    if (deltaE < 0){
        double Rnew = S[k];
        S[k] = S[j];
        S[j] = Rnew;
    }
    else if (exp(-deltaE/T) > ranf()){
        double Rnew = S[k];
        S[k] = S[j];
        S[j] = Rnew;
    }
}

// Monte Carlo Simulation
void MC(std::string out){
    int dataCounter = 0; stepCounter = 0;
    double deltaX[N], deltaY[N], deltaR2[N], R2Max = 0;
    log_obs.open(out + "obs.txt");
    log_obs << std::scientific << std::setprecision(8);
    UpdateList();
    
    for(int x = 0; x < steps; x++){
        
        // Updating NL
        if(x % 150 == 0) {//Change number?
            // every 150 steps we check if we need to update the NL
            for (int i = 0; i < N; i++){
                deltaX[i] = bcs(X[i],X0[i]);
                deltaY[i] = bcs(Y[i],Y0[i]);
                deltaR2[i] = deltaX[i]*deltaX[i] + deltaY[i]*deltaY[i];
                R2Max = std::max_element(deltaR2,deltaR2+N)[0];
            }
            if(R2Max > RUpdate){
                UpdateList();
                R2Max = 0;
                for(int j = 0; j < N; j++){
                    X0[j] = X[j];
                    Y0[j] = Y[j];
                }
            }
        }
        // Writing observables to text file

        if(Find(samplePoints, dataPoints, 1.0*x) != -1){ // checking if saving time
            if(samplePoints[dataCounter] != 0){
                // Configs
                int idx = int (samplePoints[dataCounter]);
                log_cfg.open(out + "cfg_" + std::to_string(idx) + ".xy");
                log_cfg << std::scientific << std::setprecision(8);
                for (int i = 0; i<N; i++){
                    log_cfg << S[i] << " " << X[i] << " " << Y[i] << std::endl;
                }
                log_cfg.close();
                // Fs
                double FSavg = 0;
                for(int deg = 0; deg < 90; deg++){
                    FSavg += FS(0, steps, deg);
                }
                log_obs << samplePoints[dataCounter] << " " << VTotal()/(2*N) << " " 
                        << MSD() << " " << FSavg/90 << " " << CB() << std::endl;
                // saving format: timestep Vtot MSD Fs CB
            } 
            dataCounter++;
        }
        
        // Doing the MC
        for (int i = 0; i < N; i++){
            if (ranf() > 0.2) TryDisp(i); //Displacement probability 0.8
            else TrySwap(i,floor(ranf()*N)); //Swap probability 0.2
        }

        stepCounter ++; if(stepCounter%100==0) std:: cout << stepCounter << std::endl; // Counting steps
    }
}