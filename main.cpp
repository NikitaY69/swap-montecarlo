#include "swap.h"

// Run parameters
const int tau = 5000000;
const int tw = 1;
const int cycles = 1;
const int steps = tw*(cycles-1)+tau;
const double T = 0.04; 
const int nr = 50;

// Snapshots
const int dataPoints = 100;

// Initialization of external variables
double X[N], Y[N], S[N], X0[N], Y0[N];
double Xfull[N], Yfull[N], Xref[N], Yref[N];
std::vector < std::array <double, N>> Xtw, Ytw;
std::vector < std::vector<int> > NL(N), NN(N);
std::vector < std::vector < std::vector <int>>> NN_tw;
std::vector < std::vector < std::vector <int>>> RL(N, std::vector < std::vector <int>>(nr));
//-----------------------------------------------------------------------------
//  main.cpp
int main(int argc, const char * argv[]) {
    
    // User-defined variables
    srand(time(NULL)*1.0); //Random number generator
    std::string input = motherdir + argv[1];
    std::string outdir = motherdir + argv[2] + "results/";
    fs::path out_path = outdir;
    if(!fs::is_directory(out_path)){
        // creating outdir if not existing
        fs::create_directory(outdir);
    }

    // // Do simulation with timer
    double t0 = time(NULL); // Timer
    MC(input, outdir, dataPoints); 
    std::cout << "Time taken: " << (time(NULL) - t0) << "s" << std::endl; 
    std::cout << "Done" << std::endl;
    return 0;
}