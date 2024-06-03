#include "swap.h"
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

// Run parameters
const int tau = 100000;
const int tw = 50000;
const int cycles = 1;
const int steps = tw*(cycles-1)+tau;
const double T = 0.1; 
std::string motherdir = "/home/allaglo/benchmarks/";

// Snapshots
const int dataPoints = 30;

// Initialization of external variables
double X[N], Y[N], S[N], X0[N], Y0[N];
double Xfull[N], Yfull[N], Xref[N], Yref[N];
std::vector < std::array <double, N>> Xtw, Ytw;
std::vector < std::vector<int> > NL(N), NN(N);
std::vector < std::vector < std::vector <int>>> NN_tw;
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
    
    // Read init config
    std::string line;
    std::ifstream input_file(input);
    if (input_file.is_open()){
        int i = 0; // particle index
        std::vector<std::vector<double>> cfg; // array of configurations
        while (std::getline(input_file, line)){
            double value;
            std::stringstream ss(line);

            cfg.push_back(std::vector<double>());
            while (ss >> value){
                cfg[i].push_back(value);
            }
            S[i] = cfg[i][0]; X[i] = cfg[i][1]; Y[i] = cfg[i][2];
            X0[i] = X[i]; Xfull[i] = X[i]; Xref[i] = X[i];
            Y0[i] = Y[i]; Yfull[i] = Y[i]; Yref[i] = Y[i];
            i++;}
        input_file.close();

    } else {
        std::cout << input << std::endl;
        return 0;
    }

    UpdateNL(); // First list of neighbours

    // // Do simulation with timer
    double t0 = time(NULL); // Timer
    MC(outdir, dataPoints); 
    std::cout << "Time taken: " << (time(NULL) - t0) << "s" << std::endl; 
    std::cout << "Done" << std::endl;
    return 0;
}