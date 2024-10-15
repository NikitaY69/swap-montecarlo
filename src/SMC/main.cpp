#include "swap.h"

// Default run parameters
int N = 5;
double Size = std::sqrt (N);
double T = 0.04; 
int tau = 100000;
int tw = 1;
int cycles = 1;
int steps = tw*(cycles-1)+tau;
int linPoints = 50;
int logPoints = 50;
double p_swap = 0.2;
const int nr = 50;
const int ns = 100;

// Setting arrays
double *X = nullptr, *Y = nullptr, *S = nullptr, *Sref = nullptr, *X0 = nullptr, *Y0 = nullptr;
double *Xfull = nullptr, *Yfull = nullptr, *Xref = nullptr, *Yref = nullptr;
std::vector < std::vector <double>> Xtw, Ytw;
std::vector < std::vector<int> > NL, NN;
std::vector < std::vector < std::vector <int>>> NN_tw, RL;
std::vector < std::string > allObs;

std::string input;
std::string outdir;

//-----------------------------------------------------------------------------
//  main.cpp
int main(int argc, const char * argv[]) {
    
    // Random number generator
    srand(time(NULL)*1.0);

    // Define the command-line options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("input", po::value<std::string>(&input)->required(), "set input file")
        ("outdir", po::value<std::string>(&outdir)->required(), "set out directory")
        ("N", po::value<int>(&N)->default_value(N), "set system size")
        ("T", po::value<double>(&T)->default_value(T), "set temperature")
        ("tau", po::value<int>(&tau)->default_value(tau), "set single-run time")
        ("tw", po::value<int>(&tw)->default_value(tw), "set waiting time")
        ("cycles", po::value<int>(&cycles)->default_value(cycles), "set number of cycles")
        ("lin", po::value<int>(&linPoints)->default_value(linPoints), "set number of lin-spaced snapshots")
        ("log", po::value<int>(&logPoints)->default_value(logPoints), "set number of log-spaced snapshots")
        ("p_swap", po::value<double>(&p_swap)->default_value(p_swap), "set swap-attempt probability")
        ("MSD", "Flag to compute MSD")
        ("Cb", "Flag to compute Cb")
        ("Fs", "Flag to compute Fs")
        ("U", "Flag to compute U");
    // std::string input = motherdir + argv[1];
    // std::string outdir = motherdir + argv[2] + "results/";

    // Parse the command-line arguments
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (const po::error &ex) {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    // Handle the help option
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    // Parsing the observables in order of appearance
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--MSD") {
            allObs.push_back("MSD");
        } else if (arg == "--Cb") {
            allObs.push_back("Cb");
        } else if (arg == "--Fs") {
            allObs.push_back("Fs");
        } else if (arg == "--U") {
            allObs.push_back("U");
        }
    }

    // Resizing arrays
    Size = std::sqrt (N);
    steps = tw*(cycles-1)+tau;
    X = new double[N]; Y = new double[N]; S = new double[N]; Sref = new double[N]; 
    X0 = new double[N]; Y0 = new double[N];
    Xfull = new double[N]; Yfull = new double[N]; Xref = new double[N]; Yref = new double[N];

    // Xtw.resize(cycles, std::vector<double>(N));
    // Ytw.resize(cycles, std::vector<double>(N));
    // NL.resize(N, std::vector<int>(N));
    // NN.resize(N, std::vector<int>(N));
    // NN_tw.resize(N, std::vector<std::vector<int>>(tw, std::vector<int>(N)));
    // RL.resize(N, std::vector<std::vector<int>>(nr, std::vector<int>(N)));

    // creating outdir if not existing
    fs::path out_path = outdir;
    if(!fs::is_directory(out_path)){
        
        fs::create_directory(outdir);
    }
    
     // Writing params.txt file
    std::ofstream params;
    params.open(outdir + "params.txt");
    params << "rootdir" << " " << "N" << " " << "T" << " " 
           << "tau" << " " << "tw" << " " << "cycles" << " " 
           << "logPoints" << " " << "linPoints" << " " << "p_swap" << std::endl;
    params << outdir << " " << N << " " << T << " " << tau << " " << tw << " "
           << cycles << " " << logPoints << " " << linPoints << " " << p_swap << std::endl;
    params.close();

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
            Sref[i] = S[i];
            i++;}
        input_file.close();

    } else {
        std::cout << input << std::endl;
        return 0;
    }
    UpdateNL(); // First list of neighbours

    // Do simulation with timer
    double t0 = time(NULL); // Timer
    MC(outdir, logPoints, linPoints); 
    std::cout << "Time taken: " << (time(NULL) - t0) << "s" << std::endl; 
    std::cout << "Done" << std::endl;

    // Freeing allocated memory
    delete[] X; delete[] Y; delete[] S; delete[] Sref; delete[] X0; delete[] Y0;
    delete[] Xfull; delete[] Yfull; delete[] Xref; delete[] Yref;

    return 0;
}