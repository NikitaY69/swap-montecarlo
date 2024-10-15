#include "swap.h"

// Read params file
void ReadParams(std::string outdir){
    std::string line;
    std::ifstream input_file(outdir + "params.txt");
    if (input_file.is_open()){
        int idx = 0; // line index
        while (std::getline(input_file, line)){
            if (idx==1){
                std::stringstream ss(line);
                std::string discard;
                ss >> discard >> N >> T >> tau >> tw >> cycles >> logPoints
                   >> linPoints >> p_swap;
            } idx++;
        }
        input_file.close();
        
    } else {
        std::string error = outdir + "params.txt" + " not found";
        throw error;
    }
}
// Read init config
void ReadCFG(std::string input){
    auto start = std::chrono::high_resolution_clock::now();
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
            S[i] = cfg[i][0]; Xfull[i] = cfg[i][1]; Yfull[i] = cfg[i][2];
            X[i] = Pshift(Xfull[i]); Y[i] = Pshift(Yfull[i]);
            i++;}
        input_file.close();
        
    } else {
        std::string error = input + " not found";
        throw error;
    }
}
