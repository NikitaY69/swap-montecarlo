#include "swap.h"

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
