#include "swap.h"

double dXCM, dYCM;
int dataCounter=0;
int cycle;
// double swapCount[N] = {0};

// Monte Carlo Simulation
void MC(std::string out, int n_log, int n_lin){
    int cycleCounter = 0;
    double deltaX[N], deltaY[N], deltaR2[N], R2Max = 0;
    // Building snapshots list (log-spaced)
    std::vector < std::pair <double, double>> pairs;
    std::vector <double> samplePoints, twPoints;
    double endingPoints[cycles], linPoints[n_lin];
    double exponents = log10(tau)/(n_log-1);

    for(int c=0; c<cycles; c++){
        for (int x = 0; x < n_log; x++){
            double value = tw*c + floor(pow(10,exponents*(x)));
            std::pair <double,double> p = {value, c};
            int f = std::count(pairs.begin(), pairs.end(), p);
            if(f==0){
                pairs.emplace_back(value, c);
            // this if condition is actually relevent because of the floor function
            }
        }
    }

    // Sorting
    std::sort(pairs.begin(), pairs.end());
    for (auto p: pairs){
        samplePoints.push_back(p.first); twPoints.push_back(p.second);
    }

    // Ending points
    for(int c=0;c<cycles;c++){
        endingPoints[c] = c*tw + tau;
    }
    // Linspaced points
    for (int k=1;k<=n_lin;k++){
        linPoints[k] = (tau/(n_lin))*k;
    }
    // File writing
    std::ofstream log_obs, log_cfg; 
    // log_ploc, log_p;
    // log_sigma;
    std::string out_cfg = out + "configs/";
    // std::string out_ploc = out + "micro_corr/";
    // std::string out_sigma = out + "sigma_scan/";
    log_obs.open(out + "obs.txt");
    log_obs << "t" << " " << "cycle";
    for (std::string obs: allObs){
        log_obs << " " << obs;
    } log_obs << std::endl;
    // log_p.open(out + "space_corr.txt");
    log_obs << std::scientific << std::setprecision(8);
    // log_p << std::scientific << std::setprecision(8);
    // creating outdir if not existing
    fs::create_directory(out_cfg); 
    // fs::create_directory(out_ploc);
    // fs::create_directories(out_sigma); 

    for(int t = 1; t <= steps; t++){
        // Updating NL
        if((t-1) % 150 == 0) {//Change number?
            // every 150 steps we check if we need to update the NL
            for (int i = 0; i < N; i++){
                deltaX[i] = bcs(X[i],X0[i]);
                deltaY[i] = bcs(Y[i],Y0[i]);
                deltaR2[i] = deltaX[i]*deltaX[i] + deltaY[i]*deltaY[i];
            R2Max = std::max_element(deltaR2,deltaR2+N)[0];
            }
            if(R2Max > RUpdate){
                UpdateNL();
                R2Max = 0;
                for(int j = 0; j < N; j++){
                    X0[j] = X[j];
                    Y0[j] = Y[j];
                }
            }
        }
    
        // Updating reference observables
        if((t-1)%tw == 0 && cycleCounter < cycles){
            UpdateAge(cycleCounter); cycleCounter++;
        } 

        // // Writing observables to text file
        int lin = std::count(linPoints, linPoints+n_lin, 1.0*t);
        int log = std::count(samplePoints.begin(), samplePoints.end(), 1.0*t);

        if(lin>0){ // checking if saving time
            // UpdateRL(); // updating nearest neighbours
            // Configs
            log_cfg.open(out_cfg + "cfg_" + std::to_string(t) + ".xy");
            // log_ploc.open(out_ploc + "corr_" + std::to_string(t) + ".txt");
            // log_sigma.open(out_sigma + "scan_" + std::to_string(t) + ".txt");
            log_cfg << std::scientific << std::setprecision(8);
            // log_ploc << std::scientific << std::setprecision(8);
            // log_sigma << std::scientific << std::setprecision(8);
            for (int i = 0; i<N; i++){
                // std::vector <double> disp_loc = MicroDispCorrLoc(i);
                // std::vector <double> u_sigma = SigmaScan(i);
                log_cfg << S[i] << " " << Xfull[i] << " " << Yfull[i] << std::endl;
                // for (int k=0;k<nr;k++){
                //     log_ploc << disp_loc[k] << " ";
                // } log_ploc << std::endl;
                // for (int k=0;k<ns;k++){
                //     log_sigma << u_sigma[k] << " ";
                // } log_sigma << std::endl;
            }
            log_cfg.close();
            // log_sigma.close(); log_ploc.close();
        }
        if(log>0){ // checking if saving time
            UpdateNN(); // updating nearest neighbours
            // UpdateRL(); // updating per-radius neighbour-list
            dXCM = 0; dYCM = 0;
            for (int i=0;i<N;i++){
                double dX = Xfull[i]-Xref[i], dY = Yfull[i]-Yref[i];
                dXCM += dX; dYCM += dY;
            } dXCM /= N; dYCM /= N;

            for(int s=0; s<log; s++){
                // looping different eventual tws
                cycle = twPoints[dataCounter];
                // Configs
                if(! fs::exists (out_cfg + "cfg_" + std::to_string(t) + ".xy")){
                    log_cfg.open(out_cfg + "cfg_" + std::to_string(t) + ".xy");
                    log_cfg << std::scientific << std::setprecision(8);
                    for (int i = 0; i<N; i++){
                        log_cfg << S[i] << " " << Xfull[i] << " " << Yfull[i] << std::endl;
                    }
                    log_cfg.close();
                } 

                // log_p << t << " ";
                // std::vector <double> disp = MicroDispCorr();
                // for (int k=0;k<nr;k++){
                //     log_p << disp[k] << " ";
                // } log_p << std::endl;

                log_obs << t << " " << cycle;
                for (std::string obs: allObs){
                    log_obs << " " << whichObs(obs, cycle);
                } log_obs << std::endl;

                dataCounter++;
            }  
        };
        // Doing the MC
        for (int i = 0; i < N; i++){
            // TryDisp(i);
            if (ranf() > p_swap) TryDisp(i); //Displacement probability 0.8
            else TrySwap(i,floor(ranf()*N)); //Swap probability 0.2
        }
        
        if((t-1)%100==0) std::cout << (t-1) << std::endl;; // Counting steps
    };
    log_obs.close();
    // log_p.close();
}

//  Tries displacing one particle j by vector dr = (dx, dy)
void TryDisp(int j){
    double dx = (ranf()-0.5)*deltaMax;
    double dy = (ranf()-0.5)*deltaMax;
    double Xnew = Pshift(X[j]+dx);
    double Ynew = Pshift(Y[j]+dy);
    double deltaE = V(Xnew, Ynew, S[j], j) - V(X[j], Y[j], S[j], j);
    // why is the modulus function not in deltaE ?
    if (deltaE < 0){
        // Xnew = fmod(X[j],Size);
        X[j] = Xnew; //Check modulus function
        Y[j] = Ynew;
        Xfull[j] = Xfull[j]+dx;
        Yfull[j] = Yfull[j]+dy;
    }
    else if (exp(-deltaE/T) > ranf()){
        X[j] = Xnew;
        Y[j] = Ynew;
        Xfull[j] = Xfull[j]+dx;
        Yfull[j] = Yfull[j]+dy;
    }
}

//  Tries swapping the pair of particles j, k
void TrySwap(int j, int k){
    double deltaS = std::abs (S[j]-S[k]);
    if(deltaS<=deltaSMax){
        double deltaE = V(X[j],Y[j],S[k],j)+V(X[k],Y[k],S[j],k)-V(X[j],Y[j],S[j],j)-V(X[k],Y[k],S[k],k);
        if (deltaE < 0){
            double Rnew = S[k];
            S[k] = S[j];
            S[j] = Rnew;
            // swapCount[j] += 1; swapCount[k] += 1;
        }
        else if (exp(-deltaE/T) > ranf()){
            double Rnew = S[k];
            S[k] = S[j];
            S[j] = Rnew;
            // swapCount[j] += 1; swapCount[k] += 1;
        }
    } else{
        // pass
    }
}

double whichObs(std::string obs, int cycl){
    if (obs=="MSD") return MSD();
    else if (obs=="U") return VTotal()/(2*N);
    else if (obs=="Cb") return CB(cycl);
    else if (obs=="Fs") return FS(cycl);
}