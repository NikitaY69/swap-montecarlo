#include "swap.h"

//  Calculates the potential of a pair of particles
double PairPotential(double x1, double y1, double s1, double x2, double y2, double s2){
    double sigmaij = (s1+s2)*(1-0.2*std::abs(s1-s2))/2;
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
    for (int k: NL[j]){
        total += PairPotential(xj, yj, rj, X[k], Y[k], S[k]);
    } return total;
}

//  Calculates total system energy
double VTotal(){
    double vTot = 0;
    for (int j = 0; j < N; j++)
        vTot += V(X[j], Y[j], S[j], j);
    return vTot;
}

//  Calculates avg. mean square displacements
double MSD(){
    double sum = 0, deltaX, deltaY;
        for (int i = 0; i < N; i++){
            deltaX = Xfull[i]-Xref[i];
            deltaY = Yfull[i]-Yref[i];
            deltaX -= dXCM; deltaY -= dYCM;
            sum += deltaX*deltaX + deltaY*deltaY;
    }
    return sum/N;
}

// Correlation functions

//  Calculates the self scattering function
double FS(int cycle){
    double dotProduct;
    double q = 2*pi/sigmaMax;
    double sum = 0, deltaX, deltaY;
    int ang = 360;
    
    for (int theta=0; theta<ang; theta++){
        for (int i = 0; i < N; i++){
            deltaX = Xfull[i]-Xtw[cycle][i];
            deltaY = Yfull[i]-Ytw[cycle][i];
            deltaX -= dXCM; deltaY -= dYCM;
            dotProduct = q*((cos(theta*pi/180)*deltaX)+(sin(theta*pi/180)*deltaY));
            sum += cos(dotProduct);
        }
    }
    return sum/(ang*N);
}

// Computes the bond-breaking correlation function (local)
double CBLoc(int cycle, int j){
    std::vector<int> intersect;
    std::vector<int> nn0 = NN_tw[cycle][j]; // neighbors at t=0
    std::vector<int> nn = NN[j];
    std::set_intersection(nn0.begin(), nn0.end(), nn.begin(), nn.end(),
                     std::back_inserter(intersect));

    if (nn0.size()==0){
        return 0;
    } else { 
        double frac = intersect.size()/nn0.size();
        return frac;
    } 
}

// Computes the bond-breaking correlation function (averaged)
double CB(int cycle){
    double tot = 0;
    for (int j=0; j<N; j++){
        tot += CBLoc(cycle, j);
    } return tot/N;
}

// Computes the averaged local displacement correlation over all pair of particles
double DispCorrLoc(int j){
    double sum = 0;
    double deltaXi, deltaYi, deltaXj, deltaYj;
    deltaXj = Xfull[j]-Xref[j], deltaYj = Yfull[j]-Yref[j];
    deltaXj -= dXCM; deltaYj -= dYCM;
    for (int i=0; i<N; i++){
        if (i!=j){
            deltaXi=Xfull[i]-Xref[i]; deltaYi=Yfull[i]-Yref[i];
            deltaXi -= dXCM; deltaYi -= dYCM;
            sum += deltaXi*deltaXj + deltaYi*deltaYj;
        }
    }
    return sum/(N-1);
}

// Global displacement correlation
double DispCorr(){
    double sum = 0;
    for (int i=0;i<N;i++){
        sum += DispCorrLoc(i);
    } return sum/N;
}

// Per-radius local displacement correlation
std::vector <double> MicroDispCorrLoc(int j){
    std::vector <double> sum(nr, 0);
    double deltaXi, deltaYi, deltaXj, deltaYj;
    deltaXj = Xfull[j]-Xref[j], deltaYj = Yfull[j]-Yref[j];
    deltaXj -= dXCM; deltaYj -= dYCM;
    for (int k=0; k<nr; k++){
        if (RL[j][k].size()==0) sum[k] = 0;
        else{
            for (int i: RL[j][k]){
                deltaXi=Xfull[i]-Xref[i]; deltaYi=Yfull[i]-Yref[i];
                deltaXi -= dXCM; deltaYi -= dYCM;
                sum[k] += deltaXi*deltaXj + deltaYi*deltaYj;
            } sum[k] /= RL[j][k].size();
        }
    } return sum;
}

// Per-radius global dispalcement correlation
std::vector <double> MicroDispCorr(){
    std::vector <double> sum(nr, 0);
    for (int i=0; i<N; i++){
        std::vector <double> disp_i = MicroDispCorrLoc(i);
        for (int k=0;k<nr;k++){
            sum[k] += disp_i[k];
        } 
    } return sum;
}

// Computes the sigma-dependence of the potential energy at the particle level
std::vector <double> SigmaScan(int j){
    double dS = 0.002;
    double s;
    std::vector <double> Vs(ns, 0);
    for (int k=0; k<ns; k++){
        s = S[j]+k*dS;
        Vs[k] = V(X[j], Y[j], s, j);
    } return Vs;
}

// Computes the diameter auto-correlation function
double C_sigma(){
    double sigma_m = 1.000218223;
    double C = 0, C0 = 0;
    for (int i=0; i<N; i++){
        double deltaS0 = Sref[i]-sigma_m; double deltaS = S[i]-sigma_m;
        C0 += deltaS0*deltaS0; C += deltaS*deltaS0;
    } return C/C0;
}

// Updates the reference points for the correlation functions
void UpdateAge(int cycle){
    UpdateNN(); NN_tw.push_back(NN);
    Xtw.push_back(std::vector <double>());
    Ytw.push_back(std::vector <double>());
    for (int i=0; i<N; i++){
        Xtw[cycle].push_back(Xfull[i]);
        Ytw[cycle].push_back(Yfull[i]);
    }
}

double whichObs(int idx){
    if (idx == 0){
        
    }
}