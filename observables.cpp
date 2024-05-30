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
            sum += deltaX*deltaX + deltaY*deltaY;
    }
    return sum/N;
}

// Correlation functions

//  Calculates the self scattering function
double FS(int cycle, double theta){
    double dotProduct;
    double q = 2*pi/sigmaMax;
    double sum = 0, deltaX, deltaY;
    
    for (int i = 0; i < N; i++){
        deltaX = Xfull[i]-Xtw[cycle][i];
        deltaY = Yfull[i]-Ytw[cycle][i];
        dotProduct = q*((cos(theta*pi/180)*deltaX)+(sin(theta*pi/180)*deltaY));
        sum += cos(dotProduct);
    }
    return sum/N;
}

// Computes the bond-breaking correlation function (local)
double CBLoc(int cycle, int j){
    std::vector<int> intersect;
    std::vector<int> nn0 = nn_tw[cycle][j]; // neighbors at t=0
    std::vector<int> nn = nearest_neighbours(j, x_max);
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

// Computes the averaged scalar product over all pair of particles
double DispCorrLoc(int j){
    double sum = 0;
    double deltaXi, deltaYi, deltaXj, deltaYj;
    deltaXj = Xfull[j]-Xref[j], deltaYj = Yfull[j]-Yref[j];
    for (int i=0; i<N; i++){
        if (i!=j){
            deltaXi=Xfull[i]-Xref[i]; deltaYi=Yfull[i]-Yref[i];
            sum += deltaXi*deltaXj + deltaYi*deltaYj;
        }
    }
}

double DispCorr(){
    double sum = 0;
    for (int i=0;i<N;i++){
        sum += DispCorrLoc(i);
    } return sum/(2*N);
}

double MicroDispCorrLoc(int j, double r){
    double sum = 0;
    double deltaXi, deltaYi, deltaXj, deltaYj;
    deltaXj = Xfull[j]-Xref[j], deltaYj = Yfull[j]-Yref[j];
    std::vector <int> nn_r = radius_neighbours(j, r);
    for (int i: nn_r){
        deltaXi=Xfull[i]-Xref[i]; deltaYi=Yfull[i]-Yref[i];
        sum += deltaXi*deltaXj + deltaYi*deltaYj;
    }
}

double MicroDispCorr(double r){
    double sum = 0;
    for (int i=0;i<N;i++){
        sum += MicroDispCorrLoc(i, r);
    } return sum/(2*N);
}

// Updates the reference points for the correlation functions
void UpdateAge(int cycle){
    nn_tw.push_back(std::vector < std::vector <int>>());
    for (int i=0; i<N; i++){
        nn_tw[cycle].push_back(nearest_neighbours(i, x_max));
        Xtw[cycle][i] = Xfull[i];
        Ytw[cycle][i] = Yfull[i];
    }
}