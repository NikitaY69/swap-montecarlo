#include "swap.h"

//  Calculates difference of a and b while applying periodic boundary conditions
double bcs(double a, double b) {return Size/2 - std::abs(std::abs(a-b)-Size/2);}

double Pshift(double a){
    return a - Size*floor((a+Size/2)/Size);
}

// Computes the effective neighbours of particle j
void UpdateNL(){
    NL.clear(); NL = std::vector < std::vector <int> > (N);
    for (int j=0; j<N-1; j++){
        for (int i=j+1; i<N; i++){
            double xij = bcs(X[i], X[j]); double yij = bcs(Y[i], Y[j]);
            double rij2 = (xij*xij)+(yij*yij);
            if (rij2 < rNL && i != j){
                NL[j].push_back(i);
                NL[i].push_back(j);
            }
        }
    }
}

// Computes the nearest neighbours of particle j at a given radius
std::vector<int> nearest_neighbours(int j, double x){
    std::vector<int> nn;
    for (int i=0; i<N; i++){
        double sigmaij = (S[i]+S[j])*(1-0.2*std::abs(S[i]-S[j]))/2;
        double xij = bcs(X[i], X[j]); double yij = bcs(Y[i], Y[j]);
        double rij = sqrt((xij*xij)+(yij*yij));
        if (rij < x*sigmaij && i != j){
            nn.push_back(i);
        }
    } return nn;
}

//  Creates the neighbour list for the set of particles
// void UpdateList(){
//     for (int i = 0; i < N; i++){
//         NL.push_back(effective_neighbours(i));
//     }
// }
