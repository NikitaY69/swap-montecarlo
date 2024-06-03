#include "swap.h"

//  Calculates difference of a and b while applying periodic boundary conditions
double bcs(double a, double b) {return Size/2 - std::abs(std::abs(a-b)-Size/2);}

double Pshift(double a){
    return a - Size*floor((a+Size/2)/Size);
}

// Computes the pseudo-interacting neighbours list
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

// Computes the nearest neighbours list
void UpdateNN(){
    NN.clear(); NN = std::vector < std::vector <int> > (N);
    for (int j=0; j<N-1; j++){
        for (int i=j+1; i<N; i++){
            double sigmaij = (S[i]+S[j])*(1-0.2*std::abs(S[i]-S[j]))/2;
            double xij = bcs(X[i], X[j]); double yij = bcs(Y[i], Y[j]);
            double rij = sqrt((xij*xij)+(yij*yij));
            if (rij < x_max*sigmaij && i != j){
                NN[j].push_back(i);
                NN[i].push_back(j);
            }
        } 
    }
}