#include "swap.h"

//  Calculates difference of a and b while applying periodic boundary conditions
double bcs(double a, double b) {return Size/2 - std::abs(std::abs(a-b)-Size/2);}

//  Finds index of element in array
int Find(std::vector <double> v, double seek){
    int i = 0;
    for (double v_i: v){
        if (v_i == seek) return i; 
        i++;
    }
    return -1;
}

// Computes the effective neighbours of particle j
std::vector<int> effective_neighbours(int j){
    std::vector<int> neigh;
    for (int i=0; i<N; i++){
        double xij = bcs(X[i], X[j]); double yij = bcs(Y[i], Y[j]);
        double rij2 = (xij*xij)+(yij*yij);
        if (rij2 < rNL && i != j){
            neigh.push_back(i);
        }
    }
    return neigh;
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
void UpdateList(){
    for (int i = 0; i < N; i++){
        NL.push_back(effective_neighbours(i));
    }
}
