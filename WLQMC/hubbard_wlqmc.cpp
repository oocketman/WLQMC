//
//  hubbard_wlqmc.cpp
//  WLQMC
//
//  Created by Hiren Thummar on 8/4/19.
//  Copyright © 2019 Hiren Thummar. All rights reserved.
//  World line Quantum Monte Carlo for one dimensional SU(N)
//

#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <cmath>
#include <functional>

using namespace std;

//function declaration
void Initialize();
inline int& set_elem(vector<int> &m_, size_t i_, size_t j_, size_t k_);
inline const int& get_elem(const vector<int> &m_, size_t i_, size_t j_, size_t k_);
int randInRange(int, int);
void printmat(vector<int>& v);
void printmatpart(vector<int>& v, int);
int modulus(int, int);

float S;            // spin type
int N_s;            // # of spin state
int N;              // # of spatial sites
int L;              // number of time slices
double V;           // potential energy
double Delta_tau;   // imaginary time step
int MC_steps;       // # of monte carlo steps
float t;            // hopping strength
float U;            // onsite interaction strength

float dE;

int main() {            //int main(int argc, const char * argv[])
    auto start = chrono::steady_clock::now();
    srand((unsigned int) time(NULL));
    
    
    double R, hb, random;
    int s, u, v, w, m;        // s = +/-2 (+ for L -> R; - for R -> L)
    
    cout << "World Line Quantum Monte Carlo for 1-Dimensional SU(N)\n"
    << " -----------------------------------------------------\n";
    // set simulation parameters
    cout << "Spin S = " << (S = 0.5)
    << "\n# of spatial sites N = " << (N = 10)
    << "\nHopping strength t = " << (t = 1.0)
    << "\nOnsite Interaction strength U = " << (U = 2.0)
    << "\nImaginary time step Delta_tau = " << (Delta_tau = 0.05) //0.25
    << "\nNumber of Monte Carlo steps = " << (MC_steps = 50000)
    << endl;
    Initialize();
    cout << "Number of time slices L = " << L << endl;
    
    // 2*LxNxN_s matrix filled with 0s
    vector<int> occ(2*L*N*N_s, 0);
    vector<int> nocc(2*L*N, 0);
    vector<double> E;
    // populate occ with 0's and 1's
    for (int k = 0; k < N_s; ++k) { //spin state
        for (int j = 0; j < N; ++j) { //real site (modify code for half filling every other site)
            for (int i = 0; i < 2*L; ++i) { //synthetic dimensions (imaginary time steps)
                set_elem(occ, i, (j+k)%N, k) = 1;
                nocc[i*N+(2*j+k)%N] += 1;
                dE += tanh(Delta_tau*t); //kinetic
            }}}
    
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < 2 * L; ++i) {
            if (nocc[i*N+j] >= 2) {
                dE += U*(nocc[i*N+j] - 1); //potential
            }}}
    
    E.push_back(dE);
    cout << E.back() << endl;
//    // print particle density vector
//    for (int j = 0; j < N; ++j) {
//        for (int i = 0; i < 2 * L; ++i) {
//            cout << nocc[i*N + j];
//        }
//        cout << "\n";
//    }
    
    
    ofstream data("MC.txt"); //, N, S);)
    
    data << N << " " << S << " " << N_s << " " << MC_steps << endl;
    
    // print whole tensor to file
    for (int k = 0; k < N_s; ++k) {
//        data << "N_s: " << k << endl;
        for (int i = 0; i < 2 * L; --i) {
            for (int j = 0; j < N; ++j) {
                data << get_elem(occ, i, j, k) << " ";
            }
            data << "\n";
        }
        data << "\n";
    }
    
    int lcount = 0;
    int rcount = 0;
    
    for (int step = 0; step < MC_steps; ++step) {       // loops over monte carlo steps
        dE = 0;
        for (int spin = 0; spin < N_s; ++spin) {        // loop over spin states
            for (int j = 0; j < N; ++j) {               // loop over spatial site
                for (int i = 0; i < 2*L; ++i) {         // loop over imaginary time steps
//                    cout << "Itime: " << i << " Site: " << j << " spin: " << spin << endl;
                    if ((get_elem(occ, i, j, spin) == 1 && get_elem(occ, (i+1)%(2*L), j, spin) == 1) && ((j % 2 == 0 && i % 2 == 1) || (j % 2 == 1 && i % 2 == 0)) && (get_elem(occ, i, (j+1)%N, spin) == 0 && get_elem(occ, (i+1)%(2*L), (j+1)%N, spin) == 0)) { // (even space, odd time) or (odd space, even time) World line moves L -> R
//                        cout << " L -> R" << endl;
                        s = 2;
                        u = 1 - get_elem(occ, ((2*L) + ((i-1)%(2*L)))%(2*L), (j+1)%N, spin) - get_elem(occ, (i+2)%(2*L), (j+1)%N, spin);
                        v = get_elem(occ, i, (N + ((j-1)%N))%N, spin) - get_elem(occ, i, (j+2)%N, spin);
                        int count = 0;
                        for (int u = 0; u < N_s; ++u) {
                            if (u != spin) {
                                count += get_elem(occ, i, (j+1)%N, u) + get_elem(occ, (i+1)%(2*L), (j+1)%N, u) - get_elem(occ, i, j, u) - get_elem(occ, (i+1)%(2*L), j, u);
                            }}
//                        cout << "u: " << u << " " << "v: " << v << " " << "U: " << U << endl;
                        R = pow(tanh(Delta_tau*t), s*u) * pow(cosh(Delta_tau*t), s*v) * exp(Delta_tau*V*s*v/2) * exp(-Delta_tau*U*count/2);
                        // heat bath algorithm for accepting and rejecting the step
                        hb = R/(1+R);
                        random = ((double) rand() / (RAND_MAX));
                        if (random < hb) {
//                            cout << "R: " << R << " " << "hb: " << hb << " " << "rand: " << random << endl;
                            set_elem(occ, i, j, spin) = 0;
                            set_elem(occ, (i+1)%(2*L), j, spin) = 0;
                            set_elem(occ, i, (j+1)%N, spin) = 1;
                            set_elem(occ, (i+1)%(2*L), (j+1)%N, spin) = 1;
                            nocc[i*N+j] += -1;
                            nocc[((i+1)%(2*L))*N+j] += -1;
                            nocc[i*N+(j+1)%N] += 1;
                            nocc[((i+1)%(2*L))*N+(j+1)%N] += 1;
//                            cout << "Itime: " << i << " Site: " << j << " spin: " << spin << endl;
//                            cout << "Step accepted! Left" << endl;
                            lcount += 1;
//                            printmatpart(occ, spin);
                            // calculate change in energy due to on-site interaction
                            if (nocc[i*N+j] >= 2) {
                                dE += U*(nocc[i*N+j] - 1);
                            }
                            if (nocc[((i+1)%(2*L))*N+j] >= 2) {
                                dE += U*(nocc[((i+1)%(2*L))*N+j] - 1);
                            }
                            if (nocc[i*N+(j+1)%N] >= 2) {
                                dE += U*(nocc[i*N+(j+1)%N] - 1);
                            }
                            if (nocc[((i+1)%(2*L))*N+(j+1)%N] >= 2) {
                                dE += U*(nocc[((i+1)%(2*L))*N+(j+1)%N] - 1);
                            } 
                            
                            // change in kinetic energy due to hopping
                            w = get_elem(occ, (i+2)%(2*L), j, spin) + get_elem(occ, ((2*L) + ((i-1)%(2*L)))%(2*L), j, spin) - get_elem(occ, ((2*L) + ((i-1)%(2*L)))%(2*L), (j+1)%N, spin) - get_elem(occ, (i+2)%(2*L), (j+1)%N, spin);
                            m = get_elem(occ, i, (N+((j-1)%N))%N, spin) + 2* get_elem(occ, i, (j+1)%N, spin);
                            
                            
                           
                            
                            
                            
                            
                        }}
                    
                    else if ((get_elem(occ, i, j, spin) == 1 && get_elem(occ, (i+1)%(2*L), j, spin) == 1) && ((j % 2 == 1 && i % 2 == 1) || (j % 2 == 0 && i % 2 == 0)) && (get_elem(occ, i, (N + ((j-1)%N))%N, spin) == 0 && get_elem(occ, (i+1)%(2*L), (N + ((j-1)%N))%N, spin) == 0)) { // (odd space, odd time) or (even space, even time) World line moves R -> L
//                        cout << " R -> L" << endl;
                        s = -2; // inverted equations for the sake of simplicity (R -> L)
                        u = 1 - get_elem(occ, (2*L +((i-1)%(2*L)))%(2*L), j, spin) - get_elem(occ, (i+2)%(2*L), j, spin);
                        v = get_elem(occ, i, (N + ((j-2)%N))%N, spin) - get_elem(occ, i, (j+1)%N, spin);
                        int count = 0;
                        for (int u = 0; u < N_s; ++u) {
                            if (u != spin) {
                                count += get_elem(occ, i, (N + ((j-1)%N))%N, u) + get_elem(occ, (i+1)%(2*L), (N + ((j-1)%N))%N, u) - get_elem(occ, i, j, u) - get_elem(occ, (i+1)%(2*L), j, u);
                            }}
//                        cout << "u: " << u << " " << "v: " << v << " " << "U: " << U << endl;
                        R = pow(tanh(Delta_tau*i), s*u) * pow(cosh(Delta_tau*i), s*v) * exp(Delta_tau*V*s*v/2) * exp(-Delta_tau*U*count/2);
                        // heat bath algorithm for accepting and rejecting the step
                        hb = R/(1+R);
                        random = ((double) rand() / (RAND_MAX));
                        if (random < hb) {
//                            cout << "R: " << R << " " << "hb: " << hb << " " << "rand: " << random << endl;
                            set_elem(occ, i, j, spin) = 0;
                            set_elem(occ, (i+1)%(2*L), j, spin) = 0;
                            set_elem(occ, i, (N + ((j-1)%N))%N, spin) = 1;
                            set_elem(occ, (i+1)%(2*L), (N + ((j-1)%N))%N, spin) = 1;
                            nocc[i*N+j] += -1;
                            nocc[((i+1)%(2*L))*N+j] += -1;
                            nocc[i*N+(N+((j-1)%N))%N] += 1;
                            nocc[((i+1)%(2*L))*N+(N+((j-1)%N))%N] += 1;
//                            cout << "Itime: " << i << " Site: " << j << " spin: " << spin << endl;
//                            cout << "Step accepted! Right" << endl;
                            rcount += 1;
//                            printmatpart (occ, spin);
                            // calculate change in energy due to on-site interaction
                            if (nocc[i*N+j] >= 2) {
                                dE += U*(nocc[i*N+j] - 1);
                            }
                            if (nocc[((i+1)%(2*L))*N+j] >= 2) {
                                dE += U*(nocc[((i+1)%(2*L))*N+j] - 1);
                            }
                            if (nocc[i*N+(N+((j-1)%N))%N] >= 2) {
                                dE += U*(nocc[i*N+(N+((j-1)%N))%N] - 1);
                            }
                            if (nocc[((i+1)%(2*L))*N+(N+((j-1)%N))%N] >= 2) {
                                dE += U*(nocc[((i+1)%(2*L))*N+(N+((j-1)%N))%N] - 1);
                            }
                            
                            // change in kinetic energy due to hopping
                            w = get_elem(occ, (i+2)%(2*L), j, spin) + get_elem(occ, ((2*L) + ((i-1)%(2*L)))%(2*L), j, spin) - get_elem(occ, ((2*L) + ((i-1)%(2*L)))%(2*L), (j+1)%N, spin) - get_elem(occ, (i+2)%(2*L), (j+1)%N, spin);
                            m = get_elem(occ, i, (N+((j-1)%N))%N, spin) + 2* get_elem(occ, i, (j+1)%N, spin);
                            
                            
                            
                            
                            
                        }}

                }}}
        // print whole tensor on screen
//        cout << "MC_step: " << step << endl;
        for (int k = 0; k < N_s; ++k) {
//            data << "N_s: " << k << endl;
            for (int j = 0; j < N; ++j) {
                for (int i = 0; i < 2 * L; ++i) {
                    data << get_elem(occ, i, j, k) << " ";
                }
                data << "\n";
            }
            data << "\n";
        }
    }
    data.close();
    
    cout << "lcount: " << lcount << " " << "rcount: " << rcount << endl;
    cout << "Energy: " << E.back()/(2*L*N) << endl;
 
    auto end = chrono::steady_clock::now();
    chrono::duration<double> diff = end - start;
    cout << "Time: " << diff.count() << '\n';
    
        // print particle density vector
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < 2 * L; ++i) {
            cout << nocc[i*N + j];
        }
        cout << "\n";
    }
    
//    return 0;
}


// Functions
void Initialize() {
    V = 0.0;        // considering zero potential energy case
    N_s = 2*S + 1;
    L = 1/Delta_tau;
}

inline int& set_elem(vector<int> &m_, size_t i_, size_t j_, size_t k_) {
    return m_[i_*N*N_s + j_*N_s + k_];
}

inline const int& get_elem(const vector<int> &m_, size_t i_, size_t j_, size_t k_) {
    return m_[i_*N*N_s + j_*N_s + k_];
}

// generates random number (double) in range(min, max)
int randInRange(int min, int max) {
    return (int) (rand()%max + min);
}

void printmat(vector<int>& v) {
    // print whole tensor on the screen
    for (int k = 0; k < N_s; ++k) {
        cout << "N_s: " << k << endl;
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < 2 * L; ++i) {
                cout << get_elem (v, i, j, k);
            }
            cout << "\n";
        }
        cout << "\n" << endl;
    }
}

void printmatpart(vector<int>& v, int s) {
    // print whole tensor on the screen
    cout << "N_s: " << s << endl;
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < 2 * L; ++i) {
            cout << get_elem (v, i, j, s);
        }
        cout << "\n";
    }
    cout << "\n" << endl;
}

int modulus(int i_, int n_) {
    return (int) (n_ + (i_%n_))%n_;
}
