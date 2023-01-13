//
//  hubbard_hamiltonian.cpp
//  hubbard_hamiltonian
//
//  Created by Hiren Thummar on 8/7/19.
//  Copyright © 2019 Hiren Thummar. All rights reserved.
//  Exact Diagonalization for HH
//  g++ -I /usr/local/include/eigen3/ hubbard_hamiltonian.cpp
//

#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <functional>
#include <bitset>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void Initialize();
int nCr(int, int);
long fact(long);
int randInRange(int min, int max);
vector<int> t_hop(int n_);
int u_intractn(int n_);
vector< vector<int> > hamiltonian(const int, const float);

// initializing constants
constexpr int N = 2;
constexpr float S = 1.5;
constexpr int N_s = 2*S + 1;
constexpr int np = N*N_s; // binary size NxN_s bit (or # of max. particles can be on lattice)

int nstate;
int res;
float t;
float U;
float V;

int main() {
    srand((unsigned int) time(NULL));
    Initialize();

//    ofstream data("hh_eigen_val_vec.txt"); // change cout to data 177013
    cout << "Hubbard Hamiltonian Exact Diagonalization for 1-Dimensional SU(N)\n"
    << "-----------------------------------------------------------------" << endl; // 177013
    // set simulation parameters
    cout << "# of spatial sites N = " << N
    << "\nSpin S = " << S
    << "\n# of spin flavors N_s = " << N_s
    << "\n# of possible states nstate = " << nstate
    << "\nlength of binary string = " << np
    << "\nHopping strength t = " << (t = 1.0)
    << "\nOn-site Interaction strength U = " << (U = 8.0) << endl; // 177013
    
    // hamiltonian
    vector< vector<double> > H(nstate, vector<double>(nstate, 0.0));
    vector<int> ht;
    for (int i = 1; i < nstate; ++i) {
        ht.clear();
        ht = t_hop(i);
        for (int j = 0; j < ht.size(); ++j) {
            H[i][ht[j]] += -t;
        }
        H[i][i] += u_intractn(i)*U;
    }

    int temp = nstate;
    MatrixXd matA(temp, temp);
    for (int k = 0; k < temp; ++k) {
        for (int l = 0; l < temp; ++l) {
            matA(k, l) = H[k][l];
        }
    }
    
    cout << "hamiltonian:\n" << matA << endl;
    EigenSolver<MatrixXd> es;
    es.compute(matA, true);
    MatrixXcd e_val = es.eigenvalues();
    MatrixXcd e_vec = es.eigenvectors();
    cout << "Eigenvalues:\n" << e_val << "\nEigenvectors:\n" << e_vec << endl; // 177013
//    data.close(); // 177013
    return 0;
}

// Functions
void Initialize() {
    V = 0.0;        // considering zero potential energy case
    res = 0;
    for (int i = 0; i<=N_s; ++i) {
        res += nCr(N_s, i);
    }
    nstate = pow(res, N);
}

int nCr(int n, int r) {     // calculate nCr = fact(n)/(fact(r)*fact(n-r))
    return fact(n) / (fact(r) * fact(n - r));
}

int fact(int n) { // calculate factorial of any positive integer
    return (n==1 || n==0) ? 1: n * fact(n - 1);
}

int randInRange(int min, int max) {     // generates random number (double) in range(min, max)
    return (int) (rand()%max + min);
}

vector<int> t_hop(int n_) {
    Initialize();
    vector<int> temp;
    vector<int> hop;
    int count;
    bitset<np> foo(n_);
    foo.to_string();
    for (int i = 0; i < N_s; ++i) {
        temp.clear();
        count = 0;
        for (int j = 0; j < N; ++j) {
            temp.push_back(foo.test(i+j*N_s));
            count += int(foo.test(i+j*N_s));
        }
        if (count > 0 && count < N) {
            bitset<np> state_(foo);
            vector<int> temp_(temp);
            sort(temp.begin(),temp.end());
            do {
                if (temp != temp_) {
                    for (int k = 0; k < N; ++k) {
                        state_.set(i+k*N_s, temp[k]);
                    }
                    hop.push_back((int) (state_.to_ulong()));
                }
            } while ( next_permutation(temp.begin(), temp.end()) );
        }
    }
    return hop;
}

int u_intractn(int n_) {
    Initialize();
    int count;
    int U_s = 0;
    bitset<np> foo(n_);
    foo.to_string();
    for (int i = 0; i < N; ++i) {
        count = 0;
        for (int j = 0; j < N_s; ++j) {
            count += int(foo.test(i*N_s+j));
        }
        if (count >= 2) {
            U_s += (count - 1);
        }
    }
    return U_s;
}
