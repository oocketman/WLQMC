//
//  main.cpp
//  test
//
//  Created by Hiren Thummar on 8/8/19.
//  Copyright © 2019 Hiren Thummar. All rights reserved.
//
// LAPACK test code
//compile with: g++ main.cpp -llapack -lblas -o test

//#include <iostream>
//#include <Eigen/Eigenvalues>
//using namespace std;
//using namespace Eigen;
//int main () {
//    // for loop execution
//    for( int a = -10; a < 20; a = a + 1 ) {
//        cout << "value of a: " << a << endl;
//
//        if (a == -5 || a == 5) {
//            cout << "print: " << a << endl;
//        }
//    }
//    cout <<  tanh(2.05) << endl;
//    MatrixXd A = MatrixXd::Random(6,6);
//    cout << "Here is a random 6x6 matrix, A:" << endl << A << endl << endl;
//
//    EigenSolver<MatrixXd> es(A);
//    MatrixXd D = es.pseudoEigenvalueMatrix();
//    MatrixXd V = es.pseudoEigenvectors();
//    cout << "The pseudo-eigenvalue matrix D is:" << endl << D << endl;
//    cout << "The pseudo-eigenvector matrix V is:" << endl << V << endl;
//    cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;
//    return 0;
//}


//#include <iostream>
//#include <Eigen/Dense>
//using Eigen::MatrixXd;
//int main()
//{
//    MatrixXd m(2,2);
//    m(0,0) = 3;
//    m(1,0) = 2.5;
//    m(0,1) = -1;
//    m(1,1) = m(1,0) + m(0,1);
//    std::cout << m << std::endl;
//}

#include <vector>
#include <iostream>
#include <ostream>
#include <fstream>
#include <bitset>

using namespace std;

int main() {
    int L = 10;
    int N = 1;
    int N_s = 8;
    vector<bool> occ(2*L*N*N_s);
    bitset<160> occb;
    cout << occ.size() << '\n' << occ.capacity() << '\n' << sizeof(occ) << endl;
    cout << occb.size() << '\n' << sizeof(occb) << endl;
    
//    ofstream data("mat.txt");
//    for (int k = 0; k < N_s; ++k) { //spin state
//        for (int j = 0; j < N; ++j) { //real site (modify code for half filling every other site)
//            for (int i = 0; i < 2*L; ++i) { //synthetic dimensions (imaginary time steps)
//                //                spin, rsite, synth dim, mat elem
//                data << k << "," << j << "," << i << "," << (i*N*N_s + j*N_s + k) << endl;
//            }
//            data << "\n";
//        }
//        data << "\n";
//    }
//    data.close();
    
    return 0;
}
