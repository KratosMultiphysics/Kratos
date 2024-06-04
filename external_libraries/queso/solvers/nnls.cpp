//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

//// STL includes
#include <iostream>
#include <algorithm>
//// External includes
#include "nnls/nnls_impl.h"
//// Project includes
#include "queso/solvers/nnls.h"

namespace queso {

double NNLS::solve(MatrixType& A, VectorType& B, VectorType& X) {

    // Get Dimension
    int m = B.size();
    int n = (m > 0) ? A.size()/m : 0;
    X.resize(n);

    // Pointer to doubles
    double *A_doubles = A.data();
    double *B_doubles = B.data();
    double *X_double = X.data();

    // Construct buffers
    double *W = new double[n];
    double *ZZ = new double[m];
    int *index = new int[n];

    double Rnorm = 0; // Residual norm

    int mode = 0;
    int mda = m;

    // Solve
    nnls_(A_doubles, &mda, &m, &n, B_doubles, X_double, &Rnorm, W, ZZ, index, &mode);

    if( mode == 2 ){
        std::cerr << "THE DIMENSIONS OF THE PROBLEM ARE BAD. EITHER M .LE. 0 OR N .LE. 0." << std::endl;
        Rnorm = 1e8;
    }
    else if( mode == 3){
        std::cerr << "ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS." << std::endl;
    }

    // Release buffers
    delete[] W;
    delete[] ZZ;
    delete[] index;

    return Rnorm;
}

} // End namespace queso

