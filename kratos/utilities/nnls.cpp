//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Meßmer
//                   Ricky Aristio
//

//// External includes
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "nnls/nnls_impl.h"
//// STL includes
#include <iostream>
#include <algorithm>
//// Project includes
#include "nnls.h"

namespace Kratos {

double nnls::nnls(MatrixType& A, const VectorType& b, VectorType& x){
    typedef boost::numeric::ublas::matrix<double> MatrixType;

    // Get Dimension
    int m = A.size1();
    int n = A.size2();

    // Pointer to doubles
    double *A_doubles = 0;
    double *B_doubles = 0;

    // Convert Ublas matrix and vector to a pointer to double.
    A_doubles = new double[m * n];
    B_doubles = new double[m];
    std::vector<double> tmp_vector(m*n);
    for( int i = 0; i < n; ++i){
        boost::numeric::ublas::matrix_column<MatrixType> mr(A, i);
        std::copy(mr.begin(), mr.end(), A_doubles + i*m);
    }

    std::copy(b.begin(), b.end(), B_doubles);

    // Construct buffers
    double *X = new double[n];
    double *W = new double[n];
    double *ZZ = new double[m];
    int *index = new int[n];

    double Rnorm = 0; // Residual norm

    int mode = 0;
    int mda = m;

    nnls_(A_doubles, &mda, &m, &n, B_doubles, X, &Rnorm, W, ZZ, index, &mode);

    if( mode == 2 ){
        std::cerr << "THE DIMENSIONS OF THE PROBLEM ARE BAD. EITHER M .LE. 0 OR N .LE. 0." << std::endl;
        Rnorm = 1e8;
    }
    else if( mode == 3){
        std::cerr << "ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS." << std::endl;
        std::cerr << A << std::endl;
        std::cerr << b << std::endl;
    }

    if( x.size() != static_cast<std::size_t>(n) ){
        x.resize(n);
    }
    // Copy results into x
    std::copy(X, X+n, x.begin());

    delete A_doubles;
    delete B_doubles;
    delete X;
    delete W;
    delete ZZ;
    delete index;

    return Rnorm;
}

} // End namespace Kratos

