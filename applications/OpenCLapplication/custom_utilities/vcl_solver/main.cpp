#include "../opencl_linear_solver.h"

#include "includes/matrix_market_interface.h"
#include "includes/ublas_interface.h"

#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/cg.hpp"

int main(int argc, char *argv[])
{
    int64_t T0, T1;

    if (argc < 3)
    {
        std::cout << "Not enough command line parameters." << std::endl;
        return -1;
    }

    boost::numeric::ublas::compressed_matrix <double, boost::numeric::ublas::row_major, 0, boost::numeric::ublas::unbounded_array <cl_uint>, boost::numeric::ublas::unbounded_array <double> > A;
    boost::numeric::ublas::vector <double> B;

    std::cout << "Reading matrix and vector from file..." << std::endl;

    Kratos::ReadMatrixMarketMatrix(argv[1], A);
    Kratos::ReadMatrixMarketVector(argv[2], B);

    if (A.size1() != A.size2() || A.size1() != B.size())
    {
        std::cout << "Inconsistent matrix / vector sizes." << std::endl;
        return -1;
    }

    cl_uint Size = A.size1();

    std::cout << "Copying data..." << std::endl;

    viennacl::compressed_matrix <double> VA;
    viennacl::vector <double> VB(Size);

    viennacl::copy(A, VA);
    viennacl::copy(B, VB);

    std::cout << "Solving..." << std::endl;

#ifndef MAXIT

    #define MAXIT 1000

#endif

#ifndef MAXERR

    #define MAXERR 1e-10

#endif

    viennacl::linalg::cg_tag VCLSolver(MAXERR, MAXIT);

    T0 = Kratos::OpenCL::Timer();

    viennacl::linalg::solve(VA, VB, VCLSolver);

    T1 = Kratos::OpenCL::Timer() - T0;

    std::cout << "Iterations: " << VCLSolver.iters() << ", error: " << VCLSolver.error() << std::endl;

    std::cout << "Solution took " << double(T1) / 1000000000 << "s." << std::endl;

    return 0;
}
