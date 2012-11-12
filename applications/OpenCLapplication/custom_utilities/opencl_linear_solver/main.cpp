#include "../opencl_linear_solver.h"
#include "includes/matrix_market_interface.h"
#include "includes/ublas_interface.h"

int main(int argc, char *argv[])
{
    int64_t T0, T1;

    if (argc < 3)
    {
        std::cout << "Not enough command line parameters." << std::endl;
        return -1;
    }

    std::cout << "Searching for OpenCL devices..." << std::endl;

    Kratos::OpenCL::DeviceGroup DeviceGroup(CL_DEVICE_TYPE_ALL, true);
    DeviceGroup.AddCLSearchPath("..");

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

    cl_uint A_RowIndices = DeviceGroup.CreateBuffer((Size + 1) * sizeof(cl_uint), CL_MEM_READ_ONLY);
    cl_uint A_ColumnIndices = DeviceGroup.CreateBuffer(A.nnz() * sizeof(cl_uint), CL_MEM_READ_ONLY);
    cl_uint A_Values = DeviceGroup.CreateBuffer(A.nnz() * sizeof(cl_double), CL_MEM_READ_ONLY);

    cl_uint B_Values = DeviceGroup.CreateBuffer(Size * sizeof(cl_double), CL_MEM_READ_WRITE);
    cl_uint X_Values = DeviceGroup.CreateBuffer(Size * sizeof(cl_double), CL_MEM_READ_WRITE);
    cl_uint T_Values = DeviceGroup.CreateBuffer(Size * sizeof(cl_double), CL_MEM_READ_WRITE);

    DeviceGroup.CopyBuffer(A_RowIndices, Kratos::OpenCL::HostToDevice, Kratos::OpenCL::VoidPList(1, &A.index1_data()[0]));
    DeviceGroup.CopyBuffer(A_ColumnIndices, Kratos::OpenCL::HostToDevice, Kratos::OpenCL::VoidPList(1, &A.index2_data()[0]));
    DeviceGroup.CopyBuffer(A_Values, Kratos::OpenCL::HostToDevice, Kratos::OpenCL::VoidPList(1, &A.value_data()[0]));

    DeviceGroup.CopyBuffer(B_Values, Kratos::OpenCL::HostToDevice, Kratos::OpenCL::VoidPList(1, &B[0]));
    DeviceGroup.CopyBuffer(X_Values, Kratos::OpenCL::HostToDevice, Kratos::OpenCL::VoidPList(1, &B[0]));
    DeviceGroup.CopyBuffer(T_Values, Kratos::OpenCL::HostToDevice, Kratos::OpenCL::VoidPList(1, &B[0]));

    Kratos::OpenCL::LinearSolverOptimizationParameters OptimizationParameters(DeviceGroup, Size);

    OptimizationParameters.OptimizeInnerProd(B_Values, X_Values, T_Values);
    OptimizationParameters.OptimizeSpMV(A_RowIndices, A_ColumnIndices, A_Values, B_Values, X_Values);

#ifndef MAXIT

    #define MAXIT 1000

#endif

#ifndef MAXERR

    #define MAXERR 1e-10

#endif

#ifndef CGTYPE

    #define CGTYPE 1

#endif

#if CGTYPE == 1

    std::cout << "Using original form of CG..." << std::endl;
    Kratos::OpenCL::CGSolverOriginal LinearSolver(DeviceGroup, OptimizationParameters, Size, MAXIT, MAXERR);

#elif CGTYPE == 2

    std::cout << "Using three term recurrence form of CG..." << std::endl;
    Kratos::OpenCL::CGSolverThreeTermRecurrence LinearSolver(DeviceGroup, OptimizationParameters, Size, MAXIT, MAXERR);

#elif CGTYPE == 3

    std::cout << "Using Chronopoulos form of CG..." << std::endl;
    Kratos::OpenCL::CGSolverChronopoulos LinearSolver(DeviceGroup, OptimizationParameters, Size, MAXIT, MAXERR);

#endif

    std::cout << "Solving..." << std::endl;

    T0 = Kratos::OpenCL::Timer();

    bool Result = LinearSolver.Solve(A_RowIndices, A_ColumnIndices, A_Values, B_Values, X_Values);

    T1 = Kratos::OpenCL::Timer() - T0;

    if (Result)
    {
        std::cout << "Solver converged! Est. error: " << LinearSolver.GetEstimatedError() << ", no. of iterations: " << LinearSolver.GetIterationNo() << std::endl;
    }
    else
    {
        std::cout << "Solver diverged! Est. error: " << LinearSolver.GetEstimatedError() << ", no. of iterations: " << LinearSolver.GetIterationNo() << std::endl;
    }

    std::cout << "Solution took " << double(T1) / 1000000000 << "s." << std::endl;

    Kratos::Vector X(Size);

    DeviceGroup.CopyBuffer(X_Values, Kratos::OpenCL::DeviceToHost, Kratos::OpenCL::VoidPList(1, &X[0]));

    return 0;
}
