#include "../opencl_linear_solver.h"
#include "includes/matrix_market_interface.h"

//#include "viennacl/compressed_matrix.hpp"
//#include "viennacl/linalg/cg.hpp"

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

	Kratos::CompressedMatrix TA;
	Kratos::Vector B;

	std::cout << "Reading matrix and vector from file..." << std::endl;

	Kratos::ReadMatrixMarketMatrix(argv[1], TA);
	Kratos::ReadMatrixMarketVector(argv[2], B);

	if (TA.size1() != TA.size2() || TA.size1() != B.size())
	{
		std::cout << "Inconsistent matrix / vector sizes." << std::endl;
		return -1;
	}

	cl_uint Size = TA.size1();

	std::cout << "Copying data..." << std::endl;

	boost::numeric::ublas::compressed_matrix <double, boost::numeric::ublas::row_major, 0, boost::numeric::ublas::unbounded_array<cl_uint>, boost::numeric::ublas::unbounded_array<double> > A(TA);

//	viennacl::compressed_matrix <double> VA;
//	viennacl::vector <double> VB(Size);

//	viennacl::copy(A, VA);
//	viennacl::copy(B, VB);

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

	std::cout << "Optimizing..." << std::endl;

	OptimizationParameters.OptimizeInnerProd(B_Values, X_Values, T_Values);
	OptimizationParameters.OptimizeSpMV(A_RowIndices, A_ColumnIndices, A_Values, B_Values, X_Values);

	std::cout << "Solving..." << std::endl;

	Kratos::OpenCL::CGSolver LinearSolver(DeviceGroup, OptimizationParameters, Size, 1000, 1.00e-10);

	T0 = Kratos::OpenCL::Timer();

	bool Result = LinearSolver.Solve(A_RowIndices, A_ColumnIndices, A_Values, B_Values, X_Values);

	T1 = Kratos::OpenCL::Timer() - T0;

	if (Result)
	{
		std::cout << "Solver converged! Tolerance achieved: " << LinearSolver.GetAchievedTolerance() << ", no. of iterations: " << LinearSolver.GetIterationNo() << std::endl;
	}
	else
	{
		std::cout << "Solver diverged! Tolerance achieved: " << LinearSolver.GetAchievedTolerance() << ", no. of iterations: " << LinearSolver.GetIterationNo() << std::endl;
	}

	std::cout << "Solution took " << double(T1) / 1000000000 << "s." << std::endl;

	Kratos::Vector X(Size);

	DeviceGroup.CopyBuffer(X_Values, Kratos::OpenCL::DeviceToHost, Kratos::OpenCL::VoidPList(1, &X[0]));

/*	std::cout << "Solving the system using ViennaCL..." << std::endl;

	viennacl::linalg::cg_tag VCLSolver(1.00e-10, 1000);

	T0 = Kratos::OpenCL::Timer();

	viennacl::linalg::solve(VA, VB, VCLSolver);

	T1 = Kratos::OpenCL::Timer() - T0;

	std::cout << "Iterations: " << VCLSolver.iters() << ", error: " << VCLSolver.error() << std::endl;

	std::cout << "Solution took " << double(T1) / 1000000000 << "s." << std::endl;
*/
	return 0;
}
