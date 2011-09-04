#include <iostream>
#include "opencl_interface.h"
#include "includes/matrix_market_interface.h"

#define ROWS_PER_WORKGROUP_BITS 5
#define ROWS_PER_WORKGROUP (1 << ROWS_PER_WORKGROUP_BITS)

#define WORKGROUP_SIZE_BITS 7
#define WORKGROUP_SIZE (1 << WORKGROUP_SIZE_BITS)

#define LOCAL_WORKGROUP_SIZE_BITS (WORKGROUP_SIZE_BITS - ROWS_PER_WORKGROUP_BITS)
#define LOCAL_WORKGROUP_SIZE (1 << LOCAL_WORKGROUP_SIZE_BITS)

#define N 50

int64_t Timer()
{
	struct timespec tp;

	clock_gettime(CLOCK_MONOTONIC, &tp);

	return (unsigned long long) tp.tv_sec * (1000ULL * 1000ULL * 1000ULL) + (unsigned long long) tp.tv_nsec;
}

int main(int argc, char *argv[])
{
	int64_t t1, t2, T1, T2;

	Kratos::OpenCL::DeviceGroup DeviceGroup(CL_DEVICE_TYPE_GPU, true);

	DeviceGroup.AddCLSearchPath("/home/mossaiby/kratos/applications/OpenCLapplication/custom_utilities");
	cl_uint Program = DeviceGroup.BuildProgramFromFile("opencl_spmv.cl", "-cl-fast-relaxed-math");
	cl_uint Kernel = DeviceGroup.RegisterKernel(Program, "CSR_Matrix_Vector_Multiply", WORKGROUP_SIZE);

	Kratos::CompressedMatrix A;
	Kratos::Vector X, Y1, Y2;

	Kratos::ReadMatrixMarketMatrix(argv[1], A);
	Kratos::ReadMatrixMarketVector(argv[2], X);

	Y1.resize(A.size1());
	Y2.resize(A.size1());

	cl_uint A_RowIndices = DeviceGroup.CreateBuffer((A.size1() + 1) * sizeof(cl_ulong), CL_MEM_READ_ONLY);
	cl_uint A_ColumnIndices = DeviceGroup.CreateBuffer(A.nnz() * sizeof(cl_ulong), CL_MEM_READ_ONLY);
	cl_uint A_Values = DeviceGroup.CreateBuffer(A.nnz() * sizeof(cl_double), CL_MEM_READ_ONLY);
	cl_uint X_Values = DeviceGroup.CreateBuffer(A.size1() * sizeof(cl_double), CL_MEM_READ_ONLY);
	cl_uint Y_Values = DeviceGroup.CreateBuffer(A.size1() * sizeof(cl_double), CL_MEM_WRITE_ONLY);

	DeviceGroup.CopyBuffer(A_RowIndices, Kratos::OpenCL::HostToDevice, Kratos::OpenCL::VoidPList(1, &A.index1_data()[0]));
	DeviceGroup.CopyBuffer(A_ColumnIndices, Kratos::OpenCL::HostToDevice, Kratos::OpenCL::VoidPList(1, &A.index2_data()[0]));
	DeviceGroup.CopyBuffer(A_Values, Kratos::OpenCL::HostToDevice, Kratos::OpenCL::VoidPList(1, &A.value_data()[0]));

	DeviceGroup.CopyBuffer(X_Values, Kratos::OpenCL::HostToDevice, Kratos::OpenCL::VoidPList(1, &X[0]));

	DeviceGroup.SetBufferAsKernelArg(Kernel, 0, A_RowIndices);
	DeviceGroup.SetBufferAsKernelArg(Kernel, 1, A_ColumnIndices);
	DeviceGroup.SetBufferAsKernelArg(Kernel, 2, A_Values);
	DeviceGroup.SetBufferAsKernelArg(Kernel, 3, X_Values);
	DeviceGroup.SetBufferAsKernelArg(Kernel, 4, Y_Values);
	DeviceGroup.SetKernelArg(Kernel, 5, A.size1());
	DeviceGroup.SetLocalMemAsKernelArg(Kernel, 6, (ROWS_PER_WORKGROUP + 1) * sizeof(cl_ulong));
	DeviceGroup.SetLocalMemAsKernelArg(Kernel, 7, WORKGROUP_SIZE * sizeof(cl_double));

	for (unsigned int i = 0; i < N; i++)
	{
		t1 = Timer();

		DeviceGroup.ExecuteKernel(Kernel, A.size1() * LOCAL_WORKGROUP_SIZE + 1);

		t2 = Timer();

		if (i == 0 || t2 - t1 < T1)
		{
			T1 = t2 - t1;
		}
	}

	DeviceGroup.CopyBuffer(Y_Values, Kratos::OpenCL::DeviceToHost, Kratos::OpenCL::VoidPList(1, &Y1[0]));


	for (unsigned int i = 0; i < N; i++)
	{
		t1 = Timer();

		axpy_prod(A, X, Y2);

		t2 = Timer();

		if (i == 0 || t2 - t1 < T2)
		{
			T2 = t2 - t1;
		}
	}

	for (cl_uint i = 0; i < A.size1(); i++)
	{
		if (fabs(Y1[i] - Y2[i]) > 1e-10)
		{
			std::cout << "Error in location " << i << ": " << Y1[i] << "  " << Y2[i] << std::endl;
		}
	}

	std::cout << "Norm_2 of Y1 is " << norm_2(Y1) << "." << std::endl;
	std::cout << "Norm_2 of Y2 is " << norm_2(Y2) << "." << std::endl;

	std::cout << "Test finished." << std::endl << "OpenCL SpMV:\t" << T1 / 1000000.00 << " ms" << std::endl << "uBlas:\t\t" << T2 / 1000000.00 << " ms" << std::endl;

	return 0;
}
