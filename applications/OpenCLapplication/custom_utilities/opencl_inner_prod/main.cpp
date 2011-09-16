#include <iostream>
#include "opencl_interface.h"
#include "includes/matrix_market_interface.h"

#define WORKGROUP_SIZE_BITS 8
#define WORKGROUP_SIZE (1 << WORKGROUP_SIZE_BITS)

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
	cl_uint Program = DeviceGroup.BuildProgramFromFile("opencl_inner_prod.cl", "-cl-fast-relaxed-math");
	cl_uint Kernel = DeviceGroup.RegisterKernel(Program, "Vector_Vector_Multiply", WORKGROUP_SIZE);

	Kratos::Vector X, Y, T;
	cl_double Z1, Z2;

	Kratos::ReadMatrixMarketVector(argv[1], X);
	Kratos::ReadMatrixMarketVector(argv[2], Y);

	T.resize((X.size() + WORKGROUP_SIZE - 1) / WORKGROUP_SIZE);

	cl_uint X_Values = DeviceGroup.CreateBuffer(X.size() * sizeof(cl_double), CL_MEM_READ_ONLY);
	cl_uint Y_Values = DeviceGroup.CreateBuffer(X.size() * sizeof(cl_double), CL_MEM_READ_ONLY);
	cl_uint T_Values = DeviceGroup.CreateBuffer(((X.size() + WORKGROUP_SIZE - 1) / WORKGROUP_SIZE) * sizeof(cl_double), CL_MEM_WRITE_ONLY);

	DeviceGroup.CopyBuffer(X_Values, Kratos::OpenCL::HostToDevice, Kratos::OpenCL::VoidPList(1, &X[0]));
	DeviceGroup.CopyBuffer(Y_Values, Kratos::OpenCL::HostToDevice, Kratos::OpenCL::VoidPList(1, &Y[0]));

	DeviceGroup.SetBufferAsKernelArg(Kernel, 0, X_Values);
	DeviceGroup.SetBufferAsKernelArg(Kernel, 1, Y_Values);
	DeviceGroup.SetBufferAsKernelArg(Kernel, 2, T_Values);
	DeviceGroup.SetKernelArg(Kernel, 3, X.size());
	DeviceGroup.SetLocalMemAsKernelArg(Kernel, 4, WORKGROUP_SIZE * sizeof(cl_double));


	for (unsigned int i = 0; i < N; i++)
	{
		t1 = Timer();

		DeviceGroup.ExecuteKernel(Kernel, X.size());
		DeviceGroup.CopyBuffer(T_Values, Kratos::OpenCL::DeviceToHost, Kratos::OpenCL::VoidPList(1, &T[0]));

		Z1 = 0.00;
		for (unsigned int j = 0; j < T.size(); j++)
		{
			Z1 += T[j];
		}

		t2 = Timer();

		if (i == 0 || t2 - t1 < T1)
		{
			T1 = t2 - t1;
		}
	}


	for (unsigned int i = 0; i < N; i++)
	{
		t1 = Timer();

		Z2 = inner_prod(X, Y);

		t2 = Timer();

		if (i == 0 || t2 - t1 < T2)
		{
			T2 = t2 - t1;
		}
	}

	if (fabs(Z1 - Z2) > 1e-10)
	{
		std::cout << "Error: " << Z1 << "  " << Z2 << std::endl;
	}

	std::cout << "Z1 is " << Z1 << "." << std::endl;
	std::cout << "Z2 is " << Z2 << "." << std::endl;


	std::cout << "Test finished." << std::endl << "OpenCL inner_prod:\t" << T1 / 1000000.00 << " ms" << std::endl << "uBlas:\t\t\t" << T2 / 1000000.00 << " ms" << std::endl;

	return 0;
}
