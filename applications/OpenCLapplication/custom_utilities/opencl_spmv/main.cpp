#include <iostream>
#include "opencl_interface.h"
#include "includes/matrix_market_interface.h"

#define ROWS_PER_WORKGROUP 16
#define WORKGROUP_SIZE 512

int main()
{
	Kratos::OpenCL::DeviceGroup DeviceGroup(CL_DEVICE_TYPE_GPU, true);

	DeviceGroup.AddCLSearchPath("/home/mossaiby/kratos/applications/OpenCLapplication/custom_utilities");
	cl_uint Program = DeviceGroup.BuildProgramFromFile("opencl_spmv.cl", "-cl-fast-relaxed-math");
	cl_uint Kernel = DeviceGroup.RegisterKernel(Program, "CSR_Matrix_Vector_Multiply", WORKGROUP_SIZE);

	size_t WorkgroupSize = DeviceGroup.WorkGroupSizes[Kernel][0];

	std::cout << "Kernel workgroup size: " << WorkgroupSize << std::endl;

	Kratos::CompressedMatrix A;
	Kratos::Vector X, Y1, Y2;

	Kratos::ReadMatrixMarketMatrix("/home/mossaiby/kratos/applications/OpenCLapplication/custom_utilities/opencl_spmv/A_0.mm", A);
	Kratos::ReadMatrixMarketVector("/home/mossaiby/kratos/applications/OpenCLapplication/custom_utilities/opencl_spmv/B_0.mm", X);

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
	DeviceGroup.SetLocalMemAsKernelArg(Kernel, 7, WorkgroupSize * sizeof(cl_double));

	DeviceGroup.ExecuteKernel(Kernel, A.size1() * WorkgroupSize / ROWS_PER_WORKGROUP + 1);

	DeviceGroup.CopyBuffer(Y_Values, Kratos::OpenCL::DeviceToHost, Kratos::OpenCL::VoidPList(1, &Y1[0]));

	prod(A, X, Y2);

	for (cl_uint i = 0; i < A.size1(); i++)
	{
		if (fabs(Y1[i] - Y2[i]) > 1e-10)
		{
			std::cout << "Error in location " << i << ": " << Y1[i] << "  " << Y2[i] << std::endl;
		}
	}

	/*for (int i = 0; i < 10; i++)
	{
		std::cout << A.index1_data()[i] << std::endl;
	}*/

	std::cout << "Test finished." << std::endl;

	return 0;
}
