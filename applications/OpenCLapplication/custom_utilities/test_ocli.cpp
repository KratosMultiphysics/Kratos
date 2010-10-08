#include "opencl_interface.h"
#include <cmath>

#ifdef _OPENMP
	#include <omp.h>
#else
	#include <ctime>
#endif

#ifdef _OPENMP
	#define TIMER()					omp_get_wtime()
#else
	#define TIMER()					clock()
#endif

const cl_uint DataSize = 360000;
const double Eps = 1e-10;
const double Offset = 10.00;

double OCLInput[DataSize], OCLOutput[DataSize], HostInput[DataSize], HostOutput[DataSize];

int main(int argc, char *argv[])
{
	std::cout << "Initializing OpenCL device group..." << std::endl;
	Kratos::OpenCL::DeviceGroup OCLDeviceGroup(CL_DEVICE_TYPE_ALL, false);  // Try to find all available devices

	std::cout << "Found " << OCLDeviceGroup.DeviceNo << " device(s)." << std::endl;
	for (cl_uint i = 0; i < OCLDeviceGroup.DeviceNo; i++)
	{
		std::cout << "  Device " << i << ": " << Kratos::OpenCL::DeviceTypeString(OCLDeviceGroup.DeviceTypes[i]) << std::endl;
	}

	std::cout << "Loading kernels from source files..." << std::endl;
	cl_uint TestProgram1 = OCLDeviceGroup.BuildProgramFromFile("opencl_edge_data.cl", "-cl-unsafe-math-optimizations");  // This will return 0, but we do not want to memorize it ourselves!
	cl_uint TestProgram2 = OCLDeviceGroup.BuildProgramFromFile("opencl_pure_convection_edgebased.cl", "-cl-unsafe-math-optimizations");  // This will return 1, but we do not want to memorize it ourselves!

	std::cout << "Registering kernels..." << std::endl;
	cl_uint TestKernel1 = OCLDeviceGroup.RegisterKernel(TestProgram1, "Test");  // This will return 0, but we do not want to memorize it ourselves!
	cl_uint TestKernel2 = OCLDeviceGroup.RegisterKernel(TestProgram2, "Test");  // This will return 1, but we do not want to memorize it ourselves!

	for (cl_uint i = 0; i < DataSize; i++)
	{
		OCLInput[i] = i;
		OCLOutput[i] = -1.00;
	}

	std::cout << "Creating the buffers..." << std::endl;
	cl_uint InputBuffer = OCLDeviceGroup.CreateBuffer(DataSize / OCLDeviceGroup.DeviceNo * sizeof(double), CL_MEM_READ_ONLY);  // This will return 0, but we do not want to memorize it ourselves!
	cl_uint OutputBuffer = OCLDeviceGroup.CreateBuffer(DataSize / OCLDeviceGroup.DeviceNo * sizeof(double), CL_MEM_WRITE_ONLY);  // This will return 1, but we do not want to memorize it ourselves!

	std::cout << "Copying the data to the device(s)..." << std::endl;
	for (cl_uint i = 0; i < OCLDeviceGroup.DeviceNo; i++)
	{
		OCLDeviceGroup.CopyBuffer(i, InputBuffer, Kratos::OpenCL::HostToDevice, OCLInput + i * DataSize / OCLDeviceGroup.DeviceNo);
	}

	std::cout << "Setting kernel arguments..." << std::endl;
	OCLDeviceGroup.SetBufferAsKernelArg(TestKernel1, 0, InputBuffer);
	OCLDeviceGroup.SetBufferAsKernelArg(TestKernel1, 1, OutputBuffer);
	OCLDeviceGroup.SetKernelArg(TestKernel1, 2, Offset);

	OCLDeviceGroup.SetBufferAsKernelArg(TestKernel2, 0, InputBuffer);
	OCLDeviceGroup.SetBufferAsKernelArg(TestKernel2, 1, OutputBuffer);
	OCLDeviceGroup.SetKernelArg(TestKernel2, 2, Offset);

	double OCLTimer1 = TIMER();

	std::cout << "Executing kernel 1..." << std::endl;
	OCLDeviceGroup.ExecuteKernel(TestKernel1, DataSize / OCLDeviceGroup.DeviceNo);

	OCLTimer1 = TIMER() - OCLTimer1;

	std::cout << "Copying the data from the device(s)..." << std::endl;
	for (cl_uint i = 0; i < OCLDeviceGroup.DeviceNo; i++)
	{
		OCLDeviceGroup.CopyBuffer(i, OutputBuffer, Kratos::OpenCL::DeviceToHost, OCLOutput + i * DataSize / OCLDeviceGroup.DeviceNo);
	}

	std::cout << "Verification of results..." << std::endl;
	for (cl_uint i = 0; i < DataSize; i++)
	{
		if (OCLOutput[i] - (2.00 * sin(i) * cos(i) + Offset) > Eps)
		{
			std::cout << "Error at position " << i << std::endl;
		}
	}

	std::cout << "Verification finished." << std::endl;

	double OCLTimer2 = TIMER();

	std::cout << "Executing kernel 2..." << std::endl;
	OCLDeviceGroup.ExecuteKernel(TestKernel2, DataSize / OCLDeviceGroup.DeviceNo);

	OCLTimer2 = TIMER() - OCLTimer2;

	std::cout << "Copying the data from the device(s)..." << std::endl;
	for (cl_uint i = 0; i < OCLDeviceGroup.DeviceNo; i++)
	{
		OCLDeviceGroup.CopyBuffer(i, OutputBuffer, Kratos::OpenCL::DeviceToHost, OCLOutput + i * DataSize / OCLDeviceGroup.DeviceNo);
	}

	std::cout << "Verification of results..." << std::endl;
	for (cl_uint i = 0; i < DataSize; i++)
	{
		if (OCLOutput[i] - (2.00 * sin(i) * cos(i) + Offset) > Eps)
		{
			std::cout << "Error at position " << i << std::endl;
		}
	}

	std::cout << "Verification finished." << std::endl;

	#pragma omp parallel for
	for(cl_uint i = 0; i < DataSize; i++)
	{
		HostInput[i] = i;
		// HostOutput[i] = -1;
	}

	double HostTimer = TIMER();

	#pragma omp parallel for
	for(cl_uint i = 0; i < DataSize; i++)
		HostOutput[i] = (2.00 * sin(i) * cos(i) + Offset);

	HostTimer = TIMER() - HostTimer;

#ifdef _OPENMP
		std::cout <<
			"Timing with OpenMP" << std::endl <<
			std::endl <<
			"OpenCL 1:" << OCLTimer1 << std::endl <<
			"OpenCL 2:" << OCLTimer2 << std::endl <<
			"Host:    " << HostTimer << std::endl;
#else
		std::cout <<
			"Timing without OpenMP" << std::endl <<
			std::endl <<
			"OpenCL 1:" << OCLTimer1 / CLOCKS_PER_SEC << std::endl <<
			"OpenCL 2:" << OCLTimer2 / CLOCKS_PER_SEC << std::endl <<
			"Host:    " << HostTimer / CLOCKS_PER_SEC << std::endl;
#endif

	return 0;
}
