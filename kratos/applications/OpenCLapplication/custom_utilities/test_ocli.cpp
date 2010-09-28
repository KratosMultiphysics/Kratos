#include "opencl_interface.h"
#include <cmath>

#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#else
#include <ctime>
#endif

const cl_uint DataSize = 360000;
const double Eps = 1e-10;

double Input[DataSize], Output[DataSize];

int main(int argc, char *argv[])
{
	std::cout << "Initializing OpenCL device group..." << std::endl;
	Kratos::OpenCL::DeviceGroup OCLDeviceGroup(CL_DEVICE_TYPE_CPU, false);  // Try to find all available devices

	std::cout << "Found " << OCLDeviceGroup.DeviceNo << " device(s)." << std::endl;
	for (cl_uint i = 0; i < OCLDeviceGroup.DeviceNo; i++)
	{
		std::cout << "  Device " << i << ": " << Kratos::OpenCL::DeviceTypeString(OCLDeviceGroup.DeviceTypes[i]) << std::endl;
	}

	std::cout << "Loading kernel(s) from file test_ocli.cl..." << std::endl;
	OCLDeviceGroup.BuildProgramFromFile("test_ocli.cl", "-cl-unsafe-math-optimizations");

	std::cout << "Registering kernel Test()..." << std::endl;
	cl_uint TestKernel = OCLDeviceGroup.RegisterKernel("Test");  // This will return 0, but we do not want to memorize it ourselves!

	for (cl_uint i = 0; i < DataSize; i++)
	{
		Input[i] = i;
		Output[i] = -1.00;
	}

	std::cout << "Creating the buffers..." << std::endl;
	cl_uint InputBuffer = OCLDeviceGroup.CreateBuffer(DataSize / OCLDeviceGroup.DeviceNo * sizeof(double), CL_MEM_READ_ONLY);  // This will return 0, but we do not want to memorize it ourselves!
	cl_uint OutputBuffer = OCLDeviceGroup.CreateBuffer(DataSize / OCLDeviceGroup.DeviceNo * sizeof(double), CL_MEM_WRITE_ONLY);  // This will return 1, but we do not want to memorize it ourselves!

	std::cout << "Copying the data to the device(s)..." << std::endl;
	for (cl_uint i = 0; i < OCLDeviceGroup.DeviceNo; i++)
	{
		OCLDeviceGroup.CopyBuffer(i, InputBuffer, Kratos::OpenCL::HostToDevice, Input + i * DataSize / OCLDeviceGroup.DeviceNo);
	}

	std::cout << "Setting kernel arguments..." << std::endl;
	OCLDeviceGroup.SetBufferAsKernelArg(TestKernel, 0, InputBuffer);
	OCLDeviceGroup.SetBufferAsKernelArg(TestKernel, 1, OutputBuffer);
	OCLDeviceGroup.SetKernelArg(TestKernel, 2, 10.00);

	std::cout << "Executing kernel..." << std::endl;
	OCLDeviceGroup.ExecuteKernel(TestKernel, DataSize / OCLDeviceGroup.DeviceNo);

	std::cout << "Copying the data from the device(s)..." << std::endl;
	for (cl_uint i = 0; i < OCLDeviceGroup.DeviceNo; i++)
	{
		OCLDeviceGroup.CopyBuffer(i, OutputBuffer, Kratos::OpenCL::DeviceToHost, Output + i * DataSize / OCLDeviceGroup.DeviceNo);
	}

	std::cout << "Verification of results..." << std::endl;
	for (cl_uint i = 0; i < DataSize; i++)
	{
		if (Output[i] - (2.00 * sin(i) * cos(i) + 10.00) > Eps)
		{
			std::cout << "Error at position " << i << std::endl;
		}
	}
	std::cout << "Verification finished." << std::endl;
	
	#ifndef _OPENMP
	  double start_time = clock();
	#else
	  double start_time = omp_get_wtime();
	#endif
	
	#pragma omp parallel for
	for(int i =0; i<DataSize;i++)
	  Output[i] = -1.0;
	
	#pragma omp parallel for
	for(int i =0; i<DataSize;i++)
	  Output[i] -= (2.00 * sin(i) * cos(i) + 10.00);

	
	#ifndef _OPENMP
	  double stop_time = clock();
	  std::cout << "no openmp timing " <<  (stop_time - start_time)/CLOCKS_PER_SEC <<  std::endl;
	#else
	  double stop_time = omp_get_wtime();
	  std::cout << "with openmp " <<  (stop_time - start_time) <<  std::endl;
	#endif
	
	

	return 0;
}
