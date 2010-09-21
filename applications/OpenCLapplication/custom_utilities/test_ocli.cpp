#include "opencl_interface.h"

cl_int Test()
{
	return CL_INVALID_EVENT;
}

int main(int argc, char *argv[])
{
	//KRATOS_OCL_CHECK(CL_INVALID_EVENT);
	//KRATOS_OCL_CHECKED_EXPRESSION(Test());
	Kratos::OpenCL::OpenCLManager Manager;
	std::cout << Manager.DebugData << std::endl;

	Kratos::OpenCL::OpenCLDeviceGroup DeviceGroup = Manager.CreateDeviceGroup(CL_DEVICE_TYPE_ALL);
	DeviceGroup.BuildProgramFromSource("test_ocli.cl");
	DeviceGroup.RegisterKernel("Test");
	DeviceGroup.SetKernelArg(0, 0, 1.00);
	DeviceGroup.SetKernelArg(0, 1, 2.00);

	return 0;
}
