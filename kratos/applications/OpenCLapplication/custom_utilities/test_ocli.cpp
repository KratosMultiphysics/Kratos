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

	Kratos::OpenCL::OpenCLDeviceGroup DeviceGroup = Manager.CreateDeviceGroup(CL_DEVICE_TYPE_CPU);

	return 0;
}
