/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Farshid Mossaiby
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
// opencl_interface.h
//
// Kratos OpenCL interface

#if !defined(KRATOS_OPENCL_INTERFACE_H_INCLUDED)
#define KRATOS_OPENCL_INTERFACE_H_INCLUDED

// OpenCL include path is different on Apple

#if defined(__APPLE__) || defined(__MACOSX)
	#include <OpenCL/cl.h>
#else
	#include <CL/cl.h>
#endif

// System includes

#include <iostream>
#include <sstream>
#include <vector>


// External includes 


// Project includes


namespace Kratos
{

namespace OpenCL
{

// Useful types

	typedef std::vector<cl_device_id> DeviceIDList;
	typedef std::vector<cl_command_queue> CommandQueueList;

//
// KRATOS_OCL_CHECK
//
// Used to check an OpenCL return code and abort with some debugging information

#define KRATOS_OCL_CHECK(Code)				Kratos::OpenCL::CheckError(Code, __FILE__, __FUNCTION__, __LINE__);

//
// KRATOS_OCL_CHECKED_EXPRESSION
//
// Used to evaluate an expression and check the OpenCL return code and abort with debugging information

#define KRATOS_OCL_CHECKED_EXPRESSION(Expr)	{ cl_int Err = Expr; Kratos::OpenCL::CheckError(Err, __FILE__, __FUNCTION__, __LINE__, " while executing " #Expr); }

//
// KRATOS_OCL_DEFINE_ERROR_CODE
//
// Used to define various error codes in ErrorString() function

#define KRATOS_OCL_DEFINE_ERROR_CODE(Code)	case Code: return #Code

//
// ErrorString
//
// Returns a string representation of an OpenCL error

	const char *ErrorString(cl_int Code)
	{
		switch(Code)
		{
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_SUCCESS);
		
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_DEVICE_NOT_FOUND);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_DEVICE_NOT_AVAILABLE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_COMPILER_NOT_AVAILABLE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_MEM_OBJECT_ALLOCATION_FAILURE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_OUT_OF_RESOURCES);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_OUT_OF_HOST_MEMORY);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_PROFILING_INFO_NOT_AVAILABLE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_MEM_COPY_OVERLAP);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_IMAGE_FORMAT_MISMATCH);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_IMAGE_FORMAT_NOT_SUPPORTED);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_BUILD_PROGRAM_FAILURE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_MAP_FAILURE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_MISALIGNED_SUB_BUFFER_OFFSET);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST);

			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_VALUE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_DEVICE_TYPE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_PLATFORM);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_DEVICE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_CONTEXT);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_QUEUE_PROPERTIES);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_COMMAND_QUEUE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_HOST_PTR);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_MEM_OBJECT);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_IMAGE_SIZE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_SAMPLER);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_BINARY);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_BUILD_OPTIONS);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_PROGRAM);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_PROGRAM_EXECUTABLE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_KERNEL_NAME);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_KERNEL_DEFINITION);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_KERNEL);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_ARG_INDEX);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_ARG_VALUE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_ARG_SIZE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_KERNEL_ARGS);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_WORK_DIMENSION);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_WORK_GROUP_SIZE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_WORK_ITEM_SIZE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_GLOBAL_OFFSET);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_EVENT_WAIT_LIST);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_EVENT);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_OPERATION);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_GL_OBJECT);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_BUFFER_SIZE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_MIP_LEVEL);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_GLOBAL_WORK_SIZE);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_PROPERTY);

			// Default case
			default: return "unknown";
		}	
	}

//
// RaiseError
//
// Used to raise an OpenCL error and abort
// Do not use directly; use the KRATOS_OCL_CHECK or KRATOS_OCL_CHECKED_EXPRESSION macros

	void RaiseError(cl_int Code, const char *FileName, const char *Function, unsigned int Line, const char *Expression = "")
	{
		std::cout
			<< "OpenCL reported " << ErrorString(Code) << " error"
			<< Expression
			<< " in function " << Function << ", file " << FileName << ", line " << Line << "." << std::endl
			<< "Forcing a SEGFAULT now; you may use this to find the problem in a debugger." << std::endl;

		// This will cause a SEGFAULT, so we can find it easily with a debugger
		*((int *) 0) = 0;
	}		

//
// CheckError
//
// Used to check an OpenCL return code and abort
// Do not use directly; use the KRATOS_OCL_CHECK or KRATOS_OCL_CHECKED_EXPRESSION macros

	inline void CheckError(cl_int Code, const char *FileName, const char *Function, unsigned int Line, const char *Expression = "")
	{
		if (Code != CL_SUCCESS)
			RaiseError(Code, FileName, Function, Line, Expression);
	}

// Forward declaration

	class OpenCLManager;

//
// OpenCLDeviceGroup
//
// A class to represent a group of OpenCL devices

	class OpenCLDeviceGroup
	{
		public:

			OpenCLManager &Manager;

			DeviceIDList DeviceIDs;
			CommandQueueList CommandQueues;
			
			cl_uint DeviceNo;
			cl_device_type DeviceType;
			cl_context Context;
			cl_command_queue CommandQueue;
			
			// Constructors
			OpenCLDeviceGroup(OpenCLManager &_Manager, cl_device_type _DeviceType): Manager(_Manager), DeviceType(_DeviceType)
			{
				Init();
			}

			OpenCLDeviceGroup(OpenCLManager &_Manager, DeviceIDList &_DeviceIDs): Manager(_Manager), DeviceIDs(_DeviceIDs)
			{
				Init();
			}

			// Destructor
			~OpenCLDeviceGroup()
			{
				cl_int Err;
				
				// Releasing OpenCL objects
				Err = clReleaseContext(Context);
				KRATOS_OCL_CHECK(Err);
				
				for (int i = 0; i < DeviceNo; i++)
				{
					Err = clReleaseCommandQueue(CommandQueues[i]);
					KRATOS_OCL_CHECK(Err);
				}
			}

		private:

			// Initialize OpenCL structures
			void Init();
	};

//
// OpenCLManager
//
// A helper class for OpenCL

	class OpenCLManager
	{
		public:
		
			// Used platform ID
			cl_platform_id PlatformID;

			// Debug data
			std::string DebugData;

			// Constructor
			OpenCLManager(const char *_PlatformVendor = "")
			{
				cl_int Err;

				cl_uint PlatformNo, DeviceNo;
				cl_platform_id *Platforms;

				std::string PlatformVendor(_PlatformVendor);
				std::ostringstream DebugDataStream;
				char CharData[1024];

				DebugDataStream
					<< "OpenCL debug data" << std::endl
					<< std::endl;
				
				// Query available platforms
				Err = clGetPlatformIDs(0, NULL, &PlatformNo);
				KRATOS_OCL_CHECK(Err);

				DebugDataStream
					<< "No. of platforms: " << PlatformNo << std::endl
					<< std::endl;

				Platforms = new cl_platform_id[PlatformNo];

				Err = clGetPlatformIDs(PlatformNo, Platforms, NULL);
				KRATOS_OCL_CHECK(Err);

				for (int i = 0; i < PlatformNo; i++)
				{
					DebugDataStream << "Platform " << i << ":" << std::endl;

					// Get some information about this platform; when appropriate check if this is the requested platform
					Err = clGetPlatformInfo(Platforms[i], CL_PLATFORM_NAME, sizeof(CharData), (void *) CharData, NULL);
					KRATOS_OCL_CHECK(Err);
					DebugDataStream << "  Name:             " << CharData << std::endl;

					Err = clGetPlatformInfo(Platforms[i], CL_PLATFORM_VENDOR, sizeof(CharData), (void *) CharData, NULL);
					KRATOS_OCL_CHECK(Err);
					DebugDataStream << "  Vendor:           " << CharData << std::endl;

					// Check if this is the requested platform
					if (PlatformVendor.compare(CharData) == 0)
						PlatformID = Platforms[i];

					Err = clGetPlatformInfo(Platforms[i], CL_PLATFORM_VERSION, sizeof(CharData), (void *) CharData, NULL);
					KRATOS_OCL_CHECK(Err);
					DebugDataStream << "  Version:          " << CharData << std::endl;
				
					Err = clGetPlatformInfo(Platforms[i], CL_PLATFORM_EXTENSIONS, sizeof(CharData), (void *) CharData, NULL);
					KRATOS_OCL_CHECK(Err);
					DebugDataStream << "  Extensions:       " << CharData << std::endl;

					Err = clGetPlatformInfo(Platforms[i], CL_PLATFORM_PROFILE, sizeof(CharData), (void *) CharData, NULL);
					KRATOS_OCL_CHECK(Err);
					DebugDataStream << "  Profile:          " << CharData << std::endl;
				
					// Query available devices on this platform
					Err = clGetDeviceIDs(Platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &DeviceNo);
					KRATOS_OCL_CHECK(Err);
					DebugDataStream << "  No. of devices:   " << DeviceNo << std::endl;

					// TODO: Collect debug data for each device here
				}

				// Store data
				DebugData = DebugDataStream.str();

				// If no specific platform is requested, default to first one
				if (PlatformVendor.compare("") == 0)
					PlatformID = Platforms[0];

				delete [] Platforms;
			}

			// Destructor
			~OpenCLManager()
			{
				// Nothing to do!
			}

			// Creates a device group by type
			OpenCLDeviceGroup CreateDeviceGroup(cl_device_type DeviceType)
			{
				// TODO: correct this!
				return OpenCLDeviceGroup(*this, DeviceType);
			}

			// Creates a device group by explicit device ids
			OpenCLDeviceGroup CreateDeviceGroup(DeviceIDList &DeviceIDs)
			{
				// TODO: correct this!
				return OpenCLDeviceGroup(*this, DeviceIDs);
			}

		private:
			
	};

void OpenCLDeviceGroup::Init()
{
	cl_int Err;
	cl_context_properties Properties[] = {CL_CONTEXT_PLATFORM, (cl_context_properties) Manager.PlatformID, 0};

	// Create context; try to find out which constructor has been called
	if (DeviceIDs.size() != 0)
	{
		// An explicit list of IDs are specified
		Context = clCreateContext(Properties, DeviceIDs.size(), DeviceIDs.data(), NULL, NULL, &Err);
		KRATOS_OCL_CHECK(Err);

		DeviceNo = DeviceIDs.size();
	}
	else
	{
		// A device type is specified; get all devices of that type
		Err = clGetDeviceIDs(Manager.PlatformID, DeviceType, 0, NULL, &DeviceNo);
		KRATOS_OCL_CHECK(Err);

		DeviceIDs.resize(DeviceNo);
		Err = clGetDeviceIDs(Manager.PlatformID, DeviceType, DeviceNo, DeviceIDs.data(), NULL);
		KRATOS_OCL_CHECK(Err);

		Context = clCreateContextFromType(Properties, DeviceType, NULL, NULL, &Err);
		KRATOS_OCL_CHECK(Err);
	}

	// Create command queues one per device
	CommandQueues.resize(DeviceNo);
	
	for (int i = 0; i < DeviceNo; i++)
	{
		CommandQueues[i] = clCreateCommandQueue(Context, DeviceIDs[i], 0, &Err);
		KRATOS_OCL_CHECK(Err);
	}
}

}  // namespace OpenCL 

}  // namespace Kratos

#endif // KRATOS_OPENCL_INTERFACE_H_INCLUDED
