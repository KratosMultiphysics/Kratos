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

// System includes

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>

//
// OpenCL include path is different on Apple

#if defined(__APPLE__) || defined(__MACOSX)
	#include <OpenCL/cl.h>
#else
	#include <CL/cl.h>
#endif

//
// Stringizer macro

#define KRATOS_OCL_VALUE(x)				#x
#define KRATOS_OCL_STRINGIZE(x)			KRATOS_OCL_VALUE(x)

//
// Try to find out OpenCL version

#ifdef CL_VERSION_1_1
	#define KRATOS_OCL_VERSION			110
#else
	#define KRATOS_OCL_VERSION			100
#endif

//
// OpenCL version string

#define KRATOS_OCL_VERSION_STRING		KRATOS_OCL_STRINGIZE(KRATOS_OCL_VERSION)

//
// Path handling stuff

#if defined(WINDOWS)
	#include <direct.h>
	#define GETCWD _getcwd
	#define PATH_SEPARATOR '\\'
#else
	#include <unistd.h>
	#define GETCWD getcwd
	#define PATH_SEPARATOR '/'
#endif


// External includes


// Project includes


//
// OpenCL 1.0 adjustments

#if KRATOS_OCL_VERSION < 110

// cl_double3 is the same as cl_double4 in OpenCL 1.1 and later

typedef cl_double4 cl_double3;

#endif

namespace Kratos
{

namespace OpenCL
{

// Useful vector types

	typedef enum {HostToDevice, DeviceToHost} CopyDirection;
	typedef struct _ImageDimension
	{
		size_t Sizes[3];

		_ImageDimension(size_t _Size1, size_t _Size2, size_t _Size3 = 1)
		{
			Sizes[0] = _Size1;
			Sizes[1] = _Size2;
			Sizes[2] = _Size3;
		}
	} ImageDimension;

	typedef std::vector<std::string> StringList;

	typedef std::vector<cl_device_type> DeviceTypeList;
	typedef std::vector<cl_device_id> DeviceIDList;
	typedef std::vector<cl_context> ContextList;
	typedef std::vector<cl_command_queue> CommandQueueList;

	typedef std::vector<void *> VoidPList;

	typedef std::vector<cl_program> ProgramList;
	typedef std::vector<ProgramList> ProgramList2D;

	typedef std::vector<cl_mem> MemList;
	typedef std::vector<MemList> MemList2D;

	typedef std::vector<cl_kernel> KernelList;
	typedef std::vector<KernelList> KernelList2D;

	typedef std::vector<size_t> SizeTList;
	typedef std::vector<SizeTList> SizeTList2D;

	typedef std::vector<ImageDimension> ImageDimensionList;
	typedef std::vector<ImageDimensionList> ImageDimensionList2D;

//
// KRATOS_OCL_CPU_WORK_GROUP_SIZE
//
// Size of a work group on CPU devices

#define KRATOS_OCL_CPU_WORK_GROUP_SIZE		2

//
// KRATOS_OCL_PROGRAM_BUILD_LOG_SIZE
//
// Size of program build log buffer

#define KRATOS_OCL_PROGRAM_BUILD_LOG_SIZE	102400

//
// KRATOS_OCL_PLATFORM_DATA_SIZE
//
// Size of platform vendor name buffer

#define KRATOS_OCL_PLATFORM_DATA_SIZE		256

//
// KRATOS_OCL_CHECK
//
// Used to check an OpenCL return code and abort with some debugging information

#define KRATOS_OCL_CHECK(Code)				Kratos::OpenCL::CheckError(Code, __FILE__, __FUNCTION__, __LINE__, true);

//
// KRATOS_OCL_WARN
//
// Used to check an OpenCL return code and print some debugging information

#define KRATOS_OCL_WARN(Code)				Kratos::OpenCL::CheckError(Code, __FILE__, __FUNCTION__, __LINE__, false);

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

	const char *ErrorString(cl_int _Code)
	{
		switch(_Code)
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

#if KRATOS_OCL_VERSION > 100

			KRATOS_OCL_DEFINE_ERROR_CODE(CL_MISALIGNED_SUB_BUFFER_OFFSET);
			KRATOS_OCL_DEFINE_ERROR_CODE(CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST);
#endif

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

#if KRATOS_OCL_VERSION > 100

			KRATOS_OCL_DEFINE_ERROR_CODE(CL_INVALID_PROPERTY);

#endif

			// Default case

			default:

				return "unknown";
		}
	}

//
// RaiseError
//
// Used to raise an OpenCL error and abort if requested
// Do not use directly; use the KRATOS_OCL_CHECK, KRATOS_OCL_WARN or KRATOS_OCL_CHECKED_EXPRESSION macros

	void RaiseError(cl_int _Code, const char *_FileName, const char *_Function, cl_uint _Line, bool _Abort, const char *_Expression = "")
	{
		std::cout <<
			std::endl <<
			"OpenCL reported " << ErrorString(_Code) << " error" <<
			_Expression <<
			" in function " << _Function << ", file " << _FileName << ", line " << _Line << "." << std::endl;

		if (_Abort)
		{
			std::cout <<
				"Aborting now; try setting a breakpoint on abort() to find the problem." << std::endl <<
				std::endl;

			abort();
		}
		else
		{
			std::cout <<
				"Execution will continue." << std::endl <<
				std::endl;
		}
	}

//
// CheckError
//
// Used to check an OpenCL return code and abort
// Do not use directly; use the KRATOS_OCL_CHECK or KRATOS_OCL_CHECKED_EXPRESSION macros

	inline void CheckError(cl_int _Code, const char *_FileName, const char *_Function, cl_uint _Line, bool _Abort, const char *_Expression = "")
	{
		if (_Code != CL_SUCCESS)
			RaiseError(_Code, _FileName, _Function, _Line, _Abort, _Expression);
	}

//
// DeviceTypeString
//
// Returns a string representation of an OpenCL device type

	const char *DeviceTypeString(cl_device_type _DeviceType)
	{
		switch (_DeviceType)
		{
			case CL_DEVICE_TYPE_CPU:

				return "CPU";

			case CL_DEVICE_TYPE_GPU:

				return "GPU";

			case CL_DEVICE_TYPE_ACCELERATOR:

				return "Accelerator";

			case CL_DEVICE_TYPE_DEFAULT:

				return "Default";

			case CL_DEVICE_TYPE_ALL:

				return "All";

			default:

				return "Unknown";
		}
	}

//
// DeviceGroup
//
// A class to manage a group of OpenCL devices

	class DeviceGroup
	{
		public:

			StringList CLSearchPath;

			DeviceIDList DeviceIDs;
			DeviceTypeList DeviceTypes;
			ContextList Contexts;
			CommandQueueList CommandQueues;

			ProgramList2D Programs;

			MemList2D Buffers;
			SizeTList2D BufferLengths;

			MemList2D Images;
			ImageDimensionList2D ImageDimensions;

			KernelList2D Kernels;
			SizeTList2D WorkGroupSizes;

			cl_uint DeviceNo;
			cl_device_type DeviceType;
			cl_platform_id PlatformID;

			//
			// DeviceGroup
			//
			// Constructor using a device type

			DeviceGroup(cl_device_type _DeviceType, bool _SingleDeviceOnly, const char *_PlatformVendor = ""): DeviceType(_DeviceType)
			{
				_Init(_SingleDeviceOnly, _PlatformVendor);
			}

			//
			// DeviceGroup
			//
			// Constructor using a device ID list

			DeviceGroup(DeviceIDList &_DeviceIDs, const char *_PlatformVendor = ""): DeviceIDs(_DeviceIDs)
			{
				// Check if we got an empty DeviceIDList

				if (_DeviceIDs.size() == 0)
				{
					std::cout <<
						"DeviceIDList specified cannot be empty." << std::endl <<
						"Aborting." << std::endl;

					abort();
				}

				_Init(false, _PlatformVendor);
			}

			//
			// ~DeviceGroup
			//
			// Destructor; frees OpenCL objects

			~DeviceGroup()
			{
				cl_int Err;

				// Releasing OpenCL objects

				for (cl_uint i = 0; i < Contexts.size(); i++)
				{
					Err = clReleaseContext(Contexts[i]);
					KRATOS_OCL_CHECK(Err);
				}

				for (cl_uint i = 0; i < CommandQueues.size(); i++)
				{
					Err = clReleaseCommandQueue(CommandQueues[i]);
					KRATOS_OCL_CHECK(Err);
				}

				for (cl_uint i = 0; i < Buffers.size(); i++)
				{
					for (cl_uint j = 0; j < Buffers[i].size(); j++)
					{
						Err = clReleaseMemObject(Buffers[i][j]);
						KRATOS_OCL_CHECK(Err);
					}
				}

				for (cl_uint i = 0; i < Images.size(); i++)
				{
					for (cl_uint j = 0; j < Images[i].size(); j++)
					{
						Err = clReleaseMemObject(Images[i][j]);
						KRATOS_OCL_CHECK(Err);
					}
				}

				for (cl_uint i = 0; i < Programs.size(); i++)
				{
					for (cl_uint j = 0; j < Programs[i].size(); j++)
					{
						Err = clReleaseProgram(Programs[i][j]);
						KRATOS_OCL_CHECK(Err);
					}
				}

				for (cl_uint i = 0; i < Kernels.size(); i++)
				{
					for (cl_uint j = 0; j < Kernels[i].size(); j++)
					{
						Err = clReleaseKernel(Kernels[i][j]);
						KRATOS_OCL_CHECK(Err);
					}
				}
			}

			//
			// _Init
			//
			// Initializes a device group based on given data in the constructor
			// Do not call directly

			void _Init(bool _SingleDeviceOnly, const char *_PlatformVendor)
			{
				cl_int Err;

				cl_uint PlatformNo;
				cl_platform_id *Platforms;
				std::string PlatformVendor(_PlatformVendor);
				char PlatformVendorData[KRATOS_OCL_PLATFORM_DATA_SIZE];

				// Try to find platform vendor specified, or choose a reasonabe default

				// Query available platforms

				Err = clGetPlatformIDs(0, NULL, &PlatformNo);
				KRATOS_OCL_CHECK(Err);

				Platforms = new cl_platform_id[PlatformNo];

				Err = clGetPlatformIDs(PlatformNo, Platforms, NULL);
				KRATOS_OCL_CHECK(Err);

				PlatformID = NULL;

				for (cl_uint i = 0; i < PlatformNo; i++)
				{
					PlatformID = Platforms[i];

					Err = clGetPlatformInfo(PlatformID, CL_PLATFORM_VENDOR, sizeof(PlatformVendorData), PlatformVendorData, NULL);
					KRATOS_OCL_CHECK(Err);

					// Check if this is the requested platform

					if (PlatformVendor.compare(PlatformVendorData) == 0)
					{
						break;
					}
				}

				delete [] Platforms;

				// Check if we finally found a platform

				if (PlatformID == NULL)
				{
					std::cout <<
						"No platform available." << std::endl <<
						"Aborting." << std::endl;

					abort();
				}

				// Try to find out which constructor has been called

				if (DeviceIDs.size() == 0)
				{
					// A device type is specified; check if we need only a single device

					if (_SingleDeviceOnly)
					{
						// Yes, a maximum of one device

						DeviceNo = 1;
					}
					else
					{
						// No, get all devices of that type

						Err = clGetDeviceIDs(PlatformID, DeviceType, 0, NULL, &DeviceNo);
						KRATOS_OCL_CHECK(Err);
					}

					DeviceIDs.resize(DeviceNo);

					Err = clGetDeviceIDs(PlatformID, DeviceType, DeviceNo, DeviceIDs.data(), NULL);
					KRATOS_OCL_CHECK(Err);
				}
				else
				{
					// An explicit list of device IDs is given

					DeviceNo = DeviceIDs.size();
				}

				// Get device type for each specified device

				DeviceTypes.resize(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					Err = clGetDeviceInfo(DeviceIDs[i], CL_DEVICE_TYPE, sizeof(cl_device_type), &DeviceTypes[i], NULL);
					KRATOS_OCL_CHECK(Err);
				}

				// Create contexts and command queues, one per device, so we can have fine control over everything (mostly buffer creation)

				cl_context_properties Properties[] = {CL_CONTEXT_PLATFORM, (cl_context_properties) PlatformID, 0};

				Contexts.resize(DeviceNo);
				CommandQueues.resize(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					Contexts[i] = clCreateContext(Properties, 1, &DeviceIDs[i], NULL, NULL, &Err);
					KRATOS_OCL_CHECK(Err);

					CommandQueues[i] = clCreateCommandQueue(Contexts[i], DeviceIDs[i], 0, &Err);
					KRATOS_OCL_CHECK(Err);
				}

				// Getting current working directory; Dummy is used to avoid warning

				char CurrentWorkingDirectory[FILENAME_MAX];

				char *Dummy = GETCWD(CurrentWorkingDirectory, sizeof(CurrentWorkingDirectory));
				Dummy = NULL;

				AddCLSearchPath(CurrentWorkingDirectory);
			}

			//
			// CreateBuffer
			//
			// Allocates a buffer on all devices

			cl_uint CreateBuffer(size_t _Size, cl_mem_flags _Flags)
			{
				cl_int Err;
				cl_uint BufferNo = Buffers.size();

				MemList CurrentBuffers(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					CurrentBuffers[i] = clCreateBuffer(Contexts[i], _Flags, _Size, NULL, &Err);
					KRATOS_OCL_CHECK(Err);
				}

				Buffers.push_back(CurrentBuffers);
				BufferLengths.push_back(SizeTList(DeviceNo, _Size));

				return BufferNo;
			}

			//
			// CreateBuffer
			//
			// Allocates a buffer on all devices with given sizes

			cl_uint CreateBuffer(SizeTList _Sizes, cl_mem_flags _Flags)
			{
				cl_int Err;
				cl_uint BufferNo = Buffers.size();

				MemList CurrentBuffers(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					CurrentBuffers[i] = clCreateBuffer(Contexts[i], _Flags, _Sizes[i], NULL, &Err);
					KRATOS_OCL_CHECK(Err);
				}

				Buffers.push_back(CurrentBuffers);
				BufferLengths.push_back(_Sizes);

				return BufferNo;
			}

			//
			// CreateSubBuffer
			//
			// Creates a sub-buffer on all devices
			// OpenCL 1.1 and later only

#if KRATOS_OCL_VERSION > 100

			cl_uint CreateSubBuffer(cl_uint _BufferIndex, size_t _Offset, size_t _Size)
			{
				cl_int Err;
				cl_uint BufferNo = Buffers.size();

				MemList CurrentBuffers(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					cl_mem_flags Flags;
					cl_buffer_region Region;

					// Get flags of original buffer

					Err = clGetMemObjectInfo(Buffers[_BufferIndex][i], CL_MEM_FLAGS, sizeof(Flags), &Flags, NULL);
					KRATOS_OCL_CHECK(Err);

					Region.origin = _Offset;
					Region.size = _Size;

					CurrentBuffers[i] = clCreateSubBuffer(Buffers[_BufferIndex][i], Flags, CL_BUFFER_CREATE_TYPE_REGION, &Region, &Err);
					KRATOS_OCL_CHECK(Err);
				}

				Buffers.push_back(CurrentBuffers);
				BufferLengths.push_back(SizeTList(DeviceNo, _Size));

				return BufferNo;
			}

#endif

			//
			// CreateSubBuffer
			//
			// Allocates a buffer on all devices with given offsets and sizes
			// OpenCL 1.1 and later only

#if KRATOS_OCL_VERSION > 100

			cl_uint CreateSubBuffer(cl_uint _BufferIndex, SizeTList _Offsets, SizeTList _Sizes)
			{
				cl_int Err;
				cl_uint BufferNo = Buffers.size();

				MemList CurrentBuffers(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					cl_mem_flags Flags;
					cl_buffer_region Region;

					// Get flags of original buffer

					Err = clGetMemObjectInfo(Buffers[_BufferIndex][i], CL_MEM_FLAGS, sizeof(Flags), &Flags, NULL);
					KRATOS_OCL_CHECK(Err);

					Region.origin = _Offsets[i];
					Region.size = _Sizes[i];

					CurrentBuffers[i] = clCreateSubBuffer(Buffers[_BufferIndex][i], Flags, CL_BUFFER_CREATE_TYPE_REGION, &Region, &Err);
					KRATOS_OCL_CHECK(Err);
				}

				Buffers.push_back(CurrentBuffers);
				BufferLengths.push_back(_Sizes);

				return BufferNo;
			}

#endif

			//
			// CreateBufferWithHostPointers
			//
			// Allocates a buffer on all devices

			cl_uint CreateBufferWithHostPointers(size_t _Size, cl_mem_flags _Flags, VoidPList &_HostPointers)
			{
				cl_int Err;
				cl_uint BufferNo = Buffers.size();

				MemList CurrentBuffers(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					CurrentBuffers[i] = clCreateBuffer(Contexts[i], _Flags, _Size, _HostPointers[i], &Err);
					KRATOS_OCL_CHECK(Err);
				}

				Buffers.push_back(CurrentBuffers);
				BufferLengths.push_back(SizeTList(DeviceNo, _Size));

				return BufferNo;
			}

			//
			// CreateBufferWithHostPointers
			//
			// Allocates a buffer on all devices with given sizes

			cl_uint CreateBufferWithHostPointers(SizeTList _Sizes, cl_mem_flags _Flags, VoidPList &_HostPointers)
			{
				cl_int Err;
				cl_uint BufferNo = Buffers.size();

				MemList CurrentBuffers(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					CurrentBuffers[i] = clCreateBuffer(Contexts[i], _Flags, _Sizes[i], _HostPointers[i], &Err);
					KRATOS_OCL_CHECK(Err);
				}

				Buffers.push_back(CurrentBuffers);
				BufferLengths.push_back(_Sizes);

				return BufferNo;
			}

			//
			// CopyBuffer
			//
			// Copies the content of a buffer on all devices

			void CopyBuffer(cl_uint _BufferIndex, CopyDirection _CopyDirection, VoidPList _HostPointers)
			{
				cl_int Err;

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					switch (_CopyDirection)
					{
						case HostToDevice:

							Err = clEnqueueWriteBuffer(CommandQueues[i], Buffers[_BufferIndex][i], CL_TRUE, 0, BufferLengths[_BufferIndex][i], _HostPointers[i], 0, NULL, NULL);
							KRATOS_OCL_CHECK(Err);

							break;

						case DeviceToHost:

							Err = clEnqueueReadBuffer(CommandQueues[i], Buffers[_BufferIndex][i], CL_TRUE, 0, BufferLengths[_BufferIndex][i], _HostPointers[i], 0, NULL, NULL);
							KRATOS_OCL_CHECK(Err);

							break;
					}
				}
			}

			//
			// CopyBuffer
			//
			// Copies the content of a buffer on a specific device

			void CopyBuffer(cl_uint _DeviceIndex, cl_uint _BufferIndex, CopyDirection _CopyDirection, void *_HostPtr)
			{
				cl_int Err;

				switch (_CopyDirection)
				{
					case HostToDevice:

						Err = clEnqueueWriteBuffer(CommandQueues[_DeviceIndex], Buffers[_BufferIndex][_DeviceIndex], CL_TRUE, 0, BufferLengths[_BufferIndex][_DeviceIndex], _HostPtr, 0, NULL, NULL);
						KRATOS_OCL_CHECK(Err);

						break;

					case DeviceToHost:

						Err = clEnqueueReadBuffer(CommandQueues[_DeviceIndex], Buffers[_BufferIndex][_DeviceIndex], CL_TRUE, 0, BufferLengths[_BufferIndex][_DeviceIndex], _HostPtr, 0, NULL, NULL);
						KRATOS_OCL_CHECK(Err);

						break;
				}
			}

			//
			// CopyBufferToBuffer
			//
			// Copies the content of a buffer to another buffer on all devices

			void CopyBufferToBuffer(cl_uint _SourceBufferIndex, cl_uint _DestinationBufferIndex)
			{
				cl_int Err;

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					// Check if the buffers are of the same length

					if (BufferLengths[_SourceBufferIndex][i] != BufferLengths[_DestinationBufferIndex][i])
					{
						std::cout <<
							"Buffers are not of the same size." << std::endl <<
							"Aborting." << std::endl;

						abort();
					}

					Err = clEnqueueCopyBuffer(CommandQueues[i], Buffers[_SourceBufferIndex][i], Buffers[_DestinationBufferIndex][i], 0, 0, BufferLengths[_SourceBufferIndex][i], 0, NULL, NULL);
				}
			}

			//
			// CopyBufferToBuffer
			//
			// Copies the content of a buffer to another buffer on a specific device

			void CopyBufferToBuffer(cl_uint _DeviceIndex, cl_uint _SourceBufferIndex, cl_uint _DestinationBufferIndex)
			{
				cl_int Err;

				// Check if the buffers are of the same length

				if (BufferLengths[_SourceBufferIndex][_DeviceIndex] != BufferLengths[_DestinationBufferIndex][_DeviceIndex])
				{
					std::cout <<
						"Buffers are not of the same size." << std::endl <<
						"Aborting." << std::endl;

					abort();
				}

				Err = clEnqueueCopyBuffer(CommandQueues[_DeviceIndex], Buffers[_SourceBufferIndex][_DeviceIndex], Buffers[_DestinationBufferIndex][_DeviceIndex], 0, 0, BufferLengths[_SourceBufferIndex][_DeviceIndex], 0, NULL, NULL);
			}

			//
			// MapBuffer
			//
			// Maps a buffer on all devices

			VoidPList MapBuffer(cl_uint _BufferIndex, cl_map_flags _Flags)
			{
				cl_int Err;

				VoidPList CurrentPointers(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					CurrentPointers[i] = clEnqueueMapBuffer(CommandQueues[i], Buffers[_BufferIndex][i], CL_TRUE, _Flags, 0, BufferLengths[_BufferIndex][i], 0, NULL, NULL, &Err);
					KRATOS_OCL_CHECK(Err);
				}

				return CurrentPointers;
			}

			//
			// CreateImage2D
			//
			// Creates a 2D image on all devices

			cl_uint CreateImage2D(size_t _Size1, size_t _Size2, cl_mem_flags _Flags, cl_channel_order _ChannelOrder, cl_channel_type _ChannelType)
			{
				cl_int Err;
				cl_uint ImageNo = Images.size();

				MemList CurrentImages(DeviceNo);

				cl_image_format Format = {_ChannelOrder, _ChannelType};

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					CurrentImages[i] = clCreateImage2D(Contexts[i], _Flags, &Format, _Size1, _Size2, 0, NULL, &Err);
					KRATOS_OCL_CHECK(Err);
				}

				Images.push_back(CurrentImages);
				ImageDimensions.push_back(ImageDimensionList(DeviceNo, ImageDimension(_Size1, _Size2)));

				return ImageNo;
			}

			//
			// CreateImage3D
			//
			// Creates a 3D image on all devices

			cl_uint CreateImage3D(size_t _Size1, size_t _Size2, size_t _Size3, cl_mem_flags _Flags, cl_channel_order _ChannelOrder, cl_channel_type _ChannelType)
			{
				cl_int Err;
				cl_uint ImageNo = Images.size();

				MemList CurrentImages(DeviceNo);

				cl_image_format Format = {_ChannelOrder, _ChannelType};

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					CurrentImages[i] = clCreateImage3D(Contexts[i], _Flags, &Format, _Size1, _Size2, _Size3, 0, 0, NULL, &Err);
					KRATOS_OCL_CHECK(Err);
				}

				Images.push_back(CurrentImages);
				ImageDimensions.push_back(ImageDimensionList(DeviceNo, ImageDimension(_Size1, _Size2, _Size3)));

				return ImageNo;
			}

			//
			// CreateImage2DWithHostPointers
			//
			// Creates a 2D image on all devices

			cl_uint CreateImage2DWithHostPointers(size_t _Size1, size_t _Size2, cl_mem_flags _Flags, cl_channel_order _ChannelOrder, cl_channel_type _ChannelType, VoidPList _HostPointers)
			{
				cl_int Err;
				cl_uint ImageNo = Images.size();

				MemList CurrentImages(DeviceNo);

				cl_image_format Format = {_ChannelOrder, _ChannelType};

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					CurrentImages[i] = clCreateImage2D(Contexts[i], _Flags, &Format, _Size1, _Size2, 0, _HostPointers[i], &Err);
					KRATOS_OCL_CHECK(Err);
				}

				Images.push_back(CurrentImages);
				ImageDimensions.push_back(ImageDimensionList(DeviceNo, ImageDimension(_Size1, _Size2)));

				return ImageNo;
			}

			//
			// CreateImage3DWithHostPointers
			//
			// Creates a 3D image on all devices

			cl_uint CreateImage3DWithHostPointers(size_t _Size1, size_t _Size2, size_t _Size3, cl_mem_flags _Flags, cl_channel_order _ChannelOrder, cl_channel_type _ChannelType, VoidPList _HostPointers)
			{
				cl_int Err;
				cl_uint ImageNo = Images.size();

				MemList CurrentImages(DeviceNo);

				cl_image_format Format = {_ChannelOrder, _ChannelType};

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					CurrentImages[i] = clCreateImage3D(Contexts[i], _Flags, &Format, _Size1, _Size2, _Size3, 0, 0, _HostPointers[i], &Err);
					KRATOS_OCL_CHECK(Err);
				}

				Images.push_back(CurrentImages);
				ImageDimensions.push_back(ImageDimensionList(DeviceNo, ImageDimension(_Size1, _Size2, _Size3)));

				return ImageNo;
			}

			//
			// CopyImage
			//
			// Copies the content of an image on all devices

			void CopyImage(cl_uint _ImageIndex, CopyDirection _CopyDirection, VoidPList _HostPointers)
			{
				cl_int Err;
				size_t Origin[3] = {0, 0, 0};

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					switch (_CopyDirection)
					{
						case HostToDevice:

							Err = clEnqueueWriteImage(CommandQueues[i], Images[_ImageIndex][i], CL_TRUE, Origin, ImageDimensions[_ImageIndex][i].Sizes, 0, 0, _HostPointers[i], 0, NULL, NULL);
							KRATOS_OCL_CHECK(Err);

							break;

						case DeviceToHost:

							Err = clEnqueueReadImage(CommandQueues[i], Images[_ImageIndex][i], CL_TRUE, Origin, ImageDimensions[_ImageIndex][i].Sizes, 0, 0, _HostPointers[i], 0, NULL, NULL);
							KRATOS_OCL_CHECK(Err);

							break;
					}
				}
			}

			//
			// CopyImage
			//
			// Copies the content of an image on a specific device

			void CopyImage(cl_uint _DeviceIndex, cl_uint _ImageIndex, CopyDirection _CopyDirection, VoidPList _HostPointers)
			{
				cl_int Err;
				size_t Origin[3] = {0, 0, 0};

				switch (_CopyDirection)
				{
					case HostToDevice:

						Err = clEnqueueWriteImage(CommandQueues[_DeviceIndex], Images[_ImageIndex][_DeviceIndex], CL_TRUE, Origin, ImageDimensions[_ImageIndex][_DeviceIndex].Sizes, 0, 0, _HostPointers[_DeviceIndex], 0, NULL, NULL);
						KRATOS_OCL_CHECK(Err);

						break;

					case DeviceToHost:

						Err = clEnqueueReadImage(CommandQueues[_DeviceIndex], Images[_ImageIndex][_DeviceIndex], CL_TRUE, Origin, ImageDimensions[_ImageIndex][_DeviceIndex].Sizes, 0, 0, _HostPointers[_DeviceIndex], 0, NULL, NULL);
						KRATOS_OCL_CHECK(Err);

						break;
				}
			}

			//
			// CopyImageToImage
			//
			// Copies the content of an image to another image on all devices

			void CopyImageToImage(cl_uint _SourceImageIndex, cl_uint _DestinationImageIndex)
			{
				cl_int Err;
				size_t Origin[3] = {0, 0, 0};

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					// Check if the images are of the same dimensions

					if (ImageDimensions[_SourceImageIndex][i].Sizes[0] != ImageDimensions[_DestinationImageIndex][i].Sizes[0] &&
						ImageDimensions[_SourceImageIndex][i].Sizes[1] != ImageDimensions[_DestinationImageIndex][i].Sizes[1] &&
						ImageDimensions[_SourceImageIndex][i].Sizes[2] != ImageDimensions[_DestinationImageIndex][i].Sizes[2])
					{
						std::cout <<
							"Images are not of the same dimension." << std::endl <<
							"Aborting." << std::endl;

						abort();
					}

					Err = clEnqueueCopyImage(CommandQueues[i], Images[_SourceImageIndex][i], Images[_DestinationImageIndex][i], Origin, Origin, ImageDimensions[_SourceImageIndex][i].Sizes, 0, NULL, NULL);
				}
			}

			//
			// CopyImageToImage
			//
			// Copies the content of a image to another image on a specific device

			void CopyImageToImage(cl_uint _DeviceIndex, cl_uint _SourceImageIndex, cl_uint _DestinationImageIndex)
			{
				cl_int Err;
				size_t Origin[3] = {0, 0, 0};

				// Check if the images are of the same dimensions

				if (ImageDimensions[_SourceImageIndex][_DeviceIndex].Sizes[0] != ImageDimensions[_DestinationImageIndex][_DeviceIndex].Sizes[0] &&
					ImageDimensions[_SourceImageIndex][_DeviceIndex].Sizes[1] != ImageDimensions[_DestinationImageIndex][_DeviceIndex].Sizes[1] &&
					ImageDimensions[_SourceImageIndex][_DeviceIndex].Sizes[2] != ImageDimensions[_DestinationImageIndex][_DeviceIndex].Sizes[2])
				{
					std::cout <<
						"Images are not of the same dimension." << std::endl <<
						"Aborting." << std::endl;

					abort();
				}

				Err = clEnqueueCopyImage(CommandQueues[_DeviceIndex], Images[_SourceImageIndex][_DeviceIndex], Images[_DestinationImageIndex][_DeviceIndex], Origin, Origin, ImageDimensions[_SourceImageIndex][_DeviceIndex].Sizes, 0, NULL, NULL);
			}

			//
			// CopyBufferToImage
			//
			// Copies the content of a buffer to an image on all devices

			void CopyBufferToImage(cl_uint _SourceBufferIndex, cl_uint _DestinationImageIndex)
			{
				cl_int Err;
				size_t Origin[3] = {0, 0, 0};

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					// We use the dimensions of the image; the user is responsible for the results

					Err = clEnqueueCopyBufferToImage(CommandQueues[i], Buffers[_SourceBufferIndex][i], Images[_DestinationImageIndex][i], 0, Origin, ImageDimensions[_DestinationImageIndex][i].Sizes, 0, NULL, NULL);
				}
			}

			//
			// CopyBufferToImage
			//
			// Copies the content of a buffer to an image on a specific device

			void CopyBufferToImage(cl_uint _DeviceIndex, cl_uint _SourceBufferIndex, cl_uint _DestinationImageIndex)
			{
				cl_int Err;
				size_t Origin[3] = {0, 0, 0};

				// We use the dimensions of the image; the user is responsible for the results

				Err = clEnqueueCopyBufferToImage(CommandQueues[_DeviceIndex], Buffers[_SourceBufferIndex][_DeviceIndex], Images[_DestinationImageIndex][_DeviceIndex], 0, Origin, ImageDimensions[_DestinationImageIndex][_DeviceIndex].Sizes, 0, NULL, NULL);
			}

			//
			// CopyImageToBuffer
			//
			// Copies the content of an image to a buffer on all devices

			void CopyImageToBuffer(cl_uint _DeviceIndex, cl_uint _SourceImageIndex, cl_uint _DestinationBufferIndex)
			{
				cl_int Err;
				size_t Origin[3] = {0, 0, 0};

				// We use the dimensions of the image; the user is responsible for the results

				Err = clEnqueueCopyImageToBuffer(CommandQueues[_DeviceIndex], Images[_SourceImageIndex][_DeviceIndex], Buffers[_DestinationBufferIndex][_DeviceIndex], Origin, ImageDimensions[_SourceImageIndex][_DeviceIndex].Sizes, 0, 0, NULL, NULL);
			}

			//
			// CopyImageToBuffer
			//
			// Copies the content of an image to a buffer on all devices

			void CopyImageToBuffer(cl_uint _SourceImageIndex, cl_uint _DestinationBufferIndex)
			{
				cl_int Err;
				size_t Origin[3] = {0, 0, 0};

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					// We use the dimensions of the image; the user is responsible for the results

					Err = clEnqueueCopyImageToBuffer(CommandQueues[i], Images[_SourceImageIndex][i], Buffers[_DestinationBufferIndex][i], Origin, ImageDimensions[_SourceImageIndex][i].Sizes, 0, 0, NULL, NULL);
				}
			}

			//
			// MapImage
			//
			// Maps an image on all devices

			VoidPList MapImage(cl_uint _ImageIndex, cl_map_flags _Flags)
			{
				cl_int Err;
				size_t Origin[3] = {0, 0, 0};

				VoidPList CurrentPointers(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					size_t ImageRowPitch;
					size_t ImageSlicePitch;

					CurrentPointers[i] = clEnqueueMapImage(CommandQueues[i], Images[_ImageIndex][i], CL_TRUE, _Flags, Origin, ImageDimensions[_ImageIndex][i].Sizes, &ImageRowPitch, &ImageSlicePitch, 0, NULL, NULL, &Err);
					KRATOS_OCL_CHECK(Err);
				}

				return CurrentPointers;
			}

			//
			// AddCLSearchPath
			//
			// Adds a path to search list for .cl files

			void AddCLSearchPath(const char *Path)
			{
				std::string TempPath(Path);

				// Check if we need to add a PATH_SEPARATOR
				if (TempPath[TempPath.size() - 1] != PATH_SEPARATOR)
				{
					TempPath += PATH_SEPARATOR;
				}

				CLSearchPath.push_back(TempPath);
			}

			//
			// BuildProgramFromFile
			//
			// Load a program source file and build it

			cl_uint BuildProgramFromFile(const char *_FileName, const char *_BuildOptions = "")
			{
				cl_int Err;
				cl_uint ProgramNo = Programs.size();

				std::ifstream SourceFile;
				std::stringstream Source;

				for (cl_uint i = 0; i < CLSearchPath.size(); i++)
				{
					// Try to open the file

					std::string TempFileName = CLSearchPath[i] + _FileName;
					SourceFile.open(TempFileName.c_str());

					if (SourceFile.is_open())
					{
						break;
					}
				}

				// Check if the file was found finally

				if (!SourceFile.is_open())
				{
					std::cout <<
						"An error occurred reading the kernel file. Try adding the appropriate path using AddCLSearchPath()." << std::endl <<
						"Aborting." << std::endl;

					abort();
				}

				Source << SourceFile.rdbuf();
				SourceFile.close();

				std::string SourceStr = Source.str();

				const char *SourceText = SourceStr.c_str();
				size_t SourceLen = SourceStr.size();

				// Build program for all devices

				ProgramList CurrentPrograms(DeviceNo);

				std::string Options(_BuildOptions);

				// Define KRATOS_OCL_VERSION macro inside the program

				// Add an space if not empty

				if (Options.size() != 0)
				{
					Options += " ";
				}

				Options += "-DKRATOS_OCL_VERSION=" KRATOS_OCL_VERSION_STRING;

				// Add CLSearchPath to compiler's include path, so #include's work as intended

				for (cl_uint i = 0; i < CLSearchPath.size(); i++)
				{
					// Add -I option

					Options += " -I";
					Options += CLSearchPath[i];
				}

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					CurrentPrograms[i] = clCreateProgramWithSource(Contexts[i], 1, &SourceText, &SourceLen, &Err);
					KRATOS_OCL_CHECK(Err);

					Err = clBuildProgram(CurrentPrograms[i], 0, NULL, Options.c_str(), NULL, NULL);

					if (Err == CL_BUILD_PROGRAM_FAILURE)
					{
						char *BuildLog = new char[KRATOS_OCL_PROGRAM_BUILD_LOG_SIZE];

						Err = clGetProgramBuildInfo(CurrentPrograms[i], DeviceIDs[i], CL_PROGRAM_BUILD_LOG, KRATOS_OCL_PROGRAM_BUILD_LOG_SIZE, BuildLog, NULL);
						KRATOS_OCL_CHECK(Err);

						std::cout <<
							"Build log:" << std::endl <<
							BuildLog << std::endl;

						delete [] BuildLog;

						KRATOS_OCL_CHECK(CL_BUILD_PROGRAM_FAILURE);
					}
					else
					{
						KRATOS_OCL_CHECK(Err);
					}
				}

				// Append this to the list

				Programs.push_back(CurrentPrograms);

				return ProgramNo;
			}

			//
			// RegisterKernel
			//
			// Register a kernel in the built program on all devices

			cl_uint RegisterKernel(cl_uint _ProgramIndex, const char *_KernelName)
			{
				cl_int Err;
				cl_uint KernelNo = Kernels.size();

				// Register the kernel on all devices

				KernelList CurrentKernels(DeviceNo);
				SizeTList CurrentWorkGroupSizes(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					CurrentKernels[i] = clCreateKernel(Programs[_ProgramIndex][i], _KernelName, &Err);
					KRATOS_OCL_CHECK(Err);

					size_t WorkGroupSize, MaxWorkGroupSize, PreferredWorkGroupSizeMultiple;

					Err = clGetKernelWorkGroupInfo(CurrentKernels[i], DeviceIDs[i], CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &MaxWorkGroupSize, NULL);
					KRATOS_OCL_CHECK(Err);

#if KRATOS_OCL_VERSION > 100

					Err = clGetKernelWorkGroupInfo(CurrentKernels[i], DeviceIDs[i], CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &PreferredWorkGroupSizeMultiple, NULL);
					KRATOS_OCL_CHECK(Err);

#else
					// For OpenCL 1.0 only

					PreferredWorkGroupSizeMultiple = 64;

#endif

					// Select an optimal WorkGroupSize

					if (DeviceTypes[i] == CL_DEVICE_TYPE_CPU)
					{
						// On CPU devices we use a fixed group size

						WorkGroupSize = KRATOS_OCL_CPU_WORK_GROUP_SIZE;
					}
					else
					{
						if (MaxWorkGroupSize < PreferredWorkGroupSizeMultiple)
						{
							// MaxWorkGroupSize too low; use it as WorkGroupSize

							WorkGroupSize = MaxWorkGroupSize;
						}
						else
						{
							// For best performance, we use the maximum work group size which is divisible by WorkGroupSizeMultiple

							WorkGroupSize = (MaxWorkGroupSize / PreferredWorkGroupSizeMultiple) * PreferredWorkGroupSizeMultiple;
						}
					}

					CurrentWorkGroupSizes[i] = WorkGroupSize;
				}

				// Append these to the lists

				Kernels.push_back(CurrentKernels);
				WorkGroupSizes.push_back(CurrentWorkGroupSizes);

				return KernelNo;
			}

			//
			// SetKernelArg
			//
			// Set a kernel argument on all devices

			template <typename Type> void SetKernelArg(cl_uint _KernelIndex, cl_uint _ArgIndex, Type _Value)
			{
				cl_int Err;

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					Err = clSetKernelArg(Kernels[_KernelIndex][i], _ArgIndex, sizeof(Type), &_Value);
					KRATOS_OCL_CHECK(Err);
				}
			}

			//
			// SetKernelArg
			//
			// Set a kernel argument on a specific device

			template <typename Type> void SetKernelArg(cl_uint _DeviceIndex, cl_uint _KernelIndex, cl_uint _ArgIndex, Type _Value)
			{
				cl_int Err;

				Err = clSetKernelArg(Kernels[_KernelIndex][_DeviceIndex], _ArgIndex, sizeof(Type), &_Value);
				KRATOS_OCL_CHECK(Err);
			}

			//
			// SetBufferAsKernelArg
			//
			// Sets a buffer as a kernel argument on all devices

			void SetBufferAsKernelArg(cl_uint _KernelIndex, cl_uint _ArgIndex, cl_uint _BufferIndex)
			{
				cl_int Err;

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					Err = clSetKernelArg(Kernels[_KernelIndex][i], _ArgIndex, sizeof(cl_mem), &Buffers[_BufferIndex][i]);
					KRATOS_OCL_CHECK(Err);
				}
			}

			//
			// SetBufferAsKernelArg
			//
			// Sets a buffer as a kernel argument on a specific device

			void SetBufferAsKernelArg(cl_uint _DeviceIndex, cl_uint _KernelIndex, cl_uint _ArgIndex, cl_uint _BufferIndex)
			{
				cl_int Err;

				Err = clSetKernelArg(Kernels[_KernelIndex][_DeviceIndex], _ArgIndex, sizeof(cl_mem), &Buffers[_BufferIndex][_DeviceIndex]);
				KRATOS_OCL_CHECK(Err);
			}

			//
			// SetImageAsKernelArg
			//
			// Sets an image as a kernel argument on all devices

			void SetImageAsKernelArg(cl_uint _KernelIndex, cl_uint _ArgIndex, cl_uint _ImageIndex)
			{
				cl_int Err;

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					Err = clSetKernelArg(Kernels[_KernelIndex][i], _ArgIndex, sizeof(cl_mem), &Images[_ImageIndex][i]);
					KRATOS_OCL_CHECK(Err);
				}
			}

			//
			// SetImageAsKernelArg
			//
			// Sets an image as a kernel argument on a specific device

			void SetImageAsKernelArg(cl_uint _DeviceIndex, cl_uint _KernelIndex, cl_uint _ArgIndex, cl_uint _ImageIndex)
			{
				cl_int Err;

				Err = clSetKernelArg(Kernels[_KernelIndex][_DeviceIndex], _ArgIndex, sizeof(cl_mem), &Images[_ImageIndex][_DeviceIndex]);
				KRATOS_OCL_CHECK(Err);
			}

			//
			// ExecuteKernel
			//
			// Execute a kernel on all devices

			void ExecuteKernel(cl_uint _KernelIndex, cl_uint _GlobalWorkSize)
			{
				cl_int Err;

				// Enqueue kernels

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					// We may need to use a bigger GlobalWorkSize, to keep it a multiple of preferred size

					size_t GlobalWorkSize = ((_GlobalWorkSize + WorkGroupSizes[_KernelIndex][i] - 1) / WorkGroupSizes[_KernelIndex][i]) * WorkGroupSizes[_KernelIndex][i];

					Err = clEnqueueNDRangeKernel(CommandQueues[i], Kernels[_KernelIndex][i], 1, NULL, &GlobalWorkSize, &WorkGroupSizes[_KernelIndex][i], 0, NULL, NULL);
					KRATOS_OCL_CHECK(Err);

					Err = clFlush(CommandQueues[i]);
					KRATOS_OCL_CHECK(Err);
				}

				// Wait for kernels to finish

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					Err = clFinish(CommandQueues[i]);
					KRATOS_OCL_CHECK(Err);
				}
			}

			//
			// ExecuteKernel
			//
			// Execute a kernel on a specific device

			void ExecuteKernel(cl_uint _DeviceIndex, cl_uint _KernelIndex, cl_uint _GlobalWorkSize)
			{
				cl_int Err;

				// Enqueue kernel; we may need to use a bigger GlobalWorkSize, to keep it a multiple of preferred size

				size_t GlobalWorkSize = ((_GlobalWorkSize + WorkGroupSizes[_KernelIndex][_DeviceIndex] - 1) / WorkGroupSizes[_KernelIndex][_DeviceIndex]) * WorkGroupSizes[_KernelIndex][_DeviceIndex];

				Err = clEnqueueNDRangeKernel(CommandQueues[_DeviceIndex], Kernels[_KernelIndex][_DeviceIndex], 1, NULL, &GlobalWorkSize, &WorkGroupSizes[_KernelIndex][_DeviceIndex], 0, NULL, NULL);
				KRATOS_OCL_CHECK(Err);

				Err = clFlush(CommandQueues[_DeviceIndex]);
				KRATOS_OCL_CHECK(Err);

				// Wait for kernel to finish

				Err = clFinish(CommandQueues[_DeviceIndex]);
				KRATOS_OCL_CHECK(Err);
			}
	};

}  // namespace OpenCL

}  // namespace Kratos

#endif // KRATOS_OPENCL_INTERFACE_H_INCLUDED
