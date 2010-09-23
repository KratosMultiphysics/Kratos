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
#include <fstream>
#include <sstream>
#include <vector>


// External includes


// Project includes


namespace Kratos
{

namespace OpenCL
{

// Useful vector types

	typedef enum {HostToDevice, DeviceToHost} BufferCopyDirection;

	typedef std::vector<cl_device_type> DeviceTypeList;
	typedef std::vector<cl_device_id> DeviceIDList;
	typedef std::vector<cl_context> ContextList;
	typedef std::vector<cl_command_queue> CommandQueueList;
	typedef std::vector<cl_program> ProgramList;

	typedef std::vector<void *> VoidPList;

	typedef std::vector<cl_mem> MemList;
	typedef std::vector<MemList> MemList2D;

	typedef std::vector<cl_kernel> KernelList;
	typedef std::vector<KernelList> KernelList2D;

	typedef std::vector<size_t> SizeTList;
	typedef std::vector<SizeTList> SizeTList2D;

//
// KRATOS_OCL_CPU_WORK_GROUP_SIZE
//
// Size of a work group on CPU devices

#define KRATOS_OCL_CPU_WORK_GROUP_SIZE		2

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
			default:

				return "unknown";
		}
	}

//
// RaiseError
//
// Used to raise an OpenCL error and abort if requested
// Do not use directly; use the KRATOS_OCL_CHECK, KRATOS_OCL_WARN or KRATOS_OCL_CHECKED_EXPRESSION macros

	void RaiseError(cl_int _Code, const char *_FileName, const char *_Function, unsigned int _Line, bool _Abort, const char *_Expression = "")
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

	inline void CheckError(cl_int _Code, const char *_FileName, const char *_Function, unsigned int _Line, bool _Abort, const char *_Expression = "")
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

			DeviceIDList DeviceIDs;
			DeviceTypeList DeviceTypes;
			ContextList Contexts;
			CommandQueueList CommandQueues;
			ProgramList Programs;

			MemList2D Buffers;
			SizeTList2D BufferLengths;
			KernelList2D Kernels;
			SizeTList2D WorkGroupSizes;

			cl_uint DeviceNo;
			cl_device_type DeviceType;
			cl_platform_id PlatformID;

			//
			// DeviceGroup
			//
			// Constructor using a device type

			DeviceGroup(cl_device_type _DeviceType, const char *_PlatformVendor = ""): DeviceType(_DeviceType)
			{
				_Init(_PlatformVendor);
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
					std::cout << "DeviceIDList specified cannot be empty. Aborting." << std::endl;
					abort();
				}

				_Init(_PlatformVendor);
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

				for (cl_uint i = 0; i < Programs.size(); i++)
				{
					Err = clReleaseProgram(Programs[i]);
					KRATOS_OCL_CHECK(Err);
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

			void _Init(const char *_PlatformVendor)
			{
				cl_int Err;

				cl_uint PlatformNo;
				cl_platform_id *Platforms;
				std::string PlatformVendor(_PlatformVendor);
				char CharData[1024];

				if (PlatformVendor.compare("") == 0)
				{
					// No specific platform vendor is specified, default to first one
					Err = clGetPlatformIDs(1, &PlatformID, NULL);
					KRATOS_OCL_CHECK(Err);
				}
				else
				{
					// Try to find platform vendor specified
					Err = clGetPlatformIDs(0, NULL, &PlatformNo);
					KRATOS_OCL_CHECK(Err);

					Platforms = new cl_platform_id[PlatformNo];

					Err = clGetPlatformIDs(PlatformNo, Platforms, NULL);
					KRATOS_OCL_CHECK(Err);

					// We default to first platform, in case we cannot find the requested platform vendor
					PlatformID = Platforms[0];

					for (cl_uint i = 0; i < PlatformNo; i++)
					{
						Err = clGetPlatformInfo(Platforms[i], CL_PLATFORM_VENDOR, sizeof(CharData), CharData, NULL);
						KRATOS_OCL_CHECK(Err);

						// Check if this is the requested platform
						if (PlatformVendor.compare(CharData) == 0)
							PlatformID = Platforms[i];
					}

					delete [] Platforms;
				}

				// Try to find out which constructor has been called
				if (DeviceIDs.size() == 0)
				{
					// A device type is specified; get all devices of that type
					Err = clGetDeviceIDs(PlatformID, DeviceType, 0, NULL, &DeviceNo);
					KRATOS_OCL_CHECK(Err);

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
			}

			//
			// CreateBuffer
			//
			// Allocates a buffer on all devices

			void CreateBuffer(size_t _Size, cl_mem_flags _Flags)
			{
				cl_int Err;

				MemList CurrentBuffers(DeviceNo);
				SizeTList CurrentBufferLengths(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					CurrentBufferLengths[i] = _Size;

					CurrentBuffers[i] = clCreateBuffer(Contexts[i], _Flags, _Size, NULL, &Err);
					KRATOS_OCL_CHECK(Err);
				}

				Buffers.push_back(CurrentBuffers);
				BufferLengths.push_back(CurrentBufferLengths);
			}

			//
			// CreateBufferWithHostPtrs
			//
			// Allocates a buffer on all devices

			void CreateBufferWithHostPtrs(size_t _Size, cl_mem_flags _Flags, VoidPList &_HostPtrs)
			{
				cl_int Err;

				MemList CurrentBuffers(DeviceNo);
				SizeTList CurrentBufferLengths(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					CurrentBufferLengths[i] = _Size;

					CurrentBuffers[i] = clCreateBuffer(Contexts[i], _Flags, _Size, _HostPtrs[i], &Err);
					KRATOS_OCL_CHECK(Err);
				}

				Buffers.push_back(CurrentBuffers);
				BufferLengths.push_back(CurrentBufferLengths);
			}

			//
			// CopyBuffer
			//
			// Copies the content of a buffer on all devices

			void CopyBuffer(int _BufferIndex, BufferCopyDirection _CopyDirection, VoidPList _HostPtrs)
			{
				cl_int Err;

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					switch (_CopyDirection)
					{
						case HostToDevice:

							Err = clEnqueueWriteBuffer(CommandQueues[i], Buffers[_BufferIndex][i], CL_TRUE, 0, BufferLengths[_BufferIndex][i], _HostPtrs[i], 0, NULL, NULL);
							KRATOS_OCL_CHECK(Err);

							break;

						case DeviceToHost:

							Err = clEnqueueWriteBuffer(CommandQueues[i], Buffers[_BufferIndex][i], CL_TRUE, 0, BufferLengths[_BufferIndex][i], _HostPtrs[i], 0, NULL, NULL);
							KRATOS_OCL_CHECK(Err);

							break;
					}
				}
			}

			//
			// CopyBuffer
			//
			// Copies the content of a buffer on a specific device

			void CopyBuffer(int _DeviceIndex, int _BufferIndex, BufferCopyDirection _CopyDirection, void *_HostPtr)
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
			// MapBuffer
			//
			// Maps a buffer on all devices

			VoidPList MapBuffer(int _BufferIndex, cl_map_flags _Flags)
			{
				cl_int Err;

				VoidPList CurrentPtrs(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					CurrentPtrs[i] = clEnqueueMapBuffer(CommandQueues[i], Buffers[_BufferIndex][i], CL_TRUE, _Flags, 0, BufferLengths[_BufferIndex][i], 0, NULL, NULL, &Err);
					KRATOS_OCL_CHECK(Err);
				}

				return CurrentPtrs;
			}

			//
			// BuildProgramFromFile
			//
			// Load a program source file and build it; for simplicity we do not support multiple programs

			void BuildProgramFromFile(const char *_FileName)
			{
				cl_int Err;

				std::ifstream SourceFile(_FileName);
				std::stringstream Source;

				Source << SourceFile.rdbuf();
				std::string SourceStr = Source.str();

				const char *SourceText = SourceStr.c_str();
				size_t SourceLen = SourceStr.size();

				// Build program for all devices
				Programs.resize(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					Programs[i] = clCreateProgramWithSource(Contexts[i], 1, &SourceText, &SourceLen, &Err);
					KRATOS_OCL_CHECK(Err);

					Err = clBuildProgram(Programs[i], 0, NULL, NULL, NULL, NULL);
					KRATOS_OCL_WARN(Err);

					if (Err != CL_SUCCESS)
					{
						char CharData[1024];

						Err = clGetProgramBuildInfo(Programs[i], DeviceIDs[i], CL_PROGRAM_BUILD_LOG, sizeof(CharData), CharData, NULL);
						KRATOS_OCL_CHECK(Err);
					}
				}
			}

			//
			// RegisterKernel
			//
			// Register a kernel in the built program on all devices

			void RegisterKernel(const char *_KernelName)
			{
				cl_int Err;

				// Register the kernel on all devices
				KernelList CurrentKernels(DeviceNo);
				SizeTList CurrentWorkGroupSizes(DeviceNo);

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					CurrentKernels[i] = clCreateKernel(Programs[i], _KernelName, &Err);
					KRATOS_OCL_CHECK(Err);

					size_t WorkGroupSize, MaxWorkGroupSize, PreferredWorkGroupSizeMultiple;

					Err = clGetKernelWorkGroupInfo(CurrentKernels[i], DeviceIDs[i], CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &MaxWorkGroupSize, NULL);
					KRATOS_OCL_CHECK(Err);

					Err = clGetKernelWorkGroupInfo(CurrentKernels[i], DeviceIDs[i], CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &PreferredWorkGroupSizeMultiple, NULL);
					KRATOS_OCL_CHECK(Err);

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
			}

			//
			// SetKernelArg
			//
			// Set a kernel argument on all devices

			template <typename Type> void SetKernelArg(int _KernelIndex, int _ArgIndex, Type _Value)
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

			template <typename Type> void SetKernelArg(int _DeviceIndex, int _KernelIndex, int _ArgIndex, Type _Value)
			{
				cl_int Err;

				Err = clSetKernelArg(Kernels[_KernelIndex][_DeviceIndex], _ArgIndex, sizeof(Type), &_Value);
				KRATOS_OCL_CHECK(Err);
			}

			//
			// SetBufferAsKernelArg
			//
			// Sets a buffer as a kernel argument on all devices

			void SetBufferAsKernelArg(int _KernelIndex, int _ArgIndex, int _BufferIndex)
			{
				cl_int Err;

				for (cl_uint i = 0; i < DeviceNo; i++)
				{
					Err = clSetKernelArg(Kernels[_KernelIndex][i], _ArgIndex, sizeof(cl_mem), &Buffers[_BufferIndex][i]);
					KRATOS_OCL_CHECK(Err);
				}
			}

			//
			// ExecuteKernel
			//
			// Execute a kernel on all devices

			void ExecuteKernel(int _KernelIndex, cl_uint _GlobalWorkSize)
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

			void ExecuteKernel(int _DeviceIndex, int _KernelIndex, cl_uint _GlobalWorkSize)
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
