#ifndef VIENNACL_OCL_CONTEXT_HPP_
#define VIENNACL_OCL_CONTEXT_HPP_

/* =========================================================================
   Copyright (c) 2010-2012, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at
               
   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file viennacl/ocl/context.hpp
    @brief Represents an OpenCL context within ViennaCL
*/

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <algorithm>
#include <vector>
#include <map>
#include "viennacl/ocl/forwards.h"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/program.hpp"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/command_queue.hpp"

namespace viennacl
{
  namespace ocl
  {
    class context
    {
      typedef std::vector< viennacl::ocl::program >   ProgramContainer;
      
      public:
        context() : initialized_(false),
                    device_type_(CL_DEVICE_TYPE_DEFAULT),
                    current_device_id(0),
                    default_device_num_(1) {}
        
        //////// Get and set default number of devices per context */
        /** @brief Returns the maximum number of devices to be set up for the context */
        std::size_t default_device_num() const { return default_device_num_; }
        
        /** @brief Sets the maximum number of devices to be set up for the context */
        void default_device_num(std::size_t new_num) { default_device_num_ = new_num; }
        
        ////////// get and set preferred device type /////////////////////
        /** @brief Returns the default device type for the context */
        cl_device_type default_device_type()
        {
          return device_type_;
        }
        
        /** @brief Sets the device type for this context */
        void default_device_type(cl_device_type dtype) 
        { 
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
          std::cout << "ViennaCL: Setting new device type for context " << h_ << std::endl;
          #endif
          if (!initialized_)
            device_type_ = dtype; //assume that the user provided a correct value
        }
        
        //////////////////// get devices //////////////////
        /** @brief Returns a vector with all devices in this context */
        std::vector<viennacl::ocl::device> const & devices() const
        {
          return devices_;
        }
        
        /** @brief Returns the current device */
        viennacl::ocl::device const & current_device() const
        {
          //std::cout << "Current device id in context: " << current_device_id << std::endl;
          return devices_[current_device_id];
        }
        
        /** @brief Switches the current device to the i-th device in this context */
        void switch_device(size_t i)
        {
          assert(i >= 0 && i < devices_.size());
          current_device_id = i;
        }

        /** @brief If the supplied device is used within the context, it becomes the current active device. */
        void switch_device(viennacl::ocl::device const & d)
        {
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
          std::cout << "ViennaCL: Setting new current device for context " << h_ << std::endl;
          #endif
          bool found = false;
          for (size_t i=0; i<devices_.size(); ++i)
          {
            if (devices_[i] == d)
            {
              found = true;
              current_device_id = i;
              break;
            }
          }
          if (found == false)
            std::cerr << "ViennaCL: Warning: Could not set device " << d.name() << " for context." << std::endl;
        }
        
        /** @brief Add a device to the context. Must be done before the context is initialized */
        void add_device(viennacl::ocl::device const & d)
        {
          assert(!initialized_ && "Device must be added to context before it is initialized!");
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
          std::cout << "ViennaCL: Adding new device to context " << h_ << std::endl;
          #endif
          if (std::find(devices_.begin(), devices_.end(), d) == devices_.end())
            devices_.push_back(d);
        }

        /** @brief Add a device to the context. Must be done before the context is initialized */
        void add_device(cl_device_id d)
        {
          assert(!initialized_ && "Device must be added to context before it is initialized!");
          add_device(viennacl::ocl::device(d));
        }


        /////////////////////// initialize context ///////////////////
        
        /** @brief Initializes a new context */
        void init()
        {
          init_new();
        }

        /** @brief Initializes the context from an existing, user-supplied context */
        void init(cl_context c)
        {
          init_existing(c);
        }

/*        void existing_context(cl_context context_id)
        {
          assert(!initialized_ && "ViennaCL: FATAL error: Provided a new context for an already initialized context.");
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
          std::cout << "ViennaCL: Reusing existing context " << h_ << std::endl;
          #endif
          h_ = context_id;
        }*/
        
        ////////////////////// create memory /////////////////////////////
        /** @brief Creates a memory buffer within the context
        *
        *  @param flags  OpenCL flags for the buffer creation
        *  @param size   Size of the memory buffer in bytes
        *  @param ptr    Optional pointer to CPU memory, with which the OpenCL memory should be initialized
        */
        viennacl::ocl::handle<cl_mem> create_memory(cl_mem_flags flags, unsigned int size, void * ptr = NULL)
        {
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
          std::cout << "ViennaCL: Creating memory of size " << size << " for context " << h_ << std::endl;
          #endif
          if (ptr)
            flags |= CL_MEM_COPY_HOST_PTR;
          cl_int err;
          viennacl::ocl::handle<cl_mem> mem = clCreateBuffer(h_.get(), flags, size, ptr, &err);
          VIENNACL_ERR_CHECK(err);
          return mem;
        }

        /** @brief Creates a memory buffer within the context initialized from the supplied data
        *
        *  @param flags  OpenCL flags for the buffer creation
        *  @param _buffer A vector (STL vector, ublas vector, etc.)
        */
        template < typename SCALARTYPE, typename A, template <typename, typename> class VectorType >
        viennacl::ocl::handle<cl_mem> create_memory(cl_mem_flags flags, const VectorType<SCALARTYPE, A> & _buffer)
        {
          return create_memory(flags, static_cast<cl_uint>(sizeof(SCALARTYPE) * _buffer.size()), (void*)&_buffer[0]);
        }
        
        //////////////////// create queues ////////////////////////////////
        
        /** @brief Adds an existing queue for the given device to the context */
        void add_queue(cl_device_id dev, cl_command_queue q)
        {
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
          std::cout << "ViennaCL: Adding existing queue " << q << " for device " << dev << " to context " << h_ << std::endl;
          #endif
          queues_[dev].push_back(viennacl::ocl::command_queue(q, dev));
        }
        
        /** @brief Adds a queue for the given device to the context */
        void add_queue(cl_device_id dev)
        {
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
          std::cout << "ViennaCL: Adding new queue for device " << dev << " to context " << h_ << std::endl;
          #endif
          cl_int err;
          viennacl::ocl::handle<cl_command_queue> temp = clCreateCommandQueue(h_.get(), dev, 0, &err);
          VIENNACL_ERR_CHECK(err);
          
          queues_[dev].push_back(viennacl::ocl::command_queue(temp, dev));
        }

        /** @brief Adds a queue for the given device to the context */
        void add_queue(viennacl::ocl::device d) { add_queue(d.id()); }

        //get queue for default device:
        viennacl::ocl::command_queue & get_queue()
        {
          return queues_[devices_[current_device_id].id()][0];
        }
        
        //get a particular queue:
        /** @brief Returns the queue with the provided index for the given device */
        viennacl::ocl::command_queue & get_queue(cl_device_id dev, size_t i = 0)
        {
          assert(i >= 0 && i < queues_.size() && "In class 'context': id invalid in get_queue()");
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
          std::cout << "ViennaCL: Getting queue " << i << " for device " << dev << " in context " << h_ << std::endl;
          #endif
          unsigned int device_index;
          for (device_index = 0; device_index < devices_.size(); ++device_index)
          {
            if (devices_[device_index] == dev)
              break;
          }
          
          assert(device_index < devices_.size() && "Device not within context");
          
          return queues_[devices_[device_index].id()][i];
        }
        
        /////////////////// create program ///////////////////////////////
        /** @brief Adds a program to the context
        */
        viennacl::ocl::program & add_program(cl_program p, std::string const & prog_name)
        {
          programs_.push_back(viennacl::ocl::program(p, prog_name));
          return programs_.back();
        }
        
        /** @brief Adds a new program with the provided source to the context
        */
        viennacl::ocl::program & add_program(std::string const & source, std::string const & prog_name)
        {
          const char * source_text = source.c_str();
          size_t source_size = source.size();
          cl_int err;
          
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
          std::cout << "ViennaCL: Adding program '" << prog_name << "' to context " << h_ << std::endl;
          #endif
          
          viennacl::ocl::handle<cl_program> temp = clCreateProgramWithSource(h_.get(), 1, (const char **)&source_text, &source_size, &err);
          VIENNACL_ERR_CHECK(err);
          
          err = clBuildProgram(temp.get(), 0, NULL, NULL, NULL, NULL);
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_BUILD)
            char buffer[1024];
            cl_build_status status;
            clGetProgramBuildInfo(temp, devices_[0].id(), CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &status, NULL);
            clGetProgramBuildInfo(temp, devices_[0].id(), CL_PROGRAM_BUILD_LOG, sizeof(char)*1024, &buffer, NULL);
            std::cout << "Build Scalar: Err = " << err << " Status = " << status << std::endl;
            std::cout << "Log: " << buffer << std::endl;
            //std::cout << "Sources: " << source << std::endl;
          #endif
          VIENNACL_ERR_CHECK(err);

          programs_.push_back(viennacl::ocl::program(temp, prog_name));
          
          return programs_.back();
        }
        
        /** @brief Returns the program with the provided name */
        viennacl::ocl::program & get_program(std::string const & name)
        {
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
          std::cout << "ViennaCL: Getting program '" << name << "' from context " << h_ << std::endl;
          #endif
          for (ProgramContainer::iterator it = programs_.begin();
                it != programs_.end();
                ++it)
          {
            if (it->name() == name)
              return *it;
          }
          std::cerr << "Could not find program '" << name << "'" << std::endl;
          assert(!"In class 'context': name invalid in get_program()");
          return programs_[0];  //return a defined object
        }
        
        /** @brief Returns the program with the provided id */
        viennacl::ocl::program & get_program(size_t id)
        {
          assert(id >= 0 && id < programs_.size() && "In class 'context': id invalid in get_program()");
          return programs_[id];
        }
        
        /** @brief Returns the number of programs within this context */
        size_t program_num() { return programs_.size(); }

        /** @brief Returns the number of devices within this context */
        size_t device_num() { return devices_.size(); }
        
        /** @brief Returns the context handle */
        const viennacl::ocl::handle<cl_context> & handle() const { return h_; }
        
        /** @brief Less-than comparable for compatibility with std:map  */
        bool operator<(context const & other) const
        {
          return h_.get() < other.h_.get();
        }
        
      private:
        /** @brief Initialize a new context. Reuse any previously supplied information (devices, queues) */
        void init_new()
        {
          assert(!initialized_ && "ViennaCL FATAL error: Context already created!");

          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
          std::cout << "ViennaCL: Initializing new ViennaCL context." << std::endl;
          #endif
          
          cl_int err;
          std::vector<cl_device_id> device_id_array;
          if (devices_.empty()) //get the default device if user has not yet specified a list of devices
          {
            //create an OpenCL context for the provided devices:
            #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
            std::cout << "ViennaCL: Setting all devices for context..." << std::endl;
            #endif
            
            platform pf;
            std::vector<device> devices = pf.devices(device_type_);
            #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
            std::cout << "ViennaCL: Number of devices for context: " << devices.size() << std::endl;
            #endif
            for (size_t i=0; i<devices.size(); ++i)
              devices_.push_back(devices[i]);
            
            if (devices.size() == 0)
            {
              std::cerr << "ViennaCL: FATAL ERROR: No devices of type '";
              switch (device_type_)
              {
                case CL_DEVICE_TYPE_CPU:          std::cout << "CPU"; break;
                case CL_DEVICE_TYPE_GPU:          std::cout << "GPU"; break;
                case CL_DEVICE_TYPE_ACCELERATOR:  std::cout << "ACCELERATOR"; break;
                case CL_DEVICE_TYPE_DEFAULT:      std::cout << "DEFAULT"; break;
                default:
                  std::cout << "UNKNOWN" << std::endl;
              }
              std::cout << "' found!" << std::endl;
            }
          }
          
          //extract list of device ids:
          for (std::vector< viennacl::ocl::device >::const_iterator iter = devices_.begin();
                                                                    iter != devices_.end();
                                                                  ++iter)
            device_id_array.push_back(iter->id());
            
          // deprecated cl_uint device_num = std::max(default_device_num_, device_id_array.size());

		  cl_uint device_num = (default_device_num_>device_id_array.size())? default_device_num_ :device_id_array.size();
          h_ = clCreateContext(0, 
                               device_num,
                               &(device_id_array[0]),
                               NULL, NULL, &err);
          VIENNACL_ERR_CHECK(err);
          
          initialized_ = true;
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
          std::cout << "ViennaCL: Initialization of new ViennaCL context done." << std::endl;
          #endif
        }
        
        /** @brief Reuses a supplied context. */
        void init_existing(cl_context c)
        {
          assert(!initialized_ && "ViennaCL FATAL error: Context already created!");
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
          std::cout << "ViennaCL: Initialization of ViennaCL context from existing context." << std::endl;
          #endif
          
          //set context handle:
          h_ = c;
          
          if (devices_.empty())
          {
            //get devices for context:
            cl_int err;
            cl_uint num_devices;
            size_t temp;
            //Note: The obvious
            //  err = clGetContextInfo(h_, CL_CONTEXT_NUM_DEVICES, sizeof(cl_uint), &num_devices, NULL);
            //does not work with NVIDIA OpenCL stack!
            err = clGetContextInfo(h_.get(), CL_CONTEXT_DEVICES, VIENNACL_OCL_MAX_DEVICE_NUM * sizeof(cl_device_id), NULL, &temp);
            VIENNACL_ERR_CHECK(err);
            assert(temp > 0 && "ViennaCL: FATAL error: Provided context does not contain any devices!");
            num_devices = temp / sizeof(cl_device_id);
            
            #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
            std::cout << "ViennaCL: Reusing context with " << num_devices << " devices." << std::endl;
            #endif
            
            std::vector<cl_device_id> device_ids(num_devices);
            err = clGetContextInfo(h_.get(), CL_CONTEXT_DEVICES, num_devices * sizeof(cl_device_id), &(device_ids[0]), NULL);
            VIENNACL_ERR_CHECK(err);
            
            for (size_t i=0; i<num_devices; ++i)
              devices_.push_back(viennacl::ocl::device(device_ids[i]));
          }
          current_device_id = 0;
          
          initialized_ = true;
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_CONTEXT)
          std::cout << "ViennaCL: Initialization of ViennaCL context from existing context done." << std::endl;
          #endif
        }       
        
        
        bool initialized_;
        cl_device_type device_type_;
        viennacl::ocl::handle<cl_context> h_;
        std::vector< viennacl::ocl::device > devices_;
        unsigned int current_device_id;
        std::size_t default_device_num_;
        ProgramContainer programs_;
        std::map< cl_device_id, std::vector< viennacl::ocl::command_queue> > queues_;
    }; //context
    
  }
}

#endif
