#ifndef VIENNACL_IO_KERNEL_PARAMETERS_HPP
#define VIENNACL_IO_KERNEL_PARAMETERS_HPP

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


/** @file kernel_parameters.hpp
    @brief This file holds the code necessary for reading kernel parameters from XML files using pugixml
*/

#include "viennacl/ocl/backend.hpp"
#include "pugixml/src/pugixml.hpp"

namespace viennacl
{
  namespace io 
  {
    namespace tag 
    {
      static std::string root     = "parameters";
      static std::string devices  = "devices";   
      static std::string device   = "device";   
      static std::string name     = "name";
      static std::string driver   = "driver";
      static std::string compun   = "computeunits";      
      static std::string workgrp  = "workgroupsize";            
      static std::string tests    = "tests";
      static std::string test     = "test";                     
      static std::string numeric  = "numeric";
      static std::string kernels  = "kernels";
      static std::string kernel   = "kernel";
      static std::string params   = "params";
      static std::string param    = "param";
      static std::string value    = "value";   
      static std::string alignment = "alignment";   
    } // end namespace tag

    namespace val {
      static std::string globsize = "globalsize";
      static std::string locsize  = "localsize";   
      static std::string vec      = "vector";   
      static std::string matrix   = "matrix";   
      static std::string compmat  = "compressed_matrix";
      static std::string fl       = "float";   
      static std::string dbl      = "double";      
    }

    /** @brief  A XML parameter database using PugiXML. Allows to add tests for different devices and the like */
    struct parameter_database 
    {
      parameter_database ()
      {
          root = doc.append_child();
          root.set_name(tag::root.c_str());
          last = root;
          
          devices_open = false;
          tests_open = false;      
          kernels_open = false;
          parameters_open = false;      
      }   
      
      void add_device()
      {
          pugi::xml_node dev;
          if(devices_open)
          {
            dev = devices.append_child();
            dev.set_name(tag::device.c_str());      
          }
          else
          {
            devices = last.append_child();
            devices.set_name(tag::devices.c_str());
            
            dev = devices.append_child();
            dev.set_name(tag::device.c_str());
            
            devices_open = true;
          }
          last = dev;
      }
      
      void add_test()
      {
          pugi::xml_node test;
          if(tests_open)
          {
            test = tests.append_child();
            test.set_name(tag::test.c_str());      
          }
          else
          {
            tests = last.append_child();
            tests.set_name(tag::tests.c_str());
            
            test = tests.append_child();
            test.set_name(tag::test.c_str());
            
            tests_open = true;
          }
          last = test;
          // close the current kernels section
          // so a new one is created for this new test
          kernels_open = false;      
      }   

      void add_kernel()
      {
          pugi::xml_node kern;
          if(kernels_open)
          {
            kern = kernels.append_child();
            kern.set_name(tag::kernel.c_str());      
          }
          else
          {
            kernels = last.append_child();
            kernels.set_name(tag::kernels.c_str());
            
            kern = kernels.append_child();
            kern.set_name(tag::kernel.c_str());
            
            kernels_open = true;
          }
          last = kern;
          
          // close the current parameters section
          // so a new one is created for this new kernel
          parameters_open = false;
      }      
      
      void add_parameter()
      {
          pugi::xml_node para;
          
          if(parameters_open)
          {
            para = parameters.append_child();
            para.set_name(tag::param.c_str());      
          }
          else
          {
            parameters = last.append_child();
            parameters.set_name(tag::params.c_str());
            
            para = parameters.append_child();
            para.set_name(tag::param.c_str());
            
            parameters_open = true;
          }
          last = para;
      }         
      
      template<typename ValueT>
      void add_data_node(std::string tagstr, ValueT data)
      {
          std::stringstream ss;
          ss << data;
          add_data_node(tagstr, ss.str());
      }   
      
      void add_data_node(std::string tagstr, std::string data)
      {
          pugi::xml_node node = last.append_child();
          
          if(tagstr == tag::name)
            node.set_name(tag::name.c_str());
          else if(tagstr == tag::driver)
            node.set_name(tag::driver.c_str());      
          else if(tagstr == tag::numeric)
            node.set_name(tag::numeric.c_str());      
          else if(tagstr == tag::alignment)
            node.set_name(tag::alignment.c_str());      
          else if(tagstr == tag::value)
            node.set_name(tag::value.c_str());      
          else if(tagstr == tag::compun)
            node.set_name(tag::compun.c_str());      
          else if(tagstr == tag::workgrp)
            node.set_name(tag::workgrp.c_str());                        
          else
            std::cout << "# Error adding data node: node tag not recognized .." << std::endl;
          node.append_child(pugi::node_pcdata).set_value(data.c_str());
      }

      void load(std::string filename)
      {
          doc.load_file(filename.c_str());
      }

      void dump(std::string filename)
      {
          std::ofstream outstream(filename.c_str());
          this->dump(outstream);
          outstream.close();
      }
      
      void dump(std::ostream& stream = std::cout)
      {
          doc.save(stream, "  ");
      }

      pugi::xml_document   doc;
      pugi::xml_node       root;
      pugi::xml_node       devices;
      pugi::xml_node       tests;   
      pugi::xml_node       kernels;      
      pugi::xml_node       parameters;         
      pugi::xml_node       last;   
      
      bool devices_open;
      bool tests_open;   
      bool kernels_open;      
      bool parameters_open;         

    };
    
    /** @brief Helper meta class that returns the first letter of a particular type (float or double) */
    template <typename T>
    struct first_letter_of_type
    {
      static char get(); //intentionally not implemented, class must be specialized
    };
    
    template <>
    struct first_letter_of_type <float>
    {
      static char get() { return 'f'; } 
    };
    
    template <>
    struct first_letter_of_type <double>
    {
      static char get() { return 'd'; } 
    };
    
    template <typename T>
    struct program_for_vcltype
    {
      static std::string get();  //intentionally not implemented, class must be specialized
    };
    
    template <typename T, unsigned int ALIGNMENT>
    struct program_for_vcltype < viennacl::vector<T, ALIGNMENT> >
    {
      static std::string get()
      {
        std::stringstream ss;
        ss << first_letter_of_type<T>::get() << "_vector_" << ALIGNMENT;
        return ss.str();
      } 
    };
    
    template <typename T, unsigned int ALIGNMENT>
    struct program_for_vcltype < viennacl::matrix<T, row_major, ALIGNMENT> >
    {
      static std::string get()
      {
        std::stringstream ss;
        ss << first_letter_of_type<T>::get() << "_matrix_row_" << ALIGNMENT;
        return ss.str();
      } 
    };

    template <typename T, unsigned int ALIGNMENT>
    struct program_for_vcltype < viennacl::matrix<T, column_major, ALIGNMENT> >
    {
      static std::string get()
      {
        std::stringstream ss;
        ss << first_letter_of_type<T>::get() << "_matrix_col_" << ALIGNMENT;
        return ss.str();
      } 
    };
    
    template <typename T, unsigned int ALIGNMENT>
    struct program_for_vcltype < viennacl::compressed_matrix<T,  ALIGNMENT> >
    {
      static std::string get()
      {
        std::stringstream ss;
        ss << first_letter_of_type<T>::get() << "_compressed_matrix_" << ALIGNMENT;
        return ss.str();
      } 
    };

    template<typename SCALARTYPE, unsigned int ALIGNMENT>
    void set_kernel_params(std::string program_name,
                          std::string kernel_name,
                          unsigned int glob, //total no. of threads
                          unsigned int loc)  //threads per work group
    {
      //get kernel from pool and set work sizes:
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(program_name, kernel_name);
      k.global_work_size(0, glob);
      k.local_work_size(0, loc);
      
      //std::cout << "Setting [" << glob << ", " << loc << "] for kernel " << kernel_name << std::endl;
    }

    template<typename VclBasicType>
    void tune_impl(parameter_database& paras, std::string parent)
    {
      typedef typename VclBasicType::value_type::value_type   SCALARTYPE;
      
      // create dummy vectors; the kernels have to be created ..
      VclBasicType    dummy;

      // extract the kernels for which parameters are present
      std::string          kernel_str = parent+"/kernels/kernel/name/text()";
      pugi::xpath_node_set kernel_res = paras.doc.select_nodes(kernel_str.c_str());      

      typedef std::vector<std::string>   kernels_type;
      kernels_type kernels;
      std::cout << "Retrieving kernels..." << std::endl;
      for (pugi::xpath_node_set::const_iterator it = kernel_res.begin(); it != kernel_res.end(); ++it)
      {
          std::stringstream ss;
          it->node().print(ss, "  ");
          std::string kern(ss.str());
          kern.erase(std::remove(kern.begin(), kern.end(), '\n'), kern.end()); //trim trailing linebreak
          kernels.push_back(kern);
      }
      
      // retrieve the actual parameters
      std::cout << "Retrieving actual parameters..." << std::endl;
      for(typename kernels_type::iterator iter = kernels.begin();
          iter != kernels.end(); iter++)
      {
          // retrieving the work group ..
          std::string          wg_str = parent+"/kernels/kernel[name='"+*iter+"']/params/param[name='"+val::globsize+"']/value/text()";
          pugi::xpath_node_set wg_res = paras.doc.select_nodes(wg_str.c_str());  

          unsigned int global_size(0);
          
          std::stringstream ss;
          ss << wg_res[0].node().value();
          ss >> global_size;
          
          // retrieving the local_workers ..
          std::string          lw_str = parent+"/kernels/kernel[name='"+*iter+"']/params/param[name='"+val::locsize+"']/value/text()";
          pugi::xpath_node_set lw_res = paras.doc.select_nodes(lw_str.c_str());  

          unsigned int local_workers(0);
          
          ss.clear();
          ss << lw_res[0].node().value();
          ss >> local_workers;         
          
          //std::cout << "kernel: " << *iter << " wg: " << work_group << " lw: " << local_workers << std::endl;

          // set the parameters
          set_kernel_params<SCALARTYPE,1> (program_for_vcltype<VclBasicType>::get(), *iter, global_size, local_workers);
          //set_kernel_params<SCALARTYPE,4> (*iter, work_group * local_workers, local_workers);         
          //set_kernel_params<SCALARTYPE,16>(*iter, work_group * local_workers, local_workers);                 
      }
    }

    /** @brief Helper meta-class that converts a type to a string */
    template <typename T>
    struct to_string {};

    template <>
    struct to_string<float>
    {
      static std::string get() { return "float"; }
    };

    template <>
    struct to_string<double>
    {
      static std::string get() { return "double"; }
    };

    /** @brief The interface function for reading kernel parameters
    *
    * @tparam VclBasicType  The ViennaCL type for which parameters should be read
    * @param filename       Relative filename to the XML file where the parameters are located in
    */
    template<typename VclBasicType>
    void read_kernel_parameters(std::string filename)
    {
      typedef typename VclBasicType::value_type::value_type   SCALARTYPE;
      
      parameter_database  paras;
      paras.load(filename);
      
      std::string devname   = viennacl::ocl::current_device().name();
      
      // check if tune parameters for the current device are present
      std::string          device_str = "/parameters/devices/device[name='"+devname+"']";
      pugi::xpath_node_set device_res = paras.doc.select_nodes(device_str.c_str());
      
      if(device_res.size() == 0)
      {
          std::cout << "Tuner: There are no parameters for this device present!" << std::endl;
          // evaluate the parameters for this device?
      }
      
      // check if tune parameters for float exist
      std::string          numeric_str = device_str+"/tests/test[numeric='"+to_string<SCALARTYPE>::get()+"']";
      pugi::xpath_node_set numeric_res = paras.doc.select_nodes(numeric_str.c_str());

      if(numeric_res.size() > 0)
      {
          tune_impl<VclBasicType>(paras, numeric_str);
      }
      else
      {
          std::cout << "Tuner: There are no parameters for numeric type float present!" << std::endl;   
      }

  //    // check if tune parameters for double exist
  //    std::string          double_str = device_str+"/tests/test[numeric='"+val::dbl+"']";
  //    pugi::xpath_node_set double_res = paras.doc.select_nodes(double_str.c_str());
  // 
  //    if(double_res.size() > 0)
  //    {
  //       tune_impl<double>(paras, double_str);
  //    }
  //    else
  //    {
  //       std::cout << "Tuner: There are no parameters for numeric type double present!" << std::endl;   
  //    }

    }

  } // end namespace io

} // end namespace viennacl

#endif
