#ifndef VIENNACL_OCL_DEVICE_UTILS_HPP_
#define VIENNACL_OCL_DEVICE_UTILS_HPP_

/* =========================================================================
   Copyright (c) 2010-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file viennacl/ocl/device_utils.hpp
    @brief Various utility implementations for dispatching with respect to the different devices available on the market.
*/

#define VIENNACL_OCL_MAX_DEVICE_NUM  8

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif


#include <stddef.h>
#include <map>
#include <string>

#include "viennacl/forwards.h"

namespace viennacl
{
  namespace ocl
  {

    static const cl_uint intel_id = 32902;
    static const cl_uint nvidia_id = 4318;
    static const cl_uint amd_id = 4098;
    static const cl_uint unknown_id = 0;

    //Architecture Family
    enum device_architecture_family{
      //NVidia
      Tesla,
      Fermi,
      Kepler,

      //AMD
      Evergreen,
      NorthernIslands,
      SouthernIslands,

      UNKNOWN
    };

    static device_architecture_family get_device_architecture(cl_uint vendor_id, std::string const & name){

      /*-NVidia-*/
      if(vendor_id==nvidia_id){
        //GeForce
        vcl_size_t found=0;
        if((found= name.find("GeForce",0)) != std::string::npos){
          if((found = name.find_first_of("123456789", found)) != std::string::npos){
            switch (name[found]) {
              case '2' : return Tesla;
              case '3' : return Tesla;

              case '4' : return Fermi;
              case '5' : return Fermi;

              case '6' : return Kepler;
              case '7' : return Kepler;

              default: return UNKNOWN;
            }
          }
          else
            return UNKNOWN;
        }

        //Tesla
        else if((found = name.find("Tesla",0)) != std::string::npos){
          if((found = name.find("CMK", found)) != std::string::npos){
            switch(name[found]){
              case 'C' : return Fermi;
              case 'M' : return Fermi;

              case 'K' : return Kepler;

              default : return UNKNOWN;
            }
          }
          else
            return UNKNOWN;
        }

        else
          return UNKNOWN;
      }

      /*-AMD-*/
      else if(vendor_id==amd_id){

#define VIENNACL_DEVICE_MAP(device,arch)if(name.find(device,0)!=std::string::npos) return arch;

        //Evergreen
        VIENNACL_DEVICE_MAP("Cedar",Evergreen);
        VIENNACL_DEVICE_MAP("Redwood",Evergreen);
        VIENNACL_DEVICE_MAP("Juniper",Evergreen);
        VIENNACL_DEVICE_MAP("Cypress",Evergreen);
        VIENNACL_DEVICE_MAP("Hemlock",Evergreen);

        //NorthernIslands
        VIENNACL_DEVICE_MAP("Caicos",NorthernIslands);
        VIENNACL_DEVICE_MAP("Turks",NorthernIslands);
        VIENNACL_DEVICE_MAP("Barts",NorthernIslands);
        VIENNACL_DEVICE_MAP("Cayman",NorthernIslands);
        VIENNACL_DEVICE_MAP("Antilles",NorthernIslands);

        //SouthernIslands
        VIENNACL_DEVICE_MAP("Cape",SouthernIslands);
        VIENNACL_DEVICE_MAP("Bonaire",SouthernIslands);
        VIENNACL_DEVICE_MAP("Pitcaim",SouthernIslands);
        VIENNACL_DEVICE_MAP("Tahiti",SouthernIslands);
        VIENNACL_DEVICE_MAP("Malta",SouthernIslands);

#undef VIENNACL_DEVICE_MAP

        return UNKNOWN;

      }

      /*-Other-*/
      else{
        return UNKNOWN;
      }

    }


  }
} //namespace viennacl

#endif

/*@}*/
