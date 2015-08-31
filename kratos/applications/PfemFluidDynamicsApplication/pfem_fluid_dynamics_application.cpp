//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Last modified by:    $Author:               JMCarbonell $
//   Date:                $Date:              September 2015 $
//   Revision:            $Revision:                     0.0 $
//
//


// System includes


// External includes 


// Project includes
#include "includes/define.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"

#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"

#include "geometries/triangle_3d_3.h"

#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"

#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"

#include "geometries/line_2d.h"

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"

#include "pfem_fluid_dynamics_application.h"

namespace Kratos
{
  //Create Variables


  KratosPfemFluidDynamicsApplication::KratosPfemFluidDynamicsApplication()    
  {}
  
  void KratosPfemFluidDynamicsApplication::Register()
  {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    //KratosFluidDynamicsApplication::Register();

    std::cout << "      KRATOS' _ |  __|  __| \\   |                  " << std::endl;
    std::cout << "             '__| |_   |_   |\\ /|                  " << std::endl;
    std::cout << "            _|   _|    |__ _|  _|FLUID DYNAMICS     " << std::endl;
    std::cout << "Initializing KratosPfemFluidDynamicsApplication... " << std::endl;
    
    //Register Elements


    //Register Conditions

    //Register Constitutive Laws

    //Register Flow Rules
 
    //Register Yield Criterion
 
    //Register Hardening Laws
 
    //Register Variables


  }
  
}  // namespace Kratos.


