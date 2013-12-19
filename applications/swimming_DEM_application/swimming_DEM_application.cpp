//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "swimming_DEM_application.h"
#include "geometries/point_3d.h"
#include "geometries/line_3d_2.h"
#include "../DEM_application/DEM_application.h"
#include "../FluidDynamicsApplication/fluid_dynamics_application.h"

namespace Kratos
{
        
  KRATOS_CREATE_VARIABLE(double, AUX_DOUBLE_VAR)
  KRATOS_CREATE_VARIABLE(double, POWER_LAW_N)
  KRATOS_CREATE_VARIABLE(double, POWER_LAW_K)
  KRATOS_CREATE_VARIABLE(double, GEL_STRENGTH)
  
KratosSwimmingDEMApplication::KratosSwimmingDEMApplication():
  mMonolithicDEMCoupled2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
  mMonolithicDEMCoupled3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))))
{}

void KratosSwimmingDEMApplication::Register()
{
  // calling base class register to register Kratos components
  KratosApplication::Register();
  std::cout << "Initializing KratosSwimmingDEMApplication... " << std::endl;
                

  KRATOS_REGISTER_VARIABLE(AUX_DOUBLE_VAR)
  KRATOS_REGISTER_VARIABLE(POWER_LAW_N)
  KRATOS_REGISTER_VARIABLE(POWER_LAW_K)        
  KRATOS_REGISTER_VARIABLE(GEL_STRENGTH)

      /* Define In Global variables.cpp */

  KRATOS_REGISTER_ELEMENT("MonolithicDEMCoupled2D", mMonolithicDEMCoupled2D)
  KRATOS_REGISTER_ELEMENT("MonolithicDEMCoupled3D", mMonolithicDEMCoupled3D)

  Serializer::Register( "VariablesList", mVariablesList );
 }

}  // namespace Kratos.


