//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
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
#include "geometries/line_2d_3.h"

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

  KratosPfemFluidDynamicsApplication::KratosPfemFluidDynamicsApplication():
    KratosApplication("PfemFluidDynamicsApplication"),
    mTwoStepUpdatedLagrangianVPElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mTwoStepUpdatedLagrangianVPElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node<3> >(Element::GeometryType::PointsArrayType(6)))),
    mTwoStepUpdatedLagrangianVPElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mTwoStepUpdatedLagrangianVPElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node<3> >(Element::GeometryType::PointsArrayType(10)))),
    mTwoStepUpdatedLagrangianVPSolidElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mTwoStepUpdatedLagrangianVPSolidElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node<3> >(Element::GeometryType::PointsArrayType(6)))),
    mTwoStepUpdatedLagrangianVPSolidElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mTwoStepUpdatedLagrangianVPSolidElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node<3> >(Element::GeometryType::PointsArrayType(10)))),
    mUpdatedLagrangianVSolidElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mUpdatedLagrangianVSolidElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node<3> >(Element::GeometryType::PointsArrayType(6)))),
    mUpdatedLagrangianVSolidElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mUpdatedLagrangianVSolidElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node<3> >(Element::GeometryType::PointsArrayType(10)))),
    mTwoStepUpdatedLagrangianVPFluidElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mTwoStepUpdatedLagrangianVPFluidElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node<3> >(Element::GeometryType::PointsArrayType(6)))),
    mTwoStepUpdatedLagrangianVPFluidElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mTwoStepUpdatedLagrangianVPFluidElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node<3> >(Element::GeometryType::PointsArrayType(10))))
  {
  }
  
  void KratosPfemFluidDynamicsApplication::Register()
  {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    //KratosFluidDynamicsApplication::Register();

    std::cout << "            ___  __           ___ _      _    _          " << std::endl;
    std::cout << "     KRATOS| _ \\/ _|___ _ __ | __| |_  _(_)__| |         " << std::endl;
    std::cout << "           |  _/  _/ -_) '  \\| _|| | || | / _` |         " << std::endl;
    std::cout << "           |_| |_| \\___|_|_|_|_| |_|\\_,_|_\\__,_|DYNAMICS " << std::endl;
    std::cout << "Initializing KratosPfemFluidDynamicsApplication...       " << std::endl;
       
    //Register Variables (variables created in pfem_fluid_dynamics_application_variables.cpp)

    // Material postprocess + invariants
    // KRATOS_REGISTER_VARIABLE(M_MODULUS)  
    // KRATOS_REGISTER_VARIABLE(PATCH_INDEX);
    // KRATOS_REGISTER_VARIABLE(NORMVELOCITY);
    KRATOS_REGISTER_VARIABLE(FREESURFACE);
    KRATOS_REGISTER_VARIABLE(INITIAL_DELTA_TIME);
    KRATOS_REGISTER_VARIABLE(CURRENT_DELTA_TIME);
    KRATOS_REGISTER_VARIABLE(TIME_INTERVAL_CHANGED);
    KRATOS_REGISTER_VARIABLE(BAD_VELOCITY_CONVERGENCE);
    KRATOS_REGISTER_VARIABLE(BAD_PRESSURE_CONVERGENCE);

    KRATOS_REGISTER_VARIABLE(FLOW_INDEX)
    KRATOS_REGISTER_VARIABLE(YIELD_SHEAR)
    KRATOS_REGISTER_VARIABLE(ADAPTIVE_EXPONENT)

    KRATOS_REGISTER_VARIABLE(PRESSURE_VELOCITY)
      KRATOS_REGISTER_VARIABLE(PRESSURE_ACCELERATION)

      //Register Elements
      KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPElement2D",mTwoStepUpdatedLagrangianVPElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPElement2Dquadratic",mTwoStepUpdatedLagrangianVPElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPElement3D",mTwoStepUpdatedLagrangianVPElement3D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPElement3Dquadratic",mTwoStepUpdatedLagrangianVPElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPSolidElement2D",mTwoStepUpdatedLagrangianVPSolidElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPSolidElement2Dquadratic",mTwoStepUpdatedLagrangianVPSolidElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPSolidElement3D",mTwoStepUpdatedLagrangianVPSolidElement3D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPSolidElement3Dquadratic",mTwoStepUpdatedLagrangianVPSolidElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianVSolidElement2D",mUpdatedLagrangianVSolidElement2D);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianVSolidElement2Dquadratic",mUpdatedLagrangianVSolidElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianVSolidElement3D",mUpdatedLagrangianVSolidElement3D);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianVSolidElement3Dquadratic",mUpdatedLagrangianVSolidElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidElement2D",mTwoStepUpdatedLagrangianVPFluidElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidElement2Dquadratic",mTwoStepUpdatedLagrangianVPFluidElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidElement3D",mTwoStepUpdatedLagrangianVPFluidElement3D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidElement3Dquadratic",mTwoStepUpdatedLagrangianVPFluidElement3Dquadratic);


    //Register Conditions

    //Register Constitutive Laws

    //Register Flow Rules
 
    //Register Yield Criterion
 
    //Register Hardening Laws
 
  }
  
}  // namespace Kratos.


