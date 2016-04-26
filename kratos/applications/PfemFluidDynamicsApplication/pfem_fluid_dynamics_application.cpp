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
 //  //solution
//   KRATOS_CREATE_VARIABLE(double, IMPOSED_WATER_PRESSURE )

//   //constitutive law	
//   KRATOS_CREATE_VARIABLE(double, MEAN_ERROR )

//   //material
//   // KRATOS_CREATE_VARIABLE(double, PRE_CONSOLIDATION_STRESS )
//   // KRATOS_CREATE_VARIABLE(double, OVER_CONSOLIDATION_RATIO )
//   // KRATOS_CREATE_VARIABLE(double, INITIAL_SHEAR_MODULUS )
//   KRATOS_CREATE_VARIABLE(double, WATER_BULK_MODULUS )
//   KRATOS_CREATE_VARIABLE(double, PERMEABILITY )
//   // KRATOS_CREATE_VARIABLE(double, NORMAL_COMPRESSION_SLOPE )
//   // KRATOS_CREATE_VARIABLE(double, SWELLING_SLOPE )
//   // KRATOS_CREATE_VARIABLE(double, CRITICAL_STATE_LINE )
//   // KRATOS_CREATE_VARIABLE(double, ALPHA_SHEAR )
//   KRATOS_CREATE_VARIABLE(double, INITIAL_POROSITY )

//   //element
//   KRATOS_CREATE_VARIABLE(Vector, DARCY_FLOW )
//   KRATOS_CREATE_VARIABLE(Matrix, TOTAL_CAUCHY_STRESS )

//   // transfer variables and initial
//   // KRATOS_CREATE_VARIABLE(Matrix, ELASTIC_LEFT_CAUCHY_GREEN_TENSOR )
//   // KRATOS_CREATE_VARIABLE(Vector, ELASTIC_LEFT_CAUCHY_GREEN_VECTOR )

//   // KRATOS_CREATE_VARIABLE(Matrix, KIRCHHOFF_STRESS_TENSOR )
//   KRATOS_CREATE_VARIABLE(Vector, ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS )

//    //domain definition
//   KRATOS_CREATE_VARIABLE(unsigned int, DOMAIN_LABEL )
//   KRATOS_CREATE_VARIABLE(int         , RIGID_WALL )
//   KRATOS_CREATE_VARIABLE(double      , WALL_TIP_RADIUS )
//   KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WALL_REFERENCE_POINT )
//   KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WALL_VELOCITY )

//   //contact condition
//   KRATOS_CREATE_VARIABLE(Condition::Pointer, MASTER_CONDITION )
//   KRATOS_CREATE_VARIABLE(WeakPointerVector< Element >, MASTER_ELEMENTS )
//   KRATOS_CREATE_VARIABLE(WeakPointerVector< Node<3> >, MASTER_NODES )

//   //contact 
//   KRATOS_CREATE_VARIABLE(bool,   FRICTION_ACTIVE )
//   KRATOS_CREATE_VARIABLE(double, PENALTY_PARAMETER )
//   KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_NORMAL )
//   KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_NORMAL_REACTION )
//   KRATOS_CREATE_VARIABLE(double, TAU_STAB )
//   KRATOS_CREATE_VARIABLE(double, MU_STATIC )
//   KRATOS_CREATE_VARIABLE(double, MU_DYNAMIC )


// //geometrical
//   KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( OFFSET )
//   KRATOS_CREATE_VARIABLE(Vector, BOUNDARY_NORMAL )
//   KRATOS_CREATE_VARIABLE(double, SHRINK_FACTOR )

//   //Create Variables
//   KRATOS_CREATE_VARIABLE(double, M_MODULUS )

//   KRATOS_CREATE_VARIABLE(int, PATCH_INDEX )


  KratosPfemFluidDynamicsApplication::KratosPfemFluidDynamicsApplication():   
    mTwoStepUpdatedLagrangianVPElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mTwoStepUpdatedLagrangianVPElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mTwoStepUpdatedLagrangianVPSolidElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mTwoStepUpdatedLagrangianVPSolidElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mTwoStepUpdatedLagrangianVPFluidElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mTwoStepUpdatedLagrangianVPFluidElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType(2) ) ) ),
    mWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType(3) ) ) )
  {}
  
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
    KRATOS_REGISTER_VARIABLE(M_MODULUS)  
    KRATOS_REGISTER_VARIABLE(PATCH_INDEX);

    //Register Elements
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPElement2D",mTwoStepUpdatedLagrangianVPElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPElement3D",mTwoStepUpdatedLagrangianVPElement3D);

    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPSolidElement2D",mTwoStepUpdatedLagrangianVPSolidElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPSolidElement3D",mTwoStepUpdatedLagrangianVPSolidElement3D);

    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidElement2D",mTwoStepUpdatedLagrangianVPFluidElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidElement3D",mTwoStepUpdatedLagrangianVPFluidElement3D);

    //Register Conditions
    KRATOS_REGISTER_CONDITION("WallCondition2D",mWallCondition2D);
    KRATOS_REGISTER_CONDITION("WallCondition3D",mWallCondition3D);


    //Register Constitutive Laws

    //Register Flow Rules
 
    //Register Yield Criterion
 
    //Register Hardening Laws
 
  }
  
}  // namespace Kratos.


