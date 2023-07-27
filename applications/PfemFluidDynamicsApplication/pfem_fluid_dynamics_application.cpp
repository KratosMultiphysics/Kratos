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

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"

#include "pfem_fluid_dynamics_application.h"

namespace Kratos
{

  KratosPfemFluidDynamicsApplication::KratosPfemFluidDynamicsApplication() : KratosApplication("PfemFluidDynamicsApplication"),
                                                                             mTwoStepUpdatedLagrangianVPImplicitElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node>(Element::GeometryType::PointsArrayType(6)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node>(Element::GeometryType::PointsArrayType(10)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node>(Element::GeometryType::PointsArrayType(6)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node>(Element::GeometryType::PointsArrayType(10)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitSolidElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitSolidElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node>(Element::GeometryType::PointsArrayType(6)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitSolidElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitSolidElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node>(Element::GeometryType::PointsArrayType(10)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node>(Element::GeometryType::PointsArrayType(6)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node>(Element::GeometryType::PointsArrayType(10)))),
                                                                             mUpdatedLagrangianVImplicitSolidElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mUpdatedLagrangianVImplicitSolidElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node>(Element::GeometryType::PointsArrayType(6)))),
                                                                             mUpdatedLagrangianVImplicitSolidElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mUpdatedLagrangianVImplicitSolidElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node>(Element::GeometryType::PointsArrayType(10)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitFluidFicElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitFluidFicElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node>(Element::GeometryType::PointsArrayType(6)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitFluidFicElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitFluidFicElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node>(Element::GeometryType::PointsArrayType(10)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitFluidPspgElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitFluidPspgElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node>(Element::GeometryType::PointsArrayType(6)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitFluidPspgElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitFluidPspgElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node>(Element::GeometryType::PointsArrayType(10)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node>(Element::GeometryType::PointsArrayType(6)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node>(Element::GeometryType::PointsArrayType(10)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node>(Element::GeometryType::PointsArrayType(6)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node>(Element::GeometryType::PointsArrayType(10)))),
                                                                             mUpdatedLagrangianElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mUpdatedLagrangianElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node>(Element::GeometryType::PointsArrayType(6)))),
                                                                             mUpdatedLagrangianElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mUpdatedLagrangianElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node>(Element::GeometryType::PointsArrayType(10)))),
                                                                             mTwoStepUpdatedLagrangianElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mTwoStepUpdatedLagrangianElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node>(Element::GeometryType::PointsArrayType(6)))),
                                                                             mTwoStepUpdatedLagrangianElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mTwoStepUpdatedLagrangianElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node>(Element::GeometryType::PointsArrayType(10)))),
                                                                             mThreeStepUpdatedLagrangianElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mThreeStepUpdatedLagrangianElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node>(Element::GeometryType::PointsArrayType(6)))),
                                                                             mThreeStepUpdatedLagrangianElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mThreeStepUpdatedLagrangianElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node>(Element::GeometryType::PointsArrayType(10)))),
                                                                             mThreeStepFirstOrderUpdatedLagrangianElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mThreeStepFirstOrderUpdatedLagrangianElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node>(Element::GeometryType::PointsArrayType(6)))),
                                                                             mThreeStepFirstOrderUpdatedLagrangianElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mThreeStepFirstOrderUpdatedLagrangianElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node>(Element::GeometryType::PointsArrayType(10)))),
                                                                             mThreeStepSecondOrderUpdatedLagrangianElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mThreeStepSecondOrderUpdatedLagrangianElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node>(Element::GeometryType::PointsArrayType(6)))),
                                                                             mThreeStepSecondOrderUpdatedLagrangianElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mThreeStepSecondOrderUpdatedLagrangianElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node>(Element::GeometryType::PointsArrayType(10)))),
                                                                             mThreeStepSecondOrderPspgUpdatedLagrangianElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
                                                                             mThreeStepSecondOrderPspgUpdatedLagrangianElement2Dquadratic(0, Element::GeometryType::Pointer(new Triangle2D6<Node>(Element::GeometryType::PointsArrayType(6)))),
                                                                             mThreeStepSecondOrderPspgUpdatedLagrangianElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
                                                                             mThreeStepSecondOrderPspgUpdatedLagrangianElement3Dquadratic(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node>(Element::GeometryType::PointsArrayType(10))))
  {
  }

  void KratosPfemFluidDynamicsApplication::Register()
  {
    std::cout << "            ___  __           ___ _      _    _          " << std::endl;
    std::cout << "     KRATOS| _ \\/ _|___ _ __ | __| |_  _(_)__| |         " << std::endl;
    std::cout << "           |  _/  _/ -_) '  \\| _|| | || | / _` |         " << std::endl;
    std::cout << "           |_| |_| \\___|_|_|_|_| |_|\\_,_|_\\__,_|DYNAMICS " << std::endl;
    std::cout << "Initializing KratosPfemFluidDynamicsApplication...       " << std::endl;

    // Register Variables (variables created in pfem_fluid_dynamics_application_variables.cpp)

    // Material postprocess + invariants
    // KRATOS_REGISTER_VARIABLE(M_MODULUS)
    // KRATOS_REGISTER_VARIABLE(PATCH_INDEX);
    // KRATOS_REGISTER_VARIABLE(NORMVELOCITY);
    KRATOS_REGISTER_VARIABLE(NO_MESH);
    KRATOS_REGISTER_VARIABLE(INITIAL_DELTA_TIME);
    KRATOS_REGISTER_VARIABLE(CURRENT_DELTA_TIME);
    KRATOS_REGISTER_VARIABLE(TIME_INTERVAL_CHANGED);
    KRATOS_REGISTER_VARIABLE(BAD_VELOCITY_CONVERGENCE);
    KRATOS_REGISTER_VARIABLE(BAD_PRESSURE_CONVERGENCE);
    KRATOS_REGISTER_VARIABLE(STEPS_WITH_CHANGED_DT);
    KRATOS_REGISTER_VARIABLE(MAIN_MATERIAL_PROPERTY);
    KRATOS_REGISTER_VARIABLE(ISOLATED_NODE);
    KRATOS_REGISTER_VARIABLE(NODAL_H_WALL);

    KRATOS_REGISTER_VARIABLE(MECHANICAL_DISSIPATION);

    // Papanastasiou variables
    KRATOS_REGISTER_VARIABLE(YIELDED);
    KRATOS_REGISTER_VARIABLE(FLOW_INDEX);
    KRATOS_REGISTER_VARIABLE(YIELD_SHEAR);
    KRATOS_REGISTER_VARIABLE(ADAPTIVE_EXPONENT);

    // Frictional Viscoplastic variables
    KRATOS_REGISTER_VARIABLE(COHESION);

    // mu(I)-rheology variables
    KRATOS_REGISTER_VARIABLE(STATIC_FRICTION);
    KRATOS_REGISTER_VARIABLE(DYNAMIC_FRICTION);
    KRATOS_REGISTER_VARIABLE(INERTIAL_NUMBER_ZERO);
    KRATOS_REGISTER_VARIABLE(GRAIN_DIAMETER);
    KRATOS_REGISTER_VARIABLE(GRAIN_DENSITY);
    KRATOS_REGISTER_VARIABLE(REGULARIZATION_COEFFICIENT);

    KRATOS_REGISTER_VARIABLE(PRESSURE_VELOCITY);
    KRATOS_REGISTER_VARIABLE(PRESSURE_ACCELERATION);

    KRATOS_REGISTER_VARIABLE(NODAL_ERROR_XX);

    KRATOS_REGISTER_VARIABLE(NODAL_CAUCHY_STRESS);
    KRATOS_REGISTER_VARIABLE(NODAL_DEVIATORIC_CAUCHY_STRESS);
    KRATOS_REGISTER_VARIABLE(NODAL_SFD_NEIGHBOURS);
    KRATOS_REGISTER_VARIABLE(NODAL_SFD_NEIGHBOURS_ORDER);
    KRATOS_REGISTER_VARIABLE(NODAL_DEFORMATION_GRAD);
    KRATOS_REGISTER_VARIABLE(NODAL_DEFORMATION_GRAD_VEL);
    KRATOS_REGISTER_VARIABLE(NODAL_SPATIAL_DEF_RATE);
    KRATOS_REGISTER_VARIABLE(NODAL_VOLUMETRIC_DEF_RATE);
    KRATOS_REGISTER_VARIABLE(NODAL_EQUIVALENT_STRAIN_RATE);
    KRATOS_REGISTER_VARIABLE(NODAL_MEAN_MESH_SIZE);
    KRATOS_REGISTER_VARIABLE(NODAL_TAU);
    KRATOS_REGISTER_VARIABLE(NODAL_FREESURFACE_AREA);
    KRATOS_REGISTER_VARIABLE(VOLUMETRIC_COEFFICIENT);
    KRATOS_REGISTER_VARIABLE(DEVIATORIC_COEFFICIENT);
    KRATOS_REGISTER_VARIABLE(INTERFACE_NODE);

    KRATOS_REGISTER_VARIABLE(SOLID_NODAL_VOLUME);
    KRATOS_REGISTER_VARIABLE(SOLID_NODAL_CAUCHY_STRESS);
    KRATOS_REGISTER_VARIABLE(SOLID_NODAL_DEVIATORIC_CAUCHY_STRESS);
    KRATOS_REGISTER_VARIABLE(SOLID_NODAL_SFD_NEIGHBOURS);
    KRATOS_REGISTER_VARIABLE(SOLID_NODAL_SFD_NEIGHBOURS_ORDER);
    KRATOS_REGISTER_VARIABLE(SOLID_NODAL_DEFORMATION_GRAD);
    KRATOS_REGISTER_VARIABLE(SOLID_NODAL_DEFORMATION_GRAD_VEL);
    KRATOS_REGISTER_VARIABLE(SOLID_NODAL_SPATIAL_DEF_RATE);
    KRATOS_REGISTER_VARIABLE(SOLID_NODAL_VOLUMETRIC_DEF_RATE);
    KRATOS_REGISTER_VARIABLE(SOLID_NODAL_EQUIVALENT_STRAIN_RATE);
    KRATOS_REGISTER_VARIABLE(SOLID_NODAL_MEAN_MESH_SIZE);
    KRATOS_REGISTER_VARIABLE(SOLID_DENSITY);
    KRATOS_REGISTER_VARIABLE(SOLID_NODAL_FREESURFACE_AREA);

    // Register Elements
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPImplicitElement2D", mTwoStepUpdatedLagrangianVPImplicitElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPImplicitElement2Dquadratic", mTwoStepUpdatedLagrangianVPImplicitElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPImplicitElement3D", mTwoStepUpdatedLagrangianVPImplicitElement3D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPImplicitElement3Dquadratic", mTwoStepUpdatedLagrangianVPImplicitElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement2D", mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement2Dquadratic", mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement3D", mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement3D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement3Dquadratic", mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPSolidElement2D", mTwoStepUpdatedLagrangianVPImplicitSolidElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPSolidElement2Dquadratic", mTwoStepUpdatedLagrangianVPImplicitSolidElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPSolidElement3D", mTwoStepUpdatedLagrangianVPImplicitSolidElement3D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPSolidElement3Dquadratic", mTwoStepUpdatedLagrangianVPImplicitSolidElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianVSolidElement2D", mUpdatedLagrangianVImplicitSolidElement2D);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianVSolidElement2Dquadratic", mUpdatedLagrangianVImplicitSolidElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianVSolidElement3D", mUpdatedLagrangianVImplicitSolidElement3D);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianVSolidElement3Dquadratic", mUpdatedLagrangianVImplicitSolidElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPNodallyIntegratedSolidElement2D", mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPNodallyIntegratedSolidElement2Dquadratic", mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPNodallyIntegratedSolidElement3D", mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement3D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPNodallyIntegratedSolidElement3Dquadratic", mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidElement2D", mTwoStepUpdatedLagrangianVPImplicitFluidFicElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidElement2Dquadratic", mTwoStepUpdatedLagrangianVPImplicitFluidFicElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidElement3D", mTwoStepUpdatedLagrangianVPImplicitFluidFicElement3D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidElement3Dquadratic", mTwoStepUpdatedLagrangianVPImplicitFluidFicElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidCutFemElement2D", mTwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidCutFemElement3D", mTwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement3D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidPspgElement2D", mTwoStepUpdatedLagrangianVPImplicitFluidPspgElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidPspgElement2Dquadratic", mTwoStepUpdatedLagrangianVPImplicitFluidPspgElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidPspgElement3D", mTwoStepUpdatedLagrangianVPImplicitFluidPspgElement3D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidPspgElement3Dquadratic", mTwoStepUpdatedLagrangianVPImplicitFluidPspgElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidDEMcouplingElement2D", mTwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidDEMcouplingElement2Dquadratic", mTwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidDEMcouplingElement3D", mTwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement3D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPFluidDEMcouplingElement3Dquadratic", mTwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPNodallyIntegratedFluidElement2D", mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPNodallyIntegratedFluidElement2Dquadratic", mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPNodallyIntegratedFluidElement3D", mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement3D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianVPNodallyIntegratedFluidElement3Dquadratic", mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement2D", mUpdatedLagrangianElement2D);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement2Dquadratic", mUpdatedLagrangianElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D", mUpdatedLagrangianElement3D);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3Dquadratic", mUpdatedLagrangianElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianElement2D", mTwoStepUpdatedLagrangianElement2D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianElement2Dquadratic", mTwoStepUpdatedLagrangianElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianElement3D", mTwoStepUpdatedLagrangianElement3D);
    KRATOS_REGISTER_ELEMENT("TwoStepUpdatedLagrangianElement3Dquadratic", mTwoStepUpdatedLagrangianElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("ThreeStepUpdatedLagrangianElement2D", mThreeStepUpdatedLagrangianElement2D);
    KRATOS_REGISTER_ELEMENT("ThreeStepUpdatedLagrangianElement2Dquadratic", mThreeStepUpdatedLagrangianElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("ThreeStepUpdatedLagrangianElement3D", mThreeStepUpdatedLagrangianElement3D);
    KRATOS_REGISTER_ELEMENT("ThreeStepUpdatedLagrangianElement3Dquadratic", mThreeStepUpdatedLagrangianElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("ThreeStepFirstOrderUpdatedLagrangianElement2D", mThreeStepFirstOrderUpdatedLagrangianElement2D);
    KRATOS_REGISTER_ELEMENT("ThreeStepFirstOrderUpdatedLagrangianElement2Dquadratic", mThreeStepFirstOrderUpdatedLagrangianElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("ThreeStepFirstOrderUpdatedLagrangianElement3D", mThreeStepFirstOrderUpdatedLagrangianElement3D);
    KRATOS_REGISTER_ELEMENT("ThreeStepFirstOrderUpdatedLagrangianElement3Dquadratic", mThreeStepFirstOrderUpdatedLagrangianElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("ThreeStepSecondOrderUpdatedLagrangianElement2D", mThreeStepSecondOrderUpdatedLagrangianElement2D);
    KRATOS_REGISTER_ELEMENT("ThreeStepSecondOrderUpdatedLagrangianElement2Dquadratic", mThreeStepSecondOrderUpdatedLagrangianElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("ThreeStepSecondOrderUpdatedLagrangianElement3D", mThreeStepSecondOrderUpdatedLagrangianElement3D);
    KRATOS_REGISTER_ELEMENT("ThreeStepSecondOrderUpdatedLagrangianElement3Dquadratic", mThreeStepSecondOrderUpdatedLagrangianElement3Dquadratic);
    KRATOS_REGISTER_ELEMENT("ThreeStepSecondOrderPspgUpdatedLagrangianElement2D", mThreeStepSecondOrderPspgUpdatedLagrangianElement2D);
    KRATOS_REGISTER_ELEMENT("ThreeStepSecondOrderPspgUpdatedLagrangianElement2Dquadratic", mThreeStepSecondOrderPspgUpdatedLagrangianElement2Dquadratic);
    KRATOS_REGISTER_ELEMENT("ThreeStepSecondOrderPspgUpdatedLagrangianElement3D", mThreeStepSecondOrderPspgUpdatedLagrangianElement3D);
    KRATOS_REGISTER_ELEMENT("ThreeStepSecondOrderPspgUpdatedLagrangianElement3Dquadratic", mThreeStepSecondOrderPspgUpdatedLagrangianElement3Dquadratic);

    // Register Conditions

    // Register Fluid Constitutive Laws
    KRATOS_REGISTER_CONSTITUTIVE_LAW("Bingham2DLaw", mBingham2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("Bingham3DLaw", mBingham3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("BinghamTemperatureDependent2DLaw", mBinghamTemperatureDependent2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("BinghamTemperatureDependent3DLaw", mBinghamTemperatureDependent3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HerschelBulkley2DLaw", mHerschelBulkley2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HerschelBulkley3DLaw", mHerschelBulkley3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FrictionalViscoplastic2DLaw", mFrictionalViscoplastic2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FrictionalViscoplastic3DLaw", mFrictionalViscoplastic3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FrictionalViscoplasticTemperatureDependent2DLaw", mFrictionalViscoplasticTemperatureDependent2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FrictionalViscoplasticTemperatureDependent3DLaw", mFrictionalViscoplasticTemperatureDependent3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("Newtonian2DLaw", mNewtonian2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("Newtonian3DLaw", mNewtonian3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("NewtonianTemperatureDependent2DLaw", mNewtonianTemperatureDependent2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("NewtonianTemperatureDependent3DLaw", mNewtonianTemperatureDependent3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("MuIRheology2DLaw", mMuIRheology2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("MuIRheology3DLaw", mMuIRheology3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("MuIRheologyTemperatureDependent2DLaw", mMuIRheologyTemperatureDependent2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("MuIRheologyTemperatureDependent3DLaw", mMuIRheologyTemperatureDependent3DLaw);

    // Register Solid Constitutive Laws
    KRATOS_REGISTER_CONSTITUTIVE_LAW("Hypoelastic2DLaw", mHypoelastic2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("Hypoelastic3DLaw", mHypoelastic3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HypoelasticTemperatureDependent2DLaw", mHypoelasticTemperatureDependent2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HypoelasticTemperatureDependent3DLaw", mHypoelasticTemperatureDependent3DLaw);

    // Register Flow Rules

    // Register Yield Criterion

    // Register Hardening Laws
  }

} // namespace Kratos.
