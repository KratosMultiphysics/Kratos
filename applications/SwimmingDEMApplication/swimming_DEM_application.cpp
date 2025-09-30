//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Guillermo Casas, gcasas@cimne.upc.edu $
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
#include "swimming_dem_application_variables.h"
#include "geometries/point_3d.h"
#include "geometries/line_3d_2.h"
#include "geometries/sphere_3d_1.h"
#include "custom_constitutive/hydrodynamic_interaction_law.h"
#include "custom_constitutive/power_law_hydrodynamic_interaction_law.h"
#include "custom_constitutive/buoyancy_laws/buoyancy_law.h"
#include "custom_constitutive/drag_laws/drag_law.h"
#include "custom_constitutive/virtual_mass_force_laws/virtual_mass_force_law.h"
#include "custom_constitutive/undisturbed_force_laws/undisturbed_force_law.h"
#include "custom_constitutive/history_force_laws/history_force_law.h"
#include "custom_constitutive/vorticity_induced_lift_laws/vorticity_induced_lift_law.h"
#include "custom_constitutive/rotation_induced_lift_laws/rotation_induced_lift_law.h"
#include "custom_constitutive/steady_viscous_torque_laws/steady_viscous_torque_law.h"
#include "containers/global_pointers_vector.h"

namespace Kratos
{

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VECTORIAL_ERROR)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VECTORIAL_ERROR_1)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_OLD)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AVERAGED_FLUID_VELOCITY)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(TIME_AVERAGED_BODY_FORCE)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(HYDRODYNAMIC_REACTION_PROJECTED)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(UNDISTURBED_FLOW_FORCE)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(EXACT_MATERIAL_ACCELERATION)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(EXACT_VECTOR_VARIABLE)
KRATOS_CREATE_VARIABLE(double, WEIGHTED_SUM)
KRATOS_CREATE_VARIABLE(double, RELAXATION_ALPHA)
KRATOS_CREATE_VARIABLE(double, EXACT_PRESSURE)
KRATOS_CREATE_VARIABLE(double, SCALAR_ERROR)
KRATOS_CREATE_VARIABLE(double, ERROR_X)
KRATOS_CREATE_VARIABLE(double, ERROR_Y)
KRATOS_CREATE_VARIABLE(double, ERROR_Z)
KRATOS_CREATE_VARIABLE(double, ERROR_P)
KRATOS_CREATE_VARIABLE(double, FLUID_FRACTION_OLD_2)
KRATOS_CREATE_VARIABLE(double, NODAL_DENSITY)
KRATOS_CREATE_VARIABLE(double, EXACT_SCALAR_VARIABLE)
KRATOS_CREATE_VARIABLE(double, NODAL_DENSITY_PROJECTED)
KRATOS_CREATE_VARIABLE(double, GENTLE_INITIATION_COUPLING_COEFFICIENT)
KRATOS_CREATE_VARIABLE(double, PERMEABILITY_1_DAY)
KRATOS_CREATE_VARIABLE(std::string, SDEM_HYDRODYNAMIC_INTERACTION_LAW_NAME)
KRATOS_CREATE_VARIABLE(bool, MANUFACTURED)
KRATOS_CREATE_VARIABLE(bool, IS_EDGE_NODE)
KRATOS_CREATE_VARIABLE(GlobalPointersVector<Node>, CONTIGUOUS_NODES)
KRATOS_CREATE_VARIABLE(DenseMatrix<double>, MATRIX_WEIGHTS)
KRATOS_CREATE_VARIABLE(HydrodynamicInteractionLaw::Pointer, SDEM_HYDRODYNAMIC_INTERACTION_LAW_POINTER)
KRATOS_CREATE_VARIABLE(std::string, SDEM_BUOYANCY_LAW_NAME)
KRATOS_CREATE_VARIABLE(std::string, SDEM_DRAG_LAW_NAME)
KRATOS_CREATE_VARIABLE(std::string, SDEM_VIRTUAL_MASS_FORCE_LAW_NAME)
KRATOS_CREATE_VARIABLE(std::string, SDEM_UNDISTURBED_FORCE_LAW_NAME)
KRATOS_CREATE_VARIABLE(std::string, SDEM_HISTORY_FORCE_LAW_NAME)
KRATOS_CREATE_VARIABLE(std::string, SDEM_VORTICITY_LIFT_LAW_NAME)
KRATOS_CREATE_VARIABLE(std::string, SDEM_ROTATION_LIFT_LAW_NAME)
KRATOS_CREATE_VARIABLE(std::string, SDEM_STEADY_VISCOUS_TORQUE_LAW_NAME)
KRATOS_CREATE_VARIABLE(BuoyancyLaw::Pointer, SDEM_BUOYANCY_LAW_POINTER)
KRATOS_CREATE_VARIABLE(DragLaw::Pointer, SDEM_DRAG_LAW_POINTER)
KRATOS_CREATE_VARIABLE(VirtualMassForceLaw::Pointer, SDEM_VIRTUAL_MASS_FORCE_LAW_POINTER)
KRATOS_CREATE_VARIABLE(UndisturbedForceLaw::Pointer, SDEM_UNDISTURBED_FORCE_LAW_POINTER)
KRATOS_CREATE_VARIABLE(HistoryForceLaw::Pointer, SDEM_HISTORY_FORCE_LAW_POINTER)
KRATOS_CREATE_VARIABLE(VorticityInducedLiftLaw::Pointer, SDEM_VORTICITY_INDUCED_LIFT_LAW_POINTER)
KRATOS_CREATE_VARIABLE(RotationInducedLiftLaw::Pointer, SDEM_ROTATION_INDUCED_LIFT_LAW_POINTER)
KRATOS_CREATE_VARIABLE(SteadyViscousTorqueLaw::Pointer, SDEM_STEADY_VISCOUS_TORQUE_LAW_POINTER)


KratosSwimmingDEMApplication::KratosSwimmingDEMApplication():
  KratosApplication("SwimmingDEMApplication"),
  mMonolithicDEMCoupled2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
  mMonolithicDEMCoupled3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
  mMonolithicDEMCoupledWeak2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
  mMonolithicDEMCoupledWeak3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
  mComputeLaplacianSimplex2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
  mComputeLaplacianSimplex3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
  mComputeLaplacianSimplex3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node >(Element::GeometryType::PointsArrayType(10)))),
  mComputeMaterialDerivativeSimplex2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
  mComputeMaterialDerivativeSimplex3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
  mComputeMaterialDerivativeSimplex3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node >(Element::GeometryType::PointsArrayType(10)))),
  mComputeFluidFractionGradient2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
  mComputeFluidFractionGradient2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node >(Element::GeometryType::PointsArrayType(4)))),
  mComputeFluidFractionGradient2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<Node >(Element::GeometryType::PointsArrayType(9)))),
  mComputeFluidFractionGradient3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
  mComputeFluidFractionGradient3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node >(Element::GeometryType::PointsArrayType(8)))),
  mComputeFluidFractionGradient3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<Node >(Element::GeometryType::PointsArrayType(27)))),
  mComputeComponentGradientSimplex2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
  mComputeComponentGradientSimplex3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
  mComputeComponentGradientSimplex3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node >(Element::GeometryType::PointsArrayType(10)))),
  mComputeGradientPouliot20122DEdge(0, Element::GeometryType::Pointer(new Line2D2<Node >(Element::GeometryType::PointsArrayType(2)))),
  mComputeGradientPouliot20123DEdge(0, Element::GeometryType::Pointer(new Line3D2<Node >(Element::GeometryType::PointsArrayType(2)))),
  mComputeGradientPouliot20122D(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
  mComputeGradientPouliot20123D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
  mComputeVelocityLaplacianComponentSimplex2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
  mComputeVelocityLaplacianComponentSimplex3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
  mComputeVelocityLaplacianSimplex2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
  mComputeVelocityLaplacianSimplex3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
  mComputeVelocityLaplacianSimplex3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node >(Element::GeometryType::PointsArrayType(10)))),
  mMonolithicDEMCoupledWallCondition2D(0, Element::GeometryType::Pointer(new Line2D2<Node >( Element::GeometryType::PointsArrayType(2)))),
  mMonolithicDEMCoupledWallCondition3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node >( Element::GeometryType::PointsArrayType(3)))),
  mComputeLaplacianSimplexCondition2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node >( Element::GeometryType::PointsArrayType(2)))),
  mComputeLaplacianSimplexCondition3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node >( Element::GeometryType::PointsArrayType(3)))),
  mComputeLaplacianSimplexCondition3D6N(0, Element::GeometryType::Pointer(new Triangle3D6<Node >( Element::GeometryType::PointsArrayType(6)))),
  mRigidShellElement(0, Element::GeometryType::Pointer(new Triangle3D3<Node >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
  mSphericSwimmingParticle3D(0, Element::GeometryType::Pointer(new Sphere3D1<Node >(Element::GeometryType::PointsArrayType(1)))),
  mSwimmingNanoParticle3D(0, Element::GeometryType::Pointer(new Sphere3D1<Node >(Element::GeometryType::PointsArrayType(1)))),
  mSwimmingAnalyticParticle3D(0, Element::GeometryType::Pointer(new Sphere3D1<Node >(Element::GeometryType::PointsArrayType(1))))
{}

void KratosSwimmingDEMApplication::Register()
{
  std::cout << "Initializing KratosSwimmingDEMApplication... " << std::endl;

  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VECTORIAL_ERROR)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VECTORIAL_ERROR_1)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_OLD)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(AVERAGED_FLUID_VELOCITY)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TIME_AVERAGED_BODY_FORCE)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(HYDRODYNAMIC_REACTION_PROJECTED)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(UNDISTURBED_FLOW_FORCE)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(EXACT_MATERIAL_ACCELERATION)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(EXACT_VECTOR_VARIABLE)
  KRATOS_REGISTER_VARIABLE(WEIGHTED_SUM)
  KRATOS_REGISTER_VARIABLE(RELAXATION_ALPHA)
  KRATOS_REGISTER_VARIABLE(EXACT_PRESSURE)
  KRATOS_REGISTER_VARIABLE(SCALAR_ERROR)
  KRATOS_REGISTER_VARIABLE(ERROR_X)
  KRATOS_REGISTER_VARIABLE(ERROR_Y)
  KRATOS_REGISTER_VARIABLE(ERROR_Z)
  KRATOS_REGISTER_VARIABLE(ERROR_P)
  KRATOS_REGISTER_VARIABLE(FLUID_FRACTION_OLD_2)
  KRATOS_REGISTER_VARIABLE(NODAL_DENSITY)
  KRATOS_REGISTER_VARIABLE(EXACT_SCALAR_VARIABLE)
  KRATOS_REGISTER_VARIABLE(NODAL_DENSITY_PROJECTED)
  KRATOS_REGISTER_VARIABLE(GENTLE_INITIATION_COUPLING_COEFFICIENT)
  KRATOS_REGISTER_VARIABLE(PERMEABILITY_1_DAY)
  KRATOS_REGISTER_VARIABLE(SDEM_HYDRODYNAMIC_INTERACTION_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_HYDRODYNAMIC_INTERACTION_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(SDEM_BUOYANCY_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_DRAG_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_VIRTUAL_MASS_FORCE_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_UNDISTURBED_FORCE_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_HISTORY_FORCE_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_VORTICITY_LIFT_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_ROTATION_LIFT_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_STEADY_VISCOUS_TORQUE_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_BUOYANCY_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(SDEM_DRAG_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(SDEM_VIRTUAL_MASS_FORCE_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(SDEM_UNDISTURBED_FORCE_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(SDEM_HISTORY_FORCE_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(SDEM_VORTICITY_INDUCED_LIFT_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(SDEM_ROTATION_INDUCED_LIFT_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(SDEM_STEADY_VISCOUS_TORQUE_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(MANUFACTURED)
  KRATOS_REGISTER_VARIABLE(IS_EDGE_NODE)
  KRATOS_REGISTER_VARIABLE(CONTIGUOUS_NODES)
  KRATOS_REGISTER_VARIABLE(MATRIX_WEIGHTS)

  /* Define In Global variables.cpp */

  KRATOS_REGISTER_ELEMENT("MonolithicDEMCoupled2D", mMonolithicDEMCoupled2D)
  KRATOS_REGISTER_ELEMENT("MonolithicDEMCoupled3D", mMonolithicDEMCoupled3D)
  KRATOS_REGISTER_ELEMENT("MonolithicDEMCoupledWeak2D", mMonolithicDEMCoupledWeak2D)
  KRATOS_REGISTER_ELEMENT("MonolithicDEMCoupledWeak3D", mMonolithicDEMCoupledWeak3D)
  KRATOS_REGISTER_ELEMENT("ComputeLaplacianSimplex2D3N", mComputeLaplacianSimplex2D3N)
  KRATOS_REGISTER_ELEMENT("ComputeLaplacianSimplex3D4N", mComputeLaplacianSimplex3D4N)
  KRATOS_REGISTER_ELEMENT("ComputeLaplacianSimplex3D10N", mComputeLaplacianSimplex3D10N)
  KRATOS_REGISTER_ELEMENT("RigidShellElement", mRigidShellElement)
  KRATOS_REGISTER_ELEMENT("SphericSwimmingParticle3D", mSphericSwimmingParticle3D)
  KRATOS_REGISTER_ELEMENT("SwimmingNanoParticle3D", mSwimmingNanoParticle3D)
  KRATOS_REGISTER_ELEMENT("SwimmingAnalyticParticle3D", mSwimmingAnalyticParticle3D)
  KRATOS_REGISTER_ELEMENT("ComputeMaterialDerivativeSimplex2D3N", mComputeMaterialDerivativeSimplex2D3N)
  KRATOS_REGISTER_ELEMENT("ComputeMaterialDerivativeSimplex3D4N", mComputeMaterialDerivativeSimplex3D4N)
  KRATOS_REGISTER_ELEMENT("ComputeFluidFractionGradient2D3N", mComputeFluidFractionGradient2D3N)
  KRATOS_REGISTER_ELEMENT("ComputeFluidFractionGradient2D4N", mComputeFluidFractionGradient2D4N)
  KRATOS_REGISTER_ELEMENT("ComputeFluidFractionGradient2D9N", mComputeFluidFractionGradient2D9N)
  KRATOS_REGISTER_ELEMENT("ComputeFluidFractionGradient3D4N", mComputeFluidFractionGradient3D4N)
  KRATOS_REGISTER_ELEMENT("ComputeFluidFractionGradient3D8N", mComputeFluidFractionGradient3D8N)
  KRATOS_REGISTER_ELEMENT("ComputeFluidFractionGradient3D27N", mComputeFluidFractionGradient3D27N)
  KRATOS_REGISTER_ELEMENT("ComputeMaterialDerivativeSimplex2D3N", mComputeMaterialDerivativeSimplex2D3N)
  KRATOS_REGISTER_ELEMENT("ComputeMaterialDerivativeSimplex3D4N", mComputeMaterialDerivativeSimplex3D4N)
  KRATOS_REGISTER_ELEMENT("ComputeMaterialDerivativeSimplex3D10N", mComputeMaterialDerivativeSimplex3D10N)
  KRATOS_REGISTER_ELEMENT("ComputeComponentGradientSimplex2D3N", mComputeComponentGradientSimplex2D3N)
  KRATOS_REGISTER_ELEMENT("ComputeComponentGradientSimplex3D4N", mComputeComponentGradientSimplex3D4N)
  KRATOS_REGISTER_ELEMENT("ComputeComponentGradientSimplex3D10N", mComputeComponentGradientSimplex3D10N)
  KRATOS_REGISTER_ELEMENT("ComputeGradientPouliot20122DEdge", mComputeGradientPouliot20122DEdge)
  KRATOS_REGISTER_ELEMENT("ComputeGradientPouliot20123DEdge", mComputeGradientPouliot20123DEdge)
  KRATOS_REGISTER_ELEMENT("ComputeGradientPouliot20122D", mComputeGradientPouliot20122D)
  KRATOS_REGISTER_ELEMENT("ComputeGradientPouliot20123D", mComputeGradientPouliot20123D)
  KRATOS_REGISTER_ELEMENT("ComputeVelocityLaplacianComponentSimplex2D", mComputeVelocityLaplacianComponentSimplex2D)
  KRATOS_REGISTER_ELEMENT("ComputeVelocityLaplacianComponentSimplex3D", mComputeVelocityLaplacianComponentSimplex3D)
  KRATOS_REGISTER_ELEMENT("ComputeVelocityLaplacianSimplex2D3N", mComputeVelocityLaplacianSimplex2D3N)
  KRATOS_REGISTER_ELEMENT("ComputeVelocityLaplacianSimplex3D4N", mComputeVelocityLaplacianSimplex3D4N)
  KRATOS_REGISTER_ELEMENT("ComputeVelocityLaplacianSimplex3D10N", mComputeVelocityLaplacianSimplex3D10N)
  KRATOS_REGISTER_CONDITION("MonolithicDEMCoupledWallCondition2D", mMonolithicDEMCoupledWallCondition2D)
  KRATOS_REGISTER_CONDITION("MonolithicDEMCoupledWallCondition3D",mMonolithicDEMCoupledWallCondition3D)
  KRATOS_REGISTER_CONDITION("ComputeLaplacianSimplexCondition2D2N", mComputeLaplacianSimplexCondition2D2N)
  KRATOS_REGISTER_CONDITION("ComputeLaplacianSimplexCondition3D3N", mComputeLaplacianSimplexCondition3D3N)
  KRATOS_REGISTER_CONDITION("ComputeLaplacianSimplexCondition3D6N", mComputeLaplacianSimplexCondition3D6N)
}

}  // namespace Kratos.


