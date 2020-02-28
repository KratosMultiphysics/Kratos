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
#include "custom_constitutive/inviscid_force_laws/inviscid_force_law.h"
#include "custom_constitutive/history_force_laws/history_force_law.h"
#include "custom_constitutive/vorticity_induced_lift_laws/vorticity_induced_lift_law.h"
#include "custom_constitutive/rotation_induced_lift_laws/rotation_induced_lift_law.h"
#include "custom_constitutive/steady_viscous_torque_laws/steady_viscous_torque_law.h"

namespace Kratos
{

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VECTORIAL_ERROR)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VECTORIAL_ERROR_1)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_OLD)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AVERAGED_FLUID_VELOCITY)
KRATOS_CREATE_VARIABLE(double, PERMEABILITY_1_DAY)
KRATOS_CREATE_VARIABLE(std::string, SDEM_HYDRODYNAMIC_INTERACTION_LAW_NAME)
KRATOS_CREATE_VARIABLE(HydrodynamicInteractionLaw::Pointer, SDEM_HYDRODYNAMIC_INTERACTION_LAW_POINTER)
KRATOS_CREATE_VARIABLE(std::string, SDEM_BUOYANCY_LAW_NAME)
KRATOS_CREATE_VARIABLE(std::string, SDEM_DRAG_LAW_NAME)
KRATOS_CREATE_VARIABLE(std::string, SDEM_INVISCID_FORCE_LAW_NAME)
KRATOS_CREATE_VARIABLE(std::string, SDEM_HISTORY_FORCE_LAW_NAME)
KRATOS_CREATE_VARIABLE(std::string, SDEM_VORTICITY_LIFT_LAW_NAME)
KRATOS_CREATE_VARIABLE(std::string, SDEM_ROTATION_LIFT_LAW_NAME)
KRATOS_CREATE_VARIABLE(std::string, SDEM_STEADY_VISCOUS_TORQUE_LAW_NAME)
KRATOS_CREATE_VARIABLE(BuoyancyLaw::Pointer, SDEM_BUOYANCY_LAW_POINTER)
KRATOS_CREATE_VARIABLE(DragLaw::Pointer, SDEM_DRAG_LAW_POINTER)
KRATOS_CREATE_VARIABLE(InviscidForceLaw::Pointer, SDEM_INVISCID_FORCE_LAW_POINTER)
KRATOS_CREATE_VARIABLE(HistoryForceLaw::Pointer, SDEM_HISTORY_FORCE_LAW_POINTER)
KRATOS_CREATE_VARIABLE(VorticityInducedLiftLaw::Pointer, SDEM_VORTICITY_INDUCED_LIFT_LAW_POINTER)
KRATOS_CREATE_VARIABLE(RotationInducedLiftLaw::Pointer, SDEM_ROTATION_INDUCED_LIFT_LAW_POINTER)
KRATOS_CREATE_VARIABLE(SteadyViscousTorqueLaw::Pointer, SDEM_STEADY_VISCOUS_TORQUE_LAW_POINTER)


KratosSwimmingDEMApplication::KratosSwimmingDEMApplication():
  KratosApplication("SwimmingDEMApplication"),
  mMonolithicDEMCoupled2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
  mMonolithicDEMCoupled3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
  mMonolithicDEMCoupledWeak2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
  mMonolithicDEMCoupledWeak3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
  mComputeLaplacianSimplex2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
  mComputeLaplacianSimplex3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
  mComputeMaterialDerivativeSimplex2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
  mComputeMaterialDerivativeSimplex3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
  mComputeComponentGradientSimplex2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
  mComputeComponentGradientSimplex3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
  mComputeGradientPouliot20122DEdge(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2)))),
  mComputeGradientPouliot20123DEdge(0, Element::GeometryType::Pointer(new Line3D2<Node<3> >(Element::GeometryType::PointsArrayType(2)))),
  mComputeGradientPouliot20122D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
  mComputeGradientPouliot20123D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
  mComputeVelocityLaplacianComponentSimplex2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
  mComputeVelocityLaplacianComponentSimplex3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
  mComputeVelocityLaplacianSimplex2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
  mComputeVelocityLaplacianSimplex3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
  mMonolithicDEMCoupledWallCondition2D(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType(2)))),
  mMonolithicDEMCoupledWallCondition3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType(3)))),
  mComputeLaplacianSimplexCondition2D(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType(2)))),
  mComputeLaplacianSimplexCondition3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType(3)))),
  mRigidShellElement(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
  mSphericSwimmingParticle3D(0, Element::GeometryType::Pointer(new Sphere3D1<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
  mSwimmingNanoParticle3D(0, Element::GeometryType::Pointer(new Sphere3D1<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
  mSwimmingAnalyticParticle3D(0, Element::GeometryType::Pointer(new Sphere3D1<Node<3> >(Element::GeometryType::PointsArrayType(1))))
{}

void KratosSwimmingDEMApplication::Register()
{
  // calling base class register to register Kratos components
  KratosApplication::Register();
  std::cout << "Initializing KratosSwimmingDEMApplication... " << std::endl;

  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VECTORIAL_ERROR)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VECTORIAL_ERROR_1)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_OLD)
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(AVERAGED_FLUID_VELOCITY)
  KRATOS_REGISTER_VARIABLE(PERMEABILITY_1_DAY)
  KRATOS_REGISTER_VARIABLE(SDEM_HYDRODYNAMIC_INTERACTION_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_HYDRODYNAMIC_INTERACTION_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(SDEM_BUOYANCY_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_DRAG_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_INVISCID_FORCE_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_HISTORY_FORCE_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_VORTICITY_LIFT_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_ROTATION_LIFT_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_STEADY_VISCOUS_TORQUE_LAW_NAME)
  KRATOS_REGISTER_VARIABLE(SDEM_BUOYANCY_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(SDEM_DRAG_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(SDEM_INVISCID_FORCE_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(SDEM_HISTORY_FORCE_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(SDEM_VORTICITY_INDUCED_LIFT_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(SDEM_ROTATION_INDUCED_LIFT_LAW_POINTER)
  KRATOS_REGISTER_VARIABLE(SDEM_STEADY_VISCOUS_TORQUE_LAW_POINTER)

  /* Define In Global variables.cpp */

  KRATOS_REGISTER_ELEMENT("MonolithicDEMCoupled2D", mMonolithicDEMCoupled2D)
  KRATOS_REGISTER_ELEMENT("MonolithicDEMCoupled3D", mMonolithicDEMCoupled3D)
  KRATOS_REGISTER_ELEMENT("MonolithicDEMCoupledWeak2D", mMonolithicDEMCoupledWeak2D)
  KRATOS_REGISTER_ELEMENT("MonolithicDEMCoupledWeak3D", mMonolithicDEMCoupledWeak3D)
  KRATOS_REGISTER_ELEMENT("RigidShellElement", mRigidShellElement)
  KRATOS_REGISTER_ELEMENT("SphericSwimmingParticle3D", mSphericSwimmingParticle3D)
  KRATOS_REGISTER_ELEMENT("SwimmingNanoParticle3D", mSwimmingNanoParticle3D)
  KRATOS_REGISTER_ELEMENT("SwimmingAnalyticParticle3D", mSwimmingAnalyticParticle3D)
  KRATOS_REGISTER_ELEMENT("ComputeLaplacianSimplex2D", mComputeLaplacianSimplex2D)
  KRATOS_REGISTER_ELEMENT("ComputeLaplacianSimplex3D", mComputeLaplacianSimplex3D)
  KRATOS_REGISTER_ELEMENT("ComputeMaterialDerivativeSimplex2D", mComputeMaterialDerivativeSimplex2D)
  KRATOS_REGISTER_ELEMENT("ComputeMaterialDerivativeSimplex3D", mComputeMaterialDerivativeSimplex3D)
  KRATOS_REGISTER_ELEMENT("ComputeComponentGradientSimplex2D", mComputeComponentGradientSimplex2D)
  KRATOS_REGISTER_ELEMENT("ComputeComponentGradientSimplex3D", mComputeComponentGradientSimplex3D)
  KRATOS_REGISTER_ELEMENT("ComputeGradientPouliot20122DEdge", mComputeGradientPouliot20122DEdge)
  KRATOS_REGISTER_ELEMENT("ComputeGradientPouliot20123DEdge", mComputeGradientPouliot20123DEdge)
  KRATOS_REGISTER_ELEMENT("ComputeGradientPouliot20122D", mComputeGradientPouliot20122D)
  KRATOS_REGISTER_ELEMENT("ComputeGradientPouliot20123D", mComputeGradientPouliot20123D)
  KRATOS_REGISTER_ELEMENT("ComputeVelocityLaplacianComponentSimplex2D", mComputeVelocityLaplacianComponentSimplex2D)
  KRATOS_REGISTER_ELEMENT("ComputeVelocityLaplacianComponentSimplex3D", mComputeVelocityLaplacianComponentSimplex3D)
  KRATOS_REGISTER_ELEMENT("ComputeVelocityLaplacianSimplex2D", mComputeVelocityLaplacianSimplex2D)
  KRATOS_REGISTER_ELEMENT("ComputeVelocityLaplacianSimplex3D", mComputeVelocityLaplacianSimplex3D)
  KRATOS_REGISTER_CONDITION("MonolithicDEMCoupledWallCondition2D", mMonolithicDEMCoupledWallCondition2D)
  KRATOS_REGISTER_CONDITION("MonolithicDEMCoupledWallCondition3D",mMonolithicDEMCoupledWallCondition3D)
  KRATOS_REGISTER_CONDITION("ComputeLaplacianSimplexCondition2D",  mComputeLaplacianSimplexCondition2D)
  KRATOS_REGISTER_CONDITION("ComputeLaplacianSimplexCondition3D", mComputeLaplacianSimplexCondition3D)
 }

}  // namespace Kratos.


