//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Riccardo Rossi
//


// System includes


// External includes


// Project includes
#include "includes/define.h"

#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"

#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "includes/condition.h"

//#include "includes/kernel.h"
//#include "includes/kratos_flags.h"
//#include "containers/flags.h"

#include "geometries/geometry.h"


namespace Kratos {

KratosIGAStructuralMechanicsApplication::KratosIGAStructuralMechanicsApplication() :
    KratosApplication("IGAStructuralMechanicsApplication"),
    mMeshlessElement(0, Element::GeometryType::Pointer(new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1)))),
    mTrussDiscreteElement(0, Element::GeometryType::Pointer(new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1)))),
    mMembraneDiscreteElement(0, Element::GeometryType::Pointer(new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1)))),
    mShellKLDiscreteElement(0, Element::GeometryType::Pointer(new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1)))),
    mMeshlessMembraneElement(0, Element::GeometryType::Pointer(new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1)))),
    //mMeshlessLaplaceElement(0, Element::GeometryType::Pointer(new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1)))),
    mMeshlessShellElement(0, Element::GeometryType::Pointer(new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1)))),
    mMeshlessShellKLElement(0, Element::GeometryType::Pointer(new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1)))),
    mMeshlessShellKLThickElement(0, Element::GeometryType::Pointer(new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1)))),
    //mMeshlessShellKLElement(0, Element::GeometryType::Pointer(new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1)))),
    //mContinuityConditionLagrange(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mLoadPointDiscreteCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mLoadCurveDiscreteCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mLoadSurfaceDiscreteCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mSupportPenaltyCurveDiscreteCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mSupportPenaltyPointDiscreteCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mSupportStrongDiscreteCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mMeshlessLagrangeCouplingCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mMeshlessLagrangeCouplingCondition2(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    //mContinuityConditionPenalty(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    //mLoadCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    //mSupportCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mMeshlessSupportRotationCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mMeshlessSurfaceSupportCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mMeshlessPenaltyCouplingRotationCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mMeshlessPenaltyCouplingCrackCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mMeshlessForceInterfaceCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1))))
{}

void KratosIGAStructuralMechanicsApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "     KRATOS .___  ________    _____     " << std::endl;
    std::cout << "            |   |/  _____/   /  _  \\    " << std::endl;
    std::cout << "            |   /   \\  ___  /  /_\\  \\   " << std::endl;
    std::cout << "            |   \\    \\_\\  \\/    |    \\  " << std::endl;
    std::cout << "            |___|\\______  /\\____|__  /  " << std::endl;
    std::cout << "                        \\/         \\/ STRUCTURAL MECHANICS" << std::endl;
    std::cout << "     Initializing KratosIGAStructuralMechanicsApplication... " << std::endl;

	//IGA-DEM
KRATOS_REGISTER_VARIABLE(COORDINATES)
KRATOS_REGISTER_VARIABLE(SURFACE_NORMAL)
//KRATOS_REGISTER_VARIABLE( INTEGRATION_WEIGHT)

KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_VALUES)
KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_LOCAL_DERIVATIVES)
KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES)

//KRATOS_REGISTER_VARIABLE( SHAPE_FUNCTION_LOCAL_DERIVATIVES_MASTER)
KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_VALUES_SLAVE)
KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE)
//KRATOS_REGISTER_VARIABLE(SHAPE_FUNCTION_DERIVATIVES_SLAVE)

// edge integration
KRATOS_REGISTER_VARIABLE(TANGENTS)
KRATOS_REGISTER_VARIABLE(TANGENTS_SLAVE)

// penalty factor
KRATOS_REGISTER_VARIABLE(PENALTY_FACTOR)

// coupling and support
KRATOS_REGISTER_VARIABLE(DISPLACEMENT_ROTATION_FIX)
// for load condition
KRATOS_REGISTER_VARIABLE(LOAD_TYPE)
KRATOS_REGISTER_VARIABLE(DISTRIBUTED_LOAD_FACTOR)

KRATOS_REGISTER_VARIABLE(PRESTRESS_CAUCHY)
KRATOS_REGISTER_VARIABLE(MEMBRANE_PRESTRESS_TENSOR_PK2)

KRATOS_REGISTER_VARIABLE(PRINCIPAL_STRESSES)
KRATOS_REGISTER_VARIABLE(PRINCIPAL_FORCES)

// Register the meshless element
KRATOS_REGISTER_ELEMENT("MeshlessElement", mMeshlessElement)
KRATOS_REGISTER_ELEMENT("MeshlessMembraneElement", mMeshlessMembraneElement)
//KRATOS_REGISTER_ELEMENT("MeshlessLaplaceElement", mMeshlessLaplaceElement)
KRATOS_REGISTER_ELEMENT("MeshlessShellElement", mMeshlessShellElement)
KRATOS_REGISTER_ELEMENT("MeshlessShellKLElement", mMeshlessShellKLElement)
KRATOS_REGISTER_ELEMENT("MeshlessShellKLThickElement", mMeshlessShellKLThickElement)

// Register discrete integration point elements
KRATOS_REGISTER_ELEMENT("TrussDiscreteElement", mTrussDiscreteElement)
KRATOS_REGISTER_ELEMENT("MembraneDiscreteElement", mMembraneDiscreteElement)
KRATOS_REGISTER_ELEMENT("ShellKLDiscreteElement", mShellKLDiscreteElement)

// Register discrete integration point conditions
KRATOS_REGISTER_CONDITION("LoadPointDiscreteCondition", mLoadPointDiscreteCondition)
KRATOS_REGISTER_CONDITION("LoadCurveDiscreteCondition", mLoadCurveDiscreteCondition)
KRATOS_REGISTER_CONDITION("LoadSurfaceDiscreteCondition", mLoadSurfaceDiscreteCondition)
KRATOS_REGISTER_CONDITION("SupportPenaltyCurveDiscreteCondition", mSupportPenaltyCurveDiscreteCondition)
KRATOS_REGISTER_CONDITION("SupportPenaltyPointDiscreteCondition", mSupportPenaltyPointDiscreteCondition)
KRATOS_REGISTER_CONDITION("SupportStrongDiscreteCondition", mSupportStrongDiscreteCondition)

// Register meshless condition
KRATOS_REGISTER_CONDITION("MeshlessSupportRotationCondition", mMeshlessSupportRotationCondition)
KRATOS_REGISTER_CONDITION("MeshlessSurfaceSupportCondition", mMeshlessSurfaceSupportCondition)
KRATOS_REGISTER_CONDITION("MeshlessLagrangeCouplingCondition", mMeshlessLagrangeCouplingCondition)
KRATOS_REGISTER_CONDITION("MeshlessLagrangeCouplingCondition2", mMeshlessLagrangeCouplingCondition2)
KRATOS_REGISTER_CONDITION("MeshlessPenaltyCouplingRotationCondition", mMeshlessPenaltyCouplingRotationCondition)
KRATOS_REGISTER_CONDITION("MeshlessPenaltyCouplingCrackCondition", mMeshlessPenaltyCouplingCrackCondition)
KRATOS_REGISTER_CONDITION("MeshlessForceInterfaceCondition", mMeshlessForceInterfaceCondition)

KRATOS_REGISTER_CONSTITUTIVE_LAW("PlaneStress2dKinematicallyEnrichedLaw", mPlaneStress2dKinematicallyEnrichedLaw);
KRATOS_REGISTER_CONSTITUTIVE_LAW("PlaneStress2dTCDamageLaw", mPlaneStress2dTCDamageLaw);


KRATOS_REGISTER_VARIABLE(UNIAXIAL_COMPRESSIVE_STRENGTH)
KRATOS_REGISTER_VARIABLE(UNIAXIAL_TENSILE_STRENGTH)
KRATOS_REGISTER_VARIABLE(RATE_BIAXIAL_UNIAXIAL)

KRATOS_REGISTER_VARIABLE(COMPRESSION_PARAMETER_A)
KRATOS_REGISTER_VARIABLE(COMPRESSION_PARAMETER_B)
KRATOS_REGISTER_VARIABLE(TENSION_PARAMETER_A)

KRATOS_REGISTER_VARIABLE(BETA)
KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_TENSION)
KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_COMPRESSION)

// for damage constitutive law
KRATOS_REGISTER_VARIABLE(DAMAGE_T)
KRATOS_REGISTER_VARIABLE(DAMAGE_T_INSIDE)
KRATOS_REGISTER_VARIABLE(DAMAGE_T_OUTSIDE)
KRATOS_REGISTER_VARIABLE(DAMAGE_C)
KRATOS_REGISTER_VARIABLE(DAMAGE_C_INSIDE)
KRATOS_REGISTER_VARIABLE(DAMAGE_C_OUTSIDE)


KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_T)
KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_C)
KRATOS_REGISTER_VARIABLE(YIELD_STRESS_T)
KRATOS_REGISTER_VARIABLE(YIELD_STRESS_C)
KRATOS_REGISTER_VARIABLE(DAMAGE_STRESS_T_0)
KRATOS_REGISTER_VARIABLE(DAMAGE_STRESS_C_0)
KRATOS_REGISTER_VARIABLE(DAMAGE_STRESS_C_P)
KRATOS_REGISTER_VARIABLE(DAMAGE_STRESS_C_M)
KRATOS_REGISTER_VARIABLE(DAMAGE_STRESS_C_R)
KRATOS_REGISTER_VARIABLE(DAMAGE_STRAIN_C_P)
KRATOS_REGISTER_VARIABLE(DAMAGE_STRAIN_C_M)
KRATOS_REGISTER_VARIABLE(DAMAGE_COMPRESSIVE_LAW_C1)
KRATOS_REGISTER_VARIABLE(DAMAGE_COMPRESSIVE_LAW_C2)
KRATOS_REGISTER_VARIABLE(DAMAGE_COMPRESSIVE_LAW_C3)
KRATOS_REGISTER_VARIABLE(BIAXIAL_COMPRESSION_MULTIPLIER)
KRATOS_REGISTER_VARIABLE(SHEAR_COMPRESSION_REDUCTION)
KRATOS_REGISTER_VARIABLE(DAMAGE_TENSILE_SURFACE_S1)
KRATOS_REGISTER_VARIABLE(LUBLINER_SURFACE_PARAM_KC)
KRATOS_REGISTER_VARIABLE(GENRANKINE_SURFACE_PARAM_A)
KRATOS_REGISTER_VARIABLE(GENRANKINE_SURFACE_PARAM_B)
KRATOS_REGISTER_VARIABLE(GENRANKINE_SURFACE_PARAM_C)
KRATOS_REGISTER_VARIABLE(DAMAGE_SECANT_MATRIX)
KRATOS_REGISTER_VARIABLE(DAMAGE_MODEL)
KRATOS_REGISTER_VARIABLE(DAMAGE_TENSILE_MODEL)


//// Register outdated conditions
//KRATOS_REGISTER_CONDITION("LoadCondition", mLoadCondition)
//KRATOS_REGISTER_CONDITION("SupportCondition", mSupportCondition)
//KRATOS_REGISTER_CONDITION("ContinuityConditionPenalty", mContinuityConditionPenalty)
//KRATOS_REGISTER_CONDITION("ContinuityConditionLagrange", mContinuityConditionLagrange)
  }
}  // namespace Kratos.