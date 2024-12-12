//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application

// System includes

// External includes

// Project includes
#include "iga_application.h"
#include "iga_application_variables.h"

namespace Kratos {

KratosIgaApplication::KratosIgaApplication()
    : KratosApplication("IgaApplication")
    , mTrussElement(0, Element::GeometryType::Pointer(
        new Geometry<Node>(Element::GeometryType::PointsArrayType(1))))
    , mTrussEmbeddedEdgeElement(0, Element::GeometryType::Pointer(
        new Geometry<Node>(Element::GeometryType::PointsArrayType(1))))
    , mIgaMembraneElement(0, Element::GeometryType::Pointer(
        new Geometry<Node>(Element::GeometryType::PointsArrayType(1))))
    , mShell3pElement(0, Element::GeometryType::Pointer(
        new Geometry<Node>(Element::GeometryType::PointsArrayType(1))))
    , mShell5pHierarchicElement(0, Element::GeometryType::Pointer(
        new Geometry<Node>(Element::GeometryType::PointsArrayType(1))))
    , mShell5pElement(0, Element::GeometryType::Pointer(
        new Geometry<Node>(Element::GeometryType::PointsArrayType(1))))
    , mLaplacianIGAElement(0, Element::GeometryType::Pointer(
        new Geometry<Node>(Element::GeometryType::PointsArrayType(1))))
    , mOutputCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node>(Condition::GeometryType::PointsArrayType(1))))
    , mLoadCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node>(Condition::GeometryType::PointsArrayType(1))))
    , mLoadMomentDirector5pCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node>(Condition::GeometryType::PointsArrayType(1))))
    , mCouplingPenaltyCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node>(Condition::GeometryType::PointsArrayType(1))))
    , mCouplingLagrangeCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node>(Condition::GeometryType::PointsArrayType(1))))
    , mCouplingNitscheCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node>(Condition::GeometryType::PointsArrayType(1))))
    , mSupportPenaltyCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node>(Condition::GeometryType::PointsArrayType(1))))
    , mSupportLagrangeCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node>(Condition::GeometryType::PointsArrayType(1))))
    , mSupportNitscheCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node>(Condition::GeometryType::PointsArrayType(1))))
    , mSupportLaplacianCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node>(Condition::GeometryType::PointsArrayType(1))))
    , mSBMLaplacianCondition(0, Condition::GeometryType::Pointer(
        new Geometry<Node>(Condition::GeometryType::PointsArrayType(1))))
{
}

void KratosIgaApplication::Register() {

KRATOS_INFO("") << "    KRATOS  _____ _____\n"
                << "           |_   _/ ____|   /\\\n"
                << "             | || |  __   /  \\\n"
                << "             | || | |_ | / /\\ \\\n"
                << "            _| || |__| |/ ____ \\\n"
                << "           |_____\\_____/_/    \\_\\\n"
                << "Initializing KratosIgaApplication..." << std::endl;

    // ELEMENTS
    KRATOS_REGISTER_ELEMENT("TrussElement", mTrussElement)
    KRATOS_REGISTER_ELEMENT("TrussEmbeddedEdgeElement", mTrussEmbeddedEdgeElement)
    KRATOS_REGISTER_ELEMENT("IgaMembraneElement", mIgaMembraneElement)
    KRATOS_REGISTER_ELEMENT("Shell3pElement", mShell3pElement)
    KRATOS_REGISTER_ELEMENT("Shell5pHierarchicElement", mShell5pHierarchicElement)
    KRATOS_REGISTER_ELEMENT("Shell5pElement", mShell5pElement)
    KRATOS_REGISTER_ELEMENT("LaplacianIGAElement", mLaplacianIGAElement)

    // CONDITIONS
    KRATOS_REGISTER_CONDITION("OutputCondition", mOutputCondition)
    KRATOS_REGISTER_CONDITION("LoadCondition", mLoadCondition)
    KRATOS_REGISTER_CONDITION("LoadMomentDirector5pCondition", mLoadMomentDirector5pCondition)
    KRATOS_REGISTER_CONDITION("CouplingPenaltyCondition", mCouplingPenaltyCondition)
    KRATOS_REGISTER_CONDITION("CouplingLagrangeCondition", mCouplingLagrangeCondition)
    KRATOS_REGISTER_CONDITION("CouplingNitscheCondition", mCouplingNitscheCondition)
    KRATOS_REGISTER_CONDITION("SupportPenaltyCondition", mSupportPenaltyCondition)
    KRATOS_REGISTER_CONDITION("SupportLagrangeCondition", mSupportLagrangeCondition)
    KRATOS_REGISTER_CONDITION("SupportNitscheCondition", mSupportNitscheCondition)
    KRATOS_REGISTER_CONDITION("SupportLaplacianCondition", mSupportLaplacianCondition)
    KRATOS_REGISTER_CONDITION("SBMLaplacianCondition", mSBMLaplacianCondition)

    KRATOS_REGISTER_MODELER("IgaModeler", mIgaModeler);
    KRATOS_REGISTER_MODELER("IgaThbModeler", mIgaThbModeler);
    KRATOS_REGISTER_MODELER("RefinementModeler", mRefinementModeler);
    KRATOS_REGISTER_MODELER("NurbsGeometryModeler", mNurbsGeometryModeler);

    // VARIABLES
    KRATOS_REGISTER_VARIABLE(CROSS_AREA)
      
    KRATOS_REGISTER_VARIABLE(TRUSS_PRESTRESS_CAUCHY)
    KRATOS_REGISTER_VARIABLE(TRUSS_PRESTRESS_PK2)
    KRATOS_REGISTER_VARIABLE(TRUSS_STRESS_CAUCHY)
    KRATOS_REGISTER_VARIABLE(TRUSS_STRESS_PK2)
    KRATOS_REGISTER_VARIABLE(TRUSS_FORCE)
    KRATOS_REGISTER_VARIABLE(TANGENT_MODULUS)
    KRATOS_REGISTER_VARIABLE(TRUSS_GREEN_LAGRANGE_STRAIN)

    // Structural Mechanics Application variables
    KRATOS_REGISTER_VARIABLE(NODAL_INERTIA)
    KRATOS_REGISTER_VARIABLE(PRESTRESS_CAUCHY)
    KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(PRESTRESS)
    KRATOS_REGISTER_VARIABLE(TANGENTS)

    KRATOS_REGISTER_VARIABLE(FORCE_PK2_1D)
    KRATOS_REGISTER_VARIABLE(FORCE_CAUCHY_1D)
    KRATOS_REGISTER_VARIABLE(PRINCIPAL_STRESS_1)
    KRATOS_REGISTER_VARIABLE(PRINCIPAL_STRESS_2)

    KRATOS_REGISTER_VARIABLE(LOCAL_ELEMENT_ORIENTATION)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LOCAL_PRESTRESS_AXIS_1)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LOCAL_PRESTRESS_AXIS_2)

    KRATOS_REGISTER_VARIABLE(RAYLEIGH_ALPHA)
    KRATOS_REGISTER_VARIABLE(RAYLEIGH_BETA)

    /// 5p Director Shell Variables
    KRATOS_REGISTER_VARIABLE(DIRECTOR_COMPUTED)
    KRATOS_REGISTER_VARIABLE(DIRECTOR)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DIRECTORINC)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MOMENTDIRECTORINC)
    KRATOS_REGISTER_VARIABLE(DIRECTORTANGENTSPACE)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MOMENT_LINE_LOAD)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(POINT_LOAD)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DEAD_LOAD)
    KRATOS_REGISTER_VARIABLE(PRESSURE_FOLLOWER_LOAD)

    KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(PK2_STRESS)
    KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(CAUCHY_STRESS)
    KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(CAUCHY_STRESS_TOP)
    KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(CAUCHY_STRESS_BOTTOM)
    KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(MEMBRANE_FORCE)
    KRATOS_REGISTER_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(INTERNAL_MOMENT)

    KRATOS_REGISTER_VARIABLE(SHEAR_FORCE_1)
    KRATOS_REGISTER_VARIABLE(SHEAR_FORCE_2)

    KRATOS_REGISTER_VARIABLE(INTEGRATE_CONSERVATIVE)

    KRATOS_REGISTER_VARIABLE(PENALTY_FACTOR)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VECTOR_LAGRANGE_MULTIPLIER_REACTION)

    KRATOS_REGISTER_VARIABLE(NITSCHE_STABILIZATION_FACTOR)
    KRATOS_REGISTER_VARIABLE(EIGENVALUE_NITSCHE_STABILIZATION_SIZE)
    KRATOS_REGISTER_VARIABLE(EIGENVALUE_NITSCHE_STABILIZATION_VECTOR)
    KRATOS_REGISTER_VARIABLE(BUILD_LEVEL)

    // SBM Variables 
    KRATOS_REGISTER_VARIABLE(INTEGRATION_POINTS)
    KRATOS_REGISTER_VARIABLE(INTEGRATION_WEIGHTS)
    KRATOS_REGISTER_VARIABLE(BOUNDARY_CONDITION_TYPE)
}

}  // namespace Kratos
