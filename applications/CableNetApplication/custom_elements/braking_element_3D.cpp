//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter

// System includes

// External includes

// Project includes
#include "custom_elements/braking_element_3D.hpp"
#include "structural_mechanics_application_variables.h"
#include "cable_net_application_variables.h"
#include "includes/checks.h"
#include "includes/variables.h"
#include "utilities/atomic_utilities.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"


namespace Kratos {
BrakingElement3D1N::BrakingElement3D1N(IndexType NewId,
                                   GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

BrakingElement3D1N::BrakingElement3D1N(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

Element::Pointer
BrakingElement3D1N::Create(IndexType NewId, NodesArrayType const& rThisNodes,
                         PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<BrakingElement3D1N>(NewId, rGeom.Create(rThisNodes),
            pProperties);
}

Element::Pointer
BrakingElement3D1N::Create(IndexType NewId, GeometryType::Pointer pGeom,
                         PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<BrakingElement3D1N>(NewId, pGeom,
            pProperties);
}

BrakingElement3D1N::~BrakingElement3D1N() {}

void BrakingElement3D1N::EquationIdVector(EquationIdVectorType& rResult,
                                                  const ProcessInfo& rCurrentProcessInfo) const
{
    if (rResult.size() != msLocalSize) {
        rResult.resize(msLocalSize);
    }

    rResult[0] = GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[1] = GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[2] = GetGeometry()[0].GetDof(DISPLACEMENT_Z).EquationId();

}
void BrakingElement3D1N::GetDofList(DofsVectorType& rElementalDofList,
                                            const ProcessInfo& rCurrentProcessInfo) const
{
    if (rElementalDofList.size() != msLocalSize) {
        rElementalDofList.resize(msLocalSize);
    }

    rElementalDofList[0] = GetGeometry()[0].pGetDof(DISPLACEMENT_X);
    rElementalDofList[1] = GetGeometry()[0].pGetDof(DISPLACEMENT_Y);
    rElementalDofList[2] = GetGeometry()[0].pGetDof(DISPLACEMENT_Z);

}


void BrakingElement3D1N::GetValuesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY
    if (rValues.size() != msLocalSize) {
        rValues.resize(msLocalSize, false);
    }
    rValues = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, Step);
    KRATOS_CATCH("")
}

void BrakingElement3D1N::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY
    if (rValues.size() != msLocalSize) {
        rValues.resize(msLocalSize, false);
    }
    rValues = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, Step);
    KRATOS_CATCH("")
}

void BrakingElement3D1N::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_TRY
    if (rValues.size() != msLocalSize) {
        rValues.resize(msLocalSize, false);
    }
    rValues = GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION, Step);
    KRATOS_CATCH("")
}

void BrakingElement3D1N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                      VectorType& rRightHandSideVector,
                                                      const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
    KRATOS_CATCH("")
}

void BrakingElement3D1N::CalculateRightHandSide(
    VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    rRightHandSideVector = ZeroVector(msLocalSize);
    Vector internal_forces = ZeroVector(msLocalSize);
    CalculateInternalForces(internal_forces);
    noalias(rRightHandSideVector) -= internal_forces;

    // add body forces

    KRATOS_CATCH("")
}


void BrakingElement3D1N::CalculateInternalForces(Vector& rInternalForces)
{
    KRATOS_TRY

    Vector current_displacement = ZeroVector(msLocalSize);
    GetValuesVector(current_displacement);
    const double total_disp = norm_2(current_displacement);
    const double elastic_disp = total_disp - mAccumulatedPlasticDisplacement;

    // interpret E, Ep as spring stiffnesses k, kp
    const double plastic_stiffness = GetProperties()[HARDENING_MODULUS_1D];
    const double elastic_stiffness = GetProperties()[YOUNG_MODULUS];


    double current_force = elastic_stiffness * elastic_disp;
    double temp_force(current_force);

    // check if forces are negative?
    // YIELD_STRESS : yield force
    // isnt mAccumulatedPlasticAlpha just the abs accumulated plastic deformation? Get rid of that variable?
    double trial_yield_function = std::abs(current_force) - (GetProperties()[YIELD_STRESS] + mAccumulatedPlasticAlpha * plastic_stiffness);

    if (trial_yield_function > 0.0)
    { 
        // plastic
        const double delta_gamma = trial_yield_function/(elastic_stiffness+plastic_stiffness);
        current_force = 1.00 - ((delta_gamma*elastic_stiffness)/std::abs(temp_force));
        current_force *= temp_force;

        mAccumulatedPlasticDisplacement += delta_gamma*MathUtils<double>::Sign(temp_force);
        mAccumulatedPlasticAlpha += delta_gamma;
    }

    if (total_disp > 0.0) rInternalForces = (current_displacement / total_disp) * current_force;

    KRATOS_CATCH("")
}


void BrakingElement3D1N::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY;
    // resizing the matrices + create memory for LHS
    rLeftHandSideMatrix = ZeroMatrix(msLocalSize, msLocalSize);
    // creating LHS

    rLeftHandSideMatrix(0,0) = 1.0;
    rLeftHandSideMatrix(1,1) = 1.0;
    rLeftHandSideMatrix(2,2) = 1.0;

    Vector current_displacement = ZeroVector(msLocalSize);
    GetValuesVector(current_displacement);
    const double total_disp = norm_2(current_displacement);
    const double elastic_disp = total_disp - mAccumulatedPlasticDisplacement;

    // interpret E, Ep as spring stiffnesses k, kp
    const double plastic_stiffness = GetProperties()[HARDENING_MODULUS_1D];
    const double elastic_stiffness = GetProperties()[YOUNG_MODULUS];

    double current_force = elastic_stiffness * elastic_disp;

    double k = elastic_stiffness;

    // check if forces are negative?
    // YIELD_STRESS : yield force
    // isnt mAccumulatedPlasticAlpha just the abs accumulated plastic deformation? Get rid of that variable?
    double trial_yield_function = std::abs(current_force) - (GetProperties()[YIELD_STRESS] + mAccumulatedPlasticAlpha * plastic_stiffness);
    if (trial_yield_function > 0.0) k = (plastic_stiffness*elastic_stiffness)/(plastic_stiffness+elastic_stiffness);

    rLeftHandSideMatrix *= k;

    KRATOS_CATCH("")
}

int BrakingElement3D1N::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if (dimension != msDimension ||number_of_nodes != msNumberOfNodes) {
        KRATOS_ERROR << "The element works only in 3D and with 2 nodes" << std::endl;
    }

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const NodeType& rnode = GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, rnode);

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode);
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode);
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode);
    }

    if (GetProperties().Has(YIELD_STRESS) == false ||
            GetProperties()[YIELD_STRESS] <= numerical_limit) {
        KRATOS_ERROR << "YIELD_STRESS not provided for this element " << Id()
                     << std::endl;
    }

    if (GetProperties().Has(YOUNG_MODULUS) == false ||
            GetProperties()[YOUNG_MODULUS] <= numerical_limit) {
        KRATOS_ERROR << "YOUNG_MODULUS not provided for this element " << Id()
                     << std::endl;
    }

    if (GetProperties().Has(HARDENING_MODULUS_1D) == false ||
            GetProperties()[HARDENING_MODULUS_1D] <= numerical_limit) {
        KRATOS_ERROR << "HARDENING_MODULUS_1D not provided for this element " << Id()
                     << std::endl;
    }

    if (GetProperties().Has(BRAKING_MASS) == false ||
            GetProperties()[BRAKING_MASS] <= numerical_limit) {
        KRATOS_ERROR << "BRAKING_MASS not provided for this element " << Id()
                     << std::endl;
    }

    return 0;

    KRATOS_CATCH("")
}

void BrakingElement3D1N::CalculateLumpedMassVector(
    VectorType& rLumpedMassVector,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Clear matrix
    if (rLumpedMassVector.size() != msLocalSize) {
        rLumpedMassVector.resize(msLocalSize, false);
    }

    const double total_mass = GetProperties()[BRAKING_MASS];

    for (int i = 0; i < msNumberOfNodes; ++i) {
        for (int j = 0; j < msDimension; ++j) {
            int index = i * msDimension + j;
            rLumpedMassVector[index] = total_mass * 0.50;
        }
    }

    KRATOS_CATCH("")
}


void BrakingElement3D1N::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    // Compute lumped mass matrix
    VectorType temp_vector(msLocalSize);
    CalculateLumpedMassVector(temp_vector, rCurrentProcessInfo);

    // Clear matrix
    if (rMassMatrix.size1() != msLocalSize || rMassMatrix.size2() != msLocalSize) {
        rMassMatrix.resize(msLocalSize, msLocalSize, false);
    }
    rMassMatrix = ZeroMatrix(msLocalSize, msLocalSize);

    // Fill the matrix
    for (IndexType i = 0; i < msLocalSize; ++i) {
        rMassMatrix(i, i) = temp_vector[i];
    }

    KRATOS_CATCH("")
}

void BrakingElement3D1N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        msLocalSize);
}


void BrakingElement3D1N::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<double >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    auto& r_geom = GetGeometry();

    if (rDestinationVariable == NODAL_MASS) {
        VectorType element_mass_vector(msLocalSize);
        CalculateLumpedMassVector(element_mass_vector, rCurrentProcessInfo);

        for (SizeType i = 0; i < msNumberOfNodes; ++i) {
            double& r_nodal_mass = r_geom[i].GetValue(NODAL_MASS);
            int index = i * msDimension;

            AtomicAdd(r_nodal_mass, element_mass_vector(index));
        }
    }

    KRATOS_CATCH("")
}

void BrakingElement3D1N::AddExplicitContribution(
    const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {

        BoundedVector<double, msLocalSize> damping_residual_contribution = ZeroVector(msLocalSize);
        Vector current_nodal_velocities = ZeroVector(msLocalSize);
        GetFirstDerivativesVector(current_nodal_velocities);
        Matrix damping_matrix;
        ProcessInfo temp_process_information; // cant pass const ProcessInfo
        CalculateDampingMatrix(damping_matrix, temp_process_information);
        // current residual contribution due to damping
        noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);

        for (size_t i = 0; i < msNumberOfNodes; ++i) {
            size_t index = msDimension * i;
            array_1d<double, 3>& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            for (size_t j = 0; j < msDimension; ++j) {
                AtomicAdd(r_force_residual[j], rRHSVector[index + j] - damping_residual_contribution[index + j]);
            }
        }
    } else if (rDestinationVariable == NODAL_INERTIA) {

        // Getting the vector mass
        VectorType mass_vector(msLocalSize);
        CalculateLumpedMassVector(mass_vector, rCurrentProcessInfo);

        for (int i = 0; i < msNumberOfNodes; ++i) {
            double& r_nodal_mass = GetGeometry()[i].GetValue(NODAL_MASS);
            int index = i * msDimension;

            AtomicAdd(r_nodal_mass, mass_vector[index]);
        }
    }

    KRATOS_CATCH("")
}

void BrakingElement3D1N::save(Serializer &rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

void BrakingElement3D1N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

} // namespace Kratos.
