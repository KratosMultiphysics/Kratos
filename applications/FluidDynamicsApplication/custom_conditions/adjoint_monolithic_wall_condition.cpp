//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// Project includes
#include "includes/cfd_variables.h"
#include "utilities/indirect_scalar.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "fluid_dynamics_application_variables.h"

// Include base h
#include "adjoint_monolithic_wall_condition.h"

namespace Kratos
{
template <>
void AdjointMonolithicWallCondition<2, 2>::ThisExtensions::GetFirstDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpCondition->GetGeometry()[NodeId];
    rVector.resize(3);
    std::size_t index = 0;
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_X, Step);
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_Y, Step);
    rVector[index] = IndirectScalar<double>{}; // pressure
}

template <>
void AdjointMonolithicWallCondition<3, 3>::ThisExtensions::GetFirstDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpCondition->GetGeometry()[NodeId];
    rVector.resize(4);
    std::size_t index = 0;
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_X, Step);
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_Y, Step);
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_Z, Step);
    rVector[index] = IndirectScalar<double>{}; // pressure
}

template <>
void AdjointMonolithicWallCondition<2, 2>::ThisExtensions::GetSecondDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpCondition->GetGeometry()[NodeId];
    rVector.resize(3);
    std::size_t index = 0;
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_X, Step);
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_Y, Step);
    rVector[index] = IndirectScalar<double>{}; // pressure
}

template <>
void AdjointMonolithicWallCondition<3, 3>::ThisExtensions::GetSecondDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpCondition->GetGeometry()[NodeId];
    rVector.resize(4);
    std::size_t index = 0;
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_X, Step);
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_Y, Step);
    rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_Z, Step);
    rVector[index] = IndirectScalar<double>{}; // pressure
}

template <>
void AdjointMonolithicWallCondition<2, 2>::ThisExtensions::GetAuxiliaryVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpCondition->GetGeometry()[NodeId];
    rVector.resize(3);
    std::size_t index = 0;
    rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_X, Step);
    rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_Y, Step);
    rVector[index] = IndirectScalar<double>{}; // pressure
}

template <>
void AdjointMonolithicWallCondition<3, 3>::ThisExtensions::GetAuxiliaryVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpCondition->GetGeometry()[NodeId];
    rVector.resize(4);
    std::size_t index = 0;
    rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_X, Step);
    rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_Y, Step);
    rVector[index++] = MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_Z, Step);
    rVector[index] = IndirectScalar<double>{}; // pressure
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::ThisExtensions::GetFirstDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &ADJOINT_FLUID_VECTOR_2;
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::ThisExtensions::GetSecondDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &ADJOINT_FLUID_VECTOR_3;
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::ThisExtensions::GetAuxiliaryVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &AUX_ADJOINT_FLUID_VECTOR_1;
}

template <>
void AdjointMonolithicWallCondition<2, 2>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rResult.size() != 6) {
        rResult.resize(6);
    }

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(ADJOINT_FLUID_VECTOR_1_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(ADJOINT_FLUID_SCALAR_1);

    IndexType local_index = 0;
    for (IndexType i = 0; i < 2; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        rResult[local_index++] = r_node.GetDof(ADJOINT_FLUID_VECTOR_1_X, xpos).EquationId();
        rResult[local_index++] = r_node.GetDof(ADJOINT_FLUID_VECTOR_1_Y, xpos + 1).EquationId();
        rResult[local_index++] = r_node.GetDof(ADJOINT_FLUID_SCALAR_1, ppos).EquationId();
    }
}

template <>
void AdjointMonolithicWallCondition<3, 3>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rResult.size() != 12) {
        rResult.resize(12);
    }

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(ADJOINT_FLUID_VECTOR_1_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(ADJOINT_FLUID_SCALAR_1);

    IndexType local_index = 0;
    for (IndexType i = 0; i < 3; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        rResult[local_index++] = r_node.GetDof(ADJOINT_FLUID_VECTOR_1_X, xpos).EquationId();
        rResult[local_index++] = r_node.GetDof(ADJOINT_FLUID_VECTOR_1_Y, xpos + 1).EquationId();
        rResult[local_index++] = r_node.GetDof(ADJOINT_FLUID_VECTOR_1_Z, xpos + 2).EquationId();
        rResult[local_index++] = r_node.GetDof(ADJOINT_FLUID_SCALAR_1, ppos).EquationId();
    }
}

template <>
void AdjointMonolithicWallCondition<2, 2>::GetDofList(
    DofsVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rResult.size() != 6) {
        rResult.resize(6);
    }

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(ADJOINT_FLUID_VECTOR_1_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(ADJOINT_FLUID_SCALAR_1);

    IndexType local_index = 0;
    for (IndexType i = 0; i < 2; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        rResult[local_index++] = r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_X, xpos);
        rResult[local_index++] = r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_Y, xpos + 1);
        rResult[local_index++] = r_node.pGetDof(ADJOINT_FLUID_SCALAR_1, ppos);
    }
}

template <>
void AdjointMonolithicWallCondition<3, 3>::GetDofList(
    DofsVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rResult.size() != 12) {
        rResult.resize(12);
    }

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(ADJOINT_FLUID_VECTOR_1_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(ADJOINT_FLUID_SCALAR_1);

    IndexType local_index = 0;
    for (IndexType i = 0; i < 3; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        rResult[local_index++] = r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_X, xpos);
        rResult[local_index++] = r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_Y, xpos + 1);
        rResult[local_index++] = r_node.pGetDof(ADJOINT_FLUID_VECTOR_1_Z, xpos + 2);
        rResult[local_index++] = r_node.pGetDof(ADJOINT_FLUID_SCALAR_1, ppos);
    }
}

template <>
void AdjointMonolithicWallCondition<2, 2>::GetValuesVector(
    Vector& Values,
    int Step) const
{
    if (Values.size() != 6) {
        Values.resize(6, false);
    }

    IndexType local_index = 0;
    for (IndexType i = 0; i < 2; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        Values[local_index++] = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1_X);
        Values[local_index++] = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1_Y);
        Values[local_index++] = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1);
    }
}

template <>
void AdjointMonolithicWallCondition<3, 3>::GetValuesVector(
    Vector& Values,
    int Step) const
{
    if (Values.size() != 12) {
        Values.resize(12, false);
    }

    IndexType local_index = 0;
    for (IndexType i = 0; i < 3; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        Values[local_index++] = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1_X);
        Values[local_index++] = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1_Y);
        Values[local_index++] = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1_Z);
        Values[local_index++] = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1);
    }
}

template<unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::GetFirstDerivativesVector(
    Vector& Values,
    int Step) const
{
    if (Values.size() != TFluidLocalSize) {
        Values.resize(TFluidLocalSize, false);
    }

    Values.clear();
}

template<>
void AdjointMonolithicWallCondition<2, 2>::GetSecondDerivativesVector(
    Vector& Values,
    int Step) const
{
    if (Values.size() != 6) {
        Values.resize(6, false);
    }

    IndexType local_index = 0;
    for (IndexType i = 0; i < 2; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        Values[local_index++] = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3_X);
        Values[local_index++] = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3_Y);
        Values[local_index++] = 0.0;
    }
}

template<>
void AdjointMonolithicWallCondition<3, 3>::GetSecondDerivativesVector(
    Vector& Values,
    int Step) const
{
    if (Values.size() != 12) {
        Values.resize(12, false);
    }

    IndexType local_index = 0;
    for (IndexType i = 0; i < 3; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        Values[local_index++] = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3_X);
        Values[local_index++] = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3_Y);
        Values[local_index++] = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3_Z);
        Values[local_index++] = 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
int AdjointMonolithicWallCondition<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(this->Has(NORMAL))
        << "NORMAL is not defined for " << this->Info() << ".\n";

    KRATOS_ERROR_IF(norm_2(this->GetValue(NORMAL)) == 0.0)
        << "NORMAL is not properly initialized in " << this->Info() << ".\n";

    KRATOS_ERROR_IF_NOT(this->Has(NORMAL_SHAPE_DERIVATIVE))
        << "NORMAL_SHAPE_DERIVATIVE is not defined for " << this->Info() << ".\n";

    return Condition::Check(rCurrentProcessInfo);

    KRATOS_CATCH("");
}

template<>
double AdjointMonolithicWallCondition<2, 2>::CalculateDetJ(
    const double Area) const
{
    return Area * 0.5;
}

template<>
double AdjointMonolithicWallCondition<3, 3>::CalculateDetJ(
    const double Area) const
{
    return Area * 2.0;
}

template<>
void AdjointMonolithicWallCondition<2, 2>::CalculateDetJShapeDerivatives(
    BoundedVector<double, TCoordsLocalSize>& rDetJDerivatives,
    const BoundedVector<double, TCoordsLocalSize>& rAreaDerivatives) const
{
    noalias(rDetJDerivatives) = rAreaDerivatives * 0.5;
}

template<>
void AdjointMonolithicWallCondition<3, 3>::CalculateDetJShapeDerivatives(
    BoundedVector<double, TCoordsLocalSize>& rDetJDerivatives,
    const BoundedVector<double, TCoordsLocalSize>& rAreaDerivatives) const
{
    noalias(rDetJDerivatives) = rAreaDerivatives * 2.0;
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::CalculateData(
    double& rArea,
    double& rDetJ,
    array_1d<double, TDim>& rUnitNormal) const
{
    KRATOS_TRY

    const auto& aux_normal = this->GetValue(NORMAL);
    for (IndexType i = 0; i < TDim; ++i) {
        rUnitNormal[i] = aux_normal[i];
    }
    rArea = norm_2(rUnitNormal);
    rUnitNormal /= rArea;
    rDetJ = this->CalculateDetJ(rArea);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::CalculateDataShapeDerivatives(
    BoundedVector<double, TCoordsLocalSize>& rAreaDerivatives,
    BoundedVector<double, TCoordsLocalSize>& rDetJDerivatives,
    BoundedMatrix<double, TCoordsLocalSize, TDim>& rUnitNormalDerivatives,
    const double Area,
    const double DetJ,
    const array_1d<double, TDim>& rUnitNormal) const
{
    KRATOS_TRY

    const auto& normal_shape_derivatives = this->GetValue(NORMAL_SHAPE_DERIVATIVE);
    noalias(rAreaDerivatives) = prod(normal_shape_derivatives, rUnitNormal);
    this->CalculateDetJShapeDerivatives(rDetJDerivatives, rAreaDerivatives);
    rUnitNormalDerivatives = normal_shape_derivatives;
    for (IndexType c = 0; c < TCoordsLocalSize; ++c){
        for (IndexType i = 0; i < TDim; ++i) {
            rUnitNormalDerivatives(c, i) -= rUnitNormal[i] * rAreaDerivatives[c];
        }
    }
    rUnitNormalDerivatives /= Area;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rLeftHandSideMatrix.size1() != TFluidLocalSize) {
        rLeftHandSideMatrix.resize(TFluidLocalSize, TFluidLocalSize, false);
    }

    if (rRightHandSideVector.size() != TFluidLocalSize) {
        rRightHandSideVector.resize(TFluidLocalSize, false);
    }

    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::CalculateLocalVelocityContribution(
    MatrixType& rDampMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rDampMatrix.size1() != TFluidLocalSize) {
        rDampMatrix.resize(TFluidLocalSize, TFluidLocalSize, false);
    }

    if (rRightHandSideVector.size() != TFluidLocalSize) {
        rRightHandSideVector.resize(TFluidLocalSize, false);
    }

    rDampMatrix.clear();
    rRightHandSideVector.clear();

    if (this->Is(SLIP)) {
        this->ApplyWallLaw(rDampMatrix, rRightHandSideVector, rCurrentProcessInfo);
    } else {
        this->ApplyNeumannCondition(rDampMatrix, rRightHandSideVector, rCurrentProcessInfo);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::CalculateLeftHandSide(
    Matrix& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rLeftHandSideMatrix.size1() != TFluidLocalSize) {
        rLeftHandSideMatrix.resize(TFluidLocalSize, TFluidLocalSize, false);
    }

    rLeftHandSideMatrix.clear();

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::CalculateFirstDerivativesLHS(
    Matrix& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rLeftHandSideMatrix.size1() != TFluidLocalSize) {
        rLeftHandSideMatrix.resize(TFluidLocalSize, TFluidLocalSize, false);
    }

    rLeftHandSideMatrix.clear();

    if (this->Is(SLIP)) {
        this->ApplyWallLawStateDerivatives(rLeftHandSideMatrix, rCurrentProcessInfo);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::CalculateSecondDerivativesLHS(
    Matrix& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rLeftHandSideMatrix.size1() != TFluidLocalSize) {
        rLeftHandSideMatrix.resize(TFluidLocalSize, TFluidLocalSize, false);
    }

    rLeftHandSideMatrix.clear();

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rDesignVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rDesignVariable == SHAPE_SENSITIVITY) {
        if (rOutput.size1() != TCoordsLocalSize || rOutput.size2() != TFluidLocalSize) {
            rOutput.resize(TCoordsLocalSize, TFluidLocalSize, false);
        }

        rOutput.clear();

        if (this->Is(SLIP)) {
            this->ApplyWallLawShapeDerivatives(rOutput, rCurrentProcessInfo);
        } else {
            this->ApplyNeumannConditionShapeDerivatives(rOutput, rCurrentProcessInfo);
        }

    } else {
        KRATOS_ERROR << "Partial sensitivity w.r.t. " << rDesignVariable.Name()
                     << " not supported.";
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::ApplyNeumannCondition(
    MatrixType& rLocalMatrix,
    VectorType& rLocalVector,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    double area, detJ;
    array_1d<double, TDim> unit_normal;
    this->CalculateData(area, detJ, unit_normal);

    const auto& r_geometry = this->GetGeometry();
    const auto& integration_points = r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const IndexType number_of_gauss_points = integration_points.size();
    const Matrix& N_container = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
    const IndexType block_size = rLocalVector.size() / TNumNodes;

    for (IndexType g = 0; g < number_of_gauss_points; ++g) {
        const Vector& N = row(N_container, g);
        const double weight = detJ * integration_points[g].Weight();

        double external_pressure;
        FluidCalculationUtilities::EvaluateInPoint(
            r_geometry, N, std::tie(external_pressure, EXTERNAL_PRESSURE));

        for (IndexType a = 0; a < TNumNodes; ++a) {
            const IndexType row = a * block_size;
            for (IndexType i = 0; i < TDim; ++i) {
                rLocalVector[row + i] -= weight * N[a] * external_pressure * unit_normal[i];
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::ApplyNeumannConditionShapeDerivatives(
    MatrixType& rLocalMatrix,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    double area, detJ;
    array_1d<double, TDim> unit_normal;
    this->CalculateData(area, detJ, unit_normal);

    const auto& r_geometry = this->GetGeometry();

    const auto& integration_points = r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const IndexType number_of_gauss_points = integration_points.size();
    const Matrix& N_container = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
    const IndexType block_size = rLocalMatrix.size2() / TNumNodes;

    // compute data derivatives
    BoundedVector<double, TCoordsLocalSize> area_derivatives, detJ_derivatives;
    BoundedMatrix<double, TCoordsLocalSize, TDim> unit_normal_derivatives;
    this->CalculateDataShapeDerivatives(area_derivatives, detJ_derivatives,
                                        unit_normal_derivatives, area, detJ, unit_normal);

    for (IndexType g = 0; g < number_of_gauss_points; ++g) {
        const Vector& N = row(N_container, g);

        const double integration_weight = integration_points[g].Weight();
        const double weight = detJ * integration_weight;

        double external_pressure;
        FluidCalculationUtilities::EvaluateInPoint(
            r_geometry, N, std::tie(external_pressure, EXTERNAL_PRESSURE));

        for (IndexType c = 0; c < TNumNodes; ++c) {
            for (IndexType k = 0; k < TDim; ++k) {
                const IndexType derivative_index = c * TDim + k;

                const double weight_derivative =
                    detJ_derivatives[derivative_index] * integration_weight;
                const auto& unit_normal_derivative =
                    row(unit_normal_derivatives, derivative_index);

                for (IndexType a = 0; a < TNumNodes; ++a) {
                    for (IndexType i = 0; i < TDim; ++i) {
                        rLocalMatrix(derivative_index, a * block_size + i) -=
                            weight_derivative * N[a] * external_pressure * unit_normal[i];
                        rLocalMatrix(derivative_index, a * block_size + i) -=
                            weight * N[a] * external_pressure * unit_normal_derivative[i];
                    }
                }
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::ApplyWallLaw(
    MatrixType& rLocalMatrix,
    VectorType& rLocalVector,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto& r_geometry = this->GetGeometry();

    const double area = norm_2(r_geometry.GetValue(NORMAL)) / TNumNodes;
    const double kappa = 0.41;
    const double limit_y_plus = 10.9931899;
    const double beta = 5.2;

    for (IndexType i = 0; i < TNumNodes; ++i) {
        const auto& r_node = r_geometry[i];
        const double y = r_node.GetValue(Y_WALL);

        if (y > 0.0 && r_node.Is(SLIP)) {
            array_1d<double, 3> velocity = r_node.FastGetSolutionStepValue(VELOCITY);
            velocity -= r_node.FastGetSolutionStepValue(MESH_VELOCITY);
            const double nu = r_node.FastGetSolutionStepValue(VISCOSITY);
            const double rho = r_node.FastGetSolutionStepValue(DENSITY);
            const double velocity_norm = norm_2(velocity);

            if (velocity_norm > 1e-12) {
                const double y_plus = FluidCalculationUtilities::CalculateLogarithmicYPlus(
                    velocity_norm, y, nu, kappa, beta, limit_y_plus, 100);
                const double u_tau = y_plus * nu / y;

                const double coeff = area * u_tau * u_tau * rho / velocity_norm;
                for (size_t d = 0; d < TDim; d++) {
                    size_t k = i * TBlockSize + d;
                    rLocalVector[k] -= velocity[d] * coeff;
                    rLocalMatrix(k, k) += coeff;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::ApplyWallLawStateDerivatives(
    MatrixType& rLocalMatrix,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto& r_geometry = this->GetGeometry();

    const double area = norm_2(r_geometry.GetValue(NORMAL)) / TNumNodes;
    const double kappa = 0.41;
    const double limit_y_plus = 10.9931899;
    const double beta = 5.2;

    for (IndexType c = 0; c < TNumNodes; ++c) {
        const auto& r_node = r_geometry[c];
        const double y = r_node.GetValue(Y_WALL);

        if (y > 0.0 && r_node.Is(SLIP)) {
            array_1d<double, 3> velocity = r_node.FastGetSolutionStepValue(VELOCITY);
            velocity -= r_node.FastGetSolutionStepValue(MESH_VELOCITY);
            const double nu = r_node.FastGetSolutionStepValue(VISCOSITY);
            const double rho = r_node.FastGetSolutionStepValue(DENSITY);
            const double velocity_norm = norm_2(velocity);

            if (velocity_norm > 1e-12) {
                const double y_plus = FluidCalculationUtilities::CalculateLogarithmicYPlus(
                    velocity_norm, y, nu, kappa, beta, limit_y_plus, 100);
                const double u_tau = y_plus * nu / y;
                const double u_plus = velocity_norm / u_tau;
                const double coeff = area * u_tau * u_tau * rho / velocity_norm;

                for (IndexType k = 0; k < TDim; k++) {
                    const IndexType derivative_row = c * TBlockSize + k;

                    const double velocity_norm_derivative = velocity[k] / velocity_norm;
                    double u_tau_derivative = 0.0;
                    if (y_plus > limit_y_plus) {
                        u_tau_derivative =
                            kappa * velocity_norm_derivative / (kappa * u_plus + 1);
                    } else {
                        u_tau_derivative = velocity_norm_derivative / (2 * y_plus);
                    }
                    const double coeff_derivative =
                        area * rho *
                        (2 * u_tau * u_tau_derivative / velocity_norm -
                         std::pow(u_tau / velocity_norm, 2) * velocity_norm_derivative);

                    for (IndexType i = 0; i < TDim; ++i) {
                        rLocalMatrix(derivative_row, c * TBlockSize + i) -=
                            velocity[i] * coeff_derivative;
                    }

                    rLocalMatrix(derivative_row, derivative_row) -= coeff;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointMonolithicWallCondition<TDim, TNumNodes>::ApplyWallLawShapeDerivatives(
    MatrixType& rLocalMatrix,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto& r_geometry = this->GetGeometry();

    const auto& aux_normal = this->GetValue(NORMAL);
    array_1d<double, TDim> normal;
    for (IndexType i = 0; i < TDim; ++i) {
        normal[i] = aux_normal[i];
    }
    const double area = norm_2(normal) / TNumNodes;
    const auto& normal_derivatives = r_geometry.GetValue(NORMAL_SHAPE_DERIVATIVE);
    const BoundedVector<double, TCoordsLocalSize>& area_derivatives =
        prod(normal_derivatives, normal) / (area * TNumNodes * TNumNodes);

    const double kappa = 0.41;
    const double limit_y_plus = 10.9931899;
    const double beta = 5.2;

    for (IndexType i = 0; i < TNumNodes; ++i) {
        const auto& r_node = r_geometry[i];
        const double y = r_node.GetValue(Y_WALL);

        if (y > 0.0 && r_node.Is(SLIP)) {
            array_1d<double, 3> velocity = r_node.FastGetSolutionStepValue(VELOCITY);
            velocity -= r_node.FastGetSolutionStepValue(MESH_VELOCITY);
            const double nu = r_node.FastGetSolutionStepValue(VISCOSITY);
            const double rho = r_node.FastGetSolutionStepValue(DENSITY);
            const double velocity_norm = norm_2(velocity);

            if (velocity_norm > 1e-12) {
                const double y_plus = FluidCalculationUtilities::CalculateLogarithmicYPlus(
                    velocity_norm, y, nu, kappa, beta, limit_y_plus, 100);
                const double u_tau = y_plus * nu / y;

                IndexType local_index = 0;
                for (IndexType c = 0; c < TNumNodes; ++c) {
                    for (IndexType k = 0; k < TDim; ++k) {
                        const double coeff_derivative =
                            area_derivatives[local_index] * u_tau * u_tau * rho / velocity_norm;
                        for (IndexType j = 0; j < TDim; ++j) {
                            rLocalMatrix(local_index, i * TBlockSize + j) -=
                                velocity[j] * coeff_derivative;
                        }
                        ++local_index;
                    }
                }
            }
        }
    }

    KRATOS_CATCH("");
}

template class AdjointMonolithicWallCondition<2, 2>;
template class AdjointMonolithicWallCondition<3, 3>;

} // namespace Kratos
