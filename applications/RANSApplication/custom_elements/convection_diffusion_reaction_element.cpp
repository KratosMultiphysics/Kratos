//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/element.h"

// Application includes
#include "custom_elements/data_containers/k_epsilon/k_element_data.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data.h"
#include "custom_elements/data_containers/k_omega/k_element_data.h"
#include "custom_elements/data_containers/k_omega/omega_element_data.h"
#include "custom_elements/data_containers/k_omega_sst/k_element_data.h"
#include "custom_elements/data_containers/k_omega_sst/omega_element_data.h"
#include "custom_utilities/rans_calculation_utilities.h"

// Include base h
#include "convection_diffusion_reaction_element.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& CurrentProcessInfo) const
{
    if (rResult.size() != TNumNodes) {
        rResult.resize(TNumNodes, false);
    }

    const Variable<double>& r_variable =
        TConvectionDiffusionReactionData::GetScalarVariable();

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rResult[i] = Element::GetGeometry()[i].GetDof(r_variable).EquationId();
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& CurrentProcessInfo) const
{
    if (rElementalDofList.size() != TNumNodes) {
        rElementalDofList.resize(TNumNodes);
    }

    const Variable<double>& r_variable =
        TConvectionDiffusionReactionData::GetScalarVariable();

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rElementalDofList[i] = Element::GetGeometry()[i].pGetDof(r_variable);
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::GetValuesVector(
    Vector& rValues,
    int Step) const
{
    this->GetFirstDerivativesVector(rValues, Step);
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::GetValuesArray(
    BoundedVector<double, TNumNodes>& rValues,
    const int Step) const
{
    const auto& r_geometry = this->GetGeometry();
    const Variable<double>& r_variable =
        TConvectionDiffusionReactionData::GetScalarVariable();

    IndexType LocalIndex = 0;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        rValues[LocalIndex++] =
            r_geometry[i_node].FastGetSolutionStepValue(r_variable, Step);
    }
}


template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::GetFirstDerivativesVector(
    Vector& rValues,
    int Step) const
{
    if (rValues.size() != TNumNodes) {
        rValues.resize(TNumNodes, false);
    }

    BoundedVector<double, TNumNodes> values;
    this->GetValuesArray(values, Step);

    noalias(rValues) = values;
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::GetSecondDerivativesVector(
    Vector& rValues,
    int Step) const
{
    if (rValues.size() != TNumNodes) {
        rValues.resize(TNumNodes, false);
    }

    const auto& r_geometry = this->GetGeometry();
    const Variable<double>& r_variable =
        TConvectionDiffusionReactionData::GetScalarVariable().GetTimeDerivative();

    IndexType LocalIndex = 0;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        rValues[LocalIndex++] =
            r_geometry[i_node].FastGetSolutionStepValue(r_variable, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Check sizes and initialize matrix
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes) {
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

    // Calculate RHS
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rRightHandSideVector.size() != TNumNodes) {
        rRightHandSideVector.resize(TNumNodes, false);
    }

    noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const IndexType num_gauss_points = gauss_weights.size();

    const auto& r_geometry = this->GetGeometry();
    TConvectionDiffusionReactionData r_current_data(r_geometry);

    r_current_data.CalculateConstants(rCurrentProcessInfo);

    for (IndexType g = 0; g < num_gauss_points; ++g) {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector gauss_shape_functions = row(shape_functions, g);

        r_current_data.CalculateGaussPointData(gauss_shape_functions, r_shape_derivatives);
        const double source = r_current_data.CalculateSourceTerm(
            gauss_shape_functions, r_shape_derivatives);

        noalias(rRightHandSideVector) +=
            gauss_shape_functions * (source * gauss_weights[g]);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::CalculateLocalVelocityContribution(
    MatrixType& rDampingMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateDampingMatrix(rDampingMatrix, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rMassMatrix.size1() != TNumNodes || rMassMatrix.size2() != TNumNodes) {
        rMassMatrix.resize(TNumNodes, TNumNodes, false);
    }

    rMassMatrix.clear();

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const IndexType num_gauss_points = gauss_weights.size();

    for (IndexType g = 0; g < num_gauss_points; ++g) {
        const double mass = gauss_weights[g] * (1.0 / TNumNodes);
        this->AddLumpedMassMatrix(rMassMatrix, mass);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rDampingMatrix.size1() != TNumNodes || rDampingMatrix.size2() != TNumNodes) {
        rDampingMatrix.resize(TNumNodes, TNumNodes, false);
    }

    rDampingMatrix.clear();

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const IndexType num_gauss_points = gauss_weights.size();

    const auto& r_geometry = this->GetGeometry();
    TConvectionDiffusionReactionData r_current_data(r_geometry);

    r_current_data.CalculateConstants(rCurrentProcessInfo);

    BoundedVector<double, TNumNodes> velocity_convective_terms;

    for (IndexType g = 0; g < num_gauss_points; ++g) {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector& r_shape_functions = row(shape_functions, g);

        r_current_data.CalculateGaussPointData(r_shape_functions, r_shape_derivatives);
        const array_1d<double, 3>& velocity = r_current_data.CalculateEffectiveVelocity(
            r_shape_functions, r_shape_derivatives);
        this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

        const double effective_kinematic_viscosity =
            r_current_data.CalculateEffectiveKinematicViscosity(
                r_shape_functions, r_shape_derivatives);

        const double reaction = r_current_data.CalculateReactionTerm(
            r_shape_functions, r_shape_derivatives);

        const Matrix& dNa_dNb = prod(r_shape_derivatives, trans(r_shape_derivatives));

        AddDampingMatrixGaussPointContributions(
            rDampingMatrix, reaction, effective_kinematic_viscosity,
            velocity_convective_terms, gauss_weights[g], r_shape_functions, dNa_dNb);
    }

    r_current_data.UpdateElementDataValueContainer(*this);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
int ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::Check(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int check = BaseType::Check(rCurrentProcessInfo);
    TConvectionDiffusionReactionData::Check(this->GetGeometry(), rCurrentProcessInfo);

    return check;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
GeometryData::IntegrationMethod ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
double ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::GetDivergenceOperator(
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rShapeDerivatives,
    const int Step) const
{
    double value = 0.0;
    const auto& r_geometry = this->GetGeometry();

    for (IndexType i = 0; i < TNumNodes; ++i) {
        const array_1d<double, 3>& r_value =
            r_geometry[i].FastGetSolutionStepValue(rVariable, Step);
        for (IndexType j = 0; j < TDim; ++j) {
            value += r_value[j] * rShapeDerivatives(i, j);
        }
    }

    return value;
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::GetConvectionOperator(
    BoundedVector<double, TNumNodes>& rOutput,
    const array_1d<double, 3>& rVector,
    const Matrix& rShapeDerivatives) const
{
    rOutput.clear();
    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType j = 0; j < TDim; ++j) {
            rOutput[i] += rVector[j] * rShapeDerivatives(i, j);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::CalculateGeometryData(
    Vector& rGaussWeights,
    Matrix& rNContainer,
    ShapeFunctionDerivativesArrayType& rDN_DX) const
{
    const auto& r_geometry = this->GetGeometry();

    RansCalculationUtilities::CalculateGeometryData(
        r_geometry, this->GetIntegrationMethod(),
        rGaussWeights, rNContainer, rDN_DX);
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::AddLumpedMassMatrix(
    Matrix& rMassMatrix,
    const double Mass) const
{
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        rMassMatrix(i_node, i_node) += Mass;
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::AddDampingMatrixGaussPointContributions(
    Matrix& rDampingMatrix,
    const double ReactionTerm,
    const double EffectiveKinematicViscosity,
    const Vector& rVelocityConvectiveTerms,
    const double GaussWeight,
    const Vector& rGaussShapeFunctions,
    const Matrix& rGaussdNa_dNb) const
{
    for (IndexType a = 0; a < TNumNodes; ++a) {
        for (IndexType b = 0; b < TNumNodes; ++b) {
            double value = 0.0;

            value += rGaussShapeFunctions[a] * rVelocityConvectiveTerms[b];
            value += rGaussShapeFunctions[a] * ReactionTerm * rGaussShapeFunctions[b];
            value += EffectiveKinematicViscosity * rGaussdNa_dNb(a, b);

            rDampingMatrix(a, b) += GaussWeight * value;
        }
    }
}

// template instantiations
template class ConvectionDiffusionReactionElement<2, 3, KEpsilonElementData::KElementData<2>>;
template class ConvectionDiffusionReactionElement<2, 3, KEpsilonElementData::EpsilonElementData<2>>;

template class ConvectionDiffusionReactionElement<2, 3, KOmegaElementData::KElementData<2>>;
template class ConvectionDiffusionReactionElement<2, 3, KOmegaElementData::OmegaElementData<2>>;

template class ConvectionDiffusionReactionElement<2, 3, KOmegaSSTElementData::KElementData<2>>;
template class ConvectionDiffusionReactionElement<2, 3, KOmegaSSTElementData::OmegaElementData<2>>;

template class ConvectionDiffusionReactionElement<3, 4, KEpsilonElementData::KElementData<3>>;
template class ConvectionDiffusionReactionElement<3, 4, KEpsilonElementData::EpsilonElementData<3>>;

template class ConvectionDiffusionReactionElement<3, 4, KOmegaElementData::KElementData<3>>;
template class ConvectionDiffusionReactionElement<3, 4, KOmegaElementData::OmegaElementData<3>>;

template class ConvectionDiffusionReactionElement<3, 4, KOmegaSSTElementData::KElementData<3>>;
template class ConvectionDiffusionReactionElement<3, 4, KOmegaSSTElementData::OmegaElementData<3>>;

} // namespace Kratos.
