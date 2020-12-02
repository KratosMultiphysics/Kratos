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
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/element.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"

// Include base h
#include "laplace_element.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
void LaplaceElement<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& CurrentProcessInfo) const
{
    if (rResult.size() != TNumNodes) {
        rResult.resize(TNumNodes, false);
    }

    const Variable<double>& r_variable = this->GetVariable();

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rResult[i] = Element::GetGeometry()[i].GetDof(r_variable).EquationId();
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void LaplaceElement<TDim, TNumNodes>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& CurrentProcessInfo) const
{
    if (rElementalDofList.size() != TNumNodes) {
        rElementalDofList.resize(TNumNodes);
    }

    const Variable<double>& r_variable = this->GetVariable();

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rElementalDofList[i] = Element::GetGeometry()[i].pGetDof(r_variable);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void LaplaceElement<TDim, TNumNodes>::GetValuesVector(
    VectorType& rValues,
    int Step) const
{
    if (rValues.size() != TNumNodes) {
        rValues.resize(TNumNodes, false);
    }

    BoundedVector<double, TNumNodes> values;
    this->GetValuesArray(values, Step);

    noalias(rValues) = values;
}

template <unsigned int TDim, unsigned int TNumNodes>
void LaplaceElement<TDim, TNumNodes>::GetValuesArray(
    BoundedVector<double, TNumNodes>& rValues,
    int Step) const
{
    const Variable<double>& r_variable = this->GetVariable();

    const auto& r_geometry = this->GetGeometry();
    IndexType LocalIndex = 0;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        rValues[LocalIndex++] =
            r_geometry[i_node].FastGetSolutionStepValue(r_variable, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod LaplaceElement<TDim, TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_2;
}

template <unsigned int TDim, unsigned int TNumNodes>
void LaplaceElement<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Calculate RHS
    if (rRightHandSideVector.size() != TNumNodes) {
        rRightHandSideVector.resize(TNumNodes, false);
    }

    noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

    // Calculate LHS
    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    BoundedVector<double, TNumNodes> values;
    this->GetValuesArray(values, 0);

    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, values);
}

template <unsigned int TDim, unsigned int TNumNodes>
void LaplaceElement<TDim, TNumNodes>::CalculateBoundedLeftHandSide(
    BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const IndexType num_gauss_points = gauss_weights.size();

    for (IndexType g = 0; g < num_gauss_points; ++g) {
        const Matrix& r_shape_derivatives = shape_derivatives[g];

        for (IndexType a = 0; a < TNumNodes; ++a) {
            for (IndexType b = 0; b < TNumNodes; ++b) {
                double dNa_dNb = 0.0;
                for (IndexType i = 0; i < TDim; ++i) {
                    dNa_dNb += r_shape_derivatives(a, i) * r_shape_derivatives(b, i);
                }

                rLeftHandSideMatrix(a, b) += gauss_weights[g] * dNa_dNb;
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void LaplaceElement<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes) {
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    }

    BoundedMatrix<double, TNumNodes, TNumNodes> local_matrix;
    this->CalculateBoundedLeftHandSide(local_matrix, rCurrentProcessInfo);

    noalias(rLeftHandSideMatrix) = local_matrix;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void LaplaceElement<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rRightHandSideVector.size() != TNumNodes) {
        rRightHandSideVector.resize(TNumNodes, false);
    }

    BoundedMatrix<double, TNumNodes, TNumNodes> local_matrix;
    this->CalculateBoundedLeftHandSide(local_matrix, rCurrentProcessInfo);

    BoundedVector<double, TNumNodes> values;
    this->GetValuesArray(values);

    noalias(rRightHandSideVector) = prod(local_matrix, values);
    noalias(rRightHandSideVector) = rRightHandSideVector * -1.0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
int LaplaceElement<TDim, TNumNodes>::Check(
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int check = BaseType::Check(rCurrentProcessInfo);

    const Variable<double>& r_variable = this->GetVariable();

    for (IndexType i_node = 0; i_node < this->GetGeometry().size(); ++i_node) {
        const auto& r_node = this->GetGeometry()[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_variable, r_node);
        KRATOS_CHECK_DOF_IN_NODE(r_variable, r_node);
    }

    return check;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void LaplaceElement<TDim, TNumNodes>::CalculateGeometryData(
    Vector& rGaussWeights,
    Matrix& rNContainer,
    ShapeFunctionDerivativesArrayType& rDN_DX) const
{
    const auto& r_geometry = this->GetGeometry();

    RansCalculationUtilities::CalculateGeometryData(
        r_geometry, this->GetIntegrationMethod(), rGaussWeights, rNContainer, rDN_DX);
}

// template instantiations
template class LaplaceElement<2, 3>;
template class LaplaceElement<3, 4>;

} // namespace Kratos.
