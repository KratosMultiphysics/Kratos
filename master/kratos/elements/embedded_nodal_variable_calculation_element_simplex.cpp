//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes


// External includes


// Project includes
#include "elements/embedded_nodal_variable_calculation_element_simplex.h"


namespace Kratos
{

template <class TVarType>
void EmbeddedNodalVariableCalculationElementSimplex<TVarType>::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

template <>
int EmbeddedNodalVariableCalculationElementSimplex<double>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Perform basic element checks
    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    if(ErrorCode != 0) {
        return ErrorCode;
    }

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(auto &r_node : this->GetGeometry()) {
        KRATOS_ERROR_IF(!r_node.SolutionStepsDataHas(NODAL_MAUX)) << "Missing NODAL_MAUX variable on solution step data for node " << r_node.Id() << std::endl;;
    }

    return 0;

    KRATOS_CATCH("");
}

template <>
int EmbeddedNodalVariableCalculationElementSimplex<array_1d<double,3>>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Perform basic element checks
    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    if(ErrorCode != 0) {
        return ErrorCode;
    }

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(auto &r_node : this->GetGeometry()) {
        KRATOS_ERROR_IF(!r_node.SolutionStepsDataHas(NODAL_VAUX)) << "Missing NODAL_VAUX variable on solution step data for node " << r_node.Id() << std::endl;;
    }

    return 0;

    KRATOS_CATCH("");
}

template <>
void EmbeddedNodalVariableCalculationElementSimplex<double>::CalculateLeftHandSide(
    MatrixType &rLeftHandSideMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Check size
    if (rLeftHandSideMatrix.size1() != 2 || rLeftHandSideMatrix.size2() != 2) {
        rLeftHandSideMatrix.resize(2, 2, false);
    }

    // Get the element shape function values from the normalized distance to node 0
    const auto &rN = this->GetDistanceBasedShapeFunctionValues();

    // Compute the Gramm matrix and gradient penalty term
    const double penalty = rCurrentProcessInfo[GRADIENT_PENALTY_COEFFICIENT];
    std::array<double, 2> aux_penalty{{penalty, -penalty}};
    for (unsigned int i = 0; i < 2; ++i) {
        for (unsigned int j = 0; j < 2; ++j) {
            rLeftHandSideMatrix(i, j) = rN[i] * rN[j] + aux_penalty[i] * aux_penalty[j];
        }
    }
}

template <>
void EmbeddedNodalVariableCalculationElementSimplex<array_1d<double,3>>::CalculateLeftHandSide(
    MatrixType &rLeftHandSideMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Check size
    if (rLeftHandSideMatrix.size1() != 6 || rLeftHandSideMatrix.size2() != 6) {
        rLeftHandSideMatrix.resize(6, 6, false);
    }

    // Initialize LHS. This is required since not all the entries of the matrix are iterated
    rLeftHandSideMatrix = ZeroMatrix(6,6);

    // Get the element shape function values from the normalized distance to node 0
    const auto &rN = this->GetDistanceBasedShapeFunctionValues();

    // Compute the Gramm matrix and gradient penalty term
    const double penalty = rCurrentProcessInfo[GRADIENT_PENALTY_COEFFICIENT];
    std::array<double, 2> aux_penalty{{penalty, -penalty}};
    for (unsigned int i = 0; i < 2; ++i) {
        for (unsigned int j = 0; j < 2; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {
                rLeftHandSideMatrix(i * 3 + k, j * 3 + k) = rN[i] * rN[j] + aux_penalty[i] * aux_penalty[j];
            }
        }
    }
}

template <>
void EmbeddedNodalVariableCalculationElementSimplex<double>::CalculateRightHandSide(
    VectorType &rRigthHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Check size
    if (rRigthHandSideVector.size() != 2) {
        rRigthHandSideVector.resize(2, false);
    }

    // Get the element shape function values from the normalized distance to node 0
    const auto &r_geom = this->GetGeometry();
    const double &rData = this->GetValue(NODAL_MAUX);
    const auto &rN = this->GetDistanceBasedShapeFunctionValues();

    // Compute the data and penalty Right Hand Side contributions
    const double penalty = rCurrentProcessInfo[GRADIENT_PENALTY_COEFFICIENT];
    std::array<double, 2> aux_penalty{{penalty, -penalty}};
    for (unsigned int i = 0; i < 2; ++i) {
        rRigthHandSideVector(i) = rN[i] * rData;
        for (unsigned int j = 0; j < 2; ++j) {
            rRigthHandSideVector(i) -=  (rN[i] * rN[j] + aux_penalty[i] * aux_penalty[j]) * r_geom[j].FastGetSolutionStepValue(NODAL_MAUX);
        }
    }
}

template <>
void EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>::CalculateRightHandSide(
    VectorType &rRigthHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Check size
    if (rRigthHandSideVector.size() != 6) {
        rRigthHandSideVector.resize(6, false);
    }

    // Get the element shape function values from the normalized distance to node 0
    const auto &r_geom = this->GetGeometry();
    const array_1d<double,3> &rData = this->GetValue(NODAL_VAUX);
    const auto &rN = this->GetDistanceBasedShapeFunctionValues();

    // Compute the data and penalty Right Hand Side contributions
    const double penalty = rCurrentProcessInfo[GRADIENT_PENALTY_COEFFICIENT];
    std::array<double, 2> aux_penalty{{penalty, -penalty}};
    for (unsigned int i = 0; i < 2; ++i) {
        const auto &r_aux = r_geom[i].FastGetSolutionStepValue(NODAL_VAUX);
        for (unsigned int k = 0; k < 3; ++k) {
            rRigthHandSideVector(i * 3 + k) = rN[i] * rData(k);
            for (unsigned int j = 0; j < 2; ++j) {
                rRigthHandSideVector(i * 3 + k) -= (rN[i] * rN[j] + aux_penalty[i] * aux_penalty[j]) * r_aux[k];
            }
        }
    }
}

template <>
void EmbeddedNodalVariableCalculationElementSimplex<double>::EquationIdVector(
    EquationIdVectorType &rResult,
    const ProcessInfo &rCurrentProcessInfo) const
{
    if (rResult.size() != 2) {
        rResult.resize(2, false);
    }

    const unsigned int pos = (this->GetGeometry())[0].GetDofPosition(NODAL_MAUX);
    for (unsigned int i = 0; i < 2; i++) {
        rResult[i] = (this->GetGeometry())[i].GetDof(NODAL_MAUX, pos).EquationId();
    }
}

template <>
void EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>::EquationIdVector(
    EquationIdVectorType &rResult,
    const ProcessInfo &rCurrentProcessInfo) const
{
    if (rResult.size() != 6) {
        rResult.resize(6, false);
    }

    const unsigned int x_pos = (this->GetGeometry())[0].GetDofPosition(NODAL_VAUX_X);
    for (unsigned int i = 0; i < 2; i++) {
        rResult[i * 3] = GetGeometry()[i].GetDof(NODAL_VAUX_X, x_pos).EquationId();
        rResult[i * 3 + 1] = GetGeometry()[i].GetDof(NODAL_VAUX_Y, x_pos + 1).EquationId();
        rResult[i * 3 + 2] = GetGeometry()[i].GetDof(NODAL_VAUX_Z, x_pos + 2).EquationId();
    }
}

template <>
void EmbeddedNodalVariableCalculationElementSimplex<double>::GetDofList(
    DofsVectorType &rElementalDofList,
    const ProcessInfo &rCurrentProcessInfo) const
{
    if (rElementalDofList.size() != 2) {
        rElementalDofList.resize(2);
    }

    for (unsigned int i = 0; i < 2; i++) {
        rElementalDofList[i] = (this->GetGeometry())[i].pGetDof(NODAL_MAUX);
    }
}

template <>
void EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>::GetDofList(
    DofsVectorType &rElementalDofList,
    const ProcessInfo &rCurrentProcessInfo) const
{
    if (rElementalDofList.size() != 6) {
        rElementalDofList.resize(6);
    }

    for (unsigned int i = 0; i < 2; i++) {
        rElementalDofList[i * 3] = (this->GetGeometry())[i].pGetDof(NODAL_VAUX_X);
        rElementalDofList[i * 3 + 1] = (this->GetGeometry())[i].pGetDof(NODAL_VAUX_Y);
        rElementalDofList[i * 3 + 2] = (this->GetGeometry())[i].pGetDof(NODAL_VAUX_Z);
    }
}

template <class TVarType>
const array_1d<double, 2> EmbeddedNodalVariableCalculationElementSimplex<TVarType>::GetDistanceBasedShapeFunctionValues()
{
    const double d = this->GetValue(DISTANCE);
    array_1d<double, 2> N;
    N[0] = 1.0 - d;
    N[1] = d;
    return N;
}

template class EmbeddedNodalVariableCalculationElementSimplex<double>;
template class EmbeddedNodalVariableCalculationElementSimplex<array_1d<double,3>>;

} // namespace Kratos.
