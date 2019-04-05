//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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
    ProcessInfo &rCurrentProcessInfo)
{
    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

template <>
void EmbeddedNodalVariableCalculationElementSimplex<double>::CalculateLeftHandSide(
    MatrixType &rLeftHandSideMatrix,
    ProcessInfo &rCurrentProcessInfo)
{
    // Initialize Left Hand Side matrix
    noalias(rLeftHandSideMatrix) = ZeroMatrix(2,2);

    // Get the element shape function values from the normalized distance to node 0
    const auto N = this->GetDistanceBasedShapeFunctionValues();

    // Compute the Gramm matrix
    for (unsigned int i = 0; i < 2; ++i) {
        for (unsigned int j = 0; j < 2; ++j) {
            rLeftHandSideMatrix(i, j) = N[i] * N[j];
        }
    }
}

template <>
void EmbeddedNodalVariableCalculationElementSimplex<array_1d<double,3>>::CalculateLeftHandSide(
    MatrixType &rLeftHandSideMatrix,
    ProcessInfo &rCurrentProcessInfo)
{
    // Initialize Left Hand Side matrix
    noalias(rLeftHandSideMatrix) = ZeroMatrix(6,6);

    // Get the element shape function values from the normalized distance to node 0
    const auto N = this->GetDistanceBasedShapeFunctionValues();

    // Compute the Gramm matrix
    for (unsigned int i = 0; i < 2; ++i) {
        for (unsigned int j = 0; j < 2; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {
                rLeftHandSideMatrix(i * 3 + k, j * 3 + k) = N[i] * N[j];
            }
        }
    }
}
template <>
void EmbeddedNodalVariableCalculationElementSimplex<double>::CalculateRightHandSide(
    VectorType &rRigthHandSideVector,
    ProcessInfo &rCurrentProcessInfo)
{
    // Initialize Left Hand Side matrix
    noalias(rRigthHandSideVector) = ZeroVector(2);

    // Get the element shape function values from the normalized distance to node 0
    const auto N = this->GetDistanceBasedShapeFunctionValues();

    // Compute the Gramm matrix
    for (unsigned int i = 0; i < 2; ++i) {
        rRigthHandSideVector(i) = N[i];
    }
}

template <>
void EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>::CalculateRightHandSide(
    VectorType &rRigthHandSideVector,
    ProcessInfo &rCurrentProcessInfo)
{
    // Initialize Left Hand Side matrix
    noalias(rRigthHandSideVector) = ZeroVector(6);

    // Get the element shape function values from the normalized distance to node 0
    const auto N = this->GetDistanceBasedShapeFunctionValues();

    // Compute the Gramm matrix
    for (unsigned int i = 0; i < 2; ++i) {
        for (unsigned int k = 0; k < 3; ++k) {
            rRigthHandSideVector(i * 3 + k) = N[i];
        }
    }
}

template <>
void EmbeddedNodalVariableCalculationElementSimplex<double>::EquationIdVector(
    EquationIdVectorType &rResult,
    ProcessInfo &rCurrentProcessInfo)
{
    if (rResult.size() != 2) {
        rResult.resize(2);
    }

    for (unsigned int i = 0; i < 2; i++) {
        rResult[i] = GetGeometry()[i].GetDof(NODAL_MAUX).EquationId();
    }
}

template <>
void EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>::EquationIdVector(
    EquationIdVectorType &rResult,
    ProcessInfo &rCurrentProcessInfo)
{
    if (rResult.size() != 6) {
        rResult.resize(6);
    }

    for (unsigned int i = 0; i < 2; i++) {
        rResult[i * 3] = GetGeometry()[i].GetDof(NODAL_VAUX_X).EquationId();
        rResult[i * 3 + 1] = GetGeometry()[i].GetDof(NODAL_VAUX_Y).EquationId();
        rResult[i * 3 + 2] = GetGeometry()[i].GetDof(NODAL_VAUX_Z).EquationId();
    }
}

template <>
void EmbeddedNodalVariableCalculationElementSimplex<double>::GetDofList(
    DofsVectorType &rElementalDofList,
    ProcessInfo &rCurrentProcessInfo)
{
    if (rElementalDofList.size() != 2) {
        rElementalDofList.resize(2);
    }

    for (unsigned int i = 0; i < 2; i++) {
        rElementalDofList[i] = GetGeometry()[i].pGetDof(NODAL_MAUX);
    }
}

template <>
void EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>::GetDofList(
    DofsVectorType &rElementalDofList,
    ProcessInfo &rCurrentProcessInfo)
{
    if (rElementalDofList.size() != 6) {
        rElementalDofList.resize(6);
    }

    for (unsigned int i = 0; i < 2; i++) {
        rElementalDofList[i * 3] = GetGeometry()[i].pGetDof(NODAL_VAUX_X);
        rElementalDofList[i * 3 + 1] = GetGeometry()[i].pGetDof(NODAL_VAUX_Y);
        rElementalDofList[i * 3 + 2] = GetGeometry()[i].pGetDof(NODAL_VAUX_Z);
    }
}

template <class TVarType>
const array_1d<double, 2> EmbeddedNodalVariableCalculationElementSimplex<TVarType>::GetDistanceBasedShapeFunctionValues()
{
    const auto d = this->GetValue(DISTANCE);
    array_1d<double, 2> N;
    N[0] = 1.0 - d;
    N[1] = d;
    return N;
}

template class EmbeddedNodalVariableCalculationElementSimplex<double>;
template class EmbeddedNodalVariableCalculationElementSimplex<array_1d<double,3>>;

} // namespace Kratos.
