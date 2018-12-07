//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Iñigo Lopez and Riccardo Rossi
//

#include "incompressible_potential_flow_element.h"

namespace Kratos
{

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::EquationIdVector(
    EquationIdVectorType &rResult,
    ProcessInfo &CurrentProcessInfo)
{
    if (this->IsNot(MARKER)) //normal element
    {
        if (rResult.size() != NumNodes)
            rResult.resize(NumNodes, false);

        const IncompressiblePotentialFlowElement &r_this = *this;
        const int &kutta = r_this.GetValue(KUTTA);

        if (kutta == 0)
            EquationIdVectorNormalElement(rResult);
        else
            EquationIdVectorKuttaElement(rResult);
    }
    else //wake element
    {
        if (rResult.size() != 2 * NumNodes)
            rResult.resize(2 * NumNodes, false);

        EquationIdVectorWakeElement(rResult);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetDofList(
    DofsVectorType &rElementalDofList,
    ProcessInfo &CurrentProcessInfo)
{
    if (this->IsNot(MARKER)) //normal element
    {
        if (rElementalDofList.size() != NumNodes)
            rElementalDofList.resize(NumNodes);

        const IncompressiblePotentialFlowElement &r_this = *this;
        const int &kutta = r_this.GetValue(KUTTA);

        if (kutta == 0)
            GetDofListNormalElement(rElementalDofList);
        else
            GetDofListKuttaElement(rElementalDofList);
    }
    else //wake element
    {
        if (rElementalDofList.size() != 2 * NumNodes)
            rElementalDofList.resize(2 * NumNodes);

        GetDofListWakeElement(rElementalDofList);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetWakeDistances(array_1d<double, NumNodes> &distances)
{
    noalias(distances) = GetValue(ELEMENTAL_DISTANCES);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::EquationIdVectorNormalElement(
    EquationIdVectorType &rResult)
{
    for (unsigned int i = 0; i < NumNodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(POSITIVE_FACE_PRESSURE).EquationId();
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::EquationIdVectorKuttaElement(
    EquationIdVectorType &rResult)
{
    //Kutta elements have only negative part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (!GetGeometry()[i].FastGetSolutionStepValue(TRAILING_EDGE))
            rResult[i] = GetGeometry()[i].GetDof(POSITIVE_FACE_PRESSURE).EquationId();
        else
            rResult[i] = GetGeometry()[i].GetDof(NEGATIVE_FACE_PRESSURE).EquationId();
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::EquationIdVectorWakeElement(
    EquationIdVectorType &rResult)
{
    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    //positive part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] > 0.0)
            rResult[i] = GetGeometry()[i].GetDof(POSITIVE_FACE_PRESSURE).EquationId();
        else
            rResult[i] = GetGeometry()[i].GetDof(NEGATIVE_FACE_PRESSURE, 0).EquationId();
    }

    //negative part - sign is opposite to the previous case
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] < 0.0)
            rResult[NumNodes + i] = GetGeometry()[i].GetDof(POSITIVE_FACE_PRESSURE).EquationId();
        else
            rResult[NumNodes + i] = GetGeometry()[i].GetDof(NEGATIVE_FACE_PRESSURE).EquationId();
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetDofListNormalElement(
    DofsVectorType &rElementalDofList)
{
    for (unsigned int i = 0; i < NumNodes; i++)
        rElementalDofList[i] = GetGeometry()[i].pGetDof(POSITIVE_FACE_PRESSURE);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetDofListKuttaElement(
    DofsVectorType &rElementalDofList)
{
    //Kutta elements have only negative part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (!GetGeometry()[i].FastGetSolutionStepValue(TRAILING_EDGE))
            rElementalDofList[i] = GetGeometry()[i].pGetDof(POSITIVE_FACE_PRESSURE);
        else
            rElementalDofList[i] = GetGeometry()[i].pGetDof(NEGATIVE_FACE_PRESSURE);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetDofListWakeElement(
    DofsVectorType &rElementalDofList)
{
    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    //positive part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] > 0)
            rElementalDofList[i] = GetGeometry()[i].pGetDof(POSITIVE_FACE_PRESSURE);
        else
            rElementalDofList[i] = GetGeometry()[i].pGetDof(NEGATIVE_FACE_PRESSURE);
    }

    //negative part - sign is opposite to the previous case
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] < 0)
            rElementalDofList[NumNodes + i] = GetGeometry()[i].pGetDof(POSITIVE_FACE_PRESSURE);
        else
            rElementalDofList[NumNodes + i] = GetGeometry()[i].pGetDof(NEGATIVE_FACE_PRESSURE);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class IncompressiblePotentialFlowElement<2, 3>;

} // namespace Kratos
