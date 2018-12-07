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

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        ProcessInfo &rCurrentProcessInfo)
{
    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    GetPotentialOnNormalElement(data.phis);

    if (this->IsNot(MARKER)) //normal element (non-wake) - eventually an embedded
    {
        if (rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
            rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
        if (rRightHandSideVector.size() != NumNodes)
            rRightHandSideVector.resize(NumNodes, false);
        rLeftHandSideMatrix.clear();

        ComputeLHSGaussPointContribution(data.vol, rLeftHandSideMatrix, data);

        noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, data.phis);
    }
    else //it is a wake element
    {
        GetWakeDistances(data.distances);
        //note that the lhs and rhs have double the size!!
        if (rLeftHandSideMatrix.size1() != 2 * NumNodes || rLeftHandSideMatrix.size2() != 2 * NumNodes)
            rLeftHandSideMatrix.resize(2 * NumNodes, 2 * NumNodes, false);
        if (rRightHandSideVector.size() != 2 * NumNodes)
            rRightHandSideVector.resize(2 * NumNodes, false);
        rLeftHandSideMatrix.clear();
        if (this->Is(STRUCTURE))
        {
            //subdivide the element
            constexpr unsigned int nvolumes = 3 * (Dim - 1);
            bounded_matrix<double, NumNodes, Dim> Points;
            array_1d<double, nvolumes> Volumes;
            bounded_matrix<double, nvolumes, NumNodes> GPShapeFunctionValues;
            array_1d<double, nvolumes> PartitionsSign;
            std::vector<Matrix> GradientsValue(nvolumes);
            bounded_matrix<double, nvolumes, 2> NEnriched;
            for (unsigned int i = 0; i < GradientsValue.size(); ++i)
                GradientsValue[i].resize(2, Dim, false);
            for (unsigned int i = 0; i < NumNodes; ++i)
            {
                const array_1d<double, 3> &coords = GetGeometry()[i].Coordinates();
                for (unsigned int k = 0; k < Dim; ++k)
                {
                    Points(i, k) = coords[k];
                }
            }
            const unsigned int nsubdivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(Points,
                                                                                                   data.DN_DX,
                                                                                                   data.distances,
                                                                                                   Volumes,
                                                                                                   GPShapeFunctionValues,
                                                                                                   PartitionsSign,
                                                                                                   GradientsValue,
                                                                                                   NEnriched);
            //compute the lhs and rhs that would correspond to it being divided
            Matrix lhs_positive = ZeroMatrix(NumNodes, NumNodes);
            Matrix lhs_negative = ZeroMatrix(NumNodes, NumNodes);
            for (unsigned int i = 0; i < nsubdivisions; ++i)
            {
                if (PartitionsSign[i] > 0)
                    ComputeLHSGaussPointContribution(Volumes[i], lhs_positive, data);
                else
                    ComputeLHSGaussPointContribution(Volumes[i], lhs_negative, data);
            }
            Matrix lhs_total = ZeroMatrix(NumNodes, NumNodes);
            ComputeLHSGaussPointContribution(data.vol, lhs_total, data);
            for (unsigned int i = 0; i < NumNodes; ++i)
            {
                //No contribution to the TE node //and to extra dofs
                if(GetGeometry()[i].FastGetSolutionStepValue(TRAILING_EDGE))
                {
                    for (unsigned int j = 0; j < NumNodes; ++j)
                    {
                        rLeftHandSideMatrix(i, j) = lhs_positive(i, j);
                        rLeftHandSideMatrix(i, j + NumNodes) = 0.0;
                        rLeftHandSideMatrix(i + NumNodes, j + NumNodes) = lhs_negative(i, j);
                        rLeftHandSideMatrix(i + NumNodes, j) = 0.0;
                    }
                }
                else
                {
                    for (unsigned int j = 0; j < NumNodes; ++j)
                    {
                        rLeftHandSideMatrix(i, j) = lhs_total(i, j);
                        rLeftHandSideMatrix(i, j + NumNodes) = 0.0;
                        rLeftHandSideMatrix(i + NumNodes, j + NumNodes) = lhs_total(i, j);
                        rLeftHandSideMatrix(i + NumNodes, j) = 0.0;
                    }
                    //Applying wake condition
                    if (data.distances[i] < 0.0)
                    {
                        for (unsigned int j = 0; j < NumNodes; ++j)
                            rLeftHandSideMatrix(i, j + NumNodes) = -lhs_total(i, j);
                    }
                    else if (data.distances[i] > 0.0)
                    {
                        for (unsigned int j = 0; j < NumNodes; ++j)
                            rLeftHandSideMatrix(i + NumNodes, j) = -lhs_total(i, j);
                    }
                }
            }
        }
        else
        {
            Matrix lhs_total = ZeroMatrix(NumNodes, NumNodes);
            ComputeLHSGaussPointContribution(data.vol, lhs_total, data);
            //Looping over rows
            for (unsigned int i = 0; i < NumNodes; ++i)
            { //Looping over columngs
                for (unsigned int j = 0; j < NumNodes; ++j)
                { //Filling the diagonal blocks (i.e. decoupling upper and lower fields)
                    rLeftHandSideMatrix(i, j) = lhs_total(i, j);
                    rLeftHandSideMatrix(i, j + NumNodes) = 0.0;
                    rLeftHandSideMatrix(i + NumNodes, j + NumNodes) = lhs_total(i, j);
                    rLeftHandSideMatrix(i + NumNodes, j) = 0.0;
                }
                if (data.distances[i] < 0.0 && !GetGeometry()[i].FastGetSolutionStepValue(DEACTIVATED_WAKE))
                {                            //side1  -assign constraint only on the NEGATIVE_FACE_PRESSURE dofs and not on the airfoil nodes
                    //Marking nodes where the wake constraint is applied
                    GetGeometry()[i].GetSolutionStepValue(TEMPERATURE) = 30.0;
                    for (unsigned int j = 0; j < NumNodes; ++j)
                        rLeftHandSideMatrix(i, j + NumNodes) = -lhs_total(i, j);
                }
                else if (data.distances[i] > 0.0)// && !GetGeometry()[i].FastGetSolutionStepValue(DEACTIVATED_WAKE))
                { //side2 -assign constraint only on the NEGATIVE_FACE_PRESSURE dofs and not on the airfoil nodes
                    GetGeometry()[i].GetSolutionStepValue(TEMPERATURE) = 30.0;
                    for (unsigned int j = 0; j < NumNodes; ++j)
                        rLeftHandSideMatrix(i + NumNodes, j) = -lhs_total(i, j);
                }
            }
        }
        
        Vector split_element_values(NumNodes * 2);
        GetValuesOnSplitElement(split_element_values, data.distances);
        noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, split_element_values);
    }
}

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

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::GetPotentialOnNormalElement(
    array_1d<double, NumNodes> &phis)
{
    // Gather nodal data
    for (unsigned int i = 0; i < NumNodes; i++)
        phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class IncompressiblePotentialFlowElement<2, 3>;

} // namespace Kratos
