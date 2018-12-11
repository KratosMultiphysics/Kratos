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
    if (this->IsNot(MARKER)) // Normal element (non-wake) - eventually an embedded
        CalculateLocalSystemNormalElement(rLeftHandSideMatrix, rRightHandSideVector);
    else // wake element
        CalculateLocalSystemWakeElement(rLeftHandSideMatrix, rRightHandSideVector);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::EquationIdVector(
    EquationIdVectorType &rResult,
    ProcessInfo &CurrentProcessInfo)
{
    if (this->IsNot(MARKER)) // Normal element
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
    else // Wake element
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
    if (this->IsNot(MARKER)) //Normal element
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
    else // wake element
    {
        if (rElementalDofList.size() != 2 * NumNodes)
            rElementalDofList.resize(2 * NumNodes);

        GetDofListWakeElement(rElementalDofList);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::FinalizeSolutionStep(
    ProcessInfo &rCurrentProcessInfo)
    {
        bool active = true;
        if ((this)->IsDefined(ACTIVE))
            active = (this)->Is(ACTIVE);

        if (this->Is(MARKER) && active == true)
        {
            CheckWakeCondition();
            ComputePotentialJump(rCurrentProcessInfo);
        }
        ComputeElementInternalEnergy();
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
    // Kutta elements have only negative part
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

    // Positive part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] > 0.0)
            rResult[i] = GetGeometry()[i].GetDof(POSITIVE_FACE_PRESSURE).EquationId();
        else
            rResult[i] = GetGeometry()[i].GetDof(NEGATIVE_FACE_PRESSURE, 0).EquationId();
    }

    // Negative part - sign is opposite to the previous case
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
    // Kutta elements have only negative part
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

    // Positive part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] > 0)
            rElementalDofList[i] = GetGeometry()[i].pGetDof(POSITIVE_FACE_PRESSURE);
        else
            rElementalDofList[i] = GetGeometry()[i].pGetDof(NEGATIVE_FACE_PRESSURE);
    }

    // Negative part - sign is opposite to the previous case
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
    const IncompressiblePotentialFlowElement &r_this = *this;
    const int &kutta = r_this.GetValue(KUTTA);

    if (kutta == 0)
        for (unsigned int i = 0; i < NumNodes; i++)
            phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
    else
        for (unsigned int i = 0; i < NumNodes; i++)
            if (!GetGeometry()[i].FastGetSolutionStepValue(TRAILING_EDGE))
                phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
            else
                phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::ComputeLHSGaussPointContribution(
    const double weight,
    Matrix &lhs,
    const ElementalData<NumNodes, Dim> &data)
{
    noalias(lhs) += weight * prod(data.DN_DX, trans(data.DN_DX));
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystemNormalElement(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector)
{
    if (rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
        rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
    if (rRightHandSideVector.size() != NumNodes)
        rRightHandSideVector.resize(NumNodes, false);
    rLeftHandSideMatrix.clear();

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    ComputeLHSGaussPointContribution(data.vol, rLeftHandSideMatrix, data);

    GetPotentialOnNormalElement(data.phis);
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, data.phis);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystemWakeElement(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector)
{
    // Note that the lhs and rhs have double the size
    if (rLeftHandSideMatrix.size1() != 2 * NumNodes || rLeftHandSideMatrix.size2() != 2 * NumNodes)
        rLeftHandSideMatrix.resize(2 * NumNodes, 2 * NumNodes, false);
    if (rRightHandSideVector.size() != 2 * NumNodes)
        rRightHandSideVector.resize(2 * NumNodes, false);
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    GetWakeDistances(data.distances);

    Matrix lhs_total = ZeroMatrix(NumNodes, NumNodes);

    ComputeLHSGaussPointContribution(data.vol, lhs_total, data);

    if (this->Is(STRUCTURE))
    {
        Matrix lhs_positive = ZeroMatrix(NumNodes, NumNodes);
        Matrix lhs_negative = ZeroMatrix(NumNodes, NumNodes);

        CalculateLocalSystemSubdividedElement(lhs_positive, lhs_negative);
        AssignLocalSystemSubdividedElement(rLeftHandSideMatrix, lhs_positive, lhs_negative, lhs_total, data);
    }
    else
        AssignLocalSystemWakeElement(rLeftHandSideMatrix, lhs_total, data);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystemSubdividedElement(
    Matrix &lhs_positive,
    Matrix &lhs_negative)
{
    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    GetWakeDistances(data.distances);

    // Subdivide the element
    constexpr unsigned int nvolumes = 3 * (Dim - 1);
    bounded_matrix<double, NumNodes, Dim> Points;
    array_1d<double, nvolumes> PartitionsSign;
    bounded_matrix<double, nvolumes, NumNodes> GPShapeFunctionValues;
    array_1d<double, nvolumes> Volumes;
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

    // Compute the lhs and rhs that would correspond to it being divided
    for (unsigned int i = 0; i < nsubdivisions; ++i)
    {
        if (PartitionsSign[i] > 0)
            ComputeLHSGaussPointContribution(Volumes[i], lhs_positive, data);
        else
            ComputeLHSGaussPointContribution(Volumes[i], lhs_negative, data);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::AssignLocalSystemSubdividedElement(
    MatrixType &rLeftHandSideMatrix,
    Matrix &lhs_positive,
    Matrix &lhs_negative,
    Matrix &lhs_total,
    const ElementalData<NumNodes, Dim> &data)
{
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        // The TE node takes the contribution of the subdivided element and
        // we do not apply the wake condition on the TE node
        if(GetGeometry()[i].FastGetSolutionStepValue(TRAILING_EDGE))
        {
            for (unsigned int j = 0; j < NumNodes; ++j)
            {
                rLeftHandSideMatrix(i, j) = lhs_positive(i, j);
                rLeftHandSideMatrix(i + NumNodes, j + NumNodes) = lhs_negative(i, j);
            }
        }
        else
            AssignLocalSystemWakeNode(rLeftHandSideMatrix, lhs_total, data, i);
    }
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::AssignLocalSystemWakeElement(
    MatrixType &rLeftHandSideMatrix,
    Matrix &lhs_total,
    const ElementalData<NumNodes, Dim> &data)
{
    for (unsigned int row = 0; row < NumNodes; ++row)
        AssignLocalSystemWakeNode(rLeftHandSideMatrix, lhs_total, data, row);
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::AssignLocalSystemWakeNode(
    MatrixType &rLeftHandSideMatrix,
    Matrix &lhs_total,
    const ElementalData<NumNodes, Dim> &data,
    unsigned int &row)
{
    // Filling the diagonal blocks (i.e. decoupling upper and lower dofs)
    for (unsigned int column = 0; column < NumNodes; ++column)
    {
        rLeftHandSideMatrix(row, column) = lhs_total(row, column);
        rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_total(row, column);
    }
    
    // Applying wake condition on the NEGATIVE_FACE_PRESSURE dofs
    if (data.distances[row] < 0.0)
        for (unsigned int column = 0; column < NumNodes; ++column)
            rLeftHandSideMatrix(row, column + NumNodes) = -lhs_total(row, column); //Side 1
    else if (data.distances[row] > 0.0)
        for (unsigned int column = 0; column < NumNodes; ++column)
            rLeftHandSideMatrix(row + NumNodes, column) = -lhs_total(row, column); //Side 2
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::CheckWakeCondition()
{
    array_1d<double, Dim> upper_wake_velocity;
    ComputeVelocityUpperWakeElement(upper_wake_velocity);
    const double vupnorm = inner_prod(upper_wake_velocity, upper_wake_velocity);

    array_1d<double, Dim> lower_wake_velocity;
    ComputeVelocityLowerWakeElement(lower_wake_velocity);
    const double vlownorm = inner_prod(lower_wake_velocity, lower_wake_velocity);

    if (std::abs(vupnorm - vlownorm) > 0.1)
        std::cout << "WAKE CONDITION NOT FULFILLED IN ELEMENT # " << this->Id() << std::endl;
}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::ComputePotentialJump(
    ProcessInfo &rCurrentProcessInfo)
{
    const array_1d<double, 3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
    const double vinfinity_norm = sqrt(inner_prod(vinfinity, vinfinity));

    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    for (unsigned int i = 0; i < NumNodes; i++)
        if (distances[i] > 0)
            GetGeometry()[i].GetSolutionStepValue(POTENTIAL_JUMP) = 2.0 / vinfinity_norm * (GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE) - GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE));
        else
            GetGeometry()[i].GetSolutionStepValue(POTENTIAL_JUMP) = 2.0 / vinfinity_norm * (GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE) - GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE));

}

template <int Dim, int NumNodes>
void IncompressiblePotentialFlowElement<Dim, NumNodes>::ComputeElementInternalEnergy()
{
    double internal_energy = 0.0;
    array_1d<double, Dim> velocity;

    if (this->IsNot(MARKER)) // Normal element (non-wake) - eventually an embedded
        ComputeVelocityNormalElement(velocity);
    else // Wake element
        ComputeVelocityUpperWakeElement(velocity);

    internal_energy = 0.5 * inner_prod(velocity, velocity);
    this->SetValue(INTERNAL_ENERGY, abs(internal_energy));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class IncompressiblePotentialFlowElement<2, 3>;

} // namespace Kratos
