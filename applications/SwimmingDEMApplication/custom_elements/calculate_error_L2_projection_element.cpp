//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//

// Project includes

// Aplication includes
#include "swimming_DEM_application.h"
#include "calculate_error_L2_projection_element.h"

namespace Kratos
{

template <unsigned int Dim, unsigned int NumNodes>
void CalculateErrorL2Projection<Dim, NumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int LocalSize((Dim+1) * NumNodes);

    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    for (unsigned int i=0; i<LocalSize; ++i){
        for (unsigned int j=0; j<LocalSize; ++j){
            rLeftHandSideMatrix(i, j) = 0.0;
        }
        rRightHandSideVector(i) = 0.0;
    }

    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

}

template <unsigned int Dim, unsigned int NumNodes>
void CalculateErrorL2Projection<Dim, NumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{

    const unsigned int LocalSize = (Dim+1) * NumNodes;
    unsigned int LocalIndex = 0;
    unsigned int local_position = this->GetGeometry()[0].GetDofPosition(ERROR_X);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(ERROR_X, local_position).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(ERROR_Y, local_position + 1).EquationId();
        if(Dim == 3)
        {
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(ERROR_Z, local_position + 2).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(ERROR_P, local_position + 3).EquationId();
        }
        else if(Dim == 2) rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(ERROR_P, local_position + 2).EquationId();
    }
}

template <unsigned int Dim, unsigned int NumNodes>
void CalculateErrorL2Projection<Dim, NumNodes>::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int LocalSize = (Dim+1) * NumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(ERROR_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(ERROR_Y);
        if(Dim == 3) rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(ERROR_Z);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(ERROR_P);
    }
}

template <unsigned int Dim, unsigned int NumNodes>
void CalculateErrorL2Projection<Dim, NumNodes>::CalculateLeftHandSide(
    MatrixType& rLHSMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int LocalSize = (Dim+1) * NumNodes;

    // Get the element's geometric parameters
    array_1d<double, NumNodes> N;
    rLHSMatrix = ZeroMatrix(LocalSize, LocalSize);
    MatrixType NContainer = this->GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = this->GetGeometry().IntegrationPoints(GeometryData::GI_GAUSS_2);

    const SizeType NumGauss = NContainer.size1();
    Vector DetJ = ZeroVector(NumGauss);
    this->GetGeometry().DeterminantOfJacobian(DetJ, GeometryData::GI_GAUSS_2);

    for (SizeType g = 0; g < NumGauss; g++){
            const double GaussWeight = DetJ[g] * IntegrationPoints[g].Weight();
            const ShapeFunctionsType Ng = row(NContainer, g);
            this->AddLHSMatrixContribution(rLHSMatrix, Ng, GaussWeight);
        }
}

template <unsigned int Dim, unsigned int NumNodes>
void CalculateErrorL2Projection<Dim, NumNodes>::AddLHSMatrixContribution(
    MatrixType& rLHSMatrix,
    const array_1d<double,NumNodes>& N,
    const double Weight)
{
    const unsigned int BlockSize = NumNodes;

    double Coef = Weight;
    unsigned int FirstRow(0), FirstCol(0);
    double K; // Temporary results

    // Note: Dof order is (vx,vy,[vz,]p) for each node
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        // Loop over columns
        unsigned int row = i*BlockSize;
        for (unsigned int j = 0; j < NumNodes; ++j)
        {
            unsigned int col = j*BlockSize;
            K = Coef * N[i] * N[j];
            for (unsigned int d = 0; d < NumNodes; ++d) // iterate over dimensions for velocity Dofs in this node combination
            {
                rLHSMatrix(row+d, col+d) += K;
            }
        }
    }
}

template <unsigned int Dim, unsigned int NumNodes>
void CalculateErrorL2Projection<Dim, NumNodes>::AddIntegrationPointRHSContribution(
    VectorType& rRHSVector,
    const array_1d<double,NumNodes>& N,
    const double Weight)
{
    double Coef = Weight;
    Vector NodalComponent = ZeroVector(Dim);
    double scalar_component = 0.0;

    for (unsigned int i = 0; i < NumNodes; ++i){
        int row = i * NumNodes;
        NodalComponent = this->GetGeometry()[i].FastGetSolutionStepValue(VECTORIAL_ERROR);
        scalar_component = this->GetGeometry()[i].FastGetSolutionStepValue(SCALAR_ERROR);
        for (unsigned int d = 0; d < Dim; ++d){
            rRHSVector[row+d] += Coef * N[i] * NodalComponent[d];
        }
        rRHSVector[row+Dim] += Coef * N[i] * scalar_component;
    }
}

template <unsigned int Dim, unsigned int NumNodes>
void CalculateErrorL2Projection<Dim, NumNodes>::CalculateRightHandSide(
    VectorType& rRHSVector,
    ProcessInfo& rCurrentProcessInfo)
{

    const unsigned int LocalSize = (Dim+1) * NumNodes;
    array_1d<double, NumNodes> N;
    MatrixType NContainer = this->GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
    VectorType GaussWeights;
    rRHSVector = ZeroVector(LocalSize);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = this->GetGeometry().IntegrationPoints(GeometryData::GI_GAUSS_2);

    const SizeType NumGauss = NContainer.size1();
    Vector DetJ = ZeroVector(NumGauss);
    this->GetGeometry().DeterminantOfJacobian(DetJ, GeometryData::GI_GAUSS_2);

    for (SizeType g = 0; g < NumGauss; g++){
        const double GaussWeight = DetJ[g] * IntegrationPoints[g].Weight();
        const ShapeFunctionsType& Ng = row(NContainer, g);
        this->AddIntegrationPointRHSContribution(rRHSVector, Ng, GaussWeight);
    }
}

// Explicit instantiations
template class CalculateErrorL2Projection<2, 3>;
template class CalculateErrorL2Projection<3, 4>;
} // namespace Kratos