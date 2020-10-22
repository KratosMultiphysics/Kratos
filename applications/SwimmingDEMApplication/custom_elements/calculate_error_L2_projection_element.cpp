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
    const unsigned int LocalSize(Dim * NumNodes);

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

    CalculateMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);

}

template <unsigned int Dim, unsigned int NumNodes>
void CalculateErrorL2Projection<Dim, NumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{

    const unsigned int LocalSize = Dim * NumNodes;
    unsigned int LocalIndex = 0;
    unsigned int local_position = this->GetGeometry()[0].GetDofPosition(ERROR_X);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(ERROR_X, local_position).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(ERROR_Y, local_position + 1).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(ERROR_Z, local_position + 2).EquationId();
    }
}

template <unsigned int Dim, unsigned int NumNodes>
void CalculateErrorL2Projection<Dim, NumNodes>::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int LocalSize = Dim * NumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(ERROR_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(ERROR_Y);
        if(Dim == 3) rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(ERROR_Z);
    }
}

template <unsigned int Dim, unsigned int NumNodes>
void CalculateErrorL2Projection<Dim, NumNodes>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int LocalSize = Dim * NumNodes;

    // Resize and set to zero
    if (rMassMatrix.size1() != LocalSize)
        rMassMatrix.resize(LocalSize, LocalSize, false);

    rMassMatrix = ZeroMatrix(LocalSize, LocalSize);

    // Get the element's geometric parameters
    array_1d<double, NumNodes> N;

    // Add 'consistent' mass matrix

    MatrixType NContainer = this->GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = this->GetGeometry().IntegrationPoints(GeometryData::GI_GAUSS_2);

    const SizeType NumGauss = NContainer.size1();
    Vector DetJ = ZeroVector(NumGauss);
    this->GetGeometry().DeterminantOfJacobian(DetJ, GeometryData::GI_GAUSS_2);

    for (SizeType g = 0; g < NumGauss; g++){
            const double GaussWeight = DetJ[g] * IntegrationPoints[g].Weight();
            const ShapeFunctionsType Ng = row(NContainer, g);
            this->AddConsistentMassMatrixContribution(rMassMatrix, Ng, GaussWeight);
        }

}

template <unsigned int Dim, unsigned int NumNodes>
void CalculateErrorL2Projection<Dim, NumNodes>::AddConsistentMassMatrixContribution(
    MatrixType& rLHSMatrix,
    const array_1d<double,NumNodes>& rShapeFunc,
    const double Weight)
{
    const unsigned int BlockSize = Dim;

    double Coef = Weight;
    unsigned int FirstRow(0), FirstCol(0);
    double K; // Temporary results

    // Note: Dof order is (vx,vy,[vz,]p) for each node
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        // Loop over columns
        for (unsigned int j = 0; j < NumNodes; ++j)
        {
            K = Coef * rShapeFunc[i] * rShapeFunc[j];
            for (unsigned int d = 0; d < Dim; ++d) // iterate over dimensions for velocity Dofs in this node combination
            {
                rLHSMatrix(FirstRow + d, FirstCol + d) += K;
            }
            // Update column index
            FirstCol += BlockSize;
        }
        // Update matrix indices
        FirstRow += BlockSize;
        FirstCol = 0;
    }
}

template <unsigned int Dim, unsigned int NumNodes>
void CalculateErrorL2Projection<Dim, NumNodes>::AddIntegrationPointRHSContribution(
    VectorType& rRHSVector,
    const array_1d<double,NumNodes>& rShapeFunc,
    const double Weight)
{
    double Coef = Weight;
    int LocalIndex = 0;

    for (unsigned int iNodeB = 0; iNodeB < NumNodes; ++iNodeB){
        for (unsigned int dj = 0; dj < Dim; ++dj){

            double value = 0.0;

            for (unsigned int iNodeA = 0; iNodeA < NumNodes; ++iNodeA){

                Vector NodalComponent = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VECTORIAL_ERROR);
                value += rShapeFunc[iNodeB] * NodalComponent[dj];

            }
            rRHSVector[LocalIndex++] += Coef * value;
        }
    }
}

template <unsigned int Dim, unsigned int NumNodes>
void CalculateErrorL2Projection<Dim, NumNodes>::CalculateRHS(
    VectorType& rRHSVector,
    ProcessInfo& rCurrentProcessInfo)
{

    array_1d<double, NumNodes> N;
    MatrixType NContainer = this->GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
    VectorType GaussWeights;
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