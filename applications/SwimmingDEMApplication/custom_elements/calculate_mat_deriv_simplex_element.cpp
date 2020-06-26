#include "swimming_DEM_application.h"
#include "calculate_mat_deriv_simplex_element.h"

namespace Kratos
{

/// Calculate the element's local contribution to the system for the current step.
template <unsigned int TDim, unsigned int TNumNodes>
void ComputeMaterialDerivativeSimplex<TDim, TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes(TDim+1), LocalSize(TDim * NumNodes);

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

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeMaterialDerivativeSimplex<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                              ProcessInfo& rCurrentProcessInfo)
{

    const unsigned int NumNodes(TDim+1), LocalSize(TDim * NumNodes);
    unsigned int LocalIndex = 0;
    unsigned int local_position = this->GetGeometry()[0].GetDofPosition(MATERIAL_ACCELERATION_X);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(MATERIAL_ACCELERATION_X, local_position).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(MATERIAL_ACCELERATION_Y, local_position + 1).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(MATERIAL_ACCELERATION_Z, local_position + 2).EquationId();
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeMaterialDerivativeSimplex<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                        ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes(TDim+1), LocalSize(TDim * NumNodes);

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(MATERIAL_ACCELERATION_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(MATERIAL_ACCELERATION_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(MATERIAL_ACCELERATION_Z);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
int ComputeMaterialDerivativeSimplex<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Perform basic element checks
    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    if(ErrorCode != 0) return ErrorCode;

    if(this->GetGeometry().size() != TDim+1)
        KRATOS_THROW_ERROR(std::invalid_argument,"wrong number of nodes for element",this->Id());

    if (MATERIAL_ACCELERATION.Key() == 0){
        KRATOS_THROW_ERROR(std::invalid_argument,"MATERIAL_ACCELERATION Key is 0. Check if the application was correctly registered.","");
    }

    // Checks on nodes

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
    {
        if(this->GetGeometry()[i].SolutionStepsDataHas(MATERIAL_ACCELERATION) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing MATERIAL_ACCELERATION variable on solution step data for node ",this->GetGeometry()[i].Id());
    }
    return 0;

    KRATOS_CATCH("");
}


template <unsigned int TDim, unsigned int TNumNodes>
void ComputeMaterialDerivativeSimplex<TDim, TNumNodes>::CalculateLumpedMassMatrix(MatrixType& rLHSMatrix,
                               const double Mass)
{
    unsigned int DofIndex = 0;
    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
    {
        for (unsigned int d = 0; d < TDim; ++d)
        {
            rLHSMatrix(DofIndex, DofIndex) += Mass;
            ++DofIndex;
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeMaterialDerivativeSimplex<TDim, TNumNodes>::AddConsistentMassMatrixContribution(MatrixType& rLHSMatrix,
        const array_1d<double,TNumNodes>& rShapeFunc,
        const double Weight)
{
    const unsigned int BlockSize = TDim;

    double Coef = Weight;
    unsigned int FirstRow(0), FirstCol(0);
    double K; // Temporary results

    // Note: Dof order is (vx,vy,[vz,]p) for each node
    for (unsigned int i = 0; i < TNumNodes; ++i)
    {
        // Loop over columns
        for (unsigned int j = 0; j < TNumNodes; ++j)
        {
            K = Coef * rShapeFunc[i] * rShapeFunc[j];

            for (unsigned int d = 0; d < TDim; ++d) // iterate over dimensions for velocity Dofs in this node combination
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

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeMaterialDerivativeSimplex<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int LocalSize = TDim * TNumNodes;

    // Resize and set to zero
    if (rMassMatrix.size1() != LocalSize)
        rMassMatrix.resize(LocalSize, LocalSize, false);

    rMassMatrix = ZeroMatrix(LocalSize, LocalSize);

    // Get the element's geometric parameters
    double Area;
    array_1d<double, TNumNodes> N;
    BoundedMatrix<double, TNumNodes, TDim> DN_DX;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

    // Add 'classical' mass matrix (lumped)
    if (rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX] == 1){
        double Coeff = Area / TNumNodes; //Optimize!
        this->CalculateLumpedMassMatrix(rMassMatrix, Coeff);
    }

    else {
        // Add 'consistent' mass matrix
        MatrixType NContainer;
        ShapeFunctionDerivativesArrayType DN_DXContainer;
        VectorType GaussWeights;
        this->CalculateWeights(DN_DXContainer, NContainer, GaussWeights);
        const SizeType NumGauss = NContainer.size1();

        for (SizeType g = 0; g < NumGauss; g++){
            const double GaussWeight = GaussWeights[g];
            const ShapeFunctionsType& Ng = row(NContainer, g);
            this->AddConsistentMassMatrixContribution(rMassMatrix, Ng, GaussWeight);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeMaterialDerivativeSimplex<TDim, TNumNodes>::CalculateRHS(VectorType& F, ProcessInfo& rCurrentProcessInfo)
{
    // Get the element's geometric parameters
    double Area;
    array_1d<double, TNumNodes> N;
    BoundedMatrix<double, TNumNodes, TDim> DN_DX;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

    MatrixType NContainer;
    ShapeFunctionDerivativesArrayType DN_DXContainer;
    VectorType GaussWeights;
    this->CalculateWeights(DN_DXContainer, NContainer, GaussWeights);
    const SizeType NumGauss = NContainer.size1();

    for (SizeType g = 0; g < NumGauss; g++){
        const double GaussWeight = GaussWeights[g];
        const ShapeFunctionsType& Ng = row(NContainer, g);
        this->AddIntegrationPointRHSContribution(F, Ng, DN_DX, GaussWeight);
    }

}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeMaterialDerivativeSimplex<TDim, TNumNodes>::CalculateWeights(ShapeFunctionDerivativesArrayType& rDN_DX,
        Matrix& rNContainer,
        Vector& rGaussWeights)
{
    const GeometryType& rGeom = this->GetGeometry();
    Vector DetJ;
    rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, GeometryData::GI_GAUSS_2);
    rNContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);

    rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2), false);

    for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2); g++)
        rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeMaterialDerivativeSimplex<TDim, TNumNodes>::EvaluateInPoint(array_1d< double, 3 > & rResult,
                             const Variable< array_1d< double, 3 > >& rVariable,
                             const array_1d< double, TNumNodes >& rShapeFunc)
{
    // Compute the weighted value of the nodal variable in the (Gauss) Point
    noalias(rResult) = rShapeFunc[0] * this->GetGeometry()[0].FastGetSolutionStepValue(rVariable);
    for (unsigned int iNode = 1; iNode < TNumNodes; ++iNode)
        noalias(rResult) += rShapeFunc[iNode] * this->GetGeometry()[iNode].FastGetSolutionStepValue(rVariable);
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeMaterialDerivativeSimplex<TDim, TNumNodes>::AddIntegrationPointRHSContribution(VectorType& F,
                             const array_1d<double, TNumNodes>& rShapeFunc,
                             const BoundedMatrix<double, TNumNodes, TDim>& rShapeDeriv,
                             const double Weight)
{
    double Coef = Weight;
    array_1d<double, 3 > Velocity;
    this->EvaluateInPoint(Velocity, VELOCITY, rShapeFunc);
    int LocalIndex = 0;
    for (unsigned int iNodeB = 0; iNodeB < TNumNodes; ++iNodeB){

        for (unsigned int dj = 0; dj < TDim; ++dj){
            double value = 0.0;

            for (unsigned int iNodeA = 0; iNodeA < TNumNodes; ++iNodeA){
                const array_1d<double, 3 >& NodalVelocity = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY);

                for (unsigned int di = 0; di < TDim; ++di){
                    value += rShapeFunc[iNodeB] * Velocity[di] * rShapeDeriv(iNodeA, di) * NodalVelocity[dj];
                }
            }

            F[LocalIndex++] += Coef * value;
        }
    }
}
// Explicit instantiations
template class ComputeMaterialDerivativeSimplex<2, 3>;
template class ComputeMaterialDerivativeSimplex<3, 4>;
} // namespace Kratos
