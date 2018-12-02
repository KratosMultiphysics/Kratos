#include "swimming_DEM_application.h"
#include "calculate_velocity_laplacian_component.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeVelocityLaplacianComponentSimplex<TDim, TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  ProcessInfo& rCurrentProcessInfo)
{
    const int current_component = rCurrentProcessInfo[CURRENT_COMPONENT];

    if (current_component == 0){
        mCurrentComponent = 'X';
    }

    else if (current_component == 1){
        mCurrentComponent = 'Y';
    }

    else if (current_component == 2){
        mCurrentComponent = 'Z';
    }

    else {
        KRATOS_THROW_ERROR(std::invalid_argument, "The value of CURRENT_COMPONENT passed to the ComputeVelocityLaplacianComponentSimplex element is not 0, 1 or 2, but ", current_component);
    }

    const unsigned int NumNodes(TDim+1), LocalSize(NumNodes);

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

    this->CalculateMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);

    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);

    const double volume_inv = 1.0 / this->GetGeometry().Volume();
    for (unsigned int i=0; i<LocalSize; ++i){
        for (unsigned int j=0; j<LocalSize; ++j){
            rLeftHandSideMatrix(i, j) *= volume_inv;
        }
        rRightHandSideVector(i) *= volume_inv;
    }

    if (this->Id() == 1000){
        KRATOS_WATCH(rLeftHandSideMatrix)
        KRATOS_WATCH(rRightHandSideVector)
    }
}

//template <unsigned int TDim, unsigned int TNumNodes>
//Variable<array_1d<double, 3> >  ComputeVelocityLaplacianComponentSimplex<TDim, TNumNodes>::GetCurrentVariable()
//{
//    if (mCurrentComponent == 'X'){
//        return(VELOCITY_LAPLACIAN_X);
//    }

//    else if (mCurrentComponent == 'Y'){
//        return(VELOCITY_LAPLACIAN_Y);
//    }

//    else {
//        return(VELOCITY_LAPLACIAN_Z);
//    }
//}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeVelocityLaplacianComponentSimplex<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                              ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != TNumNodes)
        rResult.resize(TNumNodes, false);

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode){
        rResult[iNode] = this->GetGeometry()[iNode].GetDof(VELOCITY_LAPLACIAN_Z).EquationId();

        if (this->Id()==1000)
        KRATOS_WATCH(rResult[iNode] )
    }

}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeVelocityLaplacianComponentSimplex<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                        ProcessInfo& rCurrentProcessInfo)
{
    if (rElementalDofList.size() != TNumNodes)
        rElementalDofList.resize(TNumNodes);

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode){
        rElementalDofList[iNode] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_Z);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
int ComputeVelocityLaplacianComponentSimplex<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Perform basic element checks
    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    if(ErrorCode != 0) return ErrorCode;

    if(this->GetGeometry().size() != TDim+1)
        KRATOS_THROW_ERROR(std::invalid_argument,"wrong number of nodes per element",this->Id());

    if(VELOCITY_LAPLACIAN_Z.Key() == 0)

        KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY_LAPLACIAN_Z Key is 0. Check if the application was correctly registered.","");

    // Checks on nodes

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i < this->GetGeometry().size(); ++i)
    {
        if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY_LAPLACIAN_Z) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY_LAPLACIAN_Z variable on solution step data for node ",this->GetGeometry()[i].Id());
    }
    return 0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeVelocityLaplacianComponentSimplex<TDim, TNumNodes>::CalculateLumpedMassMatrix(MatrixType& rLHSMatrix,
                               const double Mass)
{
    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode){
        rLHSMatrix(iNode, iNode) += Mass;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeVelocityLaplacianComponentSimplex<TDim, TNumNodes>::AddConsistentMassMatrixContribution(MatrixType& rLHSMatrix,
        const array_1d<double,TNumNodes>& rShapeFunc,
        const double Weight)
{
    for (unsigned int i = 0; i < TNumNodes; ++i){
        for (unsigned int j = 0; j < TNumNodes; ++j){
            rLHSMatrix(i, j) += Weight * rShapeFunc[i] * rShapeFunc[j];
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeVelocityLaplacianComponentSimplex<TDim, TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int LocalSize = TNumNodes;

    // Resize and set to zero
    if (rMassMatrix.size1() != LocalSize)
        rMassMatrix.resize(LocalSize, LocalSize, false);

    rMassMatrix = ZeroMatrix(LocalSize, LocalSize);

    // Get the element's geometric parameters
    double Area;
    array_1d<double, TNumNodes> N;
    BoundedMatrix<double, TNumNodes, TDim> DN_DX;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

    if (rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX] == 1){ // Add 'classical' mass matrix (lumped)
        double Coeff = Area / TNumNodes;
        this->CalculateLumpedMassMatrix(rMassMatrix, Coeff);
    }

    else { // Add 'consistent' mass matrix
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
void ComputeVelocityLaplacianComponentSimplex<TDim, TNumNodes>::AddIntegrationPointRHSContribution(VectorType& F,
                             const array_1d<double, TNumNodes>& rShapeFunc,
                             const BoundedMatrix<double, TNumNodes, TDim>& rShapeDeriv,
                             const double Weight)
{
    for (unsigned int iNodeB = 0; iNodeB < TNumNodes; ++iNodeB){
        double div = 0.0;

        for (unsigned int iNodeA = 0; iNodeA < TNumNodes; ++iNodeA){
            if (mCurrentComponent == 'X'){
                const array_1d<double, 3>& NodalComponent = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY_X_GRADIENT);
                for (unsigned int di = 0; di < TDim; ++di){
                    div += rShapeDeriv(iNodeA, di) * NodalComponent[di];
                }
            }

            else if (mCurrentComponent == 'Y'){
                const array_1d<double, 3>& NodalComponent = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY_Y_GRADIENT);
                for (unsigned int di = 0; di < TDim; ++di){
                    div += rShapeDeriv(iNodeA, di) * NodalComponent[di];
                }
            }

            else {
                const array_1d<double, 3>& NodalComponent = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY_Z_GRADIENT);
                for (unsigned int di = 0; di < TDim; ++di){
                    div += rShapeDeriv(iNodeA, di) * NodalComponent[di];
                }
            }
        }

        F[iNodeB] += Weight * div * rShapeFunc[iNodeB];
    }
}

// Explicit instantiations
template class ComputeVelocityLaplacianComponentSimplex<2, 3>;
template class ComputeVelocityLaplacianComponentSimplex<3, 4>;
} // namespace Kratos
