#include "swimming_DEM_application.h"
#include "calculate_laplacian_simplex_element.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
void ComputeLaplacianSimplex<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                              ProcessInfo& rCurrentProcessInfo)
{

    const unsigned int NumNodes(TDim+1), LocalSize(TDim * NumNodes);
    unsigned int LocalIndex = 0;
    unsigned int lappos = this->GetGeometry()[0].GetDofPosition(VELOCITY_LAPLACIAN_X);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_LAPLACIAN_X,lappos).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_LAPLACIAN_Y,lappos+1).EquationId();
        if (TDim == 3){
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_LAPLACIAN_Z,lappos+2).EquationId();
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeLaplacianSimplex<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                        ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes(TDim+1), LocalSize(TDim * NumNodes);

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_Y);
        if (TDim == 3){
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_Z);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
int ComputeLaplacianSimplex<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Perform basic element checks
    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    if(ErrorCode != 0) return ErrorCode;

    if(this->GetGeometry().size() != TDim+1)
        KRATOS_THROW_ERROR(std::invalid_argument,"wrong number of nodes for element", this->Id());

    if(VELOCITY_LAPLACIAN.Key() == 0)

        KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY_LAPLACIAN Key is 0. Check if the application was correctly registered.","");

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
    {
        if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY_LAPLACIAN) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY_LAPLACIAN variable on solution step data for node ",this->GetGeometry()[i].Id());
    }
    return 0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeLaplacianSimplex<TDim, TNumNodes>::AddIntegrationPointRHSContribution(VectorType& F,
                             const BoundedMatrix<double, TNumNodes, TDim>& rShapeDeriv,
                             const double Weight)
{
    double Coef = Weight;

    int LocalIndex = 0;
    for (unsigned int iNodeB = 0; iNodeB < TNumNodes; ++iNodeB){

        for (unsigned int dj = 0; dj < TDim; ++dj){
            double value = 0.0;

            for (unsigned int iNodeA = 0; iNodeA < TNumNodes; ++iNodeA){
                const array_1d<double, TDim >& Velocity = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY);

                for (unsigned int di = 0; di < TDim; ++di){
                    value -= rShapeDeriv(iNodeB, di) * Velocity[dj] * rShapeDeriv(iNodeA, di);
                }
            }

            F[LocalIndex++] += Coef * value;
        }
    }
}

// Explicit instantiations
template class ComputeLaplacianSimplex<2, 3>;
template class ComputeLaplacianSimplex<3, 4>;
} // namespace Kratos
