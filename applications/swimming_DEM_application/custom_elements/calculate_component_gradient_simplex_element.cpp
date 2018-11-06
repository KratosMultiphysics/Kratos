#include "swimming_DEM_application.h"
#include "calculate_component_gradient_simplex_element.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeComponentGradientSimplex<TDim, TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
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
        KRATOS_THROW_ERROR(std::invalid_argument, "The value of CURRENT_COMPONENT passed to the ComputeComponentGradientSimplex element is not 0, 1 or 2, but ", current_component);
    }

    BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeComponentGradientSimplex<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                              ProcessInfo& rCurrentProcessInfo)
{

    const unsigned int LocalSize(TDim * TNumNodes);
    unsigned int LocalIndex = 0;
    unsigned int pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_COMPONENT_GRADIENT_X);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_COMPONENT_GRADIENT_X,pos).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_COMPONENT_GRADIENT_Y,pos+1).EquationId();
        if (TDim == 3){
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_COMPONENT_GRADIENT_Z,pos+2).EquationId();
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeComponentGradientSimplex<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                        ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int LocalSize(TDim * TNumNodes);

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_COMPONENT_GRADIENT_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_COMPONENT_GRADIENT_Y);
        if (TDim == 3){
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_COMPONENT_GRADIENT_Z);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
int ComputeComponentGradientSimplex<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Perform basic element checks
    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    if(ErrorCode != 0) return ErrorCode;

    if(this->GetGeometry().size() != TDim+1)
        KRATOS_THROW_ERROR(std::invalid_argument,"wrong number of nodes for element",this->Id());

    if(VELOCITY_COMPONENT_GRADIENT.Key() == 0)

        KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY_COMPONENT_GRADIENT Key is 0. Check if the application was correctly registered.","");

    // Checks on nodes

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i < this->GetGeometry().size(); ++i)
    {
        if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY_COMPONENT_GRADIENT) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY_COMPONENT_GRADIENT variable on solution step data for node ",this->GetGeometry()[i].Id());
    }
    return 0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeComponentGradientSimplex<TDim, TNumNodes>::AddIntegrationPointRHSContribution(VectorType& F,
                             const array_1d<double, TNumNodes>& rShapeFunc,
                             const BoundedMatrix<double, TNumNodes, TDim>& rShapeDeriv,
                             const double Weight)
{
    double Coef = Weight;
    int LocalIndex = 0;

    for (unsigned int iNodeB = 0; iNodeB < TNumNodes; ++iNodeB){

        for (unsigned int di = 0; di < TDim; ++di){
            double value = 0.0;

            for (unsigned int iNodeA = 0; iNodeA < TNumNodes; ++iNodeA){
                double NodalComponent = 0.0;

                if (mCurrentComponent == 'X'){
                    NodalComponent = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY_X);
                }
                else if (mCurrentComponent == 'Y'){
                    NodalComponent = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY_Y);
                }
                else if (mCurrentComponent == 'Z'){
                    NodalComponent = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY_Z);
                }
                value += rShapeFunc[iNodeB] * rShapeDeriv(iNodeA, di) * NodalComponent;
            }

            F[LocalIndex++] += Coef * value;
        }
    }
}
// Explicit instantiations
template class ComputeComponentGradientSimplex<2, 3>;
template class ComputeComponentGradientSimplex<3, 4>;
} // namespace Kratos
