#include "swimming_DEM_application.h"
#include "calculate_gradient_Fortin_2012.h"

namespace Kratos
{
template <>
void ComputeGradientFortin2012<2, 3>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  ProcessInfo& rCurrentProcessInfo)
{
    BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
    const double epsilon = 1e-6;
    const unsigned int NumNodes(3), LocalSize(2 * NumNodes);

    for (unsigned int i=0; i<LocalSize; ++i){
        for (unsigned int j=0; j<LocalSize; ++j){
            rLeftHandSideMatrix(i, j) *= epsilon;
        }
        rRightHandSideVector(i) *= epsilon;
    }

    AddFortin2012LHS(rLeftHandSideMatrix, rCurrentProcessInfo);

    AddFortin2012RHS(rRightHandSideVector, rCurrentProcessInfo);
}

template <>
void ComputeGradientFortin2012<3, 4>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  ProcessInfo& rCurrentProcessInfo)
{
    BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
    const double epsilon = 1e-6;
    const unsigned int NumNodes(4), LocalSize(3 * NumNodes);

    for (unsigned int i=0; i<LocalSize; ++i){
        for (unsigned int j=0; j<LocalSize; ++j){
            rLeftHandSideMatrix(i, j) *= epsilon;
        }
        rRightHandSideVector(i) *= epsilon;
    }

    AddFortin2012LHS(rLeftHandSideMatrix, rCurrentProcessInfo);

    AddFortin2012RHS(rRightHandSideVector, rCurrentProcessInfo);
}


} // namespace Kratos
