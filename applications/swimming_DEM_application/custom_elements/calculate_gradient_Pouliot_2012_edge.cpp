#include "swimming_DEM_application.h"
#include "calculate_gradient_Pouliot_2012_edge.h"

namespace Kratos
{

/// Calculate the element's local contribution to the system for the current step.
template <unsigned int TDim, unsigned int TNumNodes>
void ComputeGradientPouliot2012Edge<TDim, TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int LocalSize(TDim * TNumNodes);

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

    AddPouliot2012LHS(rLeftHandSideMatrix, rCurrentProcessInfo);

    AddPouliot2012RHS(rRightHandSideVector, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeGradientPouliot2012Edge<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                              ProcessInfo& rCurrentProcessInfo)
{

    const unsigned int LocalSize(TDim * TNumNodes);
    unsigned int LocalIndex = 0;
    unsigned int pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_COMPONENT_GRADIENT_X);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode){
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_COMPONENT_GRADIENT_X,pos).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_COMPONENT_GRADIENT_Y,pos+1).EquationId();
        if (TDim == 3){
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_COMPONENT_GRADIENT_Z,pos+2).EquationId();
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeGradientPouliot2012Edge<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                        ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int LocalSize(TDim * TNumNodes);

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode){
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_COMPONENT_GRADIENT_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_COMPONENT_GRADIENT_Y);
        if (TDim == 3){
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_COMPONENT_GRADIENT_Z);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeGradientPouliot2012Edge<TDim, TNumNodes>::AddPouliot2012LHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const GeometryType& rGeom = this->GetGeometry();
    const array_1d<double, 3> le = rGeom[1].Coordinates() - rGeom[0].Coordinates(); // vector from node 0 to node 1
    const double h_edge = SWIMMING_MODULUS_3(le);
    const double h_edge_inv_2 = 1.0 / SWIMMING_INNER_PRODUCT_3(le, le);

    const double epsilon = 1e-6 * h_edge;
    for (unsigned int node_e = 0; node_e < TNumNodes; ++node_e){
        for (unsigned int i = 0; i < TDim; ++i){
            for (unsigned int node_f = 0; node_f < TNumNodes; ++node_f){
                for (unsigned int j = 0; j < TDim; ++j){
                    double stab = 0.0;
                    if (i == j){
                        stab = node_e == node_f ? epsilon : - epsilon;
                    }
                    rLeftHandSideMatrix(TDim * node_e + i, TDim * node_f + j) = h_edge_inv_2 * le[i] * le[j] + stab;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeGradientPouliot2012Edge<TDim, TNumNodes>::AddPouliot2012RHS(VectorType& F, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const array_1d<double, 3> le = rGeom[1].Coordinates() - rGeom[0].Coordinates(); // vector from node 0 to node 1
    const double h_edge_inv_2 = 1.0 / SWIMMING_INNER_PRODUCT_3(le, le);

    double vel_component_variation_along_edge;

    if (this->mCurrentComponent == 'X'){
        vel_component_variation_along_edge = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_X) - this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X);
    }

    else if (this->mCurrentComponent == 'Y'){
        vel_component_variation_along_edge = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_Y) - this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y);
    }

    else {
        vel_component_variation_along_edge = this->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_Z) - this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Z);
    }

    for (unsigned int node_e = 0; node_e < TNumNodes; ++node_e){
        for (unsigned int i = 0; i < TDim; ++i){
            F(TDim * node_e + i) = 2.0 * h_edge_inv_2 * le[i] * vel_component_variation_along_edge;
        }
    }
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
int ComputeGradientPouliot2012Edge<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Perform basic element checks
    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    if(ErrorCode != 0) return ErrorCode;

    if(this->GetGeometry().size() != TNumNodes)
        KRATOS_THROW_ERROR(std::invalid_argument, "wrong number of nodes for element",this->Id());

    if(VELOCITY_COMPONENT_GRADIENT.Key() == 0)

        KRATOS_THROW_ERROR(std::invalid_argument, "VELOCITY_COMPONENT_GRADIENT Key is 0. Check if the application was correctly registered.","");

    // Checks on nodes

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
    {
        Node<3> &rNode = this->GetGeometry()[i];
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_COMPONENT_GRADIENT_X,rNode);
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_COMPONENT_GRADIENT_Y,rNode);
        if (TDim == 3){
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_COMPONENT_GRADIENT_Z,rNode);
        }
        if(rNode.SolutionStepsDataHas(VELOCITY_COMPONENT_GRADIENT) == false)
            KRATOS_THROW_ERROR(std::invalid_argument, "missing VELOCITY_COMPONENT_GRADIENT variable on solution step data for node ",this->GetGeometry()[i].Id());
    }

    return 0;

    KRATOS_CATCH("");
}

// Explicit instantiations
template class ComputeGradientPouliot2012Edge<2, 2>;
template class ComputeGradientPouliot2012Edge<3, 2>;
} // namespace Kratos
