#include "swimming_DEM_application.h"
#include "calculate_gradient_Fortin_2012.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes>
void ComputeGradientFortin2012<TDim, TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  ProcessInfo& rCurrentProcessInfo)
{
    BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
    const double epsilon = 1e-6;
    const unsigned int LocalSize(TDim * TNumNodes);

    for (unsigned int i=0; i<LocalSize; ++i){
        for (unsigned int j=0; j<LocalSize; ++j){
            rLeftHandSideMatrix(i, j) *= epsilon;
        }
        rRightHandSideVector(i) *= epsilon;
    }

    AddFortin2012LHS(rLeftHandSideMatrix, rCurrentProcessInfo);

    AddFortin2012RHS(rRightHandSideVector, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeGradientFortin2012<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                        ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int LocalSize(TDim * TNumNodes);

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode){
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z_GRADIENT_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z_GRADIENT_Y);
        if (TDim == 3){
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z_GRADIENT_Z);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
int ComputeGradientFortin2012<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Perform basic element checks
    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    if(ErrorCode != 0) return ErrorCode;

    if(this->GetGeometry().size() != TDim+1)
        KRATOS_THROW_ERROR(std::invalid_argument,"wrong number of nodes for element",this->Id());

    if(VELOCITY_Z_GRADIENT.Key() == 0)

        KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY_Z_GRADIENT Key is 0. Check if the application was correctly registered.","");

    // Checks on nodes

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
    {
        if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY_Z_GRADIENT) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY_Z_GRADIENT variable on solution step data for node ",this->GetGeometry()[i].Id());
    }

    return 0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeGradientFortin2012<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                              ProcessInfo& rCurrentProcessInfo)
{

    const unsigned int LocalSize(TDim * TNumNodes);
    unsigned int LocalIndex = 0;
    unsigned int pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_Z_GRADIENT_X);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode){
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z_GRADIENT_X,pos).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z_GRADIENT_Y,pos+1).EquationId();
        if (TDim == 3){
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z_GRADIENT_Z,pos+2).EquationId();
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeGradientFortin2012<TDim, TNumNodes>::AddFortin2012LHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NEdges = 3 * TNumNodes - 6; // works in 2D and 3D
    unsigned int edges[NEdges][2];

    unsigned int i_edge = 0;
    for (unsigned int i = 0; i < TNumNodes - 1; ++i){
        for (unsigned int j = i + 1; j < TNumNodes; ++j){
            edges[i_edge][0]   = i;
            edges[i_edge++][1] = j;
        }
    }

    array_1d<array_1d<double, 3>, NEdges> EdgeVectors; // stores the [lx, ly(, lz)] vectors for all edges
    array_1d<double, NEdges> edge_lengths_inv;
    const GeometryType& rGeom = this->GetGeometry();

    for (unsigned int e = 0; e < NEdges; ++e){
        array_1d<double, 3>& le = EdgeVectors[e];
        noalias(le) = rGeom[edges[e][1]].Coordinates() - rGeom[edges[e][0]].Coordinates();
        const double he_inv = 1.0 / std::sqrt(le[0] * le[0] + le[1] * le[1] + le[2] * le[2]);
        edge_lengths_inv[e] = he_inv;
        le *= he_inv;
        AssembleEdgeLHSContribution(edges[e], le, rLeftHandSideMatrix);
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeGradientFortin2012<TDim, TNumNodes>::AssembleEdgeLHSContribution(const unsigned int edge[2], const array_1d<double, 3>& edge_normalized_vector, MatrixType& rLeftHandSideMatrix)
{
    for (unsigned int node_e = 0; node_e < 2; ++node_e){
        for (unsigned int i = 0; i < TDim; ++i){
            for (unsigned int node_f = 0; node_f < 2; ++node_f){
                for (unsigned int j = 0; j < TDim; ++j){
                    rLeftHandSideMatrix(TDim * edge[node_e] + i, TDim * edge[node_f] + j) += edge_normalized_vector[i] * edge_normalized_vector[j];
                }
            }
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeGradientFortin2012<TDim, TNumNodes>::AddFortin2012RHS(VectorType& F, ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NEdges = 3 * TNumNodes - 6; // works in 2D and 3D

    unsigned int edges[NEdges][2];

    unsigned int i_edge = 0;
    for (unsigned int i = 0; i < TNumNodes - 1; ++i){
        for (unsigned int j = i + 1; j < TNumNodes; ++j){
            edges[i_edge][0]   = i;
            edges[i_edge++][1] = j;
        }
    }

    array_1d<array_1d<double, 3>, NEdges> EdgeVectors; // stores the [lx, ly(, lz)] vectors for all edges
    array_1d<double, NEdges> edge_lengths_inv;
    const GeometryType& rGeom = this->GetGeometry();

    for (unsigned int e = 0; e < NEdges; ++e){
        array_1d<double, 3>& le = EdgeVectors[e];
        noalias(le) = rGeom[edges[e][1]].Coordinates() - rGeom[edges[e][0]].Coordinates();
        const double he_inv = 1.0 / std::sqrt(le[0] * le[0] + le[1] * le[1] + le[2] * le[2]);
        edge_lengths_inv[e] = he_inv;
        le *= he_inv;

        if (this->mCurrentComponent == 'X'){
            AssembleEdgeRHSContributionX(edges[e], he_inv, le, F);
        }

        else if (this->mCurrentComponent == 'Y'){
            AssembleEdgeRHSContributionY(edges[e], he_inv, le, F);
        }

        else {
            AssembleEdgeRHSContributionZ(edges[e], he_inv, le, F);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeGradientFortin2012<TDim, TNumNodes>::AssembleEdgeRHSContributionX(const unsigned int edge[2], const double h_edge_inv, const array_1d<double, 3>& edge_normalized_vector, VectorType& F)
{
    const double vel_component_variation_along_edge = this->GetGeometry()[edge[1]].FastGetSolutionStepValue(VELOCITY_X) - this->GetGeometry()[edge[0]].FastGetSolutionStepValue(VELOCITY_X);

    for (unsigned int node_e = 0; node_e < 2; ++node_e){
        for (unsigned int i = 0; i < TDim; ++i){
            F(TDim * edge[node_e] + i) += 2.0 * h_edge_inv * edge_normalized_vector[i] * vel_component_variation_along_edge;
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeGradientFortin2012<TDim, TNumNodes>::AssembleEdgeRHSContributionY(const unsigned int edge[2], const double h_edge_inv, const array_1d<double, 3>& edge_normalized_vector, VectorType& F)
{
    const double vel_component_variation_along_edge = this->GetGeometry()[edge[1]].FastGetSolutionStepValue(VELOCITY_Y) - this->GetGeometry()[edge[0]].FastGetSolutionStepValue(VELOCITY_Y);

    for (unsigned int node_e = 0; node_e < 2; ++node_e){
        for (unsigned int i = 0; i < TDim; ++i){
            F(TDim * edge[node_e] + i) += 2.0 * h_edge_inv * edge_normalized_vector[i] * vel_component_variation_along_edge;
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void ComputeGradientFortin2012<TDim, TNumNodes>::AssembleEdgeRHSContributionZ(const unsigned int edge[2], const double h_edge_inv, const array_1d<double, 3>& edge_normalized_vector, VectorType& F)
{
    const double vel_component_variation_along_edge = this->GetGeometry()[edge[1]].FastGetSolutionStepValue(VELOCITY_Z) - this->GetGeometry()[edge[0]].FastGetSolutionStepValue(VELOCITY_Z);

    for (unsigned int node_e = 0; node_e < 2; ++node_e){
        for (unsigned int i = 0; i < TDim; ++i){
            F(TDim * edge[node_e] + i) += 2.0 * h_edge_inv * edge_normalized_vector[i] * vel_component_variation_along_edge;
        }
    }
}

// Explicit instantiations
template class ComputeGradientFortin2012<2, 3>;
template class ComputeGradientFortin2012<3, 4>;
} // namespace Kratos
