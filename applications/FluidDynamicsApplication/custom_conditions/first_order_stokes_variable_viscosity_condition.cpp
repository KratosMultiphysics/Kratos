#include "first_order_stokes_variable_viscosity_condition.h"

namespace Kratos
{

///@name Specialized implementation of VMS for functions that depend on TDim
///@{

/**
 * @see FirstOrderStokesVariableViscosityCondition::EquationIdVector
 */
template <>
void FirstOrderStokesVariableViscosityCondition<2,2>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const unsigned int NumNodes = 2;
    const unsigned int LocalSize = 6;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}

/**
 * @see FirstOrderStokesVariableViscosityCondition::EquationIdVector
 */
template <>
void FirstOrderStokesVariableViscosityCondition<3,3>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 12;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}

/**
 * @see FirstOrderStokesVariableViscosityCondition::GetDofList
 */
template <>
void FirstOrderStokesVariableViscosityCondition<2,2>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 2;
    const SizeType LocalSize = 6;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

/**
 * @see FirstOrderStokesVariableViscosityCondition::GetDofList
 */
template <>
void FirstOrderStokesVariableViscosityCondition<3,3>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 12;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

template<>
void FirstOrderStokesVariableViscosityCondition<2,2>::ApplyNeumannCondition(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
    

    // const unsigned int TDim = 2;
    // const unsigned int TNumNodes = 2;

    // const unsigned int LocalSize = TDim+1;
    // const GeometryType& rGeom = this->GetGeometry();
    // const GeometryData::IntegrationMethod integration_method = static_cast<GeometryData::IntegrationMethod> (GeometryData::IntegrationMethod::GI_GAUSS_2);
    // const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(integration_method);
    // const unsigned int NumGauss = IntegrationPoints.size();

    // // To get the local N values with respect to the condition we would do:
    // MatrixType NConditionContainer = rGeom.ShapeFunctionsValues(integration_method);

    // // Calculate Jacobians at integration points (with respect to the condition)
    // Matrix J;
    // Matrix inv_Jt_J;
    // double det_Jt_J;
    // Vector det_Jt_J_vect(NumGauss);
    // std::vector<BoundedMatrix<double, TDim, TDim-1>> J_pseudo_inv_vect(NumGauss);
    // for (IndexType g = 0; g < NumGauss; ++g) {
    //     rGeom.Jacobian(J, g, integration_method); //J becomes of size TDim*(TDim-1)
    //     Matrix Jt_J = prod(trans(J), J);
    //     MathUtils<double>::InvertMatrix(Jt_J, inv_Jt_J, det_Jt_J);
    //     J_pseudo_inv_vect[g] = prod(inv_Jt_J,trans(J));
    //     det_Jt_J_vect[g] = det_Jt_J;
    // }

    // // These would be the DN matrices to use if using the Condition Line geometry
    // std::vector<BoundedMatrix<double, TNumNodes, TDim>> DNContainer(NumGauss);
    // const auto& DN_De_flat = rGeom.ShapeFunctionsLocalGradients(integration_method);
    // Matrix DN_De;
    // for (IndexType g = 0; g < NumGauss; ++g) {
    //     DN_De = ZeroMatrix(TDim,TDim);
    //     DN_De(0,0) = DN_De_flat[g](0,0);
    //     DN_De(0,1) = DN_De_flat[g](1,0);
    //     DNContainer[g] = trans(prod(trans(J_pseudo_inv_vect[g]), DN_De));
    // }

    // for (unsigned int g = 0; g < NumGauss; g++)
    // {
    //     // This would be the case for computing in reference to the condition
    //     Vector N = row(NConditionContainer,g);
    //     BoundedMatrix<double, TNumNodes, TDim> DN = DNContainer[g];
    //     double w_g = sqrt(det_Jt_J_vect[g]) * IntegrationPoints[g].Weight();

    //     // Compute Unit normal
    //     array_1d<double, 3> n_gauss;
    //     n_gauss = rGeom.Normal(IntegrationPoints[g].Coordinates());
    //     double A = norm_2(n_gauss);
    //     n_gauss /= A;

    //     // Set nodal data

    //     // If we were computing with the nodes within the condition, we would do like this:
    //     array_1d<double, TNumNodes> nu_nodes;
    //     BoundedMatrix<double, TNumNodes, TDim> u_nodes;
    //     for (IndexType i = 0; i < TNumNodes; ++i) {
    //         nu_nodes[i] = rGeom[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
    //         const auto& r_v = rGeom[i].FastGetSolutionStepValue(VELOCITY);
    //         for (IndexType d = 0; d < TDim; ++d) {
    //             u_nodes(i, d) = r_v[d];
    //         }
    //     }

    //     const double delta_gauss = rCurrentProcessInfo.GetValue(TAUONE);
        
    // }
}

template<>
void FirstOrderStokesVariableViscosityCondition<3,3>::ApplyNeumannCondition(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
    // const unsigned int TDim = 3;
    // const unsigned int TNumNodes = 3;

    // const unsigned int LocalSize = TDim+1;
    // const GeometryType& rGeom = this->GetGeometry();
    // const GeometryData::IntegrationMethod integration_method = static_cast<GeometryData::IntegrationMethod> (GeometryData::IntegrationMethod::GI_GAUSS_2);
    // const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(integration_method);
    // const unsigned int NumGauss = IntegrationPoints.size();

    // // Calculate Jacobians at integration points
    // Matrix J;
    // Matrix inv_Jt_J;
    // double det_Jt_J;
    // Vector det_Jt_J_vect(NumGauss);
    // std::vector<BoundedMatrix<double, TDim, TDim-1>> J_pseudo_inv_vect(NumGauss);
    // for (IndexType g = 0; g < NumGauss; ++g) {
    //     rGeom.Jacobian(J, g, integration_method); //J becomes of size TDim*(TDim-1)
    //     Matrix Jt_J = prod(trans(J), J);
    //     MathUtils<double>::InvertMatrix(Jt_J, inv_Jt_J, det_Jt_J);
    //     J_pseudo_inv_vect[g] = prod(inv_Jt_J,trans(J));
    //     det_Jt_J_vect[g] = det_Jt_J;
    // }

    // MatrixType NContainer = rGeom.ShapeFunctionsValues(integration_method);

    // std::vector<BoundedMatrix<double, TNumNodes, TDim>> DNContainer(NumGauss);
    // const auto& DN_De_flat = rGeom.ShapeFunctionsLocalGradients(integration_method);
    // Matrix DN_De;
    // for (IndexType g = 0; g < NumGauss; ++g) {
    //     DN_De = ZeroMatrix(TDim,TDim);
    //     DN_De(0,0) = DN_De_flat[g](0,0);
    //     DN_De(0,1) = DN_De_flat[g](1,0);
    //     DN_De(0,2) = DN_De_flat[g](2,0);
    //     DN_De(1,0) = DN_De_flat[g](0,1);
    //     DN_De(1,1) = DN_De_flat[g](1,1);
    //     DN_De(1,2) = DN_De_flat[g](2,1);
    //     DNContainer[g] = trans(prod(trans(J_pseudo_inv_vect[g]), DN_De));
    // }

    // for (unsigned int g = 0; g < NumGauss; g++)
    // {
    //     Vector N = row(NContainer,g);
    //     BoundedMatrix<double, TNumNodes, TDim> DN = DNContainer[g];
    //     double w_g = sqrt(det_Jt_J_vect[g]) * IntegrationPoints[g].Weight();

    //     // Compute Unit normal
    //     array_1d<double, 3> n_gauss;
    //     n_gauss = rGeom.Normal(IntegrationPoints[g].Coordinates());
    //     double A = norm_2(n_gauss);
    //     n_gauss /= A;

    //     // Set nodal data
    //     array_1d<double, TNumNodes> nu_nodes;
    //     BoundedMatrix<double, TNumNodes, TDim> u_nodes;
    //     for (IndexType i = 0; i < TNumNodes; ++i) {
    //         nu_nodes[i] = rGeom[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
    //         const auto& r_v = rGeom[i].FastGetSolutionStepValue(VELOCITY);
    //         for (IndexType d = 0; d < TDim; ++d) {
    //             u_nodes(i, d) = r_v[d];
    //         }
    //     }

    //     const double delta_gauss = rCurrentProcessInfo.GetValue(TAUONE);
        
    // }
}

template class FirstOrderStokesVariableViscosityCondition<2,2>;
template class FirstOrderStokesVariableViscosityCondition<3,3>;

} // namespace Kratos