// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//

// System includes


// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/element_size_calculator.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "custom_elements/symbolic_pure_convection_element.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
SymbolicPureConvectionElement<TDim, TNumNodes>::SymbolicPureConvectionElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
SymbolicPureConvectionElement<TDim, TNumNodes>::SymbolicPureConvectionElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

template<std::size_t TDim, std::size_t TNumNodes>
Element::Pointer SymbolicPureConvectionElement<TDim, TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicPureConvectionElement<TDim,TNumNodes>>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template<std::size_t TDim, std::size_t TNumNodes>
Element::Pointer SymbolicPureConvectionElement<TDim, TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SymbolicPureConvectionElement<TDim,TNumNodes>>(NewId, pGeom, pProperties);
}

template<std::size_t TDim, std::size_t TNumNodes>
SymbolicPureConvectionElement<TDim, TNumNodes>::~SymbolicPureConvectionElement()
{
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
void SymbolicPureConvectionElement<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false); //false says not to preserve existing storage!!

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false); //false says not to preserve existing storage!!

    const double delta_t = rCurrentProcessInfo[DELTA_TIME];
    const double dt_inv = 1.0 / delta_t;
    const double theta = rCurrentProcessInfo.Has(TIME_INTEGRATION_THETA) ? rCurrentProcessInfo[TIME_INTEGRATION_THETA] : 0.5;

    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();
    const auto& r_convection_var = p_settings->GetConvectionVariable();
    const auto& r_supg_tau_var = p_settings->GetFirstStabilizationVariable();
    const auto& r_supg_dyn_tau_var = p_settings->GetDynamicStabilizationVariable();
    const auto& r_crosswind_alpha_var = p_settings->GetSecondStabilizationVariable();

    //getting data for the given geometry
    BoundedMatrix<double, TNumNodes, TDim > DN_DX;
    array_1d<double, TNumNodes > N;     // It is temporary (unused).
    double volume;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, volume);
    double h = ComputeH(DN_DX, volume);

    //The current and old variables
    array_1d<double,TNumNodes> phi, phi_old; // These are basically the same for the linear problems (non-iterative solution)
    array_1d< array_1d<double,3 >, TNumNodes> v, v_old;

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        phi[i] = GetGeometry()[i].FastGetSolutionStepValue(r_unknown_var);
        phi_old[i] = GetGeometry()[i].FastGetSolutionStepValue(r_unknown_var,1);

        v[i] = GetGeometry()[i].FastGetSolutionStepValue(r_convection_var);
        v_old[i] = GetGeometry()[i].FastGetSolutionStepValue(r_convection_var,1);
    }

    const array_1d<double,TDim> grad_phi = prod(trans(DN_DX), phi);
    const double norm_grad_phi = norm_2(grad_phi);

    // These coefficients would be O(1.0) if we follow the analytical estimation (based on the dimensional analysis)
    const double supg_tau_coeff = this->GetValue(r_supg_tau_var);
    const double supg_dyn_tau_coeff = this->GetValue(r_supg_dyn_tau_var);
    const double cross_wind_alpha_coeff = 0.7*this->GetValue(r_crosswind_alpha_var);

    BoundedMatrix<double,TNumNodes, TNumNodes> mass_matrix_stabilized = ZeroMatrix(TNumNodes, TNumNodes);
    BoundedMatrix<double,TNumNodes, TNumNodes> conv_matrix_stabilized = ZeroMatrix(TNumNodes, TNumNodes);
    BoundedMatrix<double,TNumNodes, TDim> tmp;

    BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;
    GetShapeFunctionsOnGauss(Ncontainer);
    for(unsigned int igauss=0; igauss<TDim+1; igauss++)
    {
        noalias(N) = row(Ncontainer,igauss);

        //obtain the velocity in the middle of the tiem step
        array_1d<double, TDim > vel_gauss=ZeroVector(TDim);
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            for(unsigned int k=0; k<TDim; k++)
            vel_gauss[k] += N[i]*(theta*v[i][k] + (1.0-theta)*v_old[i][k]);
        }
        const double norm_vel = norm_2(vel_gauss);
        array_1d<double, TNumNodes > v_dot_grad_N = prod(DN_DX, vel_gauss);

        const double tau_denominator = std::max(2.0 * norm_vel / h,  1.0e-3);
        const double tau_estimated = 1.0 / (tau_denominator);

        //terms multiplying dphi/dt
        noalias(mass_matrix_stabilized) += outer_prod(N, N);
        noalias(mass_matrix_stabilized) += supg_dyn_tau_coeff*tau_estimated*outer_prod(v_dot_grad_N, N); // We reserve the freedom to exclude the time derivative from PG group

        //terms which multiply v.grad_phi
        noalias(conv_matrix_stabilized) += outer_prod(N, v_dot_grad_N);
        noalias(conv_matrix_stabilized) += supg_tau_coeff*tau_estimated*outer_prod(v_dot_grad_N, v_dot_grad_N); // SUPG

        //cross-wind term
        if(norm_grad_phi > 1.0e-6 && norm_vel > 1e-12)
        {
            // Estimating the (current) FEM residual at the Gauss point
            const double residual = -dt_inv*(inner_prod(N,phi)-inner_prod(N,phi_old)) - inner_prod(vel_gauss, grad_phi);

            const double discontinuity_capturing_coeff = 0.5*cross_wind_alpha_coeff*h*fabs(residual/norm_grad_phi);

            // Only if we want to remove stream-wise contribution (here, if commented, we have preserved the freedom to reduce tau while increasing alpha; this would do the magic)
            /**/ 
            BoundedMatrix<double,TDim,TDim> artificial_diffusion_matrix = discontinuity_capturing_coeff*IdentityMatrix(TDim);
            const double norm_vel_squared = norm_vel*norm_vel;
            artificial_diffusion_matrix += (std::max(discontinuity_capturing_coeff - supg_tau_coeff*tau_estimated*norm_vel_squared, 0.0)
                - discontinuity_capturing_coeff)/(norm_vel_squared) * outer_prod(vel_gauss,vel_gauss);

            noalias(tmp) = prod(DN_DX,artificial_diffusion_matrix);
            noalias(conv_matrix_stabilized) += prod(tmp,trans(DN_DX)); 
            /**/

            /* noalias(conv_matrix_stabilized) += discontinuity_capturing_coeff*prod(DN_DX, trans(DN_DX)); */
        }
    }

    noalias(rLeftHandSideMatrix)  = dt_inv*mass_matrix_stabilized;
    noalias(rRightHandSideVector) = dt_inv*prod(mass_matrix_stabilized,phi_old);

    noalias(rLeftHandSideMatrix) += theta*conv_matrix_stabilized;
    noalias(rRightHandSideVector) -= (1.0 - theta)*prod(conv_matrix_stabilized,phi_old);

    // Using an incremental solution scheme, it is necessary
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, phi);

    rRightHandSideVector *= volume/static_cast<double>(TNumNodes);
    rLeftHandSideMatrix *= volume/static_cast<double>(TNumNodes);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
void SymbolicPureConvectionElement<TDim, TNumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
void SymbolicPureConvectionElement<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
void SymbolicPureConvectionElement<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    if(rResult.size() != LocalSize) {
        rResult.resize(LocalSize);
    }

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(r_unknown_var).EquationId();
    }
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
void SymbolicPureConvectionElement<TDim, TNumNodes>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];
    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    if(rElementalDofList.size() != LocalSize) {
        rElementalDofList.resize(LocalSize);
    }

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rElementalDofList[i] = GetGeometry()[i].pGetDof(r_unknown_var);
    }
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
int SymbolicPureConvectionElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS)) << "No CONVECTION_DIFFUSION_SETTINGS defined in ProcessInfo." << std::endl;
    ConvectionDiffusionSettings::Pointer p_settings = rCurrentProcessInfo[CONVECTION_DIFFUSION_SETTINGS];

    KRATOS_ERROR_IF_NOT(p_settings->IsDefinedUnknownVariable()) << "No Unknown Variable is defined in CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(p_settings->IsDefinedFirstStabilizationVariable()) << "No Stabilization Variable (Tau) is defined in CONVECTION_DIFFUSION_SETTINGS." << std::endl;
    KRATOS_ERROR_IF_NOT(p_settings->IsDefinedSecondStabilizationVariable()) << "No Stabilization Variable (Alpha) is defined in CONVECTION_DIFFUSION_SETTINGS." << std::endl;

    const auto& r_unknown_var = p_settings->GetUnknownVariable();

    const auto& r_geom = GetGeometry();
    for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); i_node++) {
        const auto& r_node = r_geom[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_unknown_var, r_node);
        KRATOS_CHECK_DOF_IN_NODE(r_unknown_var, r_node);
    }

    return Element::Check(rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
Element::IntegrationMethod SymbolicPureConvectionElement<TDim, TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
double SymbolicPureConvectionElement<TDim, TNumNodes>::ComputeH(BoundedMatrix<double,TNumNodes, TDim>& DN_DX, const double Volume)
{
    double h=0.0;
    for(unsigned int i=0; i<TNumNodes; i++)
    {
        double h_inv = 0.0;
        for(unsigned int k=0; k<TDim; k++)
        {
            h_inv += DN_DX(i,k)*DN_DX(i,k);
        }
        h = std::max(h, 1.0 / h_inv);
    }
    h = std::sqrt(h);
    return h;
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
void SymbolicPureConvectionElement<TDim, TNumNodes>::GetShapeFunctionsOnGauss(BoundedMatrix<double,4, 4>& Ncontainer)
{
    Ncontainer(0,0) = 0.58541020; Ncontainer(0,1) = 0.13819660; Ncontainer(0,2) = 0.13819660; Ncontainer(0,3) = 0.13819660;
    Ncontainer(1,0) = 0.13819660; Ncontainer(1,1) = 0.58541020; Ncontainer(1,2) = 0.13819660; Ncontainer(1,3) = 0.13819660;
    Ncontainer(2,0) = 0.13819660; Ncontainer(2,1) = 0.13819660; Ncontainer(2,2) = 0.58541020; Ncontainer(2,3) = 0.13819660;
    Ncontainer(3,0) = 0.13819660; Ncontainer(3,1) = 0.13819660; Ncontainer(3,2) = 0.13819660; Ncontainer(3,3) = 0.58541020;
}

//************************************************************************************
//************************************************************************************
template<std::size_t TDim, std::size_t TNumNodes>
void SymbolicPureConvectionElement<TDim, TNumNodes>::GetShapeFunctionsOnGauss(BoundedMatrix<double,3,3>& Ncontainer)
{
    const double one_sixt = 1.0/6.0;
    const double two_third = 2.0/3.0;
    Ncontainer(0,0) = one_sixt; Ncontainer(0,1) = one_sixt; Ncontainer(0,2) = two_third;
    Ncontainer(1,0) = one_sixt; Ncontainer(1,1) = two_third; Ncontainer(1,2) = one_sixt;
    Ncontainer(2,0) = two_third; Ncontainer(2,1) = one_sixt; Ncontainer(2,2) = one_sixt;
}

template class SymbolicPureConvectionElement<2,3>;
template class SymbolicPureConvectionElement<3,4>;

} // Namespace Kratos
