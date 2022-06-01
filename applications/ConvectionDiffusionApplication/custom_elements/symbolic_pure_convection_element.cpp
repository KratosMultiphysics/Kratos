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
        vold[i] = GetGeometry()[i].FastGetSolutionStepValue(r_convection_var,1);
    }

    array_1d<double,TDim> grad_phi_halfstep = prod(trans(DN_DX), 0.5*(phi+phi_old));
    const double norm_grad = norm_2(grad_phi_halfstep);

    const double supg_tau = this->GetValue(r_supg_tau_var);
    const double cross_wing_alpha = this->GetValue(r_crosswind_alpha_var);

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    BoundedMatrix<double,TNumNodes, TNumNodes> aux1 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying dphi/dt
    BoundedMatrix<double,TNumNodes, TNumNodes> aux2 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying phi
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
                vel_gauss[k] += 0.5*N[i]*(v[i][k]+vold[i][k]);
        }
        const double norm_vel = norm_2(vel_gauss);
        array_1d<double, TNumNodes > a_dot_grad = prod(DN_DX, vel_gauss);

        const double tau_denom = std::max(dyn_st_beta *dt_inv + 2.0 * norm_vel / h + std::abs(/*beta**/div_v),  1e-2); //the term std::abs(div_v) is added following Pablo Becker's suggestion
        const double tau = 1.0 / (tau_denom);

        //terms multiplying dphi/dt (aux1)
        noalias(aux1) += (1.0+tau*beta*div_v)*outer_prod(N, N);
        noalias(aux1) +=  tau*outer_prod(a_dot_grad, N);

        //terms which multiply the gradient of phi
        noalias(aux2) += (1.0+tau*beta*div_v)*outer_prod(N, a_dot_grad);
        noalias(aux2) += tau*outer_prod(a_dot_grad, a_dot_grad);

        //cross-wind term
        if(norm_grad > 1e-3 && norm_vel > 1e-9)
        {
            const double C = rCurrentProcessInfo.GetValue(CROSS_WIND_STABILIZATION_FACTOR);
            const double time_derivative = dt_inv*(inner_prod(N,phi)-inner_prod(N,phi_old));
            const double res = -time_derivative -inner_prod(vel_gauss, grad_phi_halfstep);

            const double disc_capturing_coeff = 0.5*C*h*fabs(res/norm_grad);
            BoundedMatrix<double,TDim,TDim> D = disc_capturing_coeff*( IdentityMatrix(TDim));
            const double norm_vel_squared = norm_vel*norm_vel;
            D += (std::max( disc_capturing_coeff - tau*norm_vel_squared , 0.0) - disc_capturing_coeff)/(norm_vel_squared) * outer_prod(vel_gauss,vel_gauss);

            noalias(tmp) = prod(DN_DX,D);
            noalias(aux2) += prod(tmp,trans(DN_DX));
        }
    }

    //adding the second and third term in the formulation
    noalias(rLeftHandSideMatrix)  = (dt_inv + theta*beta*div_v)*aux1;
    noalias(rRightHandSideVector) = (dt_inv - (1.0 - theta)*beta*div_v)*prod(aux1,phi_old);

    //terms in aux2
    noalias(rLeftHandSideMatrix) += theta*aux2;
    noalias(rRightHandSideVector) -= (1.0 - theta)*prod(aux2,phi_old);

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
