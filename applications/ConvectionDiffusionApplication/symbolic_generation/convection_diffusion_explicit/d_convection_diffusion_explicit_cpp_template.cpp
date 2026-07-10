//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Riccardo Tosi
//

#include "d_convection_diffusion_explicit.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
DConvectionDiffusionExplicit<TDim,TNumNodes>::DConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : QSConvectionDiffusionExplicit<TDim,TNumNodes>(NewId, pGeometry) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
DConvectionDiffusionExplicit<TDim,TNumNodes>::DConvectionDiffusionExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : QSConvectionDiffusionExplicit<TDim,TNumNodes>(NewId, pGeometry, pProperties) {}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
DConvectionDiffusionExplicit<TDim,TNumNodes>::~DConvectionDiffusionExplicit() {}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer DConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<DConvectionDiffusionExplicit>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
Element::Pointer DConvectionDiffusionExplicit<TDim,TNumNodes>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<DConvectionDiffusionExplicit>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void DConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the CalculateLocalSystem() method for the explicit Convection-Diffusion element.";
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void DConvectionDiffusionExplicit<TDim,TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    KRATOS_ERROR << "Calling the CalculateRightHandSide() method for the explicit Convection-Diffusion element. Call the DCalculateRightHandSideInternal() method instead.";
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void DConvectionDiffusionExplicit<TDim,TNumNodes>::AddExplicitContribution(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    auto& r_geometry = this->GetGeometry();
    const unsigned int local_size = r_geometry.size();
    BoundedVector<double, TNumNodes> rhs;
    this->DCalculateRightHandSideInternal(rhs,rCurrentProcessInfo);
    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    const auto& reaction_variable = r_process_info[CONVECTION_DIFFUSION_SETTINGS]->GetReactionVariable();
    for (unsigned int i_node = 0; i_node < local_size; i_node++) {
        #pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(reaction_variable) += rhs[i_node];
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void DConvectionDiffusionExplicit<TDim,TNumNodes>::Initialize(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::Initialize(rCurrentProcessInfo);
    mUnknownSubScale = ZeroVector(TNumNodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void DConvectionDiffusionExplicit<TDim,TNumNodes>::FinalizeSolutionStep(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);
    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);
    // Reading integration points and local gradients
    const auto& r_geometry = this->GetGeometry();
    const auto& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    // Iterate over integration points to update subscales on the gauss point
    for (unsigned int g = 0; g < integration_points.size(); g++) {
        // Calculate N on the gauss point "g"
        rData.N = row(rData.N_gausspoint,g);
        // Compute tau
        this->DCalculateTau(rData);
        // Retrieve unknown belonging to subgrid scale space on gauss integration point g
        rData.unknown_subscale = mUnknownSubScale(g);
        // Update unknown belonging to subgrid scale space on gauss integration point g
        this->UpdateUnknownSubgridScaleGaussPoint(rData,g);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void DConvectionDiffusionExplicit<TDim,TNumNodes>::Calculate(
    const Variable<double>& rVariable,
    double& Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const ProcessInfo& r_process_info = rCurrentProcessInfo;
    ConvectionDiffusionSettings::Pointer p_settings = r_process_info[CONVECTION_DIFFUSION_SETTINGS];
    if (rVariable == p_settings->GetProjectionVariable()) {
        auto& r_geometry = this->GetGeometry();
        const unsigned int local_size = r_geometry.size();
        BoundedVector<double, TNumNodes> rhs_oss;
        this->DCalculateOrthogonalSubgridScaleRHSInternal(rhs_oss,rCurrentProcessInfo);
        for (unsigned int i_node = 0; i_node < local_size; i_node++) {
            #pragma omp atomic
            r_geometry[i_node].GetValue(rVariable) += rhs_oss[i_node];
        }
    }
    else {
        BaseType::Calculate(rVariable,Output,rCurrentProcessInfo);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void DConvectionDiffusionExplicit<2,3>::DCalculateRightHandSideInternal(
    BoundedVector<double, 3>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->DCalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& explicit_step_coefficient = rData.explicit_step_coefficient;
    const auto& v = rData.convective_velocity;
    const auto& tau = rData.tau;
    const auto& prj = rData.oss_projection;
    const auto& phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix accesses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rData.DN_DX(0, 0);
    const double& DN_DX_0_1 = rData.DN_DX(0, 1);
    const double& DN_DX_1_0 = rData.DN_DX(1, 0);
    const double& DN_DX_1_1 = rData.DN_DX(1, 1);
    const double& DN_DX_2_0 = rData.DN_DX(2, 0);
    const double& DN_DX_2_1 = rData.DN_DX(2, 1);
    // RHS
    auto& rhs = rData.rhs;

    //substitute_rhs_2D

    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 3;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void DConvectionDiffusionExplicit<3,4>::DCalculateRightHandSideInternal(
    BoundedVector<double, 4>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->DCalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& explicit_step_coefficient = rData.explicit_step_coefficient;
    const auto& v = rData.convective_velocity;
    const auto& tau = rData.tau;
    const auto& prj = rData.oss_projection;
    const auto& phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix accesses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rData.DN_DX(0,0);
    const double& DN_DX_0_1 = rData.DN_DX(0,1);
    const double& DN_DX_0_2 = rData.DN_DX(0,2);
    const double& DN_DX_1_0 = rData.DN_DX(1,0);
    const double& DN_DX_1_1 = rData.DN_DX(1,1);
    const double& DN_DX_1_2 = rData.DN_DX(1,2);
    const double& DN_DX_2_0 = rData.DN_DX(2,0);
    const double& DN_DX_2_1 = rData.DN_DX(2,1);
    const double& DN_DX_2_2 = rData.DN_DX(2,2);
    const double& DN_DX_3_0 = rData.DN_DX(3,0);
    const double& DN_DX_3_1 = rData.DN_DX(3,1);
    const double& DN_DX_3_2 = rData.DN_DX(3,2);
    // RHS
    auto& rhs = rData.rhs;

    //substitute_rhs_3D

    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 4;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void DConvectionDiffusionExplicit<2,3>::DCalculateOrthogonalSubgridScaleRHSInternal(
    BoundedVector<double, 3>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->DCalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& explicit_step_coefficient = rData.explicit_step_coefficient;
    const auto& v = rData.convective_velocity;
    const auto& phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix accesses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rData.DN_DX(0, 0);
    const double& DN_DX_0_1 = rData.DN_DX(0, 1);
    const double& DN_DX_1_0 = rData.DN_DX(1, 0);
    const double& DN_DX_1_1 = rData.DN_DX(1, 1);
    const double& DN_DX_2_0 = rData.DN_DX(2, 0);
    const double& DN_DX_2_1 = rData.DN_DX(2, 1);
    // RHS
    auto& rhs = rData.rhs;

    //substitute_oss_2D

    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 3;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void DConvectionDiffusionExplicit<3,4>::DCalculateOrthogonalSubgridScaleRHSInternal(
    BoundedVector<double, 4>& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Element data
    ElementData rData;
    this->InitializeEulerianElement(rData,rCurrentProcessInfo);

    // Compute tau
    this->DCalculateTau(rData);

    // Retrieve element data
    const auto& alpha = rData.diffusivity;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& explicit_step_coefficient = rData.explicit_step_coefficient;
    const auto& v = rData.convective_velocity;
    const auto& phi_subscale_gauss = mUnknownSubScale;
    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix accesses
    // The notation DN_i_j means shape function for node i in dimension j
    const double& DN_DX_0_0 = rData.DN_DX(0,0);
    const double& DN_DX_0_1 = rData.DN_DX(0,1);
    const double& DN_DX_0_2 = rData.DN_DX(0,2);
    const double& DN_DX_1_0 = rData.DN_DX(1,0);
    const double& DN_DX_1_1 = rData.DN_DX(1,1);
    const double& DN_DX_1_2 = rData.DN_DX(1,2);
    const double& DN_DX_2_0 = rData.DN_DX(2,0);
    const double& DN_DX_2_1 = rData.DN_DX(2,1);
    const double& DN_DX_2_2 = rData.DN_DX(2,2);
    const double& DN_DX_3_0 = rData.DN_DX(3,0);
    const double& DN_DX_3_1 = rData.DN_DX(3,1);
    const double& DN_DX_3_2 = rData.DN_DX(3,2);
    // RHS
    auto& rhs = rData.rhs;

    //substitute_oss_3D

    // All the weights of the gauss points are the same so we multiply by volume/n_nodes
    const double local_size = 4;
    noalias(rRightHandSideVector) = rhs * rData.volume/local_size;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void DConvectionDiffusionExplicit<2,3>::UpdateUnknownSubgridScaleGaussPoint(
    ElementData& rData,
    unsigned int g)
{
    KRATOS_TRY;

    // Retrieve element data
    const auto& N = rData.N;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& v = rData.convective_velocity;
    const auto& tau = rData.tau[g];
    const auto& phi_subscale_gauss = rData.unknown_subscale;
    const auto& prj = rData.oss_projection;
    double phi_subscale_gauss_new = 0;
    const double& DN_DX_0_0 = rData.DN_DX(0, 0);
    const double& DN_DX_0_1 = rData.DN_DX(0, 1);
    const double& DN_DX_1_0 = rData.DN_DX(1, 0);
    const double& DN_DX_1_1 = rData.DN_DX(1, 1);
    const double& DN_DX_2_0 = rData.DN_DX(2, 0);
    const double& DN_DX_2_1 = rData.DN_DX(2, 1);

    phi_subscale_gauss_new += N[0]*f[0] + N[1]*f[1] + N[2]*f[2]; // forcing term
    phi_subscale_gauss_new += - (N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]))/(delta_time); // mass term
    phi_subscale_gauss_new += - (DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) - (DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1)); // convective term 1
    phi_subscale_gauss_new += - (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2])*(DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1)); // convective term 2
    phi_subscale_gauss_new += N[0]*prj[0] + N[1]*prj[1] + N[2]*prj[2]; // OSS term
    phi_subscale_gauss_new *= tau;

    mUnknownSubScale(g) = (tau*phi_subscale_gauss/delta_time) + phi_subscale_gauss_new;

    KRATOS_CATCH("");
}

/***********************************************************************************/

template <>
void DConvectionDiffusionExplicit<3,4>::UpdateUnknownSubgridScaleGaussPoint(
    ElementData& rData,
    unsigned int g)
{
    KRATOS_TRY;

    // Retrieve element data
    const auto& N = rData.N;
    const auto& f = rData.forcing;
    const auto& phi = rData.unknown;
    const auto& phi_old = rData.unknown_old;
    const auto& delta_time = rData.delta_time;
    const auto& v = rData.convective_velocity;
    const auto& tau = rData.tau[g];
    const auto& phi_subscale_gauss = rData.unknown_subscale;
    const auto& prj = rData.oss_projection;
    double phi_subscale_gauss_new = 0;
    const double& DN_DX_0_0 = rData.DN_DX(0,0);
    const double& DN_DX_0_1 = rData.DN_DX(0,1);
    const double& DN_DX_0_2 = rData.DN_DX(0,2);
    const double& DN_DX_1_0 = rData.DN_DX(1,0);
    const double& DN_DX_1_1 = rData.DN_DX(1,1);
    const double& DN_DX_1_2 = rData.DN_DX(1,2);
    const double& DN_DX_2_0 = rData.DN_DX(2,0);
    const double& DN_DX_2_1 = rData.DN_DX(2,1);
    const double& DN_DX_2_2 = rData.DN_DX(2,2);
    const double& DN_DX_3_0 = rData.DN_DX(3,0);
    const double& DN_DX_3_1 = rData.DN_DX(3,1);
    const double& DN_DX_3_2 = rData.DN_DX(3,2);

    phi_subscale_gauss_new += N[0]*f[0] + N[1]*f[1] + N[2]*f[2] + N[3]*f[3]; // forcing term
    phi_subscale_gauss_new += - (N[0]*(phi[0] - phi_old[0]) + N[1]*(phi[1] - phi_old[1]) + N[2]*(phi[2] - phi_old[2]) + N[3]*(phi[3] - phi_old[3]))/(delta_time); // mass term
    phi_subscale_gauss_new += - (DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2] + DN_DX_3_0*phi[3])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - (DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2] + DN_DX_3_1*phi[3])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - (DN_DX_0_2*phi[0] + DN_DX_1_2*phi[1] + DN_DX_2_2*phi[2] + DN_DX_3_2*phi[3])*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)); // convective term 1
    phi_subscale_gauss_new += - (DN_DX_0_0*phi[0] + DN_DX_1_0*phi[1] + DN_DX_2_0*phi[2] + DN_DX_3_0*phi[3])*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) - (DN_DX_0_1*phi[0] + DN_DX_1_1*phi[1] + DN_DX_2_1*phi[2] + DN_DX_3_1*phi[3])*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) - (DN_DX_0_2*phi[0] + DN_DX_1_2*phi[1] + DN_DX_2_2*phi[2] + DN_DX_3_2*phi[3])*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)) - (N[0]*phi[0] + N[1]*phi[1] + N[2]*phi[2] + N[3]*phi[3])*(DN_DX_0_0*v(0,0) + DN_DX_0_1*v(0,1) + DN_DX_0_2*v(0,2) + DN_DX_1_0*v(1,0) + DN_DX_1_1*v(1,1) + DN_DX_1_2*v(1,2) + DN_DX_2_0*v(2,0) + DN_DX_2_1*v(2,1) + DN_DX_2_2*v(2,2) + DN_DX_3_0*v(3,0) + DN_DX_3_1*v(3,1) + DN_DX_3_2*v(3,2)); // convective term 2
    phi_subscale_gauss_new += N[0]*prj[0] + N[1]*prj[1] + N[2]*prj[2] + N[3]*prj[3]; // OSS term
    phi_subscale_gauss_new *= tau;

    mUnknownSubScale(g) = (tau*phi_subscale_gauss/delta_time) + phi_subscale_gauss_new;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< unsigned int TDim, unsigned int TNumNodes >
void DConvectionDiffusionExplicit<TDim,TNumNodes>::DCalculateTau(
    ElementData& rData)
{
    KRATOS_TRY;

    // Calculate h
    double h = this->ComputeH(rData.DN_DX);
    // Calculate tau for each gauss point
    for(unsigned int g = 0; g<TNumNodes; g++) {
	const auto& N = row(rData.N_gausspoint,g);
        // Calculate velocity and velocity divergence in the gauss point
        const array_1d<double,3> vel_gauss = prod(N, rData.convective_velocity);
        double div_vel = 0;
        for(unsigned int node_element = 0; node_element<TNumNodes; node_element++) {
            for(unsigned int dim = 0; dim < TDim; dim++) {
                div_vel += rData.DN_DX(node_element,dim)*rData.convective_velocity(node_element,dim);
            }
        }
        const double norm_velocity = norm_2(vel_gauss);
        // Estimate tau
        // Dynamic part
        double inv_tau = 1.0/rData.delta_time;
        // Convection
        inv_tau += 2.0 * norm_velocity / h;
        inv_tau += 1.0*div_vel; // unitary coefficient in front of \nabla \cdot convective_velocity term in the strong equation
        // Diffusion
        inv_tau += 4.0 * rData.diffusivity / (h*h);
        // Limiting
        inv_tau = std::max(inv_tau, 1e-2);
        rData.tau[g] = (1.0) / inv_tau;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template class DConvectionDiffusionExplicit<2,3>;
template class DConvectionDiffusionExplicit<3,4>;

/***********************************************************************************/
/***********************************************************************************/

}
