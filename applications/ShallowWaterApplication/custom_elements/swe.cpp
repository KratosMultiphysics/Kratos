//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "utilities/geometry_utilities.h"
#include "shallow_water_application_variables.h"
#include "swe.h"

namespace Kratos
{

template< size_t TNumNodes, ElementFramework TFramework >
int SWE<TNumNodes, TFramework>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    // Base class checks for positive Jacobian and Id > 0
    int err = Element::Check(rCurrentProcessInfo);
    if(err != 0) return err;

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( size_t i = 0; i < TNumNodes; i++ )
    {
        const Node<3>& node = this->GetGeometry()[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MOMENTUM, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(FREE_SURFACE_ELEVATION, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TOPOGRAPHY, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MANNING, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RAIN, node)

        KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_X, node)
        KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_Y, node)
        KRATOS_CHECK_DOF_IN_NODE(FREE_SURFACE_ELEVATION, node)
    }

    return err;
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    const size_t element_size = TNumNodes*3;
    if(rResult.size() != element_size)
        rResult.resize(element_size,false);                         // False says not to preserve existing storage!!

    const GeometryType& rGeom = GetGeometry();

    int counter=0;
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rResult[counter++] = rGeom[i].GetDof(MOMENTUM_X).EquationId();
        rResult[counter++] = rGeom[i].GetDof(MOMENTUM_Y).EquationId();
        rResult[counter++] = rGeom[i].GetDof(FREE_SURFACE_ELEVATION).EquationId();
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    const size_t element_size = TNumNodes*3;
    if(rElementalDofList.size() != element_size)
        rElementalDofList.resize(element_size);

    const GeometryType& rGeom = GetGeometry();

    int counter=0;
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rElementalDofList[counter++] = rGeom[i].pGetDof(MOMENTUM_X);
        rElementalDofList[counter++] = rGeom[i].pGetDof(MOMENTUM_Y);
        rElementalDofList[counter++] = rGeom[i].pGetDof(FREE_SURFACE_ELEVATION);
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Resize of the Left and Right Hand side
    constexpr size_t element_size = TNumNodes*3;
    if(rLeftHandSideMatrix.size1() != element_size)
        rLeftHandSideMatrix.resize(element_size,element_size,false); // False says not to preserve existing storage!!

    if(rRightHandSideVector.size() != element_size)
        rRightHandSideVector.resize(element_size,false);             // False says not to preserve existing storage!!

    const GeometryType& Geom = this->GetGeometry();

    ElementVariables variables;
    this->InitializeElementVariables(variables, rCurrentProcessInfo);

    const BoundedMatrix<double,TNumNodes, TNumNodes> NContainer = Geom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 ); // In this case, number of Gauss points and number of nodes coincides

    BoundedMatrix<double,TNumNodes, 2> DN_DX;  // Shape function gradients are constant since we are using linear functions
    array_1d<double,TNumNodes> N;
    double Area;
    this->CalculateGeometry(DN_DX,Area);

    this->GetNodalValues(variables);
    this->CalculateElementValues(DN_DX, variables);

    rLeftHandSideMatrix = ZeroMatrix(element_size, element_size);
    rRightHandSideVector = ZeroVector(element_size);

    for (size_t i_gauss = 0; i_gauss < TNumNodes; ++i_gauss)
    {
        N = row(NContainer, i_gauss);

        this->BuildAuxiliaryMatrices(N, DN_DX, variables);

        this->AddInertiaTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        this->AddConvectiveTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        this->AddWaveTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        this->AddFrictionTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        this->AddStabilizationTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        this->AddSourceTerms(rRightHandSideVector, variables);
    }

    // Substracting the Dirichlet term (since we use a residualbased approach)
    rRightHandSideVector -= prod(rLeftHandSideMatrix, variables.unknown);

    rRightHandSideVector *= Area * variables.lumping_factor;
    rLeftHandSideMatrix  *= Area * variables.lumping_factor;

    double residual = norm_1(rRightHandSideVector);
    this->SetValue(RESIDUAL_NORM, residual);
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == VEL_ART_VISC || rVariable == PR_ART_VISC || rVariable == RESIDUAL_NORM)
    {
        if (rValues.size() != TNumNodes)
            rValues.resize(TNumNodes);

        for (size_t PointNumber = 0; PointNumber < TNumNodes; PointNumber++)
            rValues[PointNumber] = this->GetValue(rVariable);
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::InitializeElementVariables(
    ElementVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    const double delta_t = rCurrentProcessInfo[DELTA_TIME];
    rVariables.epsilon = rCurrentProcessInfo[DRY_HEIGHT];
    rVariables.dt_inv = 1.0 / delta_t;
    rVariables.lumping_factor = 1.0 / static_cast<double>(TNumNodes);
    rVariables.stab_factor = rCurrentProcessInfo[STABILIZATION_FACTOR];
    rVariables.gravity = rCurrentProcessInfo[GRAVITY_Z];
    rVariables.manning2 = 0.0;
    rVariables.porosity = 0.0;
    rVariables.permeability = rCurrentProcessInfo[PERMEABILITY];
    rVariables.discharge_penalty = rCurrentProcessInfo[DRY_DISCHARGE_PENALTY];

    const GeometryType& rGeom = GetGeometry();
    for (size_t i = 0; i < TNumNodes; i++)
    {
        const double f = rGeom[i].FastGetSolutionStepValue(FREE_SURFACE_ELEVATION);
        const double z = rGeom[i].FastGetSolutionStepValue(TOPOGRAPHY);
        const double n = rGeom[i].FastGetSolutionStepValue(MANNING);
        const double h = f - z;
        if (h > rVariables.epsilon) {
            rVariables.manning2 += n;
            rVariables.porosity += 1.0;
        } else {
            const double beta = 1e4;
            rVariables.manning2 += n * (1 - beta * (h - rVariables.epsilon));
            rVariables.porosity += 0.0;
        }
    }
    rVariables.manning2 *= rVariables.lumping_factor;
    rVariables.manning2 = std::pow(rVariables.manning2, 2);

    rVariables.porosity *= rVariables.lumping_factor;
    if (rVariables.porosity < 1.0) {
        rVariables.porosity = 0.0;
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::CalculateGeometry(BoundedMatrix<double, TNumNodes, 2>& rDN_DX, double& rArea)
{
    const GeometryType& rGeom = this->GetGeometry();

    // We select GI_GAUSS_1 due to we are computing at the barycenter.
    const GeometryType::IntegrationPointsArrayType& integration_points = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_1);
    const size_t NumGPoints = integration_points.size();
    rArea = rGeom.Area();
    GeometryType::ShapeFunctionsGradientsType DN_DXContainer( NumGPoints );
    rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DXContainer, GeometryData::GI_GAUSS_1);

    rDN_DX = DN_DXContainer[0];
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::GetNodalValues(ElementVariables& rVariables)
{
    const GeometryType& rGeom = GetGeometry();
    size_t counter = 0;
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rVariables.rain[counter]  = 0;
        rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(MOMENTUM_X);
        rVariables.prev_unk[counter] = rGeom[i].FastGetSolutionStepValue(MOMENTUM_X, 1);
        counter++;

        rVariables.rain[counter]  = 0;
        rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(MOMENTUM_Y);
        rVariables.prev_unk[counter] = rGeom[i].FastGetSolutionStepValue(MOMENTUM_Y, 1);
        counter++;

        rVariables.rain[counter]  = rGeom[i].FastGetSolutionStepValue(RAIN);
        rVariables.unknown[counter]  = rGeom[i].FastGetSolutionStepValue(FREE_SURFACE_ELEVATION);
        rVariables.prev_unk[counter] = rGeom[i].FastGetSolutionStepValue(FREE_SURFACE_ELEVATION, 1);
        counter++;
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::CalculateElementValues(
    const BoundedMatrix<double,TNumNodes, 2>& rDN_DX,
    ElementVariables& rVariables)
{
    // Initialize outputs
    rVariables.projected_momentum = ZeroVector(2);
    rVariables.height = 0.0;
    rVariables.surface_grad = ZeroVector(2);
    rVariables.velocity = ZeroVector(2);
    rVariables.momentum_div = 0.0;
    rVariables.velocity_div = 0.0;

    const GeometryType& rGeom = GetGeometry();

    // integrate over the element
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rVariables.velocity += rGeom[i].FastGetSolutionStepValue(VELOCITY);
        const double z = rGeom[i].FastGetSolutionStepValue(TOPOGRAPHY);
        const double f = rGeom[i].FastGetSolutionStepValue(FREE_SURFACE_ELEVATION);
        rVariables.height += f - z;
        rVariables.surface_grad[0] += rDN_DX(i,0) * f;
        rVariables.surface_grad[1] += rDN_DX(i,1) * f;
        rVariables.momentum_div += rDN_DX(i,0) * rGeom[i].FastGetSolutionStepValue(MOMENTUM_X);
        rVariables.momentum_div += rDN_DX(i,1) * rGeom[i].FastGetSolutionStepValue(MOMENTUM_Y);
        rVariables.velocity_div += rDN_DX(i,0) * rGeom[i].FastGetSolutionStepValue(VELOCITY_X);
        rVariables.velocity_div += rDN_DX(i,1) * rGeom[i].FastGetSolutionStepValue(VELOCITY_Y);
        rVariables.projected_momentum += rGeom[i].FastGetSolutionStepValue(MOMENTUM, 1);
    }

    rVariables.velocity *= rVariables.lumping_factor;
    rVariables.height *= rVariables.lumping_factor;
    rVariables.height = std::max(rVariables.height, 0.0);
    rVariables.projected_momentum *= rVariables.lumping_factor;

    rVariables.wave_vel_2 = rVariables.gravity * rVariables.height;
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::ComputeStabilizationParameters(
    const ElementVariables& rVariables,
    double& rTauU,
    double& rTauF)
{
    // Get element values
    const double elem_size = this->GetGeometry().Length();
    const double CTau = rVariables.stab_factor;  // 0.005 ~ 0.002

    // Wave mixed form stabilization
    rTauU = CTau * elem_size * std::sqrt(rVariables.wave_vel_2);
    rTauF = CTau * elem_size / (std::sqrt(rVariables.wave_vel_2) + rVariables.epsilon);

    // Discontinuity capturing
    const double mom_div_norm = std::abs(rVariables.momentum_div);
    const double surface_grad_norm = norm_2(rVariables.surface_grad);
    rTauU += 0.5 * 0.1 * elem_size * mom_div_norm;
    rTauF += 0.5 * 0.1 * elem_size * surface_grad_norm;
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::ComputeConvectionStabilizationParameters(
    const ElementVariables& rVariables,
    double& rTau)
{
    // Get element values
    const double elem_size = this->GetGeometry().Length();
    const double CTau = rVariables.stab_factor;  // 0.005 ~ 0.002

    // Convective stabilization
    if (TFramework == Eulerian)
    {
        const double vel_modulus = norm_2(rVariables.velocity) + rVariables.epsilon;
        rTau = CTau * elem_size / vel_modulus;
    }
    else
    {
        rTau = 0.0;
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::BuildAuxiliaryMatrices(
        const array_1d<double, TNumNodes>& rN,
        const BoundedMatrix<double, TNumNodes, 2>& rDN_DX,
        ElementVariables& rVariables)
{
    rVariables.N_q     = ZeroMatrix(2, rVariables.LocalSize); // Momentum balance test functions
    rVariables.N_h     = ZeroVector(rVariables.LocalSize);   // Mass balance test functions
    rVariables.DN_DX_q = ZeroVector(rVariables.LocalSize);
    rVariables.DN_DX_h = ZeroMatrix(2, rVariables.LocalSize);
    rVariables.Grad_q1 = ZeroMatrix(2, rVariables.LocalSize);
    rVariables.Grad_q2 = ZeroMatrix(2, rVariables.LocalSize);

    // Build the shape and derivatives functions at the Gauss point
    for(size_t node = 0; node < TNumNodes; ++node)
    {
        // Momentum balance test functions
        rVariables.N_q(0,   node*3) = rN[node];
        rVariables.N_q(1, 1+node*3) = rN[node];
        // Mass balance test funtions
        rVariables.N_h[2+node*3] = rN[node];
        // Momentum balance test functions divergence
        rVariables.DN_DX_q[  node*3] = rDN_DX(node,0);
        rVariables.DN_DX_q[1+node*3] = rDN_DX(node,1);
        // Mass balance test functions gradients
        rVariables.DN_DX_h(0, 2+node*3) = rDN_DX(node,0);
        rVariables.DN_DX_h(1, 2+node*3) = rDN_DX(node,1);
        // Momentum balance test functions gradients
        rVariables.Grad_q1(0,   node*3) = rDN_DX(node,0);
        rVariables.Grad_q1(1, 1+node*3) = rDN_DX(node,0);
        rVariables.Grad_q2(0,   node*3) = rDN_DX(node,1);
        rVariables.Grad_q2(1, 1+node*3) = rDN_DX(node,1);
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::AddInertiaTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{
    LocalMatrixType mass_matrix;
    mass_matrix = prod(trans(rVariables.N_q), rVariables.N_q);
    mass_matrix += outer_prod(rVariables.N_h, rVariables.N_h);
    rLeftHandSideMatrix += rVariables.dt_inv * rVariables.porosity * mass_matrix;
    rRightHandSideVector += rVariables.dt_inv * rVariables.porosity * prod(mass_matrix, rVariables.prev_unk);
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::AddConvectiveTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{
    if (TFramework == Eulerian)
    {
        BoundedMatrix<double,2, rVariables.LocalSize> convection_operator;
        convection_operator = rVariables.velocity[0] * rVariables.Grad_q1;
        convection_operator += rVariables.velocity[1] * rVariables.Grad_q2;

        // Convective term
        rLeftHandSideMatrix += prod(trans(rVariables.N_q), convection_operator); // q * u * grad_q
        rLeftHandSideMatrix += rVariables.velocity_div * prod(trans(rVariables.N_q), rVariables.N_q); // q * div_u * q

        // Convection stabilization term
        double tau;
        this->ComputeConvectionStabilizationParameters(rVariables, tau);
        rLeftHandSideMatrix += tau * prod(trans(convection_operator), convection_operator);
        rLeftHandSideMatrix += tau * rVariables.velocity_div * prod(trans(convection_operator), rVariables.N_q);
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::AddWaveTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{
    const double p = rVariables.porosity;
    rLeftHandSideMatrix += p * rVariables.wave_vel_2 * prod(trans(rVariables.N_q), rVariables.DN_DX_h); // q * grad_h (momentum balance)
    rLeftHandSideMatrix += p * outer_prod(rVariables.N_h, rVariables.DN_DX_q); // h * div_q (mass balance)
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::AddFrictionTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{
    const double q = norm_2(rVariables.projected_momentum) + rVariables.epsilon;
    const double h73 = std::pow(rVariables.height, 2.333333333333333) + rVariables.epsilon;
    LocalMatrixType vector_mass_matrix = prod(trans(rVariables.N_q), rVariables.N_q);
    rLeftHandSideMatrix += rVariables.gravity * rVariables.manning2 * q / h73 * vector_mass_matrix;
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::AddStabilizationTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{
    double tau_u;
    double tau_f;
    this->ComputeStabilizationParameters(rVariables, tau_u, tau_f);

    // Wave term stabilization
    LocalMatrixType vector_diffusion = outer_prod(rVariables.DN_DX_q, rVariables.DN_DX_q);
    LocalMatrixType scalar_diffusion = prod(trans(rVariables.DN_DX_h), rVariables.DN_DX_h);
    const double p = rVariables.porosity;
    rLeftHandSideMatrix += p * tau_u * vector_diffusion;
    rLeftHandSideMatrix += p * tau_f * rVariables.wave_vel_2 * scalar_diffusion;

    // Dry domain stabilization
    LocalMatrixType vector_lumped_mass = ZeroMatrix(rVariables.LocalSize, rVariables.LocalSize);
    for (size_t i = 0; i < TNumNodes; ++i)
    {
        vector_lumped_mass(  3*i,  3*i) = 1.0;
        vector_lumped_mass(1+3*i,1+3*i) = 1.0;
    }
    vector_lumped_mass *= rVariables.lumping_factor;
    rLeftHandSideMatrix += (1 - p) * rVariables.discharge_penalty * vector_lumped_mass;
    rLeftHandSideMatrix += (1 - p) * rVariables.permeability * scalar_diffusion;

    bool full_subscales = false;
    if (full_subscales) {
        // Friction stabilization
        const double q = norm_2(rVariables.projected_momentum) + rVariables.epsilon;
        const double h73 = std::pow(rVariables.height, 2.333333333333333) + rVariables.epsilon;
        LocalMatrixType friction_stab = prod(trans(rVariables.DN_DX_h), rVariables.N_q);
        rLeftHandSideMatrix += p * tau_f * rVariables.gravity * rVariables.manning2 * q / h73 * friction_stab;

        // Dynamic stabilization
        const double dt_inv = rVariables.dt_inv;
        LocalMatrixType vector_dyn_stab = prod(trans(rVariables.DN_DX_h), rVariables.N_q);
        LocalMatrixType scalar_dyn_stab = outer_prod(rVariables.DN_DX_q, rVariables.N_h);
        rLeftHandSideMatrix += p * tau_f * dt_inv * vector_dyn_stab;
        rLeftHandSideMatrix += p * tau_f * dt_inv * scalar_dyn_stab;
        rRightHandSideVector += p * tau_f * dt_inv * prod(vector_dyn_stab, rVariables.prev_unk);
        rRightHandSideVector += p * tau_f * dt_inv * prod(scalar_dyn_stab, rVariables.prev_unk);
    }
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::AddSourceTerms(
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{
    LocalMatrixType scalar_mass_matrix = outer_prod(rVariables.N_h, rVariables.N_h);
    rRightHandSideVector += prod(scalar_mass_matrix, rVariables.rain);
}


template< size_t TNumNodes, ElementFramework TFramework >
void SWE<TNumNodes, TFramework>::CalculateLumpedMassMatrix(LocalMatrixType& rMassMatrix)
{
    constexpr size_t local_size = 3 * TNumNodes;
    rMassMatrix = IdentityMatrix(local_size, local_size);
    rMassMatrix /= static_cast<double>(TNumNodes);
}


template class SWE<3, Eulerian>;
template class SWE<4, Eulerian>;
template class SWE<3, PFEM2>;
template class SWE<4, PFEM2>;

} // namespace Kratos
