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
#include "includes/cfd_variables.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "shallow_water_application_variables.h"
#include "conserved_element.h"

namespace Kratos
{

template<size_t TNumNodes>
int ConservedElement<TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    // Base class checks for positive Jacobian and Id > 0
    int err = Element::Check(rCurrentProcessInfo);
    if(err != 0) return err;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(MOMENTUM)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(HEIGHT)
    KRATOS_CHECK_VARIABLE_KEY(TOPOGRAPHY)
    KRATOS_CHECK_VARIABLE_KEY(RAIN)
    KRATOS_CHECK_VARIABLE_KEY(MANNING)
    KRATOS_CHECK_VARIABLE_KEY(GRAVITY)
    KRATOS_CHECK_VARIABLE_KEY(DELTA_TIME)
    KRATOS_CHECK_VARIABLE_KEY(DYNAMIC_TAU)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( size_t i = 0; i < TNumNodes; i++ )
    {
        Node<3>& node = this->GetGeometry()[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MOMENTUM, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TOPOGRAPHY, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RAIN, node)

        KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_X, node)
        KRATOS_CHECK_DOF_IN_NODE(MOMENTUM_Y, node)
        KRATOS_CHECK_DOF_IN_NODE(HEIGHT, node)
    }

    return err;
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    const size_t element_size = TNumNodes*3;
    if(rResult.size() != element_size)
        rResult.resize(element_size, false); // False says not to preserve existing storage!!

    const GeometryType& geom = GetGeometry();

    int counter=0;
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rResult[counter++] = geom[i].GetDof(MOMENTUM_X).EquationId();
        rResult[counter++] = geom[i].GetDof(MOMENTUM_Y).EquationId();
        rResult[counter++] = geom[i].GetDof(HEIGHT).EquationId();
    }
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    const size_t element_size = TNumNodes*3;
    if(rElementalDofList.size() != element_size)
        rElementalDofList.resize(element_size);

    const GeometryType& geom = GetGeometry();

    int counter=0;
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rElementalDofList[counter++] = geom[i].pGetDof(MOMENTUM_X);
        rElementalDofList[counter++] = geom[i].pGetDof(MOMENTUM_Y);
        rElementalDofList[counter++] = geom[i].pGetDof(HEIGHT);
    }
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    // Resize of the Left and Right Hand side
    constexpr size_t element_size = TNumNodes*3;
    if(rLeftHandSideMatrix.size1() != element_size)
        rLeftHandSideMatrix.resize(element_size,element_size, false); // False says not to preserve existing storage!!

    if(rRightHandSideVector.size() != element_size)
        rRightHandSideVector.resize(element_size, false);             // False says not to preserve existing storage!!

    const GeometryType& geom = this->GetGeometry();

    ElementVariables variables;
    this->InitializeElementVariables(variables, rCurrentProcessInfo);

    const BoundedMatrix<double, TNumNodes, TNumNodes> N_container = geom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2 );
    GeometryType::ShapeFunctionsGradientsType DN_DX_container(TNumNodes);
    geom.ShapeFunctionsIntegrationPointsGradients(DN_DX_container, GeometryData::GI_GAUSS_2);

    BoundedMatrix<double, TNumNodes, 2> DN_DX;
    array_1d<double, TNumNodes> N;
    double area = geom.Area();

    this->GetNodalValues(variables);
    this->CalculateElementValues(DN_DX_container, variables);

    rLeftHandSideMatrix = ZeroMatrix(element_size, element_size);
    rRightHandSideVector = ZeroVector(element_size);

    for (size_t g = 0; g < TNumNodes; ++g)
    {
        N = row(N_container, g);
        DN_DX = DN_DX_container[g];

        this->BuildAuxiliaryMatrices(N, DN_DX, variables);

        this->AddInertiaTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        // this->AddConvectiveTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        this->AddWaveTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        this->AddFrictionTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        this->AddStabilizationTerms(rLeftHandSideMatrix, rRightHandSideVector, variables);
        this->AddSourceTerms(rRightHandSideVector, variables);
    }

    // Substracting the Dirichlet term (since we use a residualbased approach)
    rRightHandSideVector -= prod(rLeftHandSideMatrix, variables.unknown);

    rRightHandSideVector *= area * variables.lumping_factor;
    rLeftHandSideMatrix  *= area * variables.lumping_factor;

    double residual = norm_1(rRightHandSideVector);
    this->SetValue(RESIDUAL_NORM, residual);
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::GetValueOnIntegrationPoints(
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


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::InitializeElementVariables(
    ElementVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    const double delta_t = rCurrentProcessInfo[DELTA_TIME];
    rVariables.epsilon = rCurrentProcessInfo[DRY_HEIGHT];
    rVariables.dt_inv = 1.0 / delta_t;
    rVariables.lumping_factor = 1.0 / static_cast<double>(TNumNodes);
    rVariables.dyn_tau = rCurrentProcessInfo[DYNAMIC_TAU];
    rVariables.gravity = rCurrentProcessInfo[GRAVITY_Z];
    rVariables.manning2 = 0.0;

    GeometryType& geom = GetGeometry();
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rVariables.manning2 += geom[i].FastGetSolutionStepValue(MANNING);
    }
    rVariables.manning2 *= rVariables.lumping_factor;
    rVariables.manning2 = std::pow(rVariables.manning2, 2);
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::GetNodalValues(ElementVariables& rVariables)
{
    GeometryType& geom = GetGeometry();
    size_t counter = 0;
    for (size_t i = 0; i < TNumNodes; i++)
    {
        rVariables.depth[counter] = 0;
        rVariables.source[counter] = 0;
        rVariables.unknown[counter] = geom[i].FastGetSolutionStepValue(MOMENTUM_X);
        rVariables.prev_unk[counter] = geom[i].FastGetSolutionStepValue(MOMENTUM_X, 1);
        counter++;

        rVariables.depth[counter] = 0;
        rVariables.source[counter] = 0;
        rVariables.unknown[counter] = geom[i].FastGetSolutionStepValue(MOMENTUM_Y);
        rVariables.prev_unk[counter] = geom[i].FastGetSolutionStepValue(MOMENTUM_Y, 1);
        counter++;

        rVariables.depth[counter] = -geom[i].FastGetSolutionStepValue(TOPOGRAPHY);
        rVariables.source[counter] = geom[i].FastGetSolutionStepValue(RAIN);
        rVariables.unknown[counter] = geom[i].FastGetSolutionStepValue(HEIGHT);
        rVariables.prev_unk[counter] = geom[i].FastGetSolutionStepValue(HEIGHT, 1);
        counter++;
    }
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::CalculateElementValues(
    const GeometryType::ShapeFunctionsGradientsType& rDN_DXContainer,
    ElementVariables& rVariables)
{
    // Initialize outputs
    rVariables.momentum = ZeroVector(3);
    rVariables.height = 0.0;
    rVariables.velocity = ZeroVector(3);
    rVariables.momentum_div = 0.0;
    rVariables.velocity_div = 0.0;
    rVariables.height_grad = ZeroVector(2);

    GeometryType& geom = GetGeometry();

    // integrate over the element
    for (size_t i = 0; i < TNumNodes; ++i)
    {
        auto flow = geom[i].FastGetSolutionStepValue(MOMENTUM);
        auto vel = geom[i].FastGetSolutionStepValue(VELOCITY);
        auto h = geom[i].FastGetSolutionStepValue(HEIGHT);
        auto mom = geom[i].FastGetSolutionStepValue(MOMENTUM, 1);

        rVariables.velocity += vel;
        rVariables.height += h;
        rVariables.momentum += mom;

        for (size_t g = 0; g < TNumNodes; ++g)
        {
            rVariables.momentum_div += rDN_DXContainer[g](i,0) * geom[i].FastGetSolutionStepValue(MOMENTUM_X);
            rVariables.momentum_div += rDN_DXContainer[g](i,1) * geom[i].FastGetSolutionStepValue(MOMENTUM_Y);
            rVariables.velocity_div += rDN_DXContainer[g](i,0) * geom[i].FastGetSolutionStepValue(VELOCITY_X);
            rVariables.velocity_div += rDN_DXContainer[g](i,1) * geom[i].FastGetSolutionStepValue(VELOCITY_Y);
            rVariables.height_grad[0] += rDN_DXContainer[g](i,0) * geom[i].FastGetSolutionStepValue(HEIGHT);
            rVariables.height_grad[1] += rDN_DXContainer[g](i,1) * geom[i].FastGetSolutionStepValue(HEIGHT);
        }
    }

    rVariables.velocity *= rVariables.lumping_factor;
    rVariables.height *= rVariables.lumping_factor;
    rVariables.height = std::max(rVariables.height, 0.0);
    rVariables.momentum *= rVariables.lumping_factor;
    rVariables.momentum_div *= rVariables.lumping_factor;
    rVariables.velocity_div *= rVariables.lumping_factor;
    rVariables.height_grad *= rVariables.lumping_factor;
    rVariables.wave_vel_2 = rVariables.gravity * rVariables.height;
    if (rVariables.height < rVariables.epsilon) rVariables.velocity *= 0.0;
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::ComputeStabilizationParameters(
    const ElementVariables& rVariables,
    double& rTauQ,
    double& rTauH)
{
    // Get element values
    const double elem_size = this->GetGeometry().Length();
    const double CTau = rVariables.dyn_tau;

    // Wave mixed form stabilization
    double vel = norm_2(rVariables.velocity)*0;
    double c = std::sqrt(rVariables.wave_vel_2);
    rTauQ = CTau * elem_size * (c + vel);
    rTauH = CTau * elem_size / (c + rVariables.epsilon);

    // Discontinuity capturing
    // const double mom_div_norm = std::abs(rVariables.momentum_div);
    // const double surface_grad_norm = norm_2(rVariables.height_grad);
    // rTauQ += 0.5 * 0.1 * elem_size * mom_div_norm;
    // rTauH += 0.5 * 0.1 * elem_size * surface_grad_norm;
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::BuildAuxiliaryMatrices(
        const array_1d<double, TNumNodes>& rN,
        const BoundedMatrix<double, TNumNodes, 2>& rDN_DX,
        ElementVariables& rVariables)
{
    rVariables.N_v     = ZeroMatrix(2, rVariables.LocalSize); // Momentum balance test functions
    rVariables.N_f     = ZeroVector(rVariables.LocalSize);   // Mass balance test functions
    rVariables.Div_v   = ZeroVector(rVariables.LocalSize);
    rVariables.Grad_f  = ZeroMatrix(2, rVariables.LocalSize);
    rVariables.Grad_v1 = ZeroMatrix(2, rVariables.LocalSize);
    rVariables.Grad_v2 = ZeroMatrix(2, rVariables.LocalSize);

    // Build the shape and derivatives functions at the Gauss point
    for(size_t node = 0; node < TNumNodes; ++node)
    {
        // Momentum balance test functions
        rVariables.N_v(0,   node*3) = rN[node];
        rVariables.N_v(1, 1+node*3) = rN[node];
        // Mass balance test funtions
        rVariables.N_f[2+node*3] = rN[node];
        // Momentum balance test functions divergence
        rVariables.Div_v[  node*3] = rDN_DX(node,0);
        rVariables.Div_v[1+node*3] = rDN_DX(node,1);
        // Mass balance test functions gradients
        rVariables.Grad_f(0, 2+node*3) = rDN_DX(node,0);
        rVariables.Grad_f(1, 2+node*3) = rDN_DX(node,1);
        // Momentum balance test functions gradients
        rVariables.Grad_v1(0,   node*3) = rDN_DX(node,0);
        rVariables.Grad_v1(1, 1+node*3) = rDN_DX(node,0);
        rVariables.Grad_v2(0,   node*3) = rDN_DX(node,1);
        rVariables.Grad_v2(1, 1+node*3) = rDN_DX(node,1);
    }
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::AddInertiaTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{
    LocalMatrixType mass_matrix;
    mass_matrix = prod(trans(rVariables.N_v), rVariables.N_v); // <v, v>
    mass_matrix += outer_prod(rVariables.N_f, rVariables.N_f); // <f, f>
    if (rVariables.height > rVariables.epsilon) {
        rLeftHandSideMatrix += rVariables.dt_inv * mass_matrix;
        rRightHandSideVector += rVariables.dt_inv * prod(mass_matrix, rVariables.prev_unk);
    }
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::AddConvectiveTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{
    BoundedMatrix<double,2, rVariables.LocalSize> convection_operator;
    convection_operator = rVariables.velocity[0] * rVariables.Grad_v1;
    convection_operator += rVariables.velocity[1] * rVariables.Grad_v2;
    rLeftHandSideMatrix += prod(trans(rVariables.N_v), convection_operator); // <v, u * grad_v>
    rLeftHandSideMatrix += rVariables.velocity_div * prod(trans(rVariables.N_v), rVariables.N_v); // <v, div_u * v>
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::AddWaveTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{
    auto N_grad_f = prod(trans(rVariables.N_v), rVariables.Grad_f); // <v, grad_f> (momentum balance)
    rLeftHandSideMatrix += rVariables.wave_vel_2 * N_grad_f;
    rLeftHandSideMatrix += outer_prod(rVariables.N_f, rVariables.Div_v); // <f, div_v> (mass balance)
    rRightHandSideVector += rVariables.wave_vel_2 * prod(N_grad_f, rVariables.depth);
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::AddFrictionTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{
    const double q = norm_2(rVariables.momentum) + rVariables.epsilon;
    const double h73 = std::pow(rVariables.height, 2.333333333333333) + rVariables.epsilon;
    LocalMatrixType vector_mass_matrix = prod(trans(rVariables.N_v), rVariables.N_v); // <v, v>
    double friction_factor = rVariables.gravity * rVariables.manning2 * q / h73;
    if (rVariables.height < rVariables.epsilon) friction_factor = 1.0e6;
    rLeftHandSideMatrix += friction_factor * vector_mass_matrix;
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::AddStabilizationTerms(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{
    double tau_h;
    double tau_q;
    this->ComputeStabilizationParameters(rVariables, tau_q, tau_h);

    // Wave term stabilization
    LocalMatrixType f_grad_subscale = prod(trans(rVariables.Grad_f), rVariables.Grad_f);
    LocalMatrixType v_div_subscale = outer_prod(rVariables.Div_v, rVariables.Div_v);
    rLeftHandSideMatrix += tau_q * rVariables.wave_vel_2 * f_grad_subscale;
    rLeftHandSideMatrix += tau_h * v_div_subscale;
    rRightHandSideVector += tau_q * rVariables.wave_vel_2 * prod(f_grad_subscale, rVariables.depth);

    // Friction stabilization
    const double q = norm_2(rVariables.momentum) + rVariables.epsilon;
    const double h73 = std::pow(rVariables.height, 2.333333333333333) + rVariables.epsilon;
    LocalMatrixType friction_stab = prod(trans(rVariables.Grad_f), rVariables.N_v);
    rLeftHandSideMatrix += tau_q * rVariables.gravity * rVariables.manning2 * q / h73 * friction_stab;

    // Dynamic stabilization
    if (rVariables.height > rVariables.epsilon)
    {
        const double dt_inv = rVariables.dt_inv;
        LocalMatrixType v_mass_subscale = prod(trans(rVariables.Grad_f), rVariables.N_v);
        LocalMatrixType f_mass_subscale = outer_prod(rVariables.Div_v, rVariables.N_f);
        rLeftHandSideMatrix += tau_q * dt_inv * v_mass_subscale;
        rLeftHandSideMatrix += tau_h * dt_inv * f_mass_subscale;
        rRightHandSideVector += tau_q * dt_inv * prod(v_mass_subscale, rVariables.prev_unk);
        rRightHandSideVector += tau_h * dt_inv * prod(f_mass_subscale, rVariables.prev_unk);
    }
    else {
        // Stationary stabilization
        rLeftHandSideMatrix += 1.0e6 * f_grad_subscale;
        rLeftHandSideMatrix += 1.0e6 * v_div_subscale;
    }

    // // Convection stabilization
    // BoundedMatrix<double,2, rVariables.LocalSize> convection_operator;
    // convection_operator = rVariables.velocity[0] * rVariables.Grad_v1;
    // convection_operator += rVariables.velocity[1] * rVariables.Grad_v2;
    // rLeftHandSideMatrix += tau_q * prod(trans(rVariables.Grad_f), convection_operator);
    // rLeftHandSideMatrix += tau_q * rVariables.velocity_div * prod(trans(rVariables.Grad_f), rVariables.N_v);
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::AddSourceTerms(
    VectorType& rRightHandSideVector,
    ElementVariables& rVariables)
{
    LocalMatrixType scalar_mass_matrix = outer_prod(rVariables.N_f, rVariables.N_f);
    rRightHandSideVector += prod(scalar_mass_matrix, rVariables.source);
}


template<size_t TNumNodes>
void ConservedElement<TNumNodes>::CalculateLumpedMassMatrix(LocalMatrixType& rMassMatrix)
{
    constexpr size_t local_size = 3 * TNumNodes;
    rMassMatrix = IdentityMatrix(local_size, local_size);
    rMassMatrix /= static_cast<double>(TNumNodes);
}


template class ConservedElement<3>;
template class ConservedElement<4>;

} // namespace Kratos