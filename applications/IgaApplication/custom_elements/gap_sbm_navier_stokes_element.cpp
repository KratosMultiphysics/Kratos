//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Nicolo' Antonelli
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"

// Application includes
#include "custom_elements/gap_sbm_navier_stokes_element.h"

namespace Kratos
{

GapSbmNavierStokesElement::GapSbmNavierStokesElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(
        NewId,
        pGeometry)
{
}

GapSbmNavierStokesElement::GapSbmNavierStokesElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(
        NewId,
        pGeometry,
        pProperties)
{
}

Element::Pointer GapSbmNavierStokesElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GapSbmNavierStokesElement>(NewId, GetSurrogateGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer GapSbmNavierStokesElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GapSbmNavierStokesElement>(NewId, pGeom, pProperties);
}

// Deconstructor

GapSbmNavierStokesElement::~GapSbmNavierStokesElement()
{
}

void GapSbmNavierStokesElement:: Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMemberVariables();
    InitializeSbmMemberVariables();
    InitializeMaterial();
}


void GapSbmNavierStokesElement::InitializeMaterial()
{
    KRATOS_TRY
    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_surrogate_geometry = GetSurrogateGeometry();
        const Properties& r_properties = GetProperties();        
        const std::size_t number_of_control_points = r_surrogate_geometry.size();
        Vector N_sum_vec = ZeroVector(number_of_control_points);
        ComputeTaylorExpansionContribution(N_sum_vec);
    
        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLaw->InitializeMaterial( r_properties, r_surrogate_geometry, N_sum_vec);

    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH( "" );

}

void GapSbmNavierStokesElement::InitializeMemberVariables()
{
    // // Compute class memeber variables

    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const auto& r_DN_De = r_surrogate_geometry.ShapeFunctionsLocalGradients(r_surrogate_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();

    KRATOS_ERROR_IF(mDim != 2) << "GapSbmNavierStokesElement momentarily only supports 2D conditions, but the current dimension is" << mDim << std::endl;
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }
    mBasisFunctionsOrder *= 2; 

    const auto& r_geometry = GetGeometry(); // only for integration points
    // calculate the integration weight
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    const double integration_weight = r_integration_points[0].Weight();

    SetValue(INTEGRATION_WEIGHT, integration_weight);
}

void GapSbmNavierStokesElement::InitializeSbmMemberVariables()
{
    const auto& r_geometry = this->GetGeometry();
    const auto& r_surrogate_geometry = GetSurrogateGeometry();

    mDistanceVector.resize(3);
    noalias(mDistanceVector) = r_geometry.Center().Coordinates() - r_surrogate_geometry.Center().Coordinates();
}


void GapSbmNavierStokesElement::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo)
{
    const auto& r_surrogate_geometry = GetSurrogateGeometry();
    const unsigned int number_of_control_points = r_surrogate_geometry.size();
    
    // reading integration points and local gradients
    const std::size_t mat_size = number_of_control_points * mDim;
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS

    // compute Taylor expansion contribution: grad_H_sum
    Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(grad_N_sum_transposed);
    Matrix grad_N_sum = trans(grad_N_sum_transposed);

    Matrix B_sum = ZeroMatrix(mDim,mat_size);
    CalculateB(B_sum, grad_N_sum);

    // obtain the tangent constitutive matrix at the true position
    ConstitutiveLaw::Parameters values_true(r_surrogate_geometry, GetProperties(), rCurrentProcessInfo);
    Vector old_displacement_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(old_displacement_coefficient_vector);
    Vector old_strain_on_true = prod(B_sum, old_displacement_coefficient_vector);
    const std::size_t strain_size_true = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
    ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true);

    const Matrix& r_D_on_surrogate = values_true.GetConstitutiveMatrix();
    Vector& r_stress_vector_on_surrogate = values_true.GetStressVector();

    // DELETE:
    // noalias(rLeftHandSideMatrix) += integration_weight * prod(trans(B_sum), Matrix(prod(r_D_on_surrogate, B_sum))); 

    
    
    
    //____________________from NS element___ need to use B_sum and grad_N_sum
    
    // Obtain required constants
    // GeometryType& r_geometry = this->GetGeometry();
    const unsigned int number_of_points = r_surrogate_geometry.PointsNumber();

    std::size_t num_dofs_per_node = number_of_points * (mDim + 1);
    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != num_dofs_per_node)
        rLeftHandSideMatrix.resize(num_dofs_per_node,num_dofs_per_node,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(num_dofs_per_node,num_dofs_per_node); //resetting LHS
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != num_dofs_per_node)
        rRightHandSideVector.resize(num_dofs_per_node,false);
    noalias(rRightHandSideVector) = ZeroVector(num_dofs_per_node); //resetting RHS

    // compute Taylor expansion contribution: H_sum_vec
    Vector N_sum_vec = ZeroVector(number_of_points);
    ComputeTaylorExpansionContribution(N_sum_vec);

    double mu_effective;
    mpConstitutiveLaw->CalculateValue(values_true, MU, mu_effective);
    if (mu_effective == 0) {
        // if the mu_effective is zero, we use the viscosity value
        this->EvaluateInPoint(mu_effective,VISCOSITY,N_sum_vec,r_surrogate_geometry); //TODO: check
    }

    array_1d<double,3> body_force = ZeroVector(3);
    body_force = this->GetValue(BODY_FORCE); // TODO: maybe on the surrogate?
    double viscosity = mu_effective;
    const Properties& r_properties = GetProperties();
    double density = r_properties[DENSITY];

    double adv_norm;
    // Calculate the advective norm
    CalculateAdvectiveNorm(N_sum_vec, adv_norm);

    // Calculate stabilization constants
    CalculateTau(mu_effective, density, adv_norm);

    // Add velocity terms in momentum equation
    this->AddMomentumTerms(rLeftHandSideMatrix,rRightHandSideVector,body_force,mTauTwo,N_sum_vec,grad_N_sum,integration_weight, r_D_on_surrogate, r_stress_vector_on_surrogate);

    // Add velocity-pressure terms
    this->AddContinuityTerms(rLeftHandSideMatrix,rRightHandSideVector,body_force,density,mTauOne,N_sum_vec,grad_N_sum,integration_weight);

    // Add Second-Order stabilization terms from VMS
    this->AddSecondOrderStabilizationTerms(rLeftHandSideMatrix,rRightHandSideVector,mTauOne,N_sum_vec,grad_N_sum,integration_weight,r_D_on_surrogate);

    // Add Convective terms ["D" terms] 
    this->AddConvectiveTerms(rLeftHandSideMatrix,rRightHandSideVector,body_force,density,mTauOne,N_sum_vec,grad_N_sum,integration_weight,r_D_on_surrogate);

}

void GapSbmNavierStokesElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{   
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

void GapSbmNavierStokesElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}


void GapSbmNavierStokesElement::CalculateAdvectiveNorm(
    const ShapeFunctionsType &rN_sum,
    double& adv_norm)   
{
    // Correct r_geometry cause we are then using H_sum
    auto r_geometry = this->GetGeometry();
    const unsigned int number_of_points = r_geometry.PointsNumber();

    // Compute advective velocity: a = sum_j N_j * u_j
    array_1d<double, 3> advective_velocity = ZeroVector(3);
    for (unsigned int j = 0; j < number_of_points; ++j)
    {
        const array_1d<double, 3>& u_j = r_geometry[j].FastGetSolutionStepValue(VELOCITY); /// TODO:
        for (unsigned int d = 0; d < mDim; ++d)
            advective_velocity[d] += rN_sum[j] * u_j[d];
    }

    // Compute ||a||
    adv_norm = norm_2(advective_velocity);
}

void GapSbmNavierStokesElement::CalculateTau(
    const double MuEffective, 
    const double Density, 
    const double AdvectiveNorm)
{   
    // // Estimate element size
    double h = ElementSize();

    // The stabilization parameters from: 
    // https://www.researchgate.net/profile/Ramon-Codina-2/publication/222977355

    // mTauOne = std::pow(h, 2) / ( 4.0 * MuEffective ); ///mBasisFunctionsOrder; 
    // mTauTwo = MuEffective;


    // Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    // double h = std::min(mesh_size_uv[0], mesh_size_uv[1]);
    constexpr double c1 = 4.0;
    constexpr double c2 = 2.0;

    double TauDynamic = 1.0;
    double DeltaTime = 1e20;
    double Viscosity = MuEffective; // correct?

    // Compute tau1
    mTauOne = 1.0 / (
        (Density * TauDynamic / DeltaTime) +
        (c2 * Density * AdvectiveNorm / (h)) +
        (c1 * Viscosity / (h * h))
    );

    // Compute tau2
    mTauTwo = h * h / (c2 * mTauOne);

    // mTauOne = 0;
    // mTauTwo = 0;
}

double GapSbmNavierStokesElement::ElementSize()
{
    return 1.128379167 * sqrt(this->GetGeometry().DomainSize()); 
}

void GapSbmNavierStokesElement::AddMomentumTerms(MatrixType &rLHS,
        VectorType &rRHS,
        const array_1d<double,3> &rBodyForce,
        const double TauTwo,
        const ShapeFunctionsType &rN,
        const ShapeDerivativesType &rDN_DX,
        const double Weight,
        const Matrix& rD,
        Vector& rStressVector)
{
    const unsigned int number_of_nodes = this->GetSurrogateGeometry().PointsNumber();
    const unsigned int block_size = mDim+1;
    const std::size_t mat_size = number_of_nodes * mDim;
    auto r_geometry = GetGeometry(); // TODO: surrogate? Don't think so

    // compute Taylor expansion contribution: grad_H_sum
    Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_nodes);
    ComputeGradientTaylorExpansionContribution(grad_N_sum_transposed);
    Matrix grad_N_sum = trans(grad_N_sum_transposed);
    Matrix B_sum = ZeroMatrix(mDim,mat_size);
    CalculateB(B_sum, grad_N_sum);

    // Dynamic viscosity already included in rD
    Matrix diffusion_term_matrix = Weight * prod(trans(B_sum), Matrix(prod(rD, B_sum)));
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        for (IndexType j = 0; j < number_of_nodes; ++j)
        {
            // Add only to velocity DOFs (i.e., positions 0, 1, 3, 4, 6, 7, ...)
            for (IndexType dim1 = 0; dim1 < mDim; ++dim1)
            {
                for (IndexType dim2 = 0; dim2 < mDim; ++dim2)
                {
                    // diffusion term
                    rLHS(i * block_size + dim1, j * block_size + dim2) += diffusion_term_matrix(i * mDim + dim1, j * mDim + dim2);
                }
            }
        }
    }
    // RHS corresponding term (diffusion term)
    Vector internal_forces = Weight * prod(trans(B_sum), rStressVector);
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        // Add only to the velocity DOFs (i.e., positions 0, 1, 3, 4, 6, 7, ...)
        for (IndexType dim1 = 0; dim1 < mDim; ++dim1)
        {
            // viscous + pressure term
            rRHS(i * block_size + dim1) -= internal_forces(i * mDim + dim1);
        }  
    }

    unsigned int first_row = 0;
    unsigned int first_col = 0;

    for(unsigned int i = 0; i < number_of_nodes; ++i)
    {
        // Body force
        for(unsigned int d = 0; d < mDim; ++d) {
            rRHS[first_row + d] += Weight * rN[i] * rBodyForce[d];
        }

        for(unsigned int j = 0; j < number_of_nodes; ++j)
        {
            // Stabilization div(u) div(v)
            for (unsigned int m = 0; m < mDim; ++m) {
                for (unsigned int n = 0; n < mDim; ++n) {
                    rLHS(first_row+m,first_col+n) += TauTwo * Weight * rDN_DX(i,m) * rDN_DX(j,n);
                }
            }

            // Update column index
            first_col += block_size;
        }
        // RHS corresponding term: Stabilization div(u) div(v)
        double div_u_current = 0.0;
        for(unsigned int j = 0; j < number_of_nodes; ++j) {
            div_u_current += r_geometry[j].FastGetSolutionStepValue(VELOCITY_X) * rDN_DX(j, 0) + 
                             r_geometry[j].FastGetSolutionStepValue(VELOCITY_Y) * rDN_DX(j, 1) ;
        } 
        for (unsigned int m = 0; m < mDim; ++m) {
            rRHS(first_row+m) -= TauTwo * Weight * rDN_DX(i,m) * div_u_current ;
        }

        // Update matrix indices
        first_row += block_size;
        first_col = 0;
    }
}


void GapSbmNavierStokesElement::AddContinuityTerms(MatrixType &rLHS,
        VectorType &rRHS,
        const array_1d<double,3> &rBodyForce,
        const double Density,
        const double TauOne,
        const ShapeFunctionsType &rN,
        const ShapeDerivativesType &rDN_DX,
        const double Weight)
{
    const unsigned int number_of_nodes = this->GetGeometry().PointsNumber();
    const unsigned int block_size = mDim+1;
    auto r_geometry = GetGeometry();

    unsigned int first_row = 0;
    unsigned int first_col = 0;

    double pressure_previous_iteration = 0.0;
    double div_u_current = 0.0;
    Vector grad_p_previous_iteration = ZeroVector(2);
    for(unsigned int j = 0; j < number_of_nodes; ++j) {
        div_u_current += r_geometry[j].FastGetSolutionStepValue(VELOCITY_X) * rDN_DX(j, 0) +
                        r_geometry[j].FastGetSolutionStepValue(VELOCITY_Y) * rDN_DX(j, 1) ;
        grad_p_previous_iteration[0] += r_geometry[j].GetSolutionStepValue(PRESSURE) * rDN_DX(j,0) ;
        grad_p_previous_iteration[1] += r_geometry[j].GetSolutionStepValue(PRESSURE) * rDN_DX(j,1) ;
        pressure_previous_iteration += r_geometry[j].GetSolutionStepValue(PRESSURE) * rN[j];
    }

    double div_term = 0.0;

    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        double Qi = 0.0;
        for (unsigned int d = 0; d < mDim; ++d) {
            Qi += rDN_DX(i,d) * rBodyForce[d];
        }
        // Stabilization term: grad(q) * f  -> [B1]
        rRHS[first_row + mDim] += Weight * TauOne * Qi;

        for (unsigned int j = 0; j < number_of_nodes; ++j)
        {
            double l_ij = 0.0;
            for (unsigned int d = 0; d < mDim; ++d)
            {
                l_ij += rDN_DX(i,d) * rDN_DX(j,d);
                div_term = Weight * rN[i] * rDN_DX(j,d);
                // div(u) * q
                rLHS(first_row+mDim,first_col + d) += div_term; // Divergence term
                // div(v) * p
                rLHS(first_col + d,first_row+mDim) -= div_term; // Gradient term
            }
            // Stabilization (grad p, grad q) -> [B3]
            rLHS(first_row+mDim,first_col+mDim) += Weight * TauOne * l_ij;

            // Update column index
            first_col += block_size;
        }

        // --- RHS corresponding term ---
        rRHS(first_row+mDim) -= Weight * rN[i] * div_u_current ; // q div(u)
        for (unsigned int m = 0; m < mDim; ++m) {
            // div(v) * p
            rRHS(first_row+m) += Weight * rDN_DX(i,m) * pressure_previous_iteration ;
            // Stabilization (grad p, grad q) -> [B3]
            rRHS(first_row+mDim) -= Weight * TauOne * rDN_DX(i,m) * grad_p_previous_iteration[m];
        }
        // Update matrix indices
        first_col = 0;
        first_row += block_size;
    }



    // // NAVIER-STOKES Stabilization grad(q) * (a * grad(u)) -> [B5]
    first_row = 0;
    first_col = 0;
    // --- Compute advective velocity at integration point: a = sum_j N_j * u_j ---
    array_1d<double, 3> advective_velocity = ZeroVector(3);
    for (unsigned int j = 0; j < number_of_nodes; ++j) {
        const array_1d<double, 3>& u_j = r_geometry[j].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int d = 0; d < mDim; ++d) {
            advective_velocity[d] += rN(j) * u_j[d];
        }
    }
    // --- Stabilization term: grad(q) * (a · grad(u)) [B5] ---
    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        for (unsigned int j = 0; j < number_of_nodes; ++j)
        {
            const auto& grad_u_j = row(rDN_DX, j); // ∇(u_j)

            for (unsigned int d = 0; d < mDim; ++d)
            {
                // Compute (a · ∇u_j^d)
                double adv_grad_u = 0.0;
                for (unsigned int k = 0; k < mDim; ++k)
                    adv_grad_u += advective_velocity[k] * grad_u_j[k];

                // LHS contribution
                rLHS(first_row + mDim, first_col + d) += Density * Weight * TauOne * rDN_DX(i, d) * adv_grad_u;

                // RHS contribution: a · ∇(u^d) using u_j^d coefficient
                const double u_jd = r_geometry[j].FastGetSolutionStepValue(VELOCITY)[d];
                double adv_u = 0.0;
                for (unsigned int k = 0; k < mDim; ++k)
                    adv_u += advective_velocity[k] * u_jd * rDN_DX(j, k);

                rRHS[first_row + mDim] -= Density * Weight * TauOne * rDN_DX(i, d) * adv_u;
            }

            first_col += block_size;
        }
        first_col = 0;
        first_row += block_size;
    }

}


void GapSbmNavierStokesElement::AddSecondOrderStabilizationTerms(MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const double TauOne,
        const ShapeFunctionsType &rN,
        const ShapeDerivativesType &rDN_DX,
        const double GaussWeight,
        const Matrix& rD)
{
    const unsigned int number_of_points = this->GetSurrogateGeometry().PointsNumber();
    const unsigned int block_size = mDim+1;
    auto r_geometry = GetGeometry();

    // Add third order term of the stabilization
    if (mBasisFunctionsOrder > 1) {

        Matrix B = ZeroMatrix(3,number_of_points*mDim);
        CalculateB(B, rDN_DX);

        // TODO: Taylor expansio of second order derivatives?
        const Matrix& DDN_DDe = GetGeometry().ShapeFunctionDerivatives(2, 0, GetGeometry().GetDefaultIntegrationMethod());
        
        Matrix B_derivative_x;
        Matrix B_derivative_y;
        CalculateBDerivativeDx(B_derivative_x, DDN_DDe);
        CalculateBDerivativeDy(B_derivative_y, DDN_DDe);


        // TODO: maybe need to do it with Taylor expansion??

        // Computed using the Gauss Points in the same knot span
        Vector divergence_of_sigma = this->GetValue(RECOVERED_STRESS);

        // initilize the div(sigma) matrix 2x(2n)
        Vector div_sigma_1 = ZeroVector(mDim*number_of_points);
        Vector div_sigma_2 = ZeroVector(mDim*number_of_points);

        Vector r_D_0(3); Vector r_D_1(3); Vector r_D_2(3);
        for (std::size_t i = 0; i < 3; ++i) {
            r_D_0(i) = rD(0, i); 
            r_D_1(i) = rD(1, i); 
            r_D_2(i) = rD(2, i);
        }
        div_sigma_1 = prod(trans(B_derivative_x), r_D_0) + prod(trans(B_derivative_y), r_D_2);
        div_sigma_2 = prod(trans(B_derivative_y), r_D_1) + prod(trans(B_derivative_x), r_D_2);
        
        // NEW FORMULATION___________________________________________________________________________________________
        Matrix div_sigma = ZeroMatrix(2, mDim * number_of_points);
        for (std::size_t i = 0; i < mDim * number_of_points; ++i) {
            div_sigma(0, i) = div_sigma_1(i);
            div_sigma(1, i) = div_sigma_2(i);
        }

        Vector DN_DX_vector_1 = ZeroVector(number_of_points);
        Vector DN_DX_vector_2 = ZeroVector(number_of_points);
        for (std::size_t i = 0; i < number_of_points; ++i) {
            DN_DX_vector_1(i) = rDN_DX(i, 0); 
            DN_DX_vector_2(i) = rDN_DX(i, 1); 
        }
        Matrix temp_matrix_grad_q_div_tau = TauOne * GaussWeight * (outer_prod(div_sigma_1, DN_DX_vector_1) + outer_prod(div_sigma_2, DN_DX_vector_2));

        for (IndexType i = 0; i < number_of_points; ++i)
        {
            for (IndexType j = 0; j < number_of_points; ++j)
            {
                for (IndexType dim2 = 0; dim2 < mDim; ++dim2) 
                {
                    // Assemble 2nd order stabilization term -> gradQ-div(sigma) // this one
                    rLeftHandSideMatrix(j * block_size + mDim, i * block_size + dim2) -= temp_matrix_grad_q_div_tau(i * mDim + dim2, j);
                }
            }
        }

        // --- RHS corresponding term ---
        Vector velocity_previous = ZeroVector(mDim*number_of_points);
        Vector pressure_previous = ZeroVector(number_of_points);
        IndexType index = 0;
        for (IndexType i = 0; i < number_of_points; ++i) { 
            velocity_previous[index++] = r_geometry[i].GetSolutionStepValue(VELOCITY_X);
            velocity_previous[index++] = r_geometry[i].GetSolutionStepValue(VELOCITY_Y);
            pressure_previous[i] = r_geometry[i].GetSolutionStepValue(PRESSURE);
        }

        // --- RHS corresponding term ---
        for (IndexType i = 0; i < number_of_points; ++i) {   
            for (IndexType dim1 = 0; dim1 < mDim; ++dim1) {
                // Assemble 2nd order stabilization term -> gradQ-div(sigma) 
                rRightHandSideVector(i * block_size + mDim) += GaussWeight * TauOne * rDN_DX(i,dim1) * divergence_of_sigma[dim1];
            }
        }
    } 
}

void GapSbmNavierStokesElement::AddConvectiveTerms(MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const array_1d<double,3> &rBodyForce,
        const double Density,
        const double TauOne,
        const ShapeFunctionsType &rN,
        const ShapeDerivativesType &rDN_DX,
        const double GaussWeight,
        const Matrix& rD)
{
    const unsigned int number_of_points = this->GetGeometry().PointsNumber();
    const unsigned int block_size = mDim+1;
    auto r_geometry = GetGeometry();
    
    // Compute advective velocity a = sum_j N_j * u_j
    array_1d<double, 3> advective_velocity = ZeroVector(3);
    for (unsigned int j = 0; j < number_of_points; ++j)
    {
        const array_1d<double, 3>& u_j = r_geometry[j].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int d = 0; d < mDim; ++d)
            advective_velocity[d] += rN[j] * u_j[d];
    }

    // Convective term: (a · ∇u) · v = (a · ∇φ_j) φ_i
    for (unsigned int i = 0; i < number_of_points; ++i)
    {
        unsigned int row_index = i * block_size;
        for (unsigned int j = 0; j < number_of_points; ++j)
        {
            unsigned int col_index = j * block_size;
            for (unsigned int d = 0; d < mDim; ++d)
            {
                double adv_grad_phi_j = 0.0;
                for (unsigned int k = 0; k < mDim; ++k)
                {
                    adv_grad_phi_j += advective_velocity[k] * rDN_DX(j, k); // (a · ∇φ_j^d)
                }
                // Add convective term [No Stabilization]
                rLeftHandSideMatrix(row_index + d, col_index + d) += Density * GaussWeight * rN[i] * adv_grad_phi_j;
            }
        }
    }

    // Compute grad(u^d) = sum_j ∂N_j/∂x_k * u_j^d
    array_1d<double, 3> adv_grad_u = ZeroVector(3); // one per velocity component
    for (unsigned int d = 0; d < mDim; ++d) // for each component of velocity
    {
        double component_adv_grad = 0.0;
        for (unsigned int j = 0; j < number_of_points; ++j)
        {
            const double u_jd = r_geometry[j].FastGetSolutionStepValue(VELOCITY)[d];
            for (unsigned int k = 0; k < mDim; ++k)
            {
                // a_k * ∂u_j^d / ∂x_k = a_k * u_j^d * ∂N_j/∂x_k
                component_adv_grad += advective_velocity[k] * u_jd * rDN_DX(j, k);
            }
        }
        adv_grad_u[d] = component_adv_grad;
    }
    // Assemble RHS: (a · grad u) * v = adv_grad_u[d] * N_i
    for (unsigned int i = 0; i < number_of_points; ++i)
    {
        unsigned int row_index = i * block_size;

        for (unsigned int d = 0; d < mDim; ++d)
        {
            // Corresponding RHS term: (a · ∇u) · v
            rRightHandSideVector[row_index + d] -= Density * GaussWeight * rN[i] * adv_grad_u[d];
        }
    }

}

void GapSbmNavierStokesElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    // MODIFY THIS??? r_surrogate



    auto r_surrogate_geometry = GetSurrogateGeometry();
    auto r_geometry = GetGeometry();

    const unsigned int number_of_control_points = r_surrogate_geometry.PointsNumber();    
    const SizeType num_dofs_per_node = number_of_control_points * (mDim + 1);

    // compute Taylor expansion contribution: H_sum_vec
    Vector N_sum_vec = ZeroVector(number_of_control_points);
    ComputeTaylorExpansionContribution(N_sum_vec);

    // compute Taylor expansion contribution: grad_H_sum
    Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(grad_N_sum_transposed);
    Matrix grad_N_sum = trans(grad_N_sum_transposed);

    // from r_geometry get integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    if(rMassMatrix.size1() != num_dofs_per_node)
        rMassMatrix.resize(num_dofs_per_node,num_dofs_per_node,false);
    noalias(rMassMatrix) = ZeroMatrix(num_dofs_per_node, num_dofs_per_node);
    
    const Properties& r_properties = GetProperties();
    double density = r_properties[DENSITY];
    if (density == 0.0) {
        KRATOS_ERROR << "Density must be non-zero. Provided: " << density << std::endl;
    }

    // Define mass matrix and previous velocity vector
    Matrix mass_matrix = ZeroMatrix(num_dofs_per_node, num_dofs_per_node);
    
    // Build mass matrix for velocity DOFs
    for (std::size_t i = 0; i < number_of_control_points; ++i) {
        for (std::size_t j = 0; j < number_of_control_points; ++j) {
            for (std::size_t dim = 0; dim < mDim; ++dim) {
                mass_matrix(i * 3 + dim, j * 3 + dim) += density * N_sum_vec[i] * N_sum_vec[j] * integration_points[0].Weight();
            }
        }
    }
    // Add transient term to LHS
    noalias(rMassMatrix) += mass_matrix;

    // Add the VMS term with du/dt
    for (unsigned int i = 0; i < number_of_control_points; ++i)
    {
        for (unsigned int j = 0; j < number_of_control_points; ++j)
        {
            for (unsigned int d = 0; d < mDim; ++d) {
                // Stabilization-time (grad q, u^n+1)/delta_t
                rMassMatrix(i*3+mDim, j*3 + d) += density * integration_points[0].Weight() * mTauOne * grad_N_sum(i, d) * N_sum_vec[j] ;
            }
        }
    }
}



void GapSbmNavierStokesElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetSurrogateGeometry();
        const std::size_t number_of_control_points = r_geometry.size();

        if (rResult.size() != mDim * number_of_control_points)
            rResult.resize(mDim * number_of_control_points, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * mDim;
            const auto& r_node = r_geometry[i];
            rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        }
    }

void GapSbmNavierStokesElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    const auto& r_geometry = GetSurrogateGeometry();
    const std::size_t number_of_control_points = r_geometry.size();

    rElementalDofList.resize(0);
    rElementalDofList.reserve(mDim * number_of_control_points);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        const auto& r_node = r_geometry[i];
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
    }
};




int GapSbmNavierStokesElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    // Verify that the constitutive law exists
    if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
    {
        KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
    }
    else
    {
        // Verify that the constitutive law has the correct dimension
        KRATOS_ERROR_IF_NOT(this->GetProperties().Has(THICKNESS))
            << "THICKNESS not provided for element " << this->Id() << std::endl;

        // Check strain size
        KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize() == 3)
            << "Wrong constitutive law used. This is a 2D element! Expected strain size is 3 (el id = ) "
            << this->Id() << std::endl;
    }

    // Intentionally left blank: Cut-SBM element bypasses default geometry/size checks.
    return 0;
}


Element::IntegrationMethod GapSbmNavierStokesElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}



void GapSbmNavierStokesElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetSurrogateGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->FinalizeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);

    //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        //TODO: build a CalculateOnIntegrationPoints method
        //--------------------------------------------------------------------------------------------
        const auto& r_surrogate_geometry = GetSurrogateGeometry();
        const std::size_t number_of_control_points = r_surrogate_geometry.size();
        const std::size_t mat_size = number_of_control_points * 2;

        Vector old_displacement(mat_size);
        GetSolutionCoefficientVector(old_displacement);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        Matrix grad_N_sum_transposed = ZeroMatrix(3, number_of_control_points);
        ComputeGradientTaylorExpansionContribution(grad_N_sum_transposed);
        Matrix grad_N_sum = trans(grad_N_sum_transposed);

        Matrix B_sum = ZeroMatrix(mDim,mat_size);
        CalculateB(B_sum, grad_N_sum);

        // obtain the tangent constitutive matrix at the true position
        ConstitutiveLaw::Parameters values_true(GetSurrogateGeometry(), GetProperties(), rCurrentProcessInfo);

        Vector old_displacement_coefficient_vector(mat_size);
        GetSolutionCoefficientVector(old_displacement_coefficient_vector);
        Vector old_strain_on_true = prod(B_sum, old_displacement_coefficient_vector);

        const std::size_t strain_size_true = mpConstitutiveLaw->GetStrainSize();
        ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
        ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true);

        const Vector sigma = values_true.GetStressVector();
        
        SetValue(CAUCHY_STRESS_XX, sigma[0]);
        SetValue(CAUCHY_STRESS_YY, sigma[1]);
        SetValue(CAUCHY_STRESS_XY, sigma[2]);
        // //---------------------
}

void GapSbmNavierStokesElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
    ConstitutiveLaw::Parameters constitutive_law_parameters(
        GetSurrogateGeometry(), GetProperties(), rCurrentProcessInfo);

    mpConstitutiveLaw->InitializeMaterialResponse(constitutive_law_parameters, ConstitutiveLaw::StressMeasure_Cauchy);
}



void GapSbmNavierStokesElement::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints();

    if (rOutput.size() != r_integration_points.size())
    {
        rOutput.resize(r_integration_points.size());
    }

    if (mpConstitutiveLaw->Has(rVariable)) {
        mpConstitutiveLaw->GetValue(rVariable, rOutput[0]);
    } else {
        KRATOS_WATCH(rVariable);
        KRATOS_WARNING("VARIABLE PRINT STILL NOT IMPLEMENTED IN THE IGA FRAMEWORK");
    }
}

void GapSbmNavierStokesElement::CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3 >>& rVariable,
        std::vector<array_1d<double, 3 >>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();

        if (rOutput.size() != r_integration_points.size())
        {
            rOutput.resize(r_integration_points.size());
        }

    if (mpConstitutiveLaw->Has(rVariable)) {
        mpConstitutiveLaw->GetValue(rVariable, rOutput[0]);
    } else {
            KRATOS_WATCH(rVariable);
            KRATOS_WARNING("VARIABLE PRINT STILL NOT IMPLEMENTED IN THE IGA FRAMEWORK");
    }
}

void GapSbmNavierStokesElement::CalculateB(
        Matrix& rB, 
        const Matrix& r_DN_DX) const
    {
        const auto& r_surrogate_geometry = GetSurrogateGeometry();
        const std::size_t number_of_control_points = r_surrogate_geometry.size();
        const std::size_t mat_size = number_of_control_points * 2;

        if (rB.size1() != 3 || rB.size2() != mat_size)
            rB.resize(3, mat_size);
        noalias(rB) = ZeroMatrix(3, mat_size);

        for (IndexType r = 0; r < mat_size; r++)
        {
            // local node number kr and dof direction dirr
            IndexType kr = r / 2;
            IndexType dirr = r % 2;

            rB(0, r) = r_DN_DX(kr,0) * (1-dirr);
            rB(1, r) = r_DN_DX(kr,1) * dirr;
            rB(2, r) = r_DN_DX(kr,0) * (dirr) + r_DN_DX(kr,1) * (1-dirr);
        }
    }

void GapSbmNavierStokesElement::GetSolutionCoefficientVector(
        Vector& rValues) const
    {
        const auto& r_surrogate_geometry = GetSurrogateGeometry();
        const std::size_t number_of_control_points = r_surrogate_geometry.size();
        const std::size_t mat_size = number_of_control_points * 2;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3>& displacement = r_surrogate_geometry[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * 2;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
        }
    }

void GapSbmNavierStokesElement::ApplyConstitutiveLaw(std::size_t matSize, Vector& rStrain, ConstitutiveLaw::Parameters& rValues,
                                        ConstitutiveVariables& rConstitutiVariables)
{
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=rValues.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    
    rValues.SetStrainVector(rStrain);
    rValues.SetStressVector(rConstitutiVariables.StressVector);
    rValues.SetConstitutiveMatrix(rConstitutiVariables.D);

    mpConstitutiveLaw->CalculateMaterialResponse(rValues, ConstitutiveLaw::StressMeasure_Cauchy); 
}

void GapSbmNavierStokesElement::ComputeTaylorExpansionContribution(Vector& H_sum_vec)
{
    const auto& r_geometry = GetSurrogateGeometry();
    const std::size_t number_of_control_points = r_geometry.PointsNumber();
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    if (H_sum_vec.size() != number_of_control_points)
    {
        H_sum_vec = ZeroVector(number_of_control_points);
    }

    // Compute all the derivatives of the basis functions involved
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n-1] = r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }
    
    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        // Reset for each node
        double H_taylor_term = 0.0; 

        if (mDim == 2) {
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                // Retrieve the appropriate derivative for the term
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n-1];
                for (IndexType k = 0; k <= n; k++) {
                    IndexType n_k = n - k;
                    double derivative = r_shape_function_derivatives(i,k); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term += ComputeTaylorTerm(derivative, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
            }
        } else {
            // 3D
            for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
                Matrix& r_shape_function_derivatives = shape_function_derivatives[n-1];
                
                int countDerivativeId = 0;
                // Loop over blocks of derivatives in x
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    // Loop over the possible derivatives in y
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {
                        
                        // derivatives in z
                        IndexType k_z = n - k_x - k_y;
                        double derivative = r_shape_function_derivatives(i,countDerivativeId); 

                        H_taylor_term += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y, mDistanceVector[2], k_z);
                        countDerivativeId++;
                    }
                }
            }
        }
        H_sum_vec(i) = H_taylor_term + r_N(0,i);
    }
}

void GapSbmNavierStokesElement::ComputeGradientTaylorExpansionContribution(Matrix& grad_H_sum)
{
    const auto& r_geometry = GetSurrogateGeometry();
    const std::size_t number_of_control_points = r_geometry.PointsNumber();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());

    // Compute all the derivatives of the basis functions involved
    std::vector<Matrix> shape_function_derivatives(mBasisFunctionsOrder);
    for (IndexType n = 1; n <= mBasisFunctionsOrder; n++) {
        shape_function_derivatives[n-1] = r_geometry.ShapeFunctionDerivatives(n, 0, this->GetIntegrationMethod());
    }

    if (grad_H_sum.size1() != 3 || grad_H_sum.size2() != number_of_control_points)
    {
        grad_H_sum.resize(3, number_of_control_points);
    }

    // Neumann (Taylor expansion of the gradient)
    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        // Reset for each control point
        double H_taylor_term_X = 0.0; 
        double H_taylor_term_Y = 0.0; 
        double H_taylor_term_Z = 0.0; 

        if (mDim == 2) {
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                // Retrieve the appropriate derivative for the term
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n-1];
                for (IndexType k = 0; k <= n-1; k++) {
                    IndexType n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_X += ComputeTaylorTerm(derivative, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
                for (IndexType k = 0; k <= n-1; k++) {
                    IndexType n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k+1); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_Y += ComputeTaylorTerm(derivative, mDistanceVector[0], n_k, mDistanceVector[1], k);
                }
            }
        } else {
            // 3D
            for (IndexType n = 2; n <= mBasisFunctionsOrder; n++) {
                Matrix& shapeFunctionDerivatives = shape_function_derivatives[n-1];
            
                IndexType countDerivativeId = 0;
                // Loop over blocks of derivatives in x
                for (IndexType k_x = n; k_x >= 0; k_x--) {
                    // Loop over the possible derivatives in y
                    for (IndexType k_y = n - k_x; k_y >= 0; k_y--) {

                        // derivatives in z
                        IndexType k_z = n - k_x - k_y;
                        double derivative = shapeFunctionDerivatives(i,countDerivativeId); 
                        
                        if (k_x >= 1) {
                            H_taylor_term_X += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x-1, mDistanceVector[1], k_y, mDistanceVector[2], k_z);
                        }
                        if (k_y >= 1) {
                            H_taylor_term_Y += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y-1, mDistanceVector[2], k_z);
                        }
                        if (k_z >= 1) {
                            H_taylor_term_Z += ComputeTaylorTerm3D(derivative, mDistanceVector[0], k_x, mDistanceVector[1], k_y, mDistanceVector[2], k_z-1);
                        }     
                        countDerivativeId++;
                    }
                }
            }
        }
        grad_H_sum(0,i) = H_taylor_term_X + r_DN_De[0](i, 0);
        grad_H_sum(1,i) = H_taylor_term_Y + r_DN_De[0](i, 1);
        if (mDim == 3)
            grad_H_sum(2,i) = H_taylor_term_Z + r_DN_De[0](i, 2);
        else 
            grad_H_sum(2,i) = 0;
    }    
}

// Function to compute a single term in the Taylor expansion
double GapSbmNavierStokesElement::ComputeTaylorTerm(
    const double derivative, 
    const double dx, 
    const IndexType n_k, 
    const double dy, 
    const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));    
}

double GapSbmNavierStokesElement::ComputeTaylorTerm3D(
    const double derivative, 
    const double dx, 
    const IndexType k_x, 
    const double dy, 
    const IndexType k_y, 
    const double dz, 
    const IndexType k_z)
{   
    return derivative * std::pow(dx, k_x) * std::pow(dy, k_y) * std::pow(dz, k_z) / (MathUtils<double>::Factorial(k_x) * MathUtils<double>::Factorial(k_y) * MathUtils<double>::Factorial(k_z));    
}

void GapSbmNavierStokesElement::CalculateBDerivativeDx(
        Matrix& BDerivativeDx, 
        const ShapeDerivativesType& r_DDN_DDX) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * mDim; // Only 2 DOFs per node in 2D

    // Resize B matrix to 3 rows (strain vector size) and appropriate number of columns
    if (BDerivativeDx.size1() != 3 || BDerivativeDx.size2() != mat_size)
        BDerivativeDx.resize(3, mat_size);

    noalias(BDerivativeDx) = ZeroMatrix(3, mat_size);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        BDerivativeDx(0, 2 * i)     = r_DDN_DDX(i, 0); // ∂²N_i / ∂x²
        BDerivativeDx(1, 2 * i + 1) = r_DDN_DDX(i, 1); // ∂²N_i / ∂y∂x
        BDerivativeDx(2, 2 * i)     = r_DDN_DDX(i, 1); // ∂²N_i / ∂y∂x
        BDerivativeDx(2, 2 * i + 1) = r_DDN_DDX(i, 0); // ∂²N_i / ∂x²
    }
}

void GapSbmNavierStokesElement::CalculateBDerivativeDy(
        Matrix& BDerivativeDy, 
        const ShapeDerivativesType& r_DDN_DDX) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 2; // Only 2 DOFs per node in 2D

    // Resize B matrix to 3 rows (strain vector size) and appropriate number of columns
    if (BDerivativeDy.size1() != 3 || BDerivativeDy.size2() != mat_size)
        BDerivativeDy.resize(3, mat_size);

    noalias(BDerivativeDy) = ZeroMatrix(3, mat_size);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        BDerivativeDy(0, 2 * i)     = r_DDN_DDX(i, 1); // ∂²N_i / ∂x∂y
        BDerivativeDy(1, 2 * i + 1) = r_DDN_DDX(i, 2); // ∂²N_i / ∂y²
        BDerivativeDy(2, 2 * i)     = r_DDN_DDX(i, 2); // ∂²N_i / ∂y²
        BDerivativeDy(2, 2 * i + 1) = r_DDN_DDX(i, 1); // ∂²N_i / ∂x∂y
    }
}

} // Namespace Kratos
