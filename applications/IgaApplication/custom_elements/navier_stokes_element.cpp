//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"

// Application includes
#include "custom_elements/navier_stokes_element.h"

namespace Kratos
{

NavierStokesElement::NavierStokesElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(
        NewId,
        pGeometry)
{
}

NavierStokesElement::NavierStokesElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(
        NewId,
        pGeometry,
        pProperties)
{
}

Element::Pointer NavierStokesElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<NavierStokesElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer NavierStokesElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<NavierStokesElement>(NewId, pGeom, pProperties);
}

// Deconstructor

NavierStokesElement::~NavierStokesElement()
{
}

void NavierStokesElement:: Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
    GeometryType& r_geometry = this->GetGeometry();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    mDim = DN_De[0].size2();

    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
    
    // SetValue the integraion_weigth in order to do the postprocess in Python
    SetValue(INTEGRATION_WEIGHT, integration_points[0].Weight());
}


void NavierStokesElement::InitializeMaterial()
{
    KRATOS_TRY
    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLaw->InitializeMaterial( r_properties, r_geometry, row(N_values , 0 ));

    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH( "" );
}


void NavierStokesElement::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo)
{
    // Obtain required constants
    GeometryType& r_geometry = this->GetGeometry();
    const unsigned int number_of_points = r_geometry.PointsNumber();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    SizeType num_dofs_per_node = number_of_points * (mDim + 1);
    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != num_dofs_per_node)
        rLeftHandSideMatrix.resize(num_dofs_per_node,num_dofs_per_node,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(num_dofs_per_node,num_dofs_per_node); //resetting LHS
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != num_dofs_per_node)
        rRightHandSideVector.resize(num_dofs_per_node,false);
    noalias(rRightHandSideVector) = ZeroVector(num_dofs_per_node); //resetting RHS

    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(DN_De[0].size1()) - 1;
    }
    const ShapeFunctionsType& rN = row(N_gausspoint,0);
    const ShapeDerivativesType& DN_DX = DN_De[0];
    const double GaussWeight = integration_points[0].Weight();

    // Calculate the B matrix
    Matrix B = ZeroMatrix(3,number_of_points*mDim);
    CalculateB(B, DN_DX);

    // constitutive law
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables(3);
    ApplyConstitutiveLaw(B, Values, constitutive_variables);
    Vector& r_stress_vector = Values.GetStressVector();
    const Matrix& r_D = Values.GetConstitutiveMatrix();
    Matrix DB_voigt = Matrix(prod(r_D, B));

    double mu_effective;
    mpConstitutiveLaw->CalculateValue(Values, MU, mu_effective);

    if (mu_effective == 0) {
        // if the mu_effective is zero, we use the viscosity value
        this->EvaluateInPoint(mu_effective,VISCOSITY,rN,r_geometry);
    }

    array_1d<double,3> body_force = ZeroVector(3);
    body_force = this->GetValue(BODY_FORCE);
    double viscosity = mu_effective;
    const Properties& r_properties = GetProperties();
    double density = r_properties[DENSITY];

    double adv_norm;
    // Calculate the advective norm
    CalculateAdvectiveNorm(rN, adv_norm);

    // Calculate stabilization constants
    CalculateTau(mu_effective, density, adv_norm);

    // Add velocity terms in momentum equation
    this->AddMomentumTerms(rLeftHandSideMatrix,rRightHandSideVector,body_force,mTauTwo,rN,DN_DX,GaussWeight, r_D, r_stress_vector);

    // Add velocity-pressure terms
    this->AddContinuityTerms(rLeftHandSideMatrix,rRightHandSideVector,body_force,density,mTauOne,rN,DN_DX,GaussWeight);

    // Add Second-Order stabilization terms from VMS
    this->AddSecondOrderStabilizationTerms(rLeftHandSideMatrix,rRightHandSideVector,mTauOne,rN,DN_DX,GaussWeight,r_D);

    // Add Convective terms [D terms] 
    this->AddConvectiveTerms(rLeftHandSideMatrix,rRightHandSideVector,body_force,density,mTauOne,rN,DN_DX,GaussWeight,r_D);

}

void NavierStokesElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{   
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

void NavierStokesElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

void NavierStokesElement::CalculateAdvectiveNorm(
    const ShapeFunctionsType &rN,
    double& adv_norm)   
{
    auto r_geometry = this->GetGeometry();
    const unsigned int number_of_points = r_geometry.PointsNumber();

    // Compute advective velocity: a = sum_j N_j * u_j
    array_1d<double, 3> advective_velocity = ZeroVector(3);
    for (unsigned int j = 0; j < number_of_points; ++j)
    {
        const array_1d<double, 3>& u_j = r_geometry[j].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int d = 0; d < mDim; ++d)
            advective_velocity[d] += rN[j] * u_j[d];
    }

    // Compute ||a||
    adv_norm = norm_2(advective_velocity);
}

void NavierStokesElement::CalculateTau(
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

double NavierStokesElement::ElementSize()
{
    return 1.128379167 * sqrt(this->GetGeometry().DomainSize()); 
}

void NavierStokesElement::AddMomentumTerms(MatrixType &rLHS,
        VectorType &rRHS,
        const array_1d<double,3> &rBodyForce,
        const double TauTwo,
        const ShapeFunctionsType &rN,
        const ShapeDerivativesType &rDN_DX,
        const double Weight,
        const Matrix& rD,
        Vector& rStressVector)
{
    const unsigned int number_of_nodes = this->GetGeometry().PointsNumber();
    const unsigned int block_size = mDim+1;
    auto r_geometry = GetGeometry();

    Matrix B = ZeroMatrix(3,number_of_nodes*mDim);
    CalculateB(B, rDN_DX);
    // Dynamic viscosity already included in rD
    Matrix diffusion_term_matrix = Weight * prod(trans(B), Matrix(prod(rD, B)));
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
    Vector internal_forces = Weight * prod(trans(B), rStressVector);
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


void NavierStokesElement::AddContinuityTerms(MatrixType &rLHS,
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


    // // TO BE DELETED
    // // === Penalty per media-zero "morbida" ===
    // const double alpha_p = 1e2;

    // if (alpha_p > 0.0) {
    //     // contribuisce a LHS (p,p)
    //     for (unsigned int i = 0; i < number_of_nodes; ++i) {
    //         for (unsigned int j = 0; j < number_of_nodes; ++j) {
    //             rLHS(i*block_size + mDim, j*block_size + mDim) += alpha_p * Weight * rN[i] * rN[j];
    //         }
    //     }

    //     // contribuisce anche al RHS (residuo) con p alla iterazione corrente (o "previous" come fai tu)
    //     for (unsigned int i = 0; i < number_of_nodes; ++i) {
    //         rRHS(i*block_size + mDim) -= alpha_p * Weight * rN[i] * pressure_previous_iteration;
    //         // Se preferisci essere totalmente "consistenti" con l’unknown attuale, valuta p_i dei nodi direttamente:
    //         // double p_i = r_geometry[i].FastGetSolutionStepValue(PRESSURE);
    //         // rRHS(i*block_size + mDim) -= alpha_p * Weight * rN[i] * p_i;
    //     }
    // }



}


void NavierStokesElement::AddSecondOrderStabilizationTerms(MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const double TauOne,
        const ShapeFunctionsType &rN,
        const ShapeDerivativesType &rDN_DX,
        const double GaussWeight,
        const Matrix& rD)
{
    const unsigned int number_of_points = this->GetGeometry().PointsNumber();
    const unsigned int block_size = mDim+1;
    auto r_geometry = GetGeometry();

    // Add third order term of the stabilization
    if (mBasisFunctionsOrder > 1) {

        Matrix B = ZeroMatrix(3,number_of_points*mDim);
        CalculateB(B, rDN_DX);
        const Matrix& DDN_DDe = GetGeometry().ShapeFunctionDerivatives(2, 0, GetGeometry().GetDefaultIntegrationMethod());
        
        Matrix B_derivative_x;
        Matrix B_derivative_y;
        CalculateBDerivativeDx(B_derivative_x, DDN_DDe);
        CalculateBDerivativeDy(B_derivative_y, DDN_DDe);

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

void NavierStokesElement::AddConvectiveTerms(MatrixType &rLeftHandSideMatrix,
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

    // STABILIZATION TERMS [D1, D2, D3, D5]

    // // [D1] Stabilization: grad(v) · a * f => (a · ∇v_i) * f
    // for (unsigned int i = 0; i < number_of_points; ++i)
    // {
    //     unsigned int row_index = i * block_size;
    //     for (unsigned int d = 0; d < mDim; ++d) // for each velocity component
    //     {
    //         double a_dot_grad_vi = 0.0;
    //         for (unsigned int k = 0; k < mDim; ++k)
    //         {
    //             a_dot_grad_vi += advective_velocity[k] * rDN_DX(i, k);
    //         }
    //         rRightHandSideVector[row_index + d] += GaussWeight * TauOne * a_dot_grad_vi * rBodyForce[d];
    //     }
    // }
    // // [D2] Stabilization: grad(v) · a * div(sigma) → τ₁ (a · ∇v_i) · div(σ)
    // if (mBasisFunctionsOrder > 1) {
    //     Vector divergence_of_sigma = this->GetValue(RECOVERED_STRESS);
    //     for (IndexType i = 0; i < number_of_points; ++i)
    //     {
    //         unsigned int row_index = i * block_size;
    //         double a_dot_grad_vi = 0.0;
    //         for (IndexType dim1 = 0; dim1 < mDim; ++dim1)
    //         {
    //             a_dot_grad_vi += advective_velocity[dim1] * rDN_DX(i, dim1);
    //         }
    //         for (IndexType dim2 = 0; dim2 < mDim; ++dim2)
    //         {
    //             // For this term we omit the LHS contribution
    //             rRightHandSideVector(row_index + dim2) += GaussWeight * TauOne * a_dot_grad_vi * divergence_of_sigma[dim2];
    //         }
    //     }
    // }
    // // [D3] Stabilization: grad(v) · a * grad(p) → τ₁ (a · ∇v_i) · grad(p)_j
    // for (IndexType i = 0; i < number_of_points; ++i)
    // {
    //     unsigned int row_index = i * block_size;
    //     // Compute directional derivative: (a · ∇phi_i)
    //     double a_dot_grad_vi = 0.0;
    //     for (IndexType d = 0; d < mDim; ++d)
    //     {
    //         a_dot_grad_vi += advective_velocity[d] * rDN_DX(i, d);
    //     }

    //     for (IndexType j = 0; j < number_of_points; ++j)
    //     {
    //         unsigned int col_index = j * block_size;
    //         for (IndexType d = 0; d < mDim; ++d)
    //         {
    //             // Assemble into LHS
    //             rLeftHandSideMatrix(row_index + d, col_index + mDim) += GaussWeight * TauOne * a_dot_grad_vi * rDN_DX(j, d);
    //         }
    //     }
    // }
    // // Corresponding RHS term: τ₁ (a · ∇v_i) · grad(p)_j
    // array_1d<double, 3> grad_p = ZeroVector(3);
    // for (IndexType j = 0; j < number_of_points; ++j)
    // {
    //     double p_j = r_geometry[j].FastGetSolutionStepValue(PRESSURE);
    //     for (IndexType d = 0; d < mDim; ++d)
    //     {
    //         grad_p[d] += p_j * rDN_DX(j, d); // ∇p ≈ ∑_j p_j ∇φ_j
    //     }
    // }
    // for (IndexType i = 0; i < number_of_points; ++i)
    // {
    //     unsigned int row_index = i * block_size;
    //     double a_dot_grad_vi = 0.0;
    //     for (IndexType d = 0; d < mDim; ++d)
    //     {
    //         a_dot_grad_vi += advective_velocity[d] * rDN_DX(i, d); // a · ∇φ_i
    //     }

    //     for (IndexType d = 0; d < mDim; ++d)
    //     {
    //         rRightHandSideVector[row_index + d] -= GaussWeight * TauOne * a_dot_grad_vi * grad_p[d];
    //     }
    // }

    // // [D5] Stabilization: (a · grad v) * (a · grad u) → τ₁ (a · ∇φ_i) * (a · ∇φ_j)
    // for (IndexType i = 0; i < number_of_points; ++i)
    // {
    //     unsigned int row_index = i * block_size;
    //     double a_dot_grad_vi = 0.0;
    //     for (IndexType d = 0; d < mDim; ++d)
    //         a_dot_grad_vi += advective_velocity[d] * rDN_DX(i, d);
    //     for (IndexType j = 0; j < number_of_points; ++j)
    //     {
    //         unsigned int col_index = j * block_size;
    //         double a_dot_grad_phi_j = 0.0;
    //         for (IndexType d = 0; d < mDim; ++d)
    //             a_dot_grad_phi_j += advective_velocity[d] * rDN_DX(j, d);
    //         for (IndexType d = 0; d < mDim; ++d)
    //         {
    //             rLeftHandSideMatrix(row_index + d, col_index + d) += GaussWeight * TauOne * a_dot_grad_vi * a_dot_grad_vi;
    //         }
    //     }
    // }
    // // Corresponding RHS term
    // array_1d<double, 3> a_dot_grad_u = ZeroVector(3); // Compute (a · grad u) at the Gauss point
    // for (IndexType d = 0; d < mDim; ++d)
    // {
    //     double component = 0.0;
    //     for (IndexType j = 0; j < number_of_points; ++j)
    //     {
    //         double u_jd = r_geometry[j].FastGetSolutionStepValue(VELOCITY)[d];
    //         for (IndexType k = 0; k < mDim; ++k)
    //         {
    //             component += advective_velocity[k] * u_jd * rDN_DX(j, k); // a_k * ∂u_j^d/∂x_k
    //         }
    //     }
    //     a_dot_grad_u[d] = component;
    // }
    // // Assemble D5 RHS: - τ₁ * (a · ∇φ_i) * (a · ∇u)
    // for (IndexType i = 0; i < number_of_points; ++i)
    // {
    //     unsigned int row_index = i * block_size;
    //     double a_dot_grad_phi_i = 0.0;
    //     for (IndexType d = 0; d < mDim; ++d)
    //         a_dot_grad_phi_i += advective_velocity[d] * rDN_DX(i, d);

    //     for (IndexType d = 0; d < mDim; ++d)
    //     {
    //         rRightHandSideVector[row_index + d] -= GaussWeight * TauOne * a_dot_grad_phi_i * a_dot_grad_u[d];
    //     }
    // }
}

void NavierStokesElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    auto r_geometry = GetGeometry();
    const unsigned int NumNodes = r_geometry.PointsNumber();
    const unsigned int BlockSize = mDim+1;
    
    const SizeType num_dofs_per_node = NumNodes * (mDim + 1);
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    const ShapeFunctionsType& N = row(N_gausspoint,0);
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const ShapeDerivativesType& DN_DX = DN_De[0];
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
    for (std::size_t i = 0; i < NumNodes; ++i) {
        for (std::size_t j = 0; j < NumNodes; ++j) {
            for (std::size_t dim = 0; dim < mDim; ++dim) {
                mass_matrix(i * BlockSize + dim, j * BlockSize + dim) += density * N[i] * N[j] * integration_points[0].Weight();
            }
        }
    }
    // Add transient term to LHS
    noalias(rMassMatrix) += mass_matrix;

    // Add the VMS term with du/dt
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        for (unsigned int j = 0; j < NumNodes; ++j)
        {
            for (unsigned int d = 0; d < mDim; ++d) {
                // Stabilization-time (grad q, u^n+1)/delta_t
                rMassMatrix(i*BlockSize+mDim, j*BlockSize + d) += density * integration_points[0].Weight() * mTauOne * DN_DX(i, d) * N[j] ;
            }
        }
    }
}

void NavierStokesElement::GetFirstDerivativesVector(Vector &rValues, int Step) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = (mDim + 1) * NumNodes;

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(PRESSURE,Step);
    }
}

void NavierStokesElement::ApplyConstitutiveLaw(
        const Matrix& rB, 
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rConstitutiveVariables) const
{
    const SizeType number_of_nodes = GetGeometry().size();

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=rValues.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    Vector old_displacement(number_of_nodes*mDim);
    GetSolutionCoefficientVector(old_displacement);
    Vector old_strain = prod(rB,old_displacement);
    rValues.SetStrainVector(old_strain);
    rValues.SetStressVector(rConstitutiveVariables.StressVector);
    rValues.SetConstitutiveMatrix(rConstitutiveVariables.D);

    mpConstitutiveLaw->CalculateMaterialResponseCauchy(rValues);
}

void NavierStokesElement::EquationIdVector(Element::EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const int dim = 2;
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType number_of_control_points = GetGeometry().size();
    const unsigned int local_size = (dim + 1) * number_of_control_points;

    if (rResult.size() != local_size)
        rResult.resize(local_size);

    unsigned int Index = 0;
    for (unsigned int i = 0; i < number_of_control_points; i++)
    {
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
        if (dim > 2) rResult[Index++] = rGeom[i].GetDof(VELOCITY_Z).EquationId();
        rResult[Index++] = rGeom[i].GetDof(PRESSURE).EquationId();
    }
}


void NavierStokesElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    KRATOS_TRY;

    const SizeType number_of_control_points = GetGeometry().size();

    rElementalDofList.resize(0);
    rElementalDofList.reserve(3 * number_of_control_points);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(PRESSURE));
    }

    KRATOS_CATCH("")
};



int NavierStokesElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    int Result = Element::Check(rCurrentProcessInfo);

    // Checks on nodes
    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
    {
        if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
                this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
                ( this->GetGeometry().WorkingSpaceDimension() == 3 && this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false ) )
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if (this->GetGeometry().WorkingSpaceDimension() == 2)
    {
        for (unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if (this->GetGeometry()[i].Z() != 0.0)
                KRATOS_THROW_ERROR(std::invalid_argument,"Node with non-zero Z coordinate found. Id: ",this->GetGeometry()[i].Id());
        }
    }

    return Result;
}

Element::IntegrationMethod NavierStokesElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

void NavierStokesElement::CalculateB(
        Matrix& rB, 
        const ShapeDerivativesType& r_DN_DX) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 2; // Only 2 DOFs per node in 2D

    // Resize B matrix to 3 rows (strain vector size) and appropriate number of columns
    if (rB.size1() != 3 || rB.size2() != mat_size)
        rB.resize(3, mat_size);

    noalias(rB) = ZeroMatrix(3, mat_size);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        // x-derivatives of shape functions -> relates to strain component ε_11 (xx component)
        rB(0, 2 * i)     = r_DN_DX(i, 0); // ∂N_i / ∂x
        // y-derivatives of shape functions -> relates to strain component ε_22 (yy component)
        rB(1, 2 * i + 1) = r_DN_DX(i, 1); // ∂N_i / ∂y
        // Symmetric shear strain component ε_12 (xy component)
        rB(2, 2 * i)     = r_DN_DX(i, 1); // ∂N_i / ∂y
        rB(2, 2 * i + 1) = r_DN_DX(i, 0); // ∂N_i / ∂x
    }
}

void NavierStokesElement::CalculateBDerivativeDx(
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

void NavierStokesElement::CalculateBDerivativeDy(
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


void NavierStokesElement::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == CAUCHY_STRESS_VECTOR) {
        const unsigned int num_integration_points = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());

        if (rOutput.size() != num_integration_points) {
            rOutput.resize(num_integration_points);
        }

        Vector sigma_voigt(3);
        sigma_voigt = this->CalculateStressAtIntegrationPoint(rCurrentProcessInfo);
        rOutput[0] = sigma_voigt;
    }
}


Vector NavierStokesElement::CalculateStressAtIntegrationPoint(
    const ProcessInfo& rCurrentProcessInfo)
{
    Vector stress_voigt(3);
    
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_points = r_geometry.size();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const ShapeDerivativesType& DN_DX = DN_De[0];

    Matrix B = ZeroMatrix(3,number_of_points*mDim);
    CalculateB(B, DN_DX);

    // constitutive law
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables(3);
    ApplyConstitutiveLaw(B, Values, constitutive_variables);
    const Matrix& r_D = Values.GetConstitutiveMatrix();
    Matrix DB_voigt = Matrix(prod(r_D, B));
    stress_voigt = Values.GetStressVector();
    
    return stress_voigt;
}


void NavierStokesElement::GetSolutionCoefficientVector(
        Vector& rValues) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 2;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        const array_1d<double, 3 >& velocity = GetGeometry()[i].GetSolutionStepValue(VELOCITY);
        IndexType index = i * 2;

        rValues[index] = velocity[0];
        rValues[index + 1] = velocity[1];
    }
}

} // Namespace Kratos