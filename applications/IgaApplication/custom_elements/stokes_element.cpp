//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//  Main authors:    Nicolò Antonelli
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "custom_elements/stokes_element.h"

#include "utilities/math_utils.h"

namespace Kratos
{

StokesElement::StokesElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(
        NewId,
        pGeometry)
{
}

StokesElement::StokesElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(
        NewId,
        pGeometry,
        pProperties)
{
}

Element::Pointer StokesElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<StokesElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer StokesElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<StokesElement>(NewId, pGeom, pProperties);
}

// Deconstructor

StokesElement::~StokesElement()
{
}

void StokesElement:: Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
    GeometryType& r_geometry = this->GetGeometry();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    mDim = DN_De[0].size2();
}


void StokesElement::InitializeMaterial()
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


void StokesElement::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo)
{
    // Obtain required constants
    GeometryType& r_geometry = this->GetGeometry();
    const unsigned int number_of_points = r_geometry.PointsNumber();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    // mDim = DN_De[0].size2();
    const unsigned int BlockSize = mDim + 1;
    const unsigned int LocalSize = BlockSize * number_of_points;

    const unsigned int NumGauss = r_geometry.IntegrationPointsNumber(); // = 1 Gauss Point per element

    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    const int num_dofs_per_node = number_of_points * (mDim + 1);
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

    const ShapeFunctionsType& N = row(N_gausspoint,0);
    const ShapeDerivativesType& DN_DX = DN_De[0];
    const double GaussWeight = integration_points[0].Weight();

    // constitutive law
    Matrix B = ZeroMatrix(3,number_of_points*mDim);
    CalculateB(B, DN_DX);

    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    ConstitutiveVariables this_constitutive_variables(strain_size);
    Vector old_displacement(number_of_points*mDim);
    GetValuesVector(old_displacement);
    Vector old_strain = prod(B,old_displacement);
    // Values.SetStrainVector(this_constitutive_variables.StrainVector);
    Values.SetStrainVector(old_strain); // this is the input parameter
    Values.SetStressVector(this_constitutive_variables.StressVector); // this is an ouput parameter
    Values.SetConstitutiveMatrix(this_constitutive_variables.D);

    mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);

    double mu_effective;
    mpConstitutiveLaw->CalculateValue(Values, MU, mu_effective);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
    //this is ok under the hypothesis that no history dependent behavior is employed
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);

    Vector& r_stress_vector = Values.GetStressVector();
    const Matrix& r_D = Values.GetConstitutiveMatrix();

    double Density;
    double Viscosity;
    array_1d<double,3> BodyForce = ZeroVector(3);
    array_1d<double,3> Velocity = ZeroVector(3);

    // Interpolation
    this->EvaluateInPoint(Density,DENSITY,N,r_geometry);
    this->EvaluateInPoint(Viscosity,VISCOSITY,N,r_geometry);
    BodyForce = this->GetValue(BODY_FORCE);

    this->EvaluateInPoint(Velocity,VELOCITY,N,r_geometry);

    double TauOne = 0.0;
    double TauTwo = 0.0;
    double h_element_size1 = norm_2(r_geometry[0].Coordinates()-r_geometry[1].Coordinates());
    double h_element_size2 = norm_2(r_geometry[0].Coordinates()-r_geometry[mBasisFunctionsOrder+2].Coordinates());
    const double h_element_size = std::min(h_element_size1, h_element_size2);

    CalculateTau(TauOne,TauTwo,Density*Viscosity,h_element_size, mu_effective);

    // Add velocity terms in momentum equation
    this->AddMomentumTerms(rLeftHandSideMatrix,rRightHandSideVector,Density,Viscosity,BodyForce,TauTwo,N,DN_DX,GaussWeight, r_D, r_stress_vector);

    // Add velocity-pressure terms
    this->AddContinuityTerms(rLeftHandSideMatrix,rRightHandSideVector,Density,BodyForce,TauOne,N,DN_DX,GaussWeight);

    // Add Second-Order stabilization terms from VMS
    this->AddSecondOrderStabilizationTerms(rLeftHandSideMatrix,rRightHandSideVector,Density,BodyForce,TauOne,N,DN_DX,GaussWeight,r_D, r_stress_vector);

    // Vector divergence_of_sigma = this->GetValue(RECOVERED_STRESS);
    // if (r_geometry.Center().Y() < 0.1 && r_geometry.Center().Y() > 0.0) {
    //     if (r_geometry.Center().X() > 0.5) { 
    //         KRATOS_WATCH(r_geometry.Center().Y())
    //         KRATOS_WATCH(divergence_of_sigma)
    //     }
    // }

    // for (IndexType i = 0; i < number_of_points; ++i)
    // {
    //     for (IndexType j = 0; j < number_of_points; ++j)
    //     {
    //         // Add only to velocity DOFs (i.e., positions 0, 1, 3, 4, 6, 7, ...)
    //         for (IndexType dim1 = 0; dim1 < mDim; ++dim1)
    //         {
    //             rLeftHandSideMatrix(i * BlockSize + 2, j * BlockSize + dim1) = 0.0;
    //         }
    //     }
    //      rLeftHandSideMatrix(i * BlockSize + 2, i * BlockSize + 2) = 1.0;
    // }

    // // Add residual of previous iteration to RHS ---> TRAMPA
    // VectorType temp = ZeroVector(LocalSize);
    // // RHS = ExtForces - K*temp;
    // unsigned int index = 0 ;
    // for (unsigned int i = 0; i < number_of_points; i++) {
    //     temp[index++] = r_geometry[i].GetSolutionStepValue(VELOCITY_X);
    //     temp[index++] = r_geometry[i].GetSolutionStepValue(VELOCITY_Y);
    //     temp[index++] = r_geometry[i].GetSolutionStepValue(PRESSURE);
    // }
    // noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);
}

void StokesElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{   
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

void StokesElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}



void StokesElement::CalculateTau(double &TauOne, double &TauTwo, const double DynViscosity, const double h_element_size, double mu_effective)
{   
    if (mu_effective == 0) {
        KRATOS_ERROR << "mu_effective = 0" << std::endl;
    }
    // // Estimate element size
    double h = ElementSize();

    TauOne = std::pow(h, 2) / ( 4.0 * mu_effective ); 

    TauTwo = mu_effective;
}

double StokesElement::ElementSize()
{
    return 1.128379167 * sqrt(this->GetGeometry().DomainSize()); 
}

void StokesElement::AddMomentumTerms(MatrixType &rLHS,
        VectorType &rRHS,
        const double Density,
        const double Viscosity,
        const array_1d<double,3> &BodyForce,
        const double TauTwo,
        const ShapeFunctionsType &N,
        const ShapeDerivativesType &DN_DX,
        const double Weight,
        const Matrix& r_D,
        Vector& r_stress_vector)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int BlockSize = mDim+1;
    auto r_geometry = GetGeometry();

    Matrix B = ZeroMatrix(3,NumNodes*mDim);
    CalculateB(B, DN_DX);
    Matrix tempMatrix = Weight * prod(trans(B), Matrix(prod(r_D, B)));
    for (IndexType i = 0; i < NumNodes; ++i)
    {
        for (IndexType j = 0; j < NumNodes; ++j)
        {
            // Add only to velocity DOFs (i.e., positions 0, 1, 3, 4, 6, 7, ...)
            for (IndexType dim1 = 0; dim1 < mDim; ++dim1)
            {
                for (IndexType dim2 = 0; dim2 < mDim; ++dim2)
                {
                    // diffusion term
                    rLHS(i * BlockSize + dim1, j * BlockSize + dim2) += tempMatrix(i * mDim + dim1, j * mDim + dim2);
                }
            }
        }
    }
    // RHS corresponding term
    Vector internalForces = Weight * prod(trans(B), r_stress_vector);
    for (IndexType i = 0; i < NumNodes; ++i)
    {
        // Add only to the velocity DOFs (i.e., positions 0, 1, 3, 4, 6, 7, ...)
        for (IndexType dim1 = 0; dim1 < mDim; ++dim1)
        {
            // viscous + pressure term
            rRHS(i * BlockSize + dim1) -= internalForces(i * mDim + dim1);
        }  
    }

    unsigned int FirstRow = 0;
    unsigned int FirstCol = 0;
    double Term_ij = 0.0;

    for(unsigned int i = 0; i < NumNodes; ++i)
    {
        // Body force
        for(unsigned int d = 0; d < mDim; ++d) {
            rRHS[FirstRow + d] += Weight * N[i] * BodyForce[d];
        }

        for(unsigned int j = 0; j < NumNodes; ++j)
        {
            // Stabilization
            for (unsigned int m = 0; m < mDim; ++m) {
                for (unsigned int n = 0; n < mDim; ++n) {
                    rLHS(FirstRow+m,FirstCol+n) += TauTwo * Weight * DN_DX(i,m) * DN_DX(j,n);
                }
            }

            // Update column index
            FirstCol += BlockSize;
        }
        // RHS corresponding term
        double div_u_current = 0.0;
        for(unsigned int j = 0; j < NumNodes; ++j) {
            div_u_current += r_geometry[j].FastGetSolutionStepValue(VELOCITY_X) * DN_DX(j, 0) + 
                             r_geometry[j].FastGetSolutionStepValue(VELOCITY_Y) * DN_DX(j, 1) ;
        } 
        for (unsigned int m = 0; m < mDim; ++m) {
            rRHS(FirstRow+m) -= TauTwo * Weight * DN_DX(i,m) * div_u_current ;
        }

        // Update matrix indices
        FirstRow += BlockSize;
        FirstCol = 0;
    }
}


void StokesElement::AddContinuityTerms(MatrixType &rLHS,
        VectorType &rRHS,
        const double Density,
        const array_1d<double,3> &BodyForce,
        const double TauOne,
        const ShapeFunctionsType &N,
        const ShapeDerivativesType &DN_DX,
        const double Weight)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int BlockSize = mDim+1;
    auto r_geometry = GetGeometry();

    unsigned int FirstRow = 0;
    unsigned int FirstCol = 0;

    double pressure_previous_iteration = 0.0;
    double div_u_current = 0.0;
    Vector grad_p_previous_iteration = ZeroVector(2);
    for(unsigned int j = 0; j < NumNodes; ++j) {
        div_u_current += r_geometry[j].FastGetSolutionStepValue(VELOCITY_X) * DN_DX(j, 0) +
                        r_geometry[j].FastGetSolutionStepValue(VELOCITY_Y) * DN_DX(j, 1) ;
        grad_p_previous_iteration[0] += r_geometry[j].GetSolutionStepValue(PRESSURE) * DN_DX(j,0) ;
        grad_p_previous_iteration[1] += r_geometry[j].GetSolutionStepValue(PRESSURE) * DN_DX(j,1) ;
        pressure_previous_iteration += r_geometry[j].GetSolutionStepValue(PRESSURE) * N[j];
    }

    double DivTerm = 0.0;

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        double Qi = 0.0;
        for (unsigned int d = 0; d < mDim; ++d) {
            Qi += DN_DX(i,d) * BodyForce[d];
        }

        rRHS[FirstRow + mDim] += Weight * TauOne * Qi;

        for (unsigned int j = 0; j < NumNodes; ++j)
        {
            double Lij = 0.0;
            for (unsigned int d = 0; d < mDim; ++d)
            {
                Lij += DN_DX(i,d) * DN_DX(j,d);
                DivTerm = Weight * N[i] * DN_DX(j,d);
                rLHS(FirstRow+mDim,FirstCol + d) += DivTerm; // Divergence term
                rLHS(FirstCol + d,FirstRow+mDim) -= DivTerm; // Gradient term
            }
            // Stabilization (grad p, grad q)
            rLHS(FirstRow+mDim,FirstCol+mDim) += Weight * TauOne * Lij;

            // Update column index
            FirstCol += BlockSize;
        }

        // --- RHS corresponding term ---
        rRHS(FirstRow+mDim) -= Weight * N[i] * div_u_current ; // q div(u)
        for (unsigned int m = 0; m < mDim; ++m) {
            // rRHS(FirstRow+mDim) -= Weight * N[i] * div_u_previous_iteration[m] ; // ERRORE
            rRHS(FirstRow+m) += Weight * DN_DX(i,m) * pressure_previous_iteration ;
            // Stabilization
            rRHS(FirstRow+mDim) -= Weight * TauOne * DN_DX(i,m) * grad_p_previous_iteration[m];
        }

        // Update matrix indices
        FirstCol = 0;
        FirstRow += BlockSize;
    }
}


void StokesElement::AddSecondOrderStabilizationTerms(MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const double Density,
        const array_1d<double,3> &BodyForce,
        const double TauOne,
        const ShapeFunctionsType &N,
        const ShapeDerivativesType &DN_DX,
        const double GaussWeight,
        const Matrix& r_D,
        Vector& r_stress_vector)
{
    const unsigned int number_of_points = this->GetGeometry().PointsNumber();
    const unsigned int BlockSize = mDim+1;
    auto r_geometry = GetGeometry();

    // Add third order term of the stabilization
    if (mBasisFunctionsOrder > 1) {
        Matrix B = ZeroMatrix(3,number_of_points*mDim);
        CalculateB(B, DN_DX);
        const Matrix DDN_DDe = GetGeometry().ShapeFunctionDerivatives(2, 0, GetGeometry().GetDefaultIntegrationMethod());
        Matrix B_derivative_x;
        Matrix B_derivative_y;
        CalculateB_derivative_x(B_derivative_x, DDN_DDe);
        CalculateB_derivative_y(B_derivative_y, DDN_DDe);

        // Computed using the Gauss Points in the same knot span
        Vector divergence_of_sigma = this->GetValue(RECOVERED_STRESS);
        // Matrix constitutive_matrix_dx = this->GetValue(MATERIAL_STIFFNESS_MATRIX);
        // Matrix constitutive_matrix_dy = this->GetValue(GEOMETRIC_STIFFNESS_MATRIX);

        // initilize the div(sigma) matrix 2x(2n)
        Vector div_sigma_1 = ZeroVector(mDim*number_of_points);
        Vector div_sigma_2 = ZeroVector(mDim*number_of_points);

        Vector r_D_0(3); Vector r_D_1(3); Vector r_D_2(3);
        // Vector r_D_dx_0(3); Vector r_D_dx_2(3);
        // Vector r_D_dy_1(3); Vector r_D_dy_2(3);
        for (std::size_t i = 0; i < 3; ++i) {
            r_D_0(i) = r_D(0, i); r_D_1(i) = r_D(1, i); r_D_2(i) = r_D(2, i);
            // r_D_dx_0(i) = constitutive_matrix_dx(0, i);
            // r_D_dx_2(i) = constitutive_matrix_dx(2, i);
            // r_D_dy_1(i) = constitutive_matrix_dy(1, i); 
            // r_D_dy_2(i) = constitutive_matrix_dy(2, i);
        }
        div_sigma_1 = prod(trans(B_derivative_x), r_D_0) + prod(trans(B_derivative_y), r_D_2) ; // + prod(trans(B), r_D_dx_0) +  prod(trans(B), r_D_dy_2);
        div_sigma_2 = prod(trans(B_derivative_y), r_D_1) + prod(trans(B_derivative_x), r_D_2) ; // + prod(trans(B), r_D_dy_1) +  prod(trans(B), r_D_dx_2);
        
        // NEW FORMULATION___________________________________________________________________________________________
        Matrix div_sigma = ZeroMatrix(2, mDim * number_of_points);
        for (std::size_t i = 0; i < mDim * number_of_points; ++i) {
            div_sigma(0, i) = div_sigma_1(i);
            div_sigma(1, i) = div_sigma_2(i);
        }
        Matrix DN_DX_matrix = ZeroMatrix(2, number_of_points);
        for (std::size_t i = 0; i < number_of_points; ++i) {
            DN_DX_matrix(0, i) = DN_DX(i, 0);
            DN_DX_matrix(1, i) = DN_DX(i, 1);
        }
        Matrix M = ZeroMatrix(number_of_points*mDim, mDim);
        CalculateI(N, M);
        Matrix DB = prod(r_D, B);
        Matrix DBI_div_sigma = prod(Matrix(prod(DB, M)), div_sigma);
        Matrix DBI_DN_DX = prod(Matrix(prod(DB, M)), DN_DX_matrix);
        Matrix tempMatrix_sigma_sigma    = TauOne * GaussWeight * prod(trans(B), DBI_div_sigma);
        // Matrix tempMatrix_pressure_term2 = TauOne * GaussWeight * prod(trans(B), DBI_DN_DX);
        Matrix forcing = ZeroMatrix(2,1);
        forcing(0,0) = BodyForce[0];
        forcing(1,0) = BodyForce[1];
        Matrix DBI_f = prod(Matrix(prod(DB, M)), forcing);
        Matrix tempMatrix_forcing_term = TauOne * GaussWeight * prod(trans(B), DBI_f);
        Vector DN_DX_vector_1 = ZeroVector(number_of_points);
        Vector DN_DX_vector_2 = ZeroVector(number_of_points);
        for (std::size_t i = 0; i < number_of_points; ++i) {
            DN_DX_vector_1(i) = DN_DX(i, 0); 
            DN_DX_vector_2(i) = DN_DX(i, 1); 
        }
        Matrix tempMatrix_pressure_term = TauOne * GaussWeight * (outer_prod(div_sigma_1, DN_DX_vector_1) + outer_prod(div_sigma_2, DN_DX_vector_2));

        for (IndexType i = 0; i < number_of_points; ++i)
        {
            for (IndexType j = 0; j < number_of_points; ++j)
            {
                for (IndexType dim2 = 0; dim2 < mDim; ++dim2) 
                {
                    //// Assemble 3rd order stabilization term -> sigma-gradP
                    // rLeftHandSideMatrix(i * BlockSize + dim2, j * BlockSize + mDim) -= tempMatrix_pressure_term2(i * mDim + dim2, j); ///////// ATTENTION!!
                    
                    //// Assemble 2nd order stabilization term -> gradQ-sigma // this one
                    // rLeftHandSideMatrix(j * BlockSize + mDim, i * BlockSize + dim2) -= tempMatrix_pressure_term(i * mDim + dim2, j);
                    
                    for (IndexType dim1 = 0; dim1 < mDim; ++dim1)
                    {
                        //// Assemble 3rd order stabilization term -> sigma-sigma
                        // rLeftHandSideMatrix(i * BlockSize + dim1, j * BlockSize + dim2) += tempMatrix_sigma_sigma(i * mDim + dim1, j * mDim + dim2);
                    }
                }
            }
        }
        for(unsigned int i = 0; i < number_of_points; ++i)
        {
            // Assemble 3rd order stabilization term -> Body force
            for(unsigned int idim = 0; idim < mDim; ++idim) {
                // rRightHandSideVector[i * BlockSize + idim] -= tempMatrix_forcing_term(i * mDim + idim, 0);
            }
        }

        // Compute u' directly
        Vector grad_p_previous_iteration = ZeroVector(2);
        for(unsigned int j = 0; j < number_of_points; ++j) {
            grad_p_previous_iteration[0] += r_geometry[j].GetSolutionStepValue(PRESSURE) * DN_DX(j,0) ;
            grad_p_previous_iteration[1] += r_geometry[j].GetSolutionStepValue(PRESSURE) * DN_DX(j,1) ;
        }
        Vector u_tilde_n = ZeroVector(2);
        u_tilde_n += TauOne * (BodyForce + divergence_of_sigma - grad_p_previous_iteration);
        Matrix u_tilde_n_matrix = ZeroMatrix(2,1);
        u_tilde_n_matrix(0,0) = u_tilde_n[0];
        u_tilde_n_matrix(1,0) = u_tilde_n[1];
        
        Matrix DBI_u_tilde = prod(Matrix(prod(DB, M)), u_tilde_n_matrix);
        Matrix tempMatrix_tau_residuo = GaussWeight * prod(trans(B), DBI_u_tilde);
        for(unsigned int i = 0; i < number_of_points; ++i)
        {
            for(unsigned int idim = 0; idim < mDim; ++idim) {
                // rRightHandSideVector[i * BlockSize + idim] -= tempMatrix_tau_residuo(i * mDim + idim, 0);
            }
        }

        // Integrate by parts the term  < grad q, div(sigma) >
        Matrix Q;
        CalculateB_second_order(Q, DDN_DDe);
        Matrix tempMatrix_grad_q_div_sigma = TauOne * GaussWeight * prod(trans(Q), DB);
        for (IndexType i = 0; i < number_of_points; ++i)
        {
            for (IndexType j = 0; j < number_of_points; ++j)
            {
                for (IndexType dim2 = 0; dim2 < mDim; ++dim2) 
                {
                    //// Assemble 2nd order stabilization term -> gradQ-sigma
                    // rLeftHandSideMatrix(j * BlockSize + mDim, i * BlockSize + dim2) += tempMatrix_grad_q_div_sigma(j, i * mDim + dim2); 
                }
            }
        }
        // corresponding RHS term
        Vector tempMatrix_grad_q_div_sigma_RHS = TauOne * GaussWeight * prod(trans(Q), r_stress_vector);
        for (IndexType i = 0; i < number_of_points; ++i) {   
            // Assemble 2nd order stabilization term -> gradQ-sigma_u 
            // rRightHandSideVector(i * BlockSize + mDim) -= tempMatrix_grad_q_div_sigma_RHS(i);
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
        // Vector tempVector_sigma_u_sigma_v_rhs = prod(tempMatrix_sigma_sigma, velocity_previous) ;
        // Vector tempVector_sigma_v_grad_p_rhs = prod(tempMatrix_pressure_term2, pressure_previous) ; //////////// ATTENZIONE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // Vector tempVector_sigma_u_grad_q_rhs= prod(trans(tempMatrix_pressure_term), velocity_previous) ;

        // Matrix divergence_of_sigma_matrix = ZeroMatrix(2,1);
        // divergence_of_sigma_matrix(0,0) = divergence_of_sigma[0];
        // divergence_of_sigma_matrix(1,0) = divergence_of_sigma[1];
        // Matrix DBI_div_sigma_prev_iteration = prod(Matrix(prod(DB, M)), divergence_of_sigma_matrix);
        // Matrix tempMatrix_div_sigma_prev_iteration = TauOne * GaussWeight * prod(trans(B), DBI_div_sigma_prev_iteration);
        // --- RHS corresponding term ---
        for (IndexType i = 0; i < number_of_points; ++i) {   
            for (IndexType dim1 = 0; dim1 < mDim; ++dim1) {
                    // Assemble 3rd order stabilization term -> sigma-sigma
                    // rRightHandSideVector[i * BlockSize + dim1] -= tempMatrix_div_sigma_prev_iteration(i * mDim + dim1, 0);

                    // Assemble 3rd order stabilization term -> sigma_v-gradP
                    // rRightHandSideVector(i * BlockSize + dim1) += tempVector_sigma_v_grad_p_rhs(i * mDim + dim1);

                    // Assemble 2nd order stabilization term -> gradQ-sigma_u 
                    rRightHandSideVector(i * BlockSize + mDim) += GaussWeight * TauOne * DN_DX(i,dim1) * divergence_of_sigma[dim1];
                }
        }
        // END NEW FORMULATION---------------------------------------------------------------------------------------------



        //// OLD FORMULATION -> grad v : C : grad (u+u')
        // Matrix tempMatrix_sigma_sigma = TauOne * GaussWeight * (outer_prod(div_sigma_1, div_sigma_1) + outer_prod(div_sigma_2, div_sigma_2));
        // Vector DN_DX_vector_1 = ZeroVector(number_of_points);
        // Vector DN_DX_vector_2 = ZeroVector(number_of_points);
        // for (std::size_t i = 0; i < number_of_points; ++i) {
        //     DN_DX_vector_1(i) = DN_DX(i, 0); 
        //     DN_DX_vector_2(i) = DN_DX(i, 1); 
        // }
        // Matrix tempMatrix_pressure_term = TauOne * GaussWeight * (outer_prod(div_sigma_1, DN_DX_vector_1) + outer_prod(div_sigma_2, DN_DX_vector_2));

        // for (IndexType i = 0; i < number_of_points; ++i)
        // {
        //     for (IndexType j = 0; j < number_of_points; ++j)
        //     {
        //         for (IndexType dim2 = 0; dim2 < mDim; ++dim2) 
        //         {
        //             //// Assemble 3rd order stabilization term -> sigma-gradP
        //             rLeftHandSideMatrix(i * BlockSize + dim2, j * BlockSize + mDim) += tempMatrix_pressure_term(i * mDim + dim2, j);
                    
        //             //// Assemble 2nd order stabilization term -> gradQ-sigma
        //             rLeftHandSideMatrix(j * BlockSize + mDim, i * BlockSize + dim2) -= tempMatrix_pressure_term(i * mDim + dim2, j);
                    
        //             for (IndexType dim1 = 0; dim1 < mDim; ++dim1)
        //             {
        //                 // Assemble 3rd order stabilization term -> sigma-sigma
        //                 rLeftHandSideMatrix(i * BlockSize + dim1, j * BlockSize + dim2) -= tempMatrix_sigma_sigma(i * mDim + dim1, j * mDim + dim2);
        //             }
        //         }
        //     }
        // }
        // for(unsigned int i = 0; i < number_of_points; ++i)
        // {
        //     // Assemble 3rd order stabilization term -> Body force
        //     for(unsigned int idim = 0; idim < mDim; ++idim) {
        //         rRightHandSideVector[i * BlockSize + idim] += TauOne * GaussWeight * 
        //                                             (div_sigma_1[i * mDim + idim] * BodyForce[0] + div_sigma_2[i * mDim + idim] * BodyForce[1]);
        //     }
        // }

        // // --- RHS corresponding term ---
        // Vector velocity_previous = ZeroVector(mDim*number_of_points);
        // Vector pressure_previous = ZeroVector(number_of_points);
        // IndexType index = 0;
        // for (IndexType i = 0; i < number_of_points; ++i) { 
        //     velocity_previous[index++] = r_geometry[i].GetSolutionStepValue(VELOCITY_X);
        //     velocity_previous[index++] = r_geometry[i].GetSolutionStepValue(VELOCITY_Y);
        //     pressure_previous[i] = r_geometry[i].GetSolutionStepValue(PRESSURE);
        // }
        // Vector tempVector_sigma_u_sigma_v_rhs = prod(tempMatrix_sigma_sigma, velocity_previous) ;
        // Vector tempVector_sigma_v_grad_p_rhs = prod(tempMatrix_pressure_term, pressure_previous) ;
        // Vector tempVector_sigma_u_grad_q_rhs= prod(trans(tempMatrix_pressure_term), velocity_previous) ;

        //         // // for (IndexType i = 0; i < number_of_points; ++i)
        //         // {
        //         //     for (IndexType dim1 = 0; dim1 < mDim; ++dim1)
        //         //         {
        //         //             // Assemble 3rd order stabilization term -> sigma-sigma
        //         //             rRightHandSideVector(i * BlockSize + dim1) += tempVector_sigma_u_sigma_v_rhs(i * mDim + dim1);
        //         //             // Assemble 3rd order stabilization term -> sigma_v-gradP
        //         //             rRightHandSideVector(i * BlockSize + dim1) -= tempVector_sigma_v_grad_p_rhs(i * mDim + dim1);
        //         //         }
        //         //     // Assemble 2nd order stabilization term -> gradQ-sigma_u 
        //         //     rRightHandSideVector(i * BlockSize + mDim) += tempVector_sigma_u_grad_q_rhs(i);
        //         // }

        //         // // --- RHS corresponding term --- USING DIV(SIGMA)
        //         // Vector pressure_previous = ZeroVector(number_of_points);
        //         // IndexType index = 0;
        //         // for (IndexType i = 0; i < number_of_points; ++i) { 
        //         //     pressure_previous[i] = r_geometry[i].GetSolutionStepValue(PRESSURE);
        //         // }
        //         // Vector tempVector_sigma_v_grad_p_rhs = prod(tempMatrix_pressure_term, pressure_previous) ;

        // Matrix divergence_of_sigma_matrix = ZeroMatrix(2,1);
        // divergence_of_sigma_matrix(0,0) = divergence_of_sigma[0];
        // divergence_of_sigma_matrix(1,0) = divergence_of_sigma[1];
        // Matrix DBI_div_sigma_prev_iteration = prod(Matrix(prod(DB, M)), divergence_of_sigma_matrix);
        // Matrix tempMatrix_div_sigma_prev_iteration = TauOne * GaussWeight * prod(trans(B), DBI_div_sigma_prev_iteration);
        // // --- RHS corresponding term ---
        // for (IndexType i = 0; i < number_of_points; ++i) {   
        //     for (IndexType dim1 = 0; dim1 < mDim; ++dim1) {
        //             // Assemble 3rd order stabilization term -> sigma-sigma
        //             rRightHandSideVector(i * BlockSize + dim1) += GaussWeight * TauOne * (
        //                 div_sigma_1[i * mDim + dim1] * divergence_of_sigma[0] + div_sigma_2[i * mDim + dim1] * divergence_of_sigma[1]);

        //             // Assemble 3rd order stabilization term -> sigma_v-gradP
        //             rRightHandSideVector(i * BlockSize + dim1) -= tempVector_sigma_v_grad_p_rhs(i * mDim + dim1);

        //             // Assemble 2nd order stabilization term -> gradQ-sigma_u 
        //             rRightHandSideVector(i * BlockSize + mDim) += GaussWeight * TauOne * DN_DX(i,dim1) * divergence_of_sigma[dim1];
        //         }
        // }
    } // end if

}


void StokesElement::EquationIdVector(Element::EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const int dim = 2;
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType number_of_control_points = GetGeometry().size();
    const unsigned int LocalSize = (dim + 1) * number_of_control_points;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < number_of_control_points; i++)
    {
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
        if (dim > 2) rResult[Index++] = rGeom[i].GetDof(VELOCITY_Z).EquationId();
        rResult[Index++] = rGeom[i].GetDof(PRESSURE).EquationId();
    }
}


void StokesElement::GetDofList(
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



int StokesElement::Check(const ProcessInfo& rCurrentProcessInfo) const
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

Element::IntegrationMethod StokesElement::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

void StokesElement::CalculateB(
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


void StokesElement::CalculateI(
        const ShapeFunctionsType &N,
        Matrix& rI) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 2; // Only 2 DOFs per node in 2D

    // Resize B matrix to 3 rows (strain vector size) and appropriate number of columns
    if (rI.size1() != mat_size || rI.size2() != 2)
        rI.resize(mat_size, 2);

    noalias(rI) = ZeroMatrix(mat_size, 2);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        rI(2 * i, 0)     = 1;//;N[i];1/N[i]
        rI(2 * i, 1)     = 0.0;
        rI(2 * i + 1, 0) = 0.0;
        rI(2 * i + 1, 1) = 1;//N[i];1/N[i]
    }
}

void StokesElement::CalculateB_second_order(
        Matrix& B_second_order, 
        const ShapeDerivativesType& r_DDN_DDX) const
{
    const SizeType number_of_control_points = GetGeometry().size();

    // Resize B_second_order
    if (B_second_order.size1() != 3 || B_second_order.size2() != number_of_control_points)
        B_second_order.resize(3, number_of_control_points);

    noalias(B_second_order) = ZeroMatrix(3, number_of_control_points);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        // x-derivatives of shape functions
        B_second_order(0, i) = r_DDN_DDX(i, 0); // ∂²N_i / ∂²x
        // y-derivatives of shape functions 
        B_second_order(1, i) = r_DDN_DDX(i, 2); // ∂²N_i / ∂y²
        // ε_12 (xy component)
        B_second_order(2, i) = r_DDN_DDX(i, 1); // ∂²N_i / ∂y∂x
    }
}

    // // Add the second order term of the stabilization
    // if (mBasisFunctionsOrder > 1) {
    //     const Matrix DDN_DDe = GetGeometry().ShapeFunctionDerivatives(2, 0, GetGeometry().GetDefaultIntegrationMethod());
    //     Matrix B_second_order = ZeroMatrix(3,number_of_points);
    //     CalculateB_second_order(B_second_order, DDN_DDe);

    //     // Compute the product B_second_order^T D B
    //     Matrix tempMatrix = TauOne * GaussWeight * prod(trans(B_second_order), Matrix(prod(r_D, B)));
    //     // KRATOS_WATCH(tempMatrix)
    //     // exit(0);
    //     // Assemble the controbution
    //     for (IndexType i = 0; i < number_of_points; ++i)
    //     {
    //         for (IndexType j = 0; j < number_of_points; ++j)
    //         {
    //             for (IndexType dim2 = 0; dim2 < mDim; ++dim2)  
    //             {
    //                 // Assemble stabilization term
    //                 rLeftHandSideMatrix(i * BlockSize + mDim, j * BlockSize + dim2) += tempMatrix(i, j * mDim + dim2);
    //             }
    //         }
    //     }
    // }


    // Matrix tempMatrix = GaussWeight * prod(trans(B), Matrix(prod(r_D, B)));
    // for (IndexType i = 0; i < number_of_points; ++i)
    // {
    //     for (IndexType j = 0; j < number_of_points; ++j)
    //     {
    //         // Add only to velocity DOFs (i.e., positions 0, 1, 3, 4, 6, 7, ...)
    //         for (IndexType dim1 = 0; dim1 < mDim; ++dim1)  // Velocity components (vx, vy)
    //         {
    //             for (IndexType dim2 = 0; dim2 < mDim; ++dim2)  // Velocity components (vx, vy)
    //             {
    //                 // diffusion term
    //                 rLeftHandSideMatrix(i * BlockSize + dim1, j * BlockSize + dim2) += tempMatrix(i * mDim + dim1, j * mDim + dim2);
    //             }
    //         }
    //     }
    // }

void StokesElement::CalculateB_derivative_x(
        Matrix& B_derivative_x, 
        const ShapeDerivativesType& r_DDN_DDX) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * mDim; // Only 2 DOFs per node in 2D

    // Resize B matrix to 3 rows (strain vector size) and appropriate number of columns
    if (B_derivative_x.size1() != 3 || B_derivative_x.size2() != mat_size)
        B_derivative_x.resize(3, mat_size);

    noalias(B_derivative_x) = ZeroMatrix(3, mat_size);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        B_derivative_x(0, 2 * i)     = r_DDN_DDX(i, 0); // ∂²N_i / ∂x²
        B_derivative_x(1, 2 * i + 1) = r_DDN_DDX(i, 1); // ∂²N_i / ∂y∂x
        B_derivative_x(2, 2 * i)     = r_DDN_DDX(i, 1); // ∂²N_i / ∂y∂x
        B_derivative_x(2, 2 * i + 1) = r_DDN_DDX(i, 0); // ∂²N_i / ∂x²
    }
}

void StokesElement::CalculateB_derivative_y(
        Matrix& B_derivative_y, 
        const ShapeDerivativesType& r_DDN_DDX) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 2; // Only 2 DOFs per node in 2D

    // Resize B matrix to 3 rows (strain vector size) and appropriate number of columns
    if (B_derivative_y.size1() != 3 || B_derivative_y.size2() != mat_size)
        B_derivative_y.resize(3, mat_size);

    noalias(B_derivative_y) = ZeroMatrix(3, mat_size);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        B_derivative_y(0, 2 * i)     = r_DDN_DDX(i, 1); // ∂²N_i / ∂x∂y
        B_derivative_y(1, 2 * i + 1) = r_DDN_DDX(i, 2); // ∂²N_i / ∂y²
        B_derivative_y(2, 2 * i)     = r_DDN_DDX(i, 2); // ∂²N_i / ∂y²
        B_derivative_y(2, 2 * i + 1) = r_DDN_DDX(i, 1); // ∂²N_i / ∂x∂y
    }
}


void StokesElement::CalculateOnIntegrationPoints(
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

void StokesElement::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& D_constitutive_matrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == CONSTITUTIVE_MATRIX) {
        const unsigned int num_integration_points = this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());

        if (D_constitutive_matrix.size() != 1) {
            D_constitutive_matrix.resize(1);
        }

        // Resize each matrix in the vector to 3x3
        D_constitutive_matrix[0].resize(3, 3, false);
        noalias(D_constitutive_matrix[0]) = ZeroMatrix(3, 3); // Optionally initialize each matrix to zero

        D_constitutive_matrix[0] = this->CalculateConstitutiveMatrixAtIntegrationPoint(rCurrentProcessInfo);
    }
}

Vector StokesElement::CalculateStressAtIntegrationPoint(
    const ProcessInfo& rCurrentProcessInfo)
{
    Vector stress_voigt(3);
    
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_points = r_geometry.size();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const ShapeDerivativesType& DN_DX = DN_De[0];
    // constitutive law
    Matrix B = ZeroMatrix(3,number_of_points*mDim);
    CalculateB(B, DN_DX);
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);
    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    ConstitutiveVariables this_constitutive_variables(strain_size);
    Vector old_displacement(number_of_points*mDim);
    GetValuesVector(old_displacement);
    Vector old_strain = prod(B,old_displacement);
    // Values.SetStrainVector(this_constitutive_variables.StrainVector);
    Values.SetStrainVector(old_strain); // this is the input parameter
    Values.SetStressVector(this_constitutive_variables.StressVector); // this is an ouput parameter
    Values.SetConstitutiveMatrix(this_constitutive_variables.D);
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);
    stress_voigt = Values.GetStressVector();
    
    return stress_voigt;
}


Matrix StokesElement::CalculateConstitutiveMatrixAtIntegrationPoint(
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_points = r_geometry.size();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const ShapeDerivativesType& DN_DX = DN_De[0];
    // constitutive law
    Matrix B = ZeroMatrix(3,number_of_points*mDim);
    CalculateB(B, DN_DX);
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);
    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    ConstitutiveVariables this_constitutive_variables(strain_size);
    Vector old_displacement(number_of_points*mDim);
    GetValuesVector(old_displacement);
    Vector old_strain = prod(B,old_displacement);
    // Values.SetStrainVector(this_constitutive_variables.StrainVector);
    Values.SetStrainVector(old_strain); // this is the input parameter
    Values.SetStressVector(this_constitutive_variables.StressVector); // this is an ouput parameter
    Values.SetConstitutiveMatrix(this_constitutive_variables.D);
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);
    const Matrix& r_D = Values.GetConstitutiveMatrix();
    
    return r_D;
}


void StokesElement::GetValuesVector(
        Vector& rValues) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 2;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        const array_1d<double, 2 >& velocity = GetGeometry()[i].GetSolutionStepValue(VELOCITY);
        IndexType index = i * 2;

        rValues[index] = velocity[0];
        rValues[index + 1] = velocity[1];
    }
}


void StokesElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
    const auto& r_geometry = GetGeometry();
    const SizeType nb_nodes = r_geometry.size();

    // Integration Points
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
    // Shape function values
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

    const ProcessInfo& r_process_info = rCurrentProcessInfo;

    double rOutput_vel_x = 0;   double rOutput_vel_y = 0;   double rOutput_pressure = 0;
    for (IndexType i = 0; i < nb_nodes; ++i)
    {
        double output_solution_step_value_vel_x = r_geometry[i].GetSolutionStepValue(VELOCITY_X);
        double output_solution_step_value_vel_y = r_geometry[i].GetSolutionStepValue(VELOCITY_Y);
        double output_solution_step_value_pressure = r_geometry[i].GetSolutionStepValue(PRESSURE);
        rOutput_vel_x += r_N(0, i) * output_solution_step_value_vel_x;
        rOutput_vel_y += r_N(0, i) * output_solution_step_value_vel_y;
        rOutput_pressure += r_N(0, i) * output_solution_step_value_pressure;
    } 

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    const ShapeDerivativesType& DN_DX = DN_De[0];

    // constitutive law
    Matrix B = ZeroMatrix(3,nb_nodes*mDim);
    CalculateB(B, DN_DX);
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);
    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    ConstitutiveVariables this_constitutive_variables(strain_size);
    Vector old_displacement(nb_nodes*mDim);
    GetValuesVector(old_displacement);
    Vector old_strain = prod(B,old_displacement);
    Values.SetStrainVector(old_strain); // this is the input parameter
    Values.SetStressVector(this_constitutive_variables.StressVector); // this is an ouput parameter
    Values.SetConstitutiveMatrix(this_constitutive_variables.D);
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);

    // double yielded_state = 0.0;
    // const double sigma_y = 0;

    // const double sigma_y = 500;

    // const double sigma_y = 395;
    // const double sigma_y = 790;
    const double sigma_y = 987.5;
    // const double sigma_y = 1086;
    // const double sigma_y = 1145.5;
    Vector& StressVector = this_constitutive_variables.StressVector;

    // const double tau = std::sqrt(StressVector[0] * StressVector[0] + 
    //                              StressVector[1] * StressVector[1] + 
    //                              0.5 * StressVector[2] * StressVector[2]);
    // KRATOS_WATCH(tau)

    const double trS = StressVector[0] + StressVector[1];
    const double hydrostatic_stress = trS / 3.0;

    // Calculate the deviatoric part of the stress tensor
    Vector deviatoric_stress(3);
    deviatoric_stress[0] = StressVector[0] - hydrostatic_stress;
    deviatoric_stress[1] = StressVector[1] - hydrostatic_stress;
    deviatoric_stress[2] = StressVector[2]; // Shear stress remains the same in Voigt notation

    // Calculate the equivalent deviatoric stress (tau) using von Mises criterion
    const double tau = std::sqrt(
        0.5*deviatoric_stress[0] * deviatoric_stress[0] + 
        0.5*deviatoric_stress[1] * deviatoric_stress[1] + 
        1.0 * deviatoric_stress[2] * deviatoric_stress[2]
    );

    double yielded_state = (tau < sigma_y) ? 1.0 : 0.0;
    // mpConstitutiveLaw->CalculateValue(Values, STRAIN_ENERGY, yielded_state); // Usa STRAIN_ENERGY come proxy
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    #pragma omp critical
    {
        std::ofstream output_file("txt_files/output_results_GPs.txt", std::ios::app);
        if (output_file.is_open()) {
        output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
        output_file << rOutput_vel_x << " " << rOutput_vel_y << " " << rOutput_pressure << " " 
                    << r_geometry.Center().X() << " " << r_geometry.Center().Y() << " " << r_geometry.Center().Z() << " " << integration_points[0].Weight() << std::endl;
        output_file.close();
        }   

        // Same but saving the yielded stress value
        std::ofstream yield_output_file("txt_files/yielded_states_GPs.txt", std::ios::app);
        if (yield_output_file.is_open()) {
            yield_output_file << std::scientific << std::setprecision(14); // Precisione 10^-14
            yield_output_file << yielded_state << " "
                              << r_geometry.Center().X() << " " << r_geometry.Center().Y() << std::endl;
            yield_output_file.close();
        }

    }  
    // SetValue the integraion_weigth in order to do the postprocess in Python
    SetValue(INTEGRATION_WEIGHT, integration_points[0].Weight() );
}

void StokesElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
}

/// Reads in a json formatted file and returns its KratosParameters instance.
Parameters StokesElement::ReadParamatersFile(
    const std::string& rDataFileName) const
{
    std::ifstream infile(rDataFileName);

    std::stringstream buffer;
    buffer << infile.rdbuf();

    return Parameters(buffer.str());
};



} // Namespace Kratos