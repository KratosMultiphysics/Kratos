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
#include "custom_elements/stokes_element_transient.h"

#include "utilities/math_utils.h"

namespace Kratos
{

StokesElementTransient::StokesElementTransient(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(
        NewId,
        pGeometry)
{
}

StokesElementTransient::StokesElementTransient(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(
        NewId,
        pGeometry,
        pProperties)
{
}

Element::Pointer StokesElementTransient::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<StokesElementTransient>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer StokesElementTransient::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<StokesElementTransient>(NewId, pGeom, pProperties);
}

// Deconstructor

StokesElementTransient::~StokesElementTransient()
{
}

void StokesElementTransient:: Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
    GeometryType& r_geometry = this->GetGeometry();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    mDim = DN_De[0].size2();
}


void StokesElementTransient::InitializeMaterial()
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


void StokesElementTransient::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo)
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
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points,number_of_points); //resetting LHS
    
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
    //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
    //this is ok under the hypothesis that no history dependent behavior is employed
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);

    Vector& r_stress_vector = Values.GetStressVector();
    const Matrix& r_D = Values.GetConstitutiveMatrix();

    double Density;
    double Viscosity;
    array_1d<double,3> BodyForce = ZeroVector(3);
    array_1d<double,3> BodyForce_old = ZeroVector(3);
    array_1d<double,3> Velocity = ZeroVector(3);

    // Interpolation
    this->EvaluateInPoint(Density,DENSITY,N,r_geometry);
    this->EvaluateInPoint(Viscosity,VISCOSITY,N,r_geometry);
    BodyForce = this->GetValue(BODY_FORCE);
    BodyForce_old = this->GetValue(FORCE);

    this->EvaluateInPoint(Velocity,VELOCITY,N ,r_geometry);

    mTauOne = 0.0;
    mTauTwo = 0.0;
    CalculateTau(Density*Viscosity);

    // Add velocity terms in momentum equation
    this->AddMomentumTerms(rLeftHandSideMatrix,rRightHandSideVector,Density,Viscosity,BodyForce,BodyForce_old,mTauTwo,N,DN_DX,GaussWeight, r_D, r_stress_vector);

    // Add velocity-pressure terms
    this->AddContinuityTerms(rLeftHandSideMatrix,rRightHandSideVector,Density,BodyForce,BodyForce_old,mTauOne,N,DN_DX,GaussWeight);

    // Add Second-Order stabilization terms from VMS
    this->AddSecondOrderStabilizationTerms(rLeftHandSideMatrix,rRightHandSideVector,Density,mTauOne,N,DN_DX,GaussWeight,r_D, r_stress_vector);

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
    // }
}

void StokesElementTransient::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{   
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

void StokesElementTransient::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

void StokesElementTransient::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int BlockSize = mDim+1;
    auto r_geometry = GetGeometry();
    const int num_dofs_per_node = NumNodes * (mDim + 1);
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    const ShapeFunctionsType& N = row(N_gausspoint,0);
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const ShapeDerivativesType& DN_DX = DN_De[0];
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    if(rMassMatrix.size1() != num_dofs_per_node)
        rMassMatrix.resize(num_dofs_per_node,num_dofs_per_node,false);
    noalias(rMassMatrix) = ZeroMatrix(NumNodes,NumNodes); //resetting MassMatrix

    // Define mass matrix and previous velocity vector
    Matrix mass_matrix;
    Vector coeff_previous_velocity;
    Vector coeff_current_iteration_velocity;
    // Resize mass matrix if necessary
    if (mass_matrix.size1() != num_dofs_per_node) {
        mass_matrix.resize(num_dofs_per_node, num_dofs_per_node, false);
    }
    noalias(mass_matrix) = ZeroMatrix(num_dofs_per_node, num_dofs_per_node);
    
    // Build mass matrix for velocity DOFs
    for (std::size_t i = 0; i < NumNodes; ++i) {
        for (std::size_t j = 0; j < NumNodes; ++j) {
            for (std::size_t dim = 0; dim < mDim; ++dim) {
                mass_matrix(i * BlockSize + dim, j * BlockSize + dim) += N[i] * N[j] * integration_points[0].Weight();
            }
        }
    }

    // Add transient term to LHS
    noalias(rMassMatrix) += mass_matrix; //(1 / delta_t) * mass_matrix; 

    // Add the VMS term with du/dt
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        for (unsigned int j = 0; j < NumNodes; ++j)
        {
            for (unsigned int d = 0; d < mDim; ++d) {
                // Stabilization-time (grad q, u^n+1)/delta_t
                rMassMatrix(i*BlockSize+mDim, j*BlockSize + d) += integration_points[0].Weight() * mTauOne * DN_DX(i, d) * N[j] ;
            }
        }
    }
}

void StokesElementTransient::CalculateTau(const double DynViscosity)
{   
    // // Estimate element size
    double h = ElementSize();
    // double h = h_element_size;

    // const Parameters refinements_parameters = ReadParamatersFile("refinements.iga.json");
    // int insertions = refinements_parameters["refinements"][0]["parameters"]["insert_nb_per_span_u"].GetInt();
    // h = 2.0/(insertions+1) ;

    mTauOne = std::pow(h, 2) / ( 4.0 * DynViscosity ); //  /  pow(mBasisFunctionsOrder, 4);
    mTauTwo = DynViscosity;
}

double StokesElementTransient::ElementSize()
{
    return 1.128379167 * sqrt(this->GetGeometry().DomainSize()); 
}

void StokesElementTransient::AddMomentumTerms(MatrixType &rLHS,
        VectorType &rRHS,
        const double Density,
        const double Viscosity,
        const array_1d<double,3> &BodyForce,
        const array_1d<double,3> &BodyForce_old,
        const double mTauTwo,
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
    // Without constitutive law
    // Matrix tempMatrix = Weight * 2* prod(trans(B), B);

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
    Vector u_old(NumNodes * mDim);
    for (IndexType i = 0; i < NumNodes; ++i)
    {
        u_old[i * mDim]     = r_geometry[i].FastGetSolutionStepValue(VELOCITY_X); // X velocity at previous time step
        u_old[i * mDim+1]   = r_geometry[i].FastGetSolutionStepValue(VELOCITY_Y); // Y velocity at previous time step
    }
    Vector rhs_previous = prod(tempMatrix, u_old);

    Vector internalForces = Weight * prod(trans(B), r_stress_vector);
    for (IndexType i = 0; i < NumNodes; ++i)
    {
        // Add only to the velocity DOFs (i.e., positions 0, 1, 3, 4, 6, 7, ...)
        for (IndexType dim1 = 0; dim1 < mDim; ++dim1)
        {
            rRHS(i * BlockSize +  dim1) -= internalForces(i * mDim + dim1); // BDBu at previous iteration
            // Without constitutive law
            // rRHS(i * BlockSize + dim1) -= rhs_previous(i * mDim + dim1); // (B^T B u) at previous iteration
        }  
    }
    

    unsigned int FirstRow = 0;
    unsigned int FirstCol = 0;

    double div_u_current = 0.0;
    for(unsigned int j = 0; j < NumNodes; ++j) {
        div_u_current += r_geometry[j].FastGetSolutionStepValue(VELOCITY_X) * DN_DX(j, 0) + 
                         r_geometry[j].FastGetSolutionStepValue(VELOCITY_Y) * DN_DX(j, 1) ;
    } 

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
                    rLHS(FirstRow+m,FirstCol+n) += mTauTwo * Weight * DN_DX(i ,m) * DN_DX(j, n);
                }
            }

            // Update column index
            FirstCol += BlockSize;
        }
        // RHS corresponding term
        for (unsigned int m = 0; m < mDim; ++m) {
            rRHS(FirstRow + m) -= mTauTwo * Weight * DN_DX(i, m) * div_u_current ;
        }

        // Update matrix indices
        FirstRow += BlockSize;
        FirstCol = 0;
    }
}


void StokesElementTransient::AddContinuityTerms(MatrixType &rLHS,
        VectorType &rRHS,
        const double Density,
        const array_1d<double,3> &BodyForce,
        const array_1d<double,3> &BodyForce_old,
        const double mTauOne,
        const ShapeFunctionsType &N,
        const ShapeDerivativesType &DN_DX,
        const double Weight)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int BlockSize = mDim+1;
    auto r_geometry = GetGeometry();

    unsigned int FirstRow = 0;
    unsigned int FirstCol = 0;

    // Compute the current and previous velocities
    double pressure_current = 0.0;
    double div_u_current = 0.0;
    Vector grad_p_current = ZeroVector(mDim);

    for(unsigned int j = 0; j < NumNodes; ++j) {
        div_u_current += r_geometry[j].FastGetSolutionStepValue(VELOCITY_X) * DN_DX(j, 0) +
                         r_geometry[j].FastGetSolutionStepValue(VELOCITY_Y) * DN_DX(j, 1) ;
        grad_p_current[0] += r_geometry[j].FastGetSolutionStepValue(PRESSURE) * DN_DX(j, 0) ;
        grad_p_current[1] += r_geometry[j].FastGetSolutionStepValue(PRESSURE) * DN_DX(j, 1) ;
        pressure_current += r_geometry[j].FastGetSolutionStepValue(PRESSURE) * N[j];
    }

    double DivTerm = 0.0;
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        double Qi = 0.0;
        double Qi_old = 0.0;
        for (unsigned int d = 0; d < mDim; ++d) {
            Qi += DN_DX(i,d) * BodyForce[d];
            Qi_old += DN_DX(i,d) * BodyForce_old[d];
        }

        rRHS[FirstRow + mDim] += Weight * mTauOne * Qi;

        for (unsigned int j = 0; j < NumNodes; ++j)
        {
            double Lij = 0.0;
            for (unsigned int d = 0; d < mDim; ++d)
            {
                Lij += DN_DX(i,d) * DN_DX(j,d);
                DivTerm = Weight * N[i] * DN_DX(j,d);
                rLHS(FirstRow+mDim,FirstCol + d) += DivTerm; // Divergence term -> q div(u)
                rLHS(FirstCol + d,FirstRow+mDim) -= DivTerm; // Gradient term -> p div(v)
            }
            // Stabilization (grad q, grad p)
            rLHS(FirstRow+mDim,FirstCol+mDim) += Weight * mTauOne * Lij;

            // Update column index
            FirstCol += BlockSize;
        }

        // --- RHS corresponding term ---
        rRHS(FirstRow + mDim) -= Weight * N[i] * div_u_current ; // q div(u)

        for (unsigned int m = 0; m < mDim; ++m) {
            rRHS(FirstRow + m) += Weight * DN_DX(i, m) * pressure_current ; // p div(v)

            // Stabilization (grad q, grad p)
            rRHS(FirstRow + mDim) -= Weight * mTauOne * DN_DX(i, m) * grad_p_current[m];
        }

        // Update matrix indices
        FirstCol = 0;
        FirstRow += BlockSize;
    }
}

void StokesElementTransient::AddTimeTerms(MatrixType &rLHS,
        VectorType &rRHS,
        const double Density,
        const double mTauOne,
        const ShapeFunctionsType &N,
        const ShapeDerivativesType &DN_DX,
        const double Weight,
        const double delta_t)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int BlockSize = mDim+1;
    auto r_geometry = GetGeometry();
    const int num_dofs_per_node = NumNodes * (mDim + 1);

    // Compute the current and previous velocities
    double pressure_previous_iteration = 0.0;
    Vector velocity_current = ZeroVector(mDim);
    Vector previous_velocity = r_geometry.GetValue(VELOCITY);
    for(unsigned int j = 0; j < NumNodes; ++j) {
        velocity_current[0] += r_geometry[j].FastGetSolutionStepValue(VELOCITY_X) * N[j];
        velocity_current[1] += r_geometry[j].FastGetSolutionStepValue(VELOCITY_Y) * N[j];
    }

    // Define mass matrix and previous velocity vector
    Matrix mass_matrix;
    Vector coeff_previous_velocity;
    Vector coeff_current_iteration_velocity;
    // Resize mass matrix if necessary
    if (mass_matrix.size1() != num_dofs_per_node) {
        mass_matrix.resize(num_dofs_per_node, num_dofs_per_node, false);
    }
    noalias(mass_matrix) = ZeroMatrix(num_dofs_per_node, num_dofs_per_node);
    
    // Build mass matrix for velocity DOFs
    for (std::size_t i = 0; i < NumNodes; ++i) {
        for (std::size_t j = 0; j < NumNodes; ++j) {
            for (std::size_t dim = 0; dim < mDim; ++dim) {
                mass_matrix(i * BlockSize + dim, j * BlockSize + dim) += N[i] * N[j] * Weight;
            }
        }
    }

    // Add transient term to LHS
    noalias(rLHS) += (1 / delta_t) * mass_matrix; 

    // Get previous velocity values for RHS
    coeff_current_iteration_velocity.resize(num_dofs_per_node, false);
    for (unsigned int i = 0; i < NumNodes; i++) {
        coeff_current_iteration_velocity[i * BlockSize]     = r_geometry[i].FastGetSolutionStepValue(VELOCITY_X); // X velocity at current time step
        coeff_current_iteration_velocity[i * BlockSize + 1] = r_geometry[i].FastGetSolutionStepValue(VELOCITY_Y); // Y velocity at current time step
    }

    // // Update RHS with transient term from the previous time step
    // noalias(rRHS) += prod(mass_matrix, coeff_previous_velocity) / delta_t;
    
    // Loop over nodes to compute transient term contribution to RHS
    for (std::size_t i = 0; i < NumNodes; ++i) {
        for (std::size_t dim = 0; dim < mDim; ++dim) {
            // Add transient term contribution to RHS
            rRHS[i * BlockSize + dim] += Weight * N[i] * previous_velocity[dim] / delta_t;
        }
    }
    // Corresonding RHS -> transient term from the current time step
    noalias(rRHS) -= prod(mass_matrix, coeff_current_iteration_velocity) / delta_t;


    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        for (unsigned int j = 0; j < NumNodes; ++j)
        {
            for (unsigned int d = 0; d < mDim; ++d) {
                // Stabilization-time (grad q, u^n+1)/delta_t
                rLHS(i*BlockSize+mDim, j*BlockSize + d) += Weight * mTauOne * DN_DX(i, d) * N[j] / delta_t ;
            }
        }
        
        // Stabilization-time RHS part with u_old
        for (unsigned int d = 0; d < mDim; ++d) {
            rRHS(i*BlockSize+mDim) += Weight * mTauOne * DN_DX(i, d) * previous_velocity[d] / delta_t;

            rRHS(i*BlockSize+mDim) -= Weight * mTauOne * DN_DX(i, d) * velocity_current[d] / delta_t ;
        }
        
    }
}


void StokesElementTransient::AddSecondOrderStabilizationTerms(MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const double Density,
        const double mTauOne,
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
        const Matrix DDN_DDe = GetGeometry().ShapeFunctionDerivatives(2, 0, GetGeometry().GetDefaultIntegrationMethod());

        Matrix B = ZeroMatrix(3,number_of_points*mDim);
        CalculateB(B, DN_DX);
        Matrix B_derivative_x;
        Matrix B_derivative_y;
        CalculateB_derivative_x(B_derivative_x, DDN_DDe);
        CalculateB_derivative_y(B_derivative_y, DDN_DDe);

        // Computed using the Gauss Points in the same knot span
        Vector divergence_of_sigma = this->GetValue(RECOVERED_STRESS);

        // initilize the div(sigma) matrix 2x(2n)
        Vector div_sigma_1 = ZeroVector(mDim*number_of_points);
        Vector div_sigma_2 = ZeroVector(mDim*number_of_points);

        Vector r_D_0(3); Vector r_D_1(3); Vector r_D_2(3);
        Vector r_D_dx_0(3); Vector r_D_dx_2(3);
        Vector r_D_dy_1(3); Vector r_D_dy_2(3);
        for (std::size_t i = 0; i < 3; ++i) {
            r_D_0(i) = r_D(0, i); r_D_1(i) = r_D(1, i); r_D_2(i) = r_D(2, i);
        }
        div_sigma_1 = prod(trans(B_derivative_x), r_D_0) + prod(trans(B_derivative_y), r_D_2); 
        div_sigma_2 = prod(trans(B_derivative_y), r_D_1) + prod(trans(B_derivative_x), r_D_2);

        Vector DN_DX_vector_1 = ZeroVector(number_of_points);
        Vector DN_DX_vector_2 = ZeroVector(number_of_points);
        for (std::size_t i = 0; i < number_of_points; ++i) {
            DN_DX_vector_1(i) = DN_DX(i, 0); 
            DN_DX_vector_2(i) = DN_DX(i, 1); 
        }
        Matrix tempMatrix_pressure_term = mTauOne * GaussWeight * (outer_prod(div_sigma_1, DN_DX_vector_1) + outer_prod(div_sigma_2, DN_DX_vector_2));

        for (IndexType i = 0; i < number_of_points; ++i)
        {
            for (IndexType j = 0; j < number_of_points; ++j)
            {
                for (IndexType dim2 = 0; dim2 < mDim; ++dim2) 
                {
                    //// Assemble 2nd order stabilization term -> gradQ-div(sigma)
                    rLeftHandSideMatrix(j * BlockSize + mDim, i * BlockSize + dim2) -= tempMatrix_pressure_term(i * mDim + dim2, j);    
                }
            }
        }
        // --- RHS corresponding term ---
        Matrix divergence_of_sigma_matrix = ZeroMatrix(2,1);
        divergence_of_sigma_matrix(0,0) = divergence_of_sigma[0];
        divergence_of_sigma_matrix(1,0) = divergence_of_sigma[1];
        // --- RHS corresponding term ---
        for (IndexType i = 0; i < number_of_points; ++i) {   
            for (IndexType dim1 = 0; dim1 < mDim; ++dim1) {
                    // Assemble 2nd order stabilization term -> gradQ-div(sigma) 
                    rRightHandSideVector(i * BlockSize + mDim) += GaussWeight * mTauOne * DN_DX(i, dim1) * divergence_of_sigma[dim1];
                }
        }




        // // Without Constitutive law
        // Matrix DivSymmetrycGradient = ZeroMatrix(2,number_of_points*mDim);
        // CalculateDivSymmetrycGradient(DivSymmetrycGradient, DDN_DDe);
        // Matrix GradientQ = ZeroMatrix(number_of_points, 2);
        // CalculateGradientQ(GradientQ, DN_DX);
        // Matrix tempMatrix_grad_q_dot_DivSymmetrycGradient = mTauOne * 2.0 * GaussWeight * prod(GradientQ, DivSymmetrycGradient);
        // // RHS corresponding term using "la trampa"
        // Vector u_old(number_of_points * mDim);
        // for (IndexType i = 0; i < number_of_points; ++i)
        // {
        //     u_old[i * mDim]     = r_geometry[i].FastGetSolutionStepValue(VELOCITY_X); // X velocity at previous time step
        //     u_old[i * mDim+1]   = r_geometry[i].FastGetSolutionStepValue(VELOCITY_Y); // Y velocity at previous time step
        // }
        // Vector rhs_temp = prod(tempMatrix_grad_q_dot_DivSymmetrycGradient, u_old);
        // // Assemble the matrix
        // for (IndexType i = 0; i < number_of_points; ++i)
        // {
        //     for (IndexType j = 0; j < number_of_points; ++j)
        //     {
        //         for (IndexType dim2 = 0; dim2 < mDim; ++dim2) 
        //         {
        //             //// Assemble 2nd order stabilization term -> gradQ-div(sigma)
        //             rLeftHandSideMatrix(i * BlockSize + mDim, j * BlockSize + dim2) -= tempMatrix_grad_q_dot_DivSymmetrycGradient(i, j * mDim + dim2);    
        //         }
        //     }
        //     rRightHandSideVector(i * BlockSize + mDim) += rhs_temp(i);
        // }
        

    }


}


void StokesElementTransient::EquationIdVector(Element::EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
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


void StokesElementTransient::GetDofList(
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



int StokesElementTransient::Check(const ProcessInfo& rCurrentProcessInfo) const
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

Element::IntegrationMethod StokesElementTransient::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

void StokesElementTransient::CalculateB(
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


void StokesElementTransient::CalculateB_second_order(
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


void StokesElementTransient::CalculateB_derivative_x(
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

void StokesElementTransient::CalculateB_derivative_y(
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


void StokesElementTransient::CalculateDivSymmetrycGradient(
        Matrix& DivSymmetrycGradient, 
        const ShapeDerivativesType& r_DDN_DDX) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 2; // Only 2 DOFs per node in 2D

    noalias(DivSymmetrycGradient) = ZeroMatrix(2, mat_size);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        DivSymmetrycGradient(0, 2 * i)     = r_DDN_DDX(i, 0) + 1.0 * r_DDN_DDX(i, 2); 
        DivSymmetrycGradient(0, 2 * i + 1) = 1.0 * r_DDN_DDX(i, 1); 
        DivSymmetrycGradient(1, 2 * i)     = 1.0 * r_DDN_DDX(i, 1); 
        DivSymmetrycGradient(1, 2 * i + 1) = 1.0 * r_DDN_DDX(i, 0) + r_DDN_DDX(i, 2); 
    }
}

void StokesElementTransient::CalculateGradientQ(
        Matrix& GradientQ, 
        const ShapeDerivativesType& r_DN_DX) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 1; // Only 1 DOFs per node in 2D

    noalias(GradientQ) = ZeroMatrix(mat_size, 2);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        GradientQ(i, 0)     = r_DN_DX(i, 0); // ∂N_i / ∂x
        GradientQ(i, 1)     = r_DN_DX(i, 1); // ∂N_i / ∂y
    }
}



void StokesElementTransient::CalculateOnIntegrationPoints(
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

void StokesElementTransient::CalculateOnIntegrationPoints(
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

Vector StokesElementTransient::CalculateStressAtIntegrationPoint(
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


Matrix StokesElementTransient::CalculateConstitutiveMatrixAtIntegrationPoint(
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


void StokesElementTransient::GetValuesVector(
        Vector& rValues) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 2;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        const array_1d<double, 2 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        IndexType index = i * 2;

        rValues[index] = velocity[0];
        rValues[index + 1] = velocity[1];
    }
}


void StokesElementTransient::GetFirstDerivativesVector(Vector &rValues, int Step) const
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


void StokesElementTransient::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
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

    // SetValue the integraion_weigth in order to do the postprocess in Python
    SetValue(INTEGRATION_WEIGHT, integration_points[0].Weight() );

    array_1d<double,3> BodyForce = ZeroVector(3);
    BodyForce = this->GetValue(BODY_FORCE);
    // Set value for the BodyForce_old for next time-step
    SetValue(FORCE, BodyForce);

    // Set value to velocity and pressure for next time-step
    SetValue(VELOCITY_X, rOutput_vel_x);
    SetValue(VELOCITY_Y, rOutput_vel_y);
    SetValue(PRESSURE, rOutput_pressure);

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

    double yielded_state = 0.0;
    mpConstitutiveLaw->CalculateValue(Values, STRAIN_ENERGY, yielded_state); // Usa STRAIN_ENERGY come proxy
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    #pragma omp critical
    {
        std::ofstream output_file("txt_files/output_results_GPs.txt", std::ios::app);
        if (output_file.is_open()) {
        output_file << std::scientific << std::setprecision(14); // Set precision to 10^-14
        output_file << rOutput_vel_x << " " << rOutput_vel_y << " " << rOutput_pressure << " " 
                    << r_geometry.Center().X() << " " << r_geometry.Center().Y() << " " << r_geometry.Center().Z() << " " << integration_points[0].Weight()  << std::endl;
        output_file.close();
        }   

        // Same but saving the yielded stress value
        std::ofstream yield_output_file("txt_files/yielded_states_GPs.txt", std::ios::app);
        if (yield_output_file.is_open()) {
            yield_output_file << std::scientific << std::setprecision(14); // Precisione 10^-14
            yield_output_file << yielded_state << " "
                              << r_geometry.Center().X() << " " << r_geometry.Center().Y() << " " << r_geometry.Center().Z() << std::endl;
            yield_output_file.close();
        }

    }  
    
}

void StokesElementTransient::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
}

/// Reads in a json formatted file and returns its KratosParameters instance.
Parameters StokesElementTransient::ReadParamatersFile(
    const std::string& rDataFileName) const
{
    std::ifstream infile(rDataFileName);

    std::stringstream buffer;
    buffer << infile.rdbuf();

    return Parameters(buffer.str());
};



} // Namespace Kratos



    // for (IndexType i = 0; i < number_of_points; ++i)
    // {
    //     for (IndexType j = 0; j < number_of_points; ++j)
    //     {
    //         // diffusion term
    //         rLeftHandSideMatrix(i * BlockSize + 0, j * BlockSize + 2) = 0.0; 
    //         rLeftHandSideMatrix(i * BlockSize + 1, j * BlockSize + 2) = 0.0; 

    //         for (IndexType dim1 = 0; dim1 < 3; ++dim1)
    //         {
    //             rLeftHandSideMatrix(i * BlockSize + 2, j * BlockSize + dim1) = 0.0; 
    //         }
            
    //     }
    //     rRightHandSideVector(i * BlockSize + 2) = 0.0; 
    //     rLeftHandSideMatrix(i * BlockSize + 2, i * BlockSize + 2) = 1.0; 
    // }