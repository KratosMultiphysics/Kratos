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
#include "custom_elements/stokes_element.h"

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

    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
    
    // SetValue the integraion_weigth in order to do the postprocess in Python
    SetValue(INTEGRATION_WEIGHT, integration_points[0].Weight());
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
    const ShapeFunctionsType& N = row(N_gausspoint,0);
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
        this->EvaluateInPoint(mu_effective,VISCOSITY,N,r_geometry);
    }

    array_1d<double,3> BodyForce = ZeroVector(3);
    BodyForce = this->GetValue(BODY_FORCE);

    const auto& coords = r_geometry[0].Coordinates();
    const double h_element_size = std::min(
        norm_2(coords - r_geometry[1].Coordinates()),
        norm_2(coords - r_geometry[mBasisFunctionsOrder + 2].Coordinates())
    );

    // Calculate stabilization constants
    double TauOne = 0.0;
    double TauTwo = 0.0;
    CalculateTau(TauOne,TauTwo,h_element_size, mu_effective);

    // Add velocity terms in momentum equation
    this->AddMomentumTerms(rLeftHandSideMatrix,rRightHandSideVector,BodyForce,TauTwo,N,DN_DX,GaussWeight, r_D, r_stress_vector);

    // Add velocity-pressure terms
    this->AddContinuityTerms(rLeftHandSideMatrix,rRightHandSideVector,BodyForce,TauOne,N,DN_DX,GaussWeight);

    // Add Second-Order stabilization terms from VMS
    this->AddSecondOrderStabilizationTerms(rLeftHandSideMatrix,rRightHandSideVector,BodyForce,TauOne,N,DN_DX,GaussWeight,r_D, r_stress_vector);

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

void StokesElement::CalculateTau(double &TauOne, double &TauTwo, const double h_element_size, double mu_effective)
{   
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
        const array_1d<double,3> &BodyForce,
        const double TauTwo,
        const ShapeFunctionsType &N,
        const ShapeDerivativesType &DN_DX,
        const double Weight,
        const Matrix& r_D,
        Vector& r_stress_vector)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int block_size = mDim+1;
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
                    rLHS(i * block_size + dim1, j * block_size + dim2) += tempMatrix(i * mDim + dim1, j * mDim + dim2);
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
            rRHS(i * block_size + dim1) -= internalForces(i * mDim + dim1);
        }  
    }

    unsigned int first_row = 0;
    unsigned int first_col = 0;

    for(unsigned int i = 0; i < NumNodes; ++i)
    {
        // Body force
        for(unsigned int d = 0; d < mDim; ++d) {
            rRHS[first_row + d] += Weight * N[i] * BodyForce[d];
        }

        for(unsigned int j = 0; j < NumNodes; ++j)
        {
            // Stabilization
            for (unsigned int m = 0; m < mDim; ++m) {
                for (unsigned int n = 0; n < mDim; ++n) {
                    rLHS(first_row+m,first_col+n) += TauTwo * Weight * DN_DX(i,m) * DN_DX(j,n);
                }
            }

            // Update column index
            first_col += block_size;
        }
        // RHS corresponding term
        double div_u_current = 0.0;
        for(unsigned int j = 0; j < NumNodes; ++j) {
            div_u_current += r_geometry[j].FastGetSolutionStepValue(VELOCITY_X) * DN_DX(j, 0) + 
                             r_geometry[j].FastGetSolutionStepValue(VELOCITY_Y) * DN_DX(j, 1) ;
        } 
        for (unsigned int m = 0; m < mDim; ++m) {
            rRHS(first_row+m) -= TauTwo * Weight * DN_DX(i,m) * div_u_current ;
        }

        // Update matrix indices
        first_row += block_size;
        first_col = 0;
    }
}


void StokesElement::AddContinuityTerms(MatrixType &rLHS,
        VectorType &rRHS,
        const array_1d<double,3> &BodyForce,
        const double TauOne,
        const ShapeFunctionsType &N,
        const ShapeDerivativesType &DN_DX,
        const double Weight)
{
    const unsigned int NumNodes = this->GetGeometry().PointsNumber();
    const unsigned int block_size = mDim+1;
    auto r_geometry = GetGeometry();

    unsigned int first_row = 0;
    unsigned int first_col = 0;

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

    double div_term = 0.0;

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        double Qi = 0.0;
        for (unsigned int d = 0; d < mDim; ++d) {
            Qi += DN_DX(i,d) * BodyForce[d];
        }

        rRHS[first_row + mDim] += Weight * TauOne * Qi;

        for (unsigned int j = 0; j < NumNodes; ++j)
        {
            double l_ij = 0.0;
            for (unsigned int d = 0; d < mDim; ++d)
            {
                l_ij += DN_DX(i,d) * DN_DX(j,d);
                div_term = Weight * N[i] * DN_DX(j,d);
                rLHS(first_row+mDim,first_col + d) += div_term; // Divergence term
                rLHS(first_col + d,first_row+mDim) -= div_term; // Gradient term
            }
            // Stabilization (grad p, grad q)
            rLHS(first_row+mDim,first_col+mDim) += Weight * TauOne * l_ij;

            // Update column index
            first_col += block_size;
        }

        // --- RHS corresponding term ---
        rRHS(first_row+mDim) -= Weight * N[i] * div_u_current ; // q div(u)
        for (unsigned int m = 0; m < mDim; ++m) {
            rRHS(first_row+m) += Weight * DN_DX(i,m) * pressure_previous_iteration ;
            // Stabilization
            rRHS(first_row+mDim) -= Weight * TauOne * DN_DX(i,m) * grad_p_previous_iteration[m];
        }

        // Update matrix indices
        first_col = 0;
        first_row += block_size;
    }
}


void StokesElement::AddSecondOrderStabilizationTerms(MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const array_1d<double,3> &BodyForce,
        const double TauOne,
        const ShapeFunctionsType &N,
        const ShapeDerivativesType &DN_DX,
        const double GaussWeight,
        const Matrix& r_D,
        Vector& r_stress_vector)
{
    const unsigned int number_of_points = this->GetGeometry().PointsNumber();
    const unsigned int block_size = mDim+1;
    auto r_geometry = GetGeometry();

    // Add third order term of the stabilization
    if (mBasisFunctionsOrder > 1) {

        Matrix B = ZeroMatrix(3,number_of_points*mDim);
        CalculateB(B, DN_DX);
        const Matrix DDN_DDe = GetGeometry().ShapeFunctionDerivatives(2, 0, GetGeometry().GetDefaultIntegrationMethod());
        
        Matrix B_derivative_x;
        Matrix B_derivative_y;
        CalculateBDerivativeDx(B_derivative_x, DDN_DDe);
        CalculateBDerivativeDy(B_derivative_y, DDN_DDe);

        // Computed using the Gauss Points in the same knot span
        Vector divergence_of_sigma = this->GetValue(DIVERGENCE_STRESS);

        // initilize the div(sigma) matrix 2x(2n)
        Vector div_sigma_1 = ZeroVector(mDim*number_of_points);
        Vector div_sigma_2 = ZeroVector(mDim*number_of_points);

        Vector r_D_0(3); Vector r_D_1(3); Vector r_D_2(3);
        for (std::size_t i = 0; i < 3; ++i) {
            r_D_0(i) = r_D(0, i); 
            r_D_1(i) = r_D(1, i); 
            r_D_2(i) = r_D(2, i);
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
                    // Assemble 2nd order stabilization term -> gradQ-sigma // this one
                    rLeftHandSideMatrix(j * block_size + mDim, i * block_size + dim2) -= tempMatrix_pressure_term(i * mDim + dim2, j);
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
                // Assemble 2nd order stabilization term -> gradQ-sigma_u 
                rRightHandSideVector(i * block_size + mDim) += GaussWeight * TauOne * DN_DX(i,dim1) * divergence_of_sigma[dim1];
            }
        }
    } 
}

void StokesElement::ApplyConstitutiveLaw(
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

    //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
    //this is ok under the hypothesis that no history dependent behavior is employed
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(rValues);
}

void StokesElement::EquationIdVector(Element::EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
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

void StokesElement::CalculateBDerivativeDx(
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

void StokesElement::CalculateBDerivativeDy(
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


Vector StokesElement::CalculateStressAtIntegrationPoint(
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


void StokesElement::GetSolutionCoefficientVector(
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