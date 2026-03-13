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
    const ShapeFunctionsType& rN = row(N_gausspoint,0);
    const ShapeDerivativesType& DN_DX = DN_De[0];
    const double GaussWeight = integration_points[0].Weight();

    const SizeType strain_size = (mDim == 3) ? 6 : 3;
    // Calculate the B matrix
    Matrix B = ZeroMatrix(strain_size, number_of_points * mDim);
    CalculateB(B, DN_DX);

    // constitutive law
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables(strain_size);
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

    // Calculate stabilization constants
    CalculateTau(mu_effective);

    // Add velocity terms in momentum equation
    this->AddMomentumTerms(rLeftHandSideMatrix,rRightHandSideVector,body_force,mTauTwo,rN,DN_DX,GaussWeight, r_D, r_stress_vector);

    // Add velocity-pressure terms
    this->AddContinuityTerms(rLeftHandSideMatrix,rRightHandSideVector,body_force,mTauOne,rN,DN_DX,GaussWeight);

    // Add Second-Order stabilization terms from VMS
    this->AddSecondOrderStabilizationTerms(rLeftHandSideMatrix,rRightHandSideVector,body_force,mTauOne,rN,DN_DX,GaussWeight,r_D);

    // // Add corresponding RHS: -K * [u;p] (previous iteration values)
    // {
    //     VectorType current_values;
    //     GetFirstDerivativesVector(current_values, 0); // VELOCITY_X, VELOCITY_Y, PRESSURE per node
    //     noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, current_values);
    // }

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

void StokesElement::CalculateTau(double MuEffective)
{   
    // // Estimate element size
    double h = ElementSize();

    // The stabilization parameters from: 
    // https://www.researchgate.net/profile/Ramon-Codina-2/publication/222977355

    mTauOne = std::pow(h, 2) / ( 4.0 * MuEffective ); 

    mTauTwo = MuEffective;
}

double StokesElement::ElementSize()
{
    const double domain_size = this->GetGeometry().DomainSize();
    if (mDim == 3) {
        return 1.240700982 * std::cbrt(domain_size);
    }
    return 1.128379167 * std::sqrt(domain_size);
}

void StokesElement::AddMomentumTerms(MatrixType &rLHS,
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
    auto &r_geometry = GetGeometry();

    const SizeType strain_size = (mDim == 3) ? 6 : 3;
    Matrix B = ZeroMatrix(strain_size, number_of_nodes * mDim);
    CalculateB(B, rDN_DX);
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
    // RHS corresponding term
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
            // Stabilization
            for (unsigned int m = 0; m < mDim; ++m) {
                for (unsigned int n = 0; n < mDim; ++n) {
                    rLHS(first_row+m,first_col+n) += TauTwo * Weight * rDN_DX(i,m) * rDN_DX(j,n);
                }
            }

            // Update column index
            first_col += block_size;
        }
        // RHS corresponding term
        double div_u_current = 0.0;
        for(unsigned int j = 0; j < number_of_nodes; ++j) {
            const auto& r_velocity = r_geometry[j].FastGetSolutionStepValue(VELOCITY);
            for (unsigned int d = 0; d < mDim; ++d) {
                div_u_current += r_velocity[d] * rDN_DX(j, d);
            }
        }
        for (unsigned int m = 0; m < mDim; ++m) {
            rRHS(first_row+m) -= TauTwo * Weight * rDN_DX(i,m) * div_u_current ;
        }

        // Update matrix indices
        first_row += block_size;
        first_col = 0;
    }
}


void StokesElement::AddContinuityTerms(MatrixType &rLHS,
        VectorType &rRHS,
        const array_1d<double,3> &rBodyForce,
        const double TauOne,
        const ShapeFunctionsType &rN,
        const ShapeDerivativesType &rDN_DX,
        const double Weight)
{
    const unsigned int number_of_nodes = this->GetGeometry().PointsNumber();
    const unsigned int block_size = mDim+1;
    auto &r_geometry = GetGeometry();

    unsigned int first_row = 0;
    unsigned int first_col = 0;

    double pressure_previous_iteration = 0.0;
    double div_u_current = 0.0;
    Vector grad_p_previous_iteration = ZeroVector(mDim);
    for(unsigned int j = 0; j < number_of_nodes; ++j) {
        const auto& r_velocity = r_geometry[j].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int d = 0; d < mDim; ++d) {
            div_u_current += r_velocity[d] * rDN_DX(j, d);
            grad_p_previous_iteration[d] += r_geometry[j].GetSolutionStepValue(PRESSURE) * rDN_DX(j, d);
        }
        pressure_previous_iteration += r_geometry[j].GetSolutionStepValue(PRESSURE) * rN[j];
    }

    double div_term = 0.0;

    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        double Qi = 0.0;
        for (unsigned int d = 0; d < mDim; ++d) {
            Qi += rDN_DX(i,d) * rBodyForce[d];
        }

        rRHS[first_row + mDim] += Weight * TauOne * Qi;

        for (unsigned int j = 0; j < number_of_nodes; ++j)
        {
            double l_ij = 0.0;
            for (unsigned int d = 0; d < mDim; ++d)
            {
                l_ij += rDN_DX(i,d) * rDN_DX(j,d);
                div_term = Weight * rN[i] * rDN_DX(j,d);
                rLHS(first_row+mDim,first_col + d) += div_term; // Divergence term
                rLHS(first_col + d,first_row+mDim) -= div_term; // Gradient term
            }
            // Stabilization (grad p, grad q)
            rLHS(first_row+mDim,first_col+mDim) += Weight * TauOne * l_ij;

            // Update column index
            first_col += block_size;
        }

        // --- RHS corresponding term ---
        rRHS(first_row+mDim) -= Weight * rN[i] * div_u_current ; // q div(u)
        for (unsigned int m = 0; m < mDim; ++m) {
            rRHS(first_row+m) += Weight * rDN_DX(i,m) * pressure_previous_iteration ;
            // Stabilization
            rRHS(first_row+mDim) -= Weight * TauOne * rDN_DX(i,m) * grad_p_previous_iteration[m];
        }

        // Update matrix indices
        first_col = 0;
        first_row += block_size;
    }
}


// void StokesElement::AddSecondOrderStabilizationTerms(MatrixType &rLeftHandSideMatrix,
//         VectorType &rRightHandSideVector,
//         const array_1d<double,3> &rBodyForce,
//         const double TauOne,
//         const ShapeFunctionsType &rN,
//         const ShapeDerivativesType &rDN_DX,
//         const double GaussWeight,
//         const Matrix& rD)
// {
//     const unsigned int number_of_points = this->GetGeometry().PointsNumber();
//     const unsigned int block_size = mDim+1;
//     auto &r_geometry = GetGeometry();

//     // Add third order term of the stabilization
//     // Second-order stabilization is only implemented for 2D.
//     if (mBasisFunctionsOrder > 1 && mDim == 2) {

//         Matrix B = ZeroMatrix(3,number_of_points*mDim);
//         CalculateB(B, rDN_DX);
//         const Matrix& DDN_DDe = GetGeometry().ShapeFunctionDerivatives(2, 0, GetGeometry().GetDefaultIntegrationMethod());
        
//         Matrix B_derivative_x;
//         Matrix B_derivative_y;
//         CalculateBDerivativeDx(B_derivative_x, DDN_DDe);
//         CalculateBDerivativeDy(B_derivative_y, DDN_DDe);

//         // Computed using the Gauss Points in the same knot span
//         Vector divergence_of_sigma = this->GetValue(RECOVERED_STRESS);

//         // initilize the div(sigma) matrix 2x(2n)
//         Vector div_sigma_1 = ZeroVector(mDim*number_of_points);
//         Vector div_sigma_2 = ZeroVector(mDim*number_of_points);

//         Vector r_D_0(3); Vector r_D_1(3); Vector r_D_2(3);
//         for (std::size_t i = 0; i < 3; ++i) {
//             r_D_0(i) = rD(0, i); 
//             r_D_1(i) = rD(1, i); 
//             r_D_2(i) = rD(2, i);
//         }
//         div_sigma_1 = prod(trans(B_derivative_x), r_D_0) + prod(trans(B_derivative_y), r_D_2);
//         div_sigma_2 = prod(trans(B_derivative_y), r_D_1) + prod(trans(B_derivative_x), r_D_2);
        
//         // NEW FORMULATION___________________________________________________________________________________________
//         Matrix div_sigma = ZeroMatrix(2, mDim * number_of_points);
//         for (std::size_t i = 0; i < mDim * number_of_points; ++i) {
//             div_sigma(0, i) = div_sigma_1(i);
//             div_sigma(1, i) = div_sigma_2(i);
//         }

//         Vector DN_DX_vector_1 = ZeroVector(number_of_points);
//         Vector DN_DX_vector_2 = ZeroVector(number_of_points);
//         for (std::size_t i = 0; i < number_of_points; ++i) {
//             DN_DX_vector_1(i) = rDN_DX(i, 0); 
//             DN_DX_vector_2(i) = rDN_DX(i, 1); 
//         }
//         Matrix tempMatrix_pressure_term = TauOne * GaussWeight * (outer_prod(div_sigma_1, DN_DX_vector_1) + outer_prod(div_sigma_2, DN_DX_vector_2));

//         for (IndexType i = 0; i < number_of_points; ++i)
//         {
//             for (IndexType j = 0; j < number_of_points; ++j)
//             {
//                 for (IndexType dim2 = 0; dim2 < mDim; ++dim2) 
//                 {
//                     // Assemble 2nd order stabilization term -> gradQ-sigma // this one
//                     rLeftHandSideMatrix(j * block_size + mDim, i * block_size + dim2) -= tempMatrix_pressure_term(i * mDim + dim2, j);
//                 }
//             }
//         }

//         // --- RHS corresponding term ---
//         Vector velocity_previous = ZeroVector(mDim*number_of_points);
//         Vector pressure_previous = ZeroVector(number_of_points);
//         IndexType index = 0;
//         for (IndexType i = 0; i < number_of_points; ++i) { 
//             const auto& r_velocity = r_geometry[i].GetSolutionStepValue(VELOCITY);
//             for (IndexType d = 0; d < mDim; ++d) {
//                 velocity_previous[index++] = r_velocity[d];
//             }
//             pressure_previous[i] = r_geometry[i].GetSolutionStepValue(PRESSURE);
//         }

//         // --- RHS corresponding term ---
//         for (IndexType i = 0; i < number_of_points; ++i) {   
//             for (IndexType dim1 = 0; dim1 < mDim; ++dim1) {
//                 // Assemble 2nd order stabilization term -> gradQ-sigma_u 
//                 rRightHandSideVector(i * block_size + mDim) += GaussWeight * TauOne * rDN_DX(i,dim1) * divergence_of_sigma[dim1];
//             }
//         }
//     } 
// }

void StokesElement::AddSecondOrderStabilizationTerms(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const array_1d<double,3> &rBodyForce,
    const double TauOne,
    const ShapeFunctionsType &rN,
    const ShapeDerivativesType &rDN_DX,
    const double GaussWeight,
    const Matrix& rD)
{
    const unsigned int number_of_points = this->GetGeometry().PointsNumber();
    const unsigned int block_size = mDim + 1;
    auto &r_geometry = GetGeometry();

    // Second-order stabilization: implement for 2D and 3D when basis order > 1
    if (mBasisFunctionsOrder > 1) {
        // ---------------------------------------------------------------------
        // RECOVERED_STRESS stores the recovered divergence of viscous stress:
        // 2D: size 2 -> [div_tau_x, div_tau_y]
        // 3D: size 3 -> [div_tau_x, div_tau_y, div_tau_z]

        // ---------------------------------------------------------------------
        const Vector divergence_of_sigma = this->GetValue(RECOVERED_STRESS);
        // KRATOS_WATCH(divergence_of_sigma)

        KRATOS_ERROR_IF(divergence_of_sigma.size() != mDim)
            << "RECOVERED_STRESS must store div(tau) of size " << mDim
            << " but has size " << divergence_of_sigma.size() << std::endl;

        // Get 2nd derivatives of shape functions in physical space (as provided by the geometry)
        // For 2D, it is expected to contain 3 components per shape function: [xx, xy, yy] (or equivalent)
        // For 3D, it is expected to contain 6 components per shape function: [xx, xy, xz, yy, yz, zz] (or equivalent)
        const Matrix& DDN_DDe = GetGeometry().ShapeFunctionDerivatives(
            2, 0, GetGeometry().GetDefaultIntegrationMethod());

        // Build derivative B-matrices: B_derivative_x/y/z are of size (strain_size x (mDim*n))
        Matrix B_derivative_x;
        Matrix B_derivative_y;
        Matrix B_derivative_z;

        if (mDim == 2) {
            CalculateBDerivativeDx(B_derivative_x, DDN_DDe);
            CalculateBDerivativeDy(B_derivative_y, DDN_DDe);
        } else { // mDim == 3
            CalculateBDerivativeDx3D(B_derivative_x, DDN_DDe);
            CalculateBDerivativeDy3D(B_derivative_y, DDN_DDe);
            CalculateBDerivativeDz3D(B_derivative_z, DDN_DDe);
        }

        // ---------------------------------------------------------------------
        // Build div_sigma_k vectors mapping velocity dofs -> (div tau)_k at the GP
        // Each div_sigma_k has size (mDim * number_of_points).
        // ---------------------------------------------------------------------
        Vector div_sigma_1 = ZeroVector(mDim * number_of_points);
        Vector div_sigma_2 = ZeroVector(mDim * number_of_points);
        Vector div_sigma_3; // only used in 3D

        const unsigned int strain_size = (mDim == 3) ? 6 : 3;

        if (mDim == 2) {
            // 2D Voigt: [xx, yy, xy]
            Vector r_D_0(3), r_D_1(3), r_D_2(3);
            for (std::size_t i = 0; i < 3; ++i) {
                r_D_0(i) = rD(0, i);
                r_D_1(i) = rD(1, i);
                r_D_2(i) = rD(2, i);
            }

            // div_tau_x = d/dx(tau_xx) + d/dy(tau_xy)
            // div_tau_y = d/dx(tau_xy) + d/dy(tau_yy)
            div_sigma_1 = prod(trans(B_derivative_x), r_D_0) + prod(trans(B_derivative_y), r_D_2);
            div_sigma_2 = prod(trans(B_derivative_x), r_D_2) + prod(trans(B_derivative_y), r_D_1);

        } else {
            // 3D Voigt: [xx, yy, zz, xy, yz, xz]
            // div_tau_x = d/dx(tau_xx) + d/dy(tau_xy) + d/dz(tau_xz)
            // div_tau_y = d/dx(tau_xy) + d/dy(tau_yy) + d/dz(tau_yz)
            // div_tau_z = d/dx(tau_xz) + d/dy(tau_yz) + d/dz(tau_zz)

            div_sigma_3 = ZeroVector(mDim * number_of_points);

            // Extract the rows of D corresponding to each stress component
            Vector r_D_xx(6), r_D_yy(6), r_D_zz(6), r_D_xy(6), r_D_yz(6), r_D_xz(6);
            for (std::size_t i = 0; i < 6; ++i) {
                r_D_xx(i) = rD(0, i);
                r_D_yy(i) = rD(1, i);
                r_D_zz(i) = rD(2, i);
                r_D_xy(i) = rD(3, i);
                r_D_yz(i) = rD(4, i);
                r_D_xz(i) = rD(5, i);
            }

            // Note: B_derivative_x/y/z map u_dofs -> d/dx(eps_voigt), d/dy(eps_voigt), d/dz(eps_voigt)
            // Then (d/dx tau) = (B'_x)^T * D_row, etc. assembled per div component.

            // div_tau_x
            div_sigma_1 =
                prod(trans(B_derivative_x), r_D_xx) +
                prod(trans(B_derivative_y), r_D_xy) +
                prod(trans(B_derivative_z), r_D_xz);

            // div_tau_y
            div_sigma_2 =
                prod(trans(B_derivative_x), r_D_xy) +
                prod(trans(B_derivative_y), r_D_yy) +
                prod(trans(B_derivative_z), r_D_yz);

            // div_tau_z
            div_sigma_3 =
                prod(trans(B_derivative_x), r_D_xz) +
                prod(trans(B_derivative_y), r_D_yz) +
                prod(trans(B_derivative_z), r_D_zz);
        }

        // ---------------------------------------------------------------------
        // Build DN_DX vectors (one per spatial direction)
        // ---------------------------------------------------------------------
        Vector DN_DX_vector_1 = ZeroVector(number_of_points);
        Vector DN_DX_vector_2 = ZeroVector(number_of_points);
        Vector DN_DX_vector_3; // 3D only

        for (std::size_t i = 0; i < number_of_points; ++i) {
            DN_DX_vector_1(i) = rDN_DX(i, 0);
            DN_DX_vector_2(i) = rDN_DX(i, 1);
        }
        if (mDim == 3) {
            DN_DX_vector_3 = ZeroVector(number_of_points);
            for (std::size_t i = 0; i < number_of_points; ++i) {
                DN_DX_vector_3(i) = rDN_DX(i, 2);
            }
        }

        // tempMatrix_pressure_term has size (mDim*n) x n
        Matrix tempMatrix_pressure_term;
        if (mDim == 2) {
            tempMatrix_pressure_term =
                TauOne * GaussWeight *
                (outer_prod(div_sigma_1, DN_DX_vector_1) +
                 outer_prod(div_sigma_2, DN_DX_vector_2));
        } else {
            tempMatrix_pressure_term =
                TauOne * GaussWeight *
                (outer_prod(div_sigma_1, DN_DX_vector_1) +
                 outer_prod(div_sigma_2, DN_DX_vector_2) +
                 outer_prod(div_sigma_3, DN_DX_vector_3));
        }

        // ---------------------------------------------------------------------
        // Assemble LHS: q-row (pressure dof) vs u-columns (velocity dofs)
        // ---------------------------------------------------------------------
        for (IndexType i = 0; i < number_of_points; ++i) {
            for (IndexType j = 0; j < number_of_points; ++j) {
                for (IndexType dim2 = 0; dim2 < mDim; ++dim2) {
                    rLeftHandSideMatrix(j * block_size + mDim, i * block_size + dim2) -=
                        tempMatrix_pressure_term(i * mDim + dim2, j);
                }
            }
        }

        // ---------------------------------------------------------------------
        // Assemble RHS: recovered divergence term
        // rRHS(p_i) += w * TauOne * grad(N_i)·div_tau_rec
        // ---------------------------------------------------------------------
        for (IndexType i = 0; i < number_of_points; ++i) {
            double contrib = 0.0;
            for (IndexType dim1 = 0; dim1 < mDim; ++dim1) {
                contrib += rDN_DX(i, dim1) * divergence_of_sigma[dim1];
            }
            rRightHandSideVector(i * block_size + mDim) += GaussWeight * TauOne * contrib;
        }
    }
}


void StokesElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    auto r_geometry = GetGeometry();
    const unsigned int NumNodes = r_geometry.PointsNumber();
    const unsigned int BlockSize = mDim+1;
    
    const int num_dofs_per_node = NumNodes * (mDim + 1);
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

void StokesElement::GetFirstDerivativesVector(Vector &rValues, int Step) const
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
        if (mDim == 3) {
            rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Z,Step);
        }
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(PRESSURE,Step);
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

    mpConstitutiveLaw->CalculateMaterialResponseCauchy(rValues);
}

void StokesElement::EquationIdVector(Element::EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType number_of_control_points = GetGeometry().size();
    const unsigned int local_size = (mDim + 1) * number_of_control_points;

    if (rResult.size() != local_size)
        rResult.resize(local_size);

    unsigned int Index = 0;
    for (unsigned int i = 0; i < number_of_control_points; i++)
    {
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
        if (mDim > 2) rResult[Index++] = rGeom[i].GetDof(VELOCITY_Z).EquationId();
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
    rElementalDofList.reserve((mDim + 1) * number_of_control_points);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
        if (mDim == 3) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
        }
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
                ( mDim == 3 && this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false ) )
            KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if (mDim == 2)
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
    const SizeType mat_size = number_of_control_points * mDim;
    const SizeType strain_size = (mDim == 3) ? 6 : 3;

    // Resize B matrix to Voigt strain size and appropriate number of columns.
    if (rB.size1() != strain_size || rB.size2() != mat_size)
        rB.resize(strain_size, mat_size);

    noalias(rB) = ZeroMatrix(strain_size, mat_size);

    if (mDim == 2) {
        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            rB(0, 2 * i)     = r_DN_DX(i, 0);
            rB(1, 2 * i + 1) = r_DN_DX(i, 1);
            rB(2, 2 * i)     = r_DN_DX(i, 1);
            rB(2, 2 * i + 1) = r_DN_DX(i, 0);
        }
    } else {
        // 3D small-strain Voigt order: xx, yy, zz, xy, yz, xz.
        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            rB(0, 3 * i)     = r_DN_DX(i, 0);
            rB(1, 3 * i + 1) = r_DN_DX(i, 1);
            rB(2, 3 * i + 2) = r_DN_DX(i, 2);
            rB(3, 3 * i)     = r_DN_DX(i, 1);
            rB(3, 3 * i + 1) = r_DN_DX(i, 0);
            rB(4, 3 * i + 1) = r_DN_DX(i, 2);
            rB(4, 3 * i + 2) = r_DN_DX(i, 1);
            rB(5, 3 * i)     = r_DN_DX(i, 2);
            rB(5, 3 * i + 2) = r_DN_DX(i, 0);
        }
    }
}

void StokesElement::CalculateBDerivativeDx(
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

void StokesElement::CalculateBDerivativeDy(
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

void StokesElement::CalculateBDerivativeDx3D(
    Matrix& BDerivativeDx,
    const ShapeDerivativesType& r_DDN_DDX) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 3; // 3 dofs per node

    // Voigt size = 6 in 3D
    if (BDerivativeDx.size1() != 6 || BDerivativeDx.size2() != mat_size)
        BDerivativeDx.resize(6, mat_size);

    noalias(BDerivativeDx) = ZeroMatrix(6, mat_size);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        const double dxx = r_DDN_DDX(i,0); // d2N/dx2
        const double dxy = r_DDN_DDX(i,1); // d2N/dxdy
        const double dxz = r_DDN_DDX(i,2); // d2N/dxdz

        const IndexType col = 3 * i;

        // ε_xx,x
        BDerivativeDx(0, col    ) = dxx;

        // ε_yy,x
        BDerivativeDx(1, col + 1) = dxy;

        // ε_zz,x
        BDerivativeDx(2, col + 2) = dxz;

        // ε_xy,x
        BDerivativeDx(3, col    ) = dxy;
        BDerivativeDx(3, col + 1) = dxx;

        // ε_yz,x
        BDerivativeDx(4, col + 1) = dxz;
        BDerivativeDx(4, col + 2) = dxy;

        // ε_xz,x
        BDerivativeDx(5, col    ) = dxz;
        BDerivativeDx(5, col + 2) = dxx;
    }
}


void StokesElement::CalculateBDerivativeDy3D(
    Matrix& BDerivativeDy,
    const ShapeDerivativesType& r_DDN_DDX) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 3;

    if (BDerivativeDy.size1() != 6 || BDerivativeDy.size2() != mat_size)
        BDerivativeDy.resize(6, mat_size);

    noalias(BDerivativeDy) = ZeroMatrix(6, mat_size);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        const double dxy = r_DDN_DDX(i,1); // d2N/dxdy
        const double dyy = r_DDN_DDX(i,3); // d2N/dy2
        const double dyz = r_DDN_DDX(i,4); // d2N/dydz

        const IndexType col = 3 * i;

        // ε_xx,y
        BDerivativeDy(0, col    ) = dxy;

        // ε_yy,y
        BDerivativeDy(1, col + 1) = dyy;

        // ε_zz,y
        BDerivativeDy(2, col + 2) = dyz;

        // ε_xy,y
        BDerivativeDy(3, col    ) = dyy;
        BDerivativeDy(3, col + 1) = dxy;

        // ε_yz,y
        BDerivativeDy(4, col + 1) = dyz;
        BDerivativeDy(4, col + 2) = dyy;

        // ε_xz,y
        BDerivativeDy(5, col    ) = dyz;
        BDerivativeDy(5, col + 2) = dxy;
    }
}

void StokesElement::CalculateBDerivativeDz3D(
    Matrix& BDerivativeDz,
    const ShapeDerivativesType& r_DDN_DDX) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 3;

    if (BDerivativeDz.size1() != 6 || BDerivativeDz.size2() != mat_size)
        BDerivativeDz.resize(6, mat_size);

    noalias(BDerivativeDz) = ZeroMatrix(6, mat_size);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        const double dxz = r_DDN_DDX(i,2); // d2N/dxdz
        const double dyz = r_DDN_DDX(i,4); // d2N/dydz
        const double dzz = r_DDN_DDX(i,5); // d2N/dz2

        const IndexType col = 3 * i;

        // ε_xx,z
        BDerivativeDz(0, col    ) = dxz;

        // ε_yy,z
        BDerivativeDz(1, col + 1) = dyz;

        // ε_zz,z
        BDerivativeDz(2, col + 2) = dzz;

        // ε_xy,z
        BDerivativeDz(3, col    ) = dyz;
        BDerivativeDz(3, col + 1) = dxz;

        // ε_yz,z
        BDerivativeDz(4, col + 1) = dzz;
        BDerivativeDz(4, col + 2) = dyz;

        // ε_xz,z
        BDerivativeDz(5, col    ) = dzz;
        BDerivativeDz(5, col + 2) = dxz;
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

        Vector sigma_voigt((mDim == 3) ? 6 : 3);
        sigma_voigt = this->CalculateStressAtIntegrationPoint(rCurrentProcessInfo);
        rOutput[0] = sigma_voigt;
    }
}


Vector StokesElement::CalculateStressAtIntegrationPoint(
    const ProcessInfo& rCurrentProcessInfo)
{
    Vector stress_voigt((mDim == 3) ? 6 : 3);
    
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_points = r_geometry.size();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const ShapeDerivativesType& DN_DX = DN_De[0];

    Matrix B = ZeroMatrix(stress_voigt.size(), number_of_points * mDim);
    CalculateB(B, DN_DX);

    // constitutive law
    ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables constitutive_variables(stress_voigt.size());
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
    const SizeType mat_size = number_of_control_points * mDim;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        const array_1d<double, 3 >& velocity = GetGeometry()[i].GetSolutionStepValue(VELOCITY);
        IndexType index = i * mDim;

        for (IndexType d = 0; d < mDim; ++d) {
            rValues[index + d] = velocity[d];
        }
    }
}

} // Namespace Kratos
