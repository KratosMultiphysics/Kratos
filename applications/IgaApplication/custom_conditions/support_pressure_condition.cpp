//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//                  
//

// System includes

// External includes

// Project includes
#include "custom_conditions/support_pressure_condition.h"

namespace Kratos
{

void SupportPressureCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMemberVariables();
}

void SupportPressureCondition::CalculateLocalSystem(
MatrixType& rLeftHandSideMatrix,
VectorType& rRightHandSideVector,
const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType mat_size = number_of_nodes * (mDim+1);

    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size);
    noalias(rRightHandSideVector) = ZeroVector(mat_size);

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void SupportPressureCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType mat_size = GetGeometry().size() * (mDim+1);

    if (rLeftHandSideMatrix.size1() != mat_size && rLeftHandSideMatrix.size2())
        rLeftHandSideMatrix.resize(mat_size, mat_size);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);
}

void SupportPressureCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    Matrix DN_DX(number_of_nodes,mDim);

    const SizeType mat_size = number_of_nodes * (mDim+1);
    
    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS
    
    // Compute the normals
    array_1d<double, 3> normal_parameter_space = -r_geometry.Normal(0, GetIntegrationMethod());
    if (mDim == 3) {
        r_geometry.Calculate(NORMAL, normal_parameter_space);
    }
    normal_parameter_space = normal_parameter_space / MathUtils<double>::Norm(normal_parameter_space);

    // Matrix DB = prod(r_D,B);
    const Matrix& N = r_geometry.ShapeFunctionsValues();
    // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    noalias(DN_DX) = r_DN_De[0]; // prod(DN_De[point_number],InvJ0);

    double p_D = this->GetValue(PRESSURE);

    double current_time = rCurrentProcessInfo[TIME];

    for (IndexType j = 0; j < number_of_nodes; j++) {
        
        //// integration by parts PRESSURE < v dot n , p >
        for (IndexType idim = 0; idim < mDim; idim++) {
            rRightHandSideVector((mDim+1)*j+idim) -= p_D * ( N(0,j) * normal_parameter_space[idim] ) * mIntegrationWeight;
        }
        
        Vector t_N = ZeroVector(mDim);
        double x = r_geometry.Center().X(); 
        double y = r_geometry.Center().Y();
        double z = r_geometry.Center().Z();

        Matrix grad_u(mDim, mDim); // grad_u(i,j) = du_i/dx_j

        if (mDim == 2) {
            grad_u(0, 0) = sinh(x) * sinh(y);     // ∂u_x / ∂x
            grad_u(0, 1) = cosh(x) * cosh(y);     // ∂u_x / ∂y
            grad_u(1, 0) = -cosh(x) * cosh(y);    // ∂u_y / ∂x
            grad_u(1, 1) = -sinh(x) * sinh(y);    // ∂u_y / ∂y
        } else {
            // 3D analytical gradient
            noalias(grad_u) = ZeroMatrix(mDim, mDim);

            // grad_u(0,0) =  2.0;  // dux/dx
            // grad_u(0,1) = -1.0;  // dux/dy
            // grad_u(0,2) =  3.0;  // dux/dz
            // grad_u(1,0) =  3.0;  // duy/dx
            // grad_u(1,1) =  2.0;  // duy/dy
            // grad_u(1,2) = -1.0;  // duy/dz
            // grad_u(2,0) = -5.0;  // duz/dx
            // grad_u(2,1) = -4.0;  // duz/dy
            // grad_u(2,2) = -4.0;  // duz/dz

            // // quadratic
            // grad_u(0,0) = -2.0*y;  grad_u(0,1) = -2.0*x;  grad_u(0,2) = 0.0;
            // grad_u(1,0) = -2.0*x;  grad_u(1,1) =  2.0*y + z;  grad_u(1,2) = y - z;
            // grad_u(2,0) =  0.0;    grad_u(2,1) =  0.0;    grad_u(2,2) = -z;

            // // cubic
            // grad_u(0,0) = 2.0*y*z - 2.0*x*y;   grad_u(0,1) = 2.0*x*z - x*x;      grad_u(0,2) = 2.0*x*y;
            // grad_u(1,0) = 0.0;                 grad_u(1,1) = -2.0*y*z;           grad_u(1,2) = -y*y;
            // grad_u(2,0) = 2.0*y*z;             grad_u(2,1) = 2.0*x*z;            grad_u(2,2) = 2.0*x*y;

            // ux = cosh(x)cosh(y)cosh(z)          // uy = sinh(x)*sinh(z-y)          // uz = -sinh(x)sinh(y)cosh(z)
            grad_u(0,0) = sinh(x)*cosh(y)*cosh(z); grad_u(1,0) = cosh(x)*sinh(z-y);   grad_u(2,0) = -cosh(x)*sinh(y)*cosh(z);
            grad_u(0,1) = cosh(x)*sinh(y)*cosh(z); grad_u(1,1) = -sinh(x)*cosh(z-y);  grad_u(2,1) = -sinh(x)*cosh(y)*cosh(z);
            grad_u(0,2) = cosh(x)*cosh(y)*sinh(z); grad_u(1,2) =  sinh(x)*cosh(z-y);  grad_u(2,2) = -sinh(x)*sinh(y)*sinh(z);
        }

        // grad_u(0, 0) = 3.0 * x * x;      // ∂u_x / ∂x
        // grad_u(0, 1) = 0.0;              // ∂u_x / ∂y
        // grad_u(1, 0) = -6.0 * x * y;     // ∂u_y / ∂x
        // grad_u(1, 1) = -3.0 * x * x;     // ∂u_y / ∂y

        // grad_u(0, 0) = 2.0 * x;      // ∂u_x / ∂x
        // grad_u(0, 1) = 0.0;              // ∂u_x / ∂y
        // grad_u(1, 0) = -2.0 * y;     // ∂u_y / ∂x
        // grad_u(1, 1) = -2.0 * x;     // ∂u_y / ∂y

        // grad_u(0, 0) = sinh(x) * sinh(y) * std::exp(-current_time)*current_time*current_time;    // ∂u_x / ∂x
        // grad_u(0, 1) = cosh(x) * cosh(y) * std::exp(-current_time)*current_time*current_time;    // ∂u_x / ∂y
        // grad_u(1, 0) = -cosh(x) * cosh(y) * std::exp(-current_time)*current_time*current_time;    // ∂u_y / ∂x
        // grad_u(1, 1) = -sinh(x) * sinh(y) * std::exp(-current_time)*current_time*current_time;     // ∂u_y / ∂y

        Matrix sym_grad_u(mDim, mDim); // ε(u) = 0.5*(∇u + ∇u^T)
        for (IndexType i = 0; i < mDim; ++i) {
            for (IndexType j = 0; j < mDim; ++j) {
                sym_grad_u(i,j) = 0.5 * (grad_u(i,j) + grad_u(j,i));
            }
        }
        // Now compute stress vector: σ·n = 2ν ε(u)·n
        for (IndexType i = 0; i < mDim; ++i) {
            for (IndexType j = 0; j < mDim; ++j) {
                t_N[i] += 2.0 * sym_grad_u(i,j) * normal_parameter_space[j];
            }
        }

        // -------------------------------------------------------CHANNEL-------------------------------------------------------
        t_N = ZeroVector(mDim); 
        // -------------------------------------------------------CHANNEL-------------------------------------------------------

        for (IndexType idim = 0; idim < mDim; idim++) {
            rRightHandSideVector((mDim+1)*j+idim) += N(0,j) * t_N[idim] * mIntegrationWeight;
        }
        
    }
    KRATOS_CATCH("")
}


void SupportPressureCondition::InitializeMemberVariables()
{
    // Compute class memeber variables
    const auto& r_geometry = this->GetGeometry();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();
    
    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
    double h = std::min(mesh_size_uv[0], mesh_size_uv[1]);
    if (mDim == 3) {h = std::min(h,  mesh_size_uv[2]);}
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }

    // Integration
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints();
    GeometryType::JacobiansType J0;
    Matrix InvJ0(mDim,mDim);
    r_geometry.Jacobian(J0,r_geometry.GetDefaultIntegrationMethod());
    Matrix Jacobian = ZeroMatrix(3,3);
    double DetJ0;
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);
    Jacobian(2,2) = 1.0; // 2D case
    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);
    array_1d<double, 3> tangent_parameter_space;
    r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
    Vector add_factor = prod(Jacobian, tangent_parameter_space);
    add_factor[2] = 0.0; // 2D case
    DetJ0 = norm_2(add_factor);

    if (mDim == 3) {
        Matrix tangent_matrix;
        r_geometry.Calculate(LOCAL_TANGENT_MATRIX, tangent_matrix);  // 3x2

        array_1d<double,3> t1, t2;
        for (std::size_t i = 0; i < 3; ++i) {
            t1[i] = tangent_matrix(i, 0);
            t2[i] = tangent_matrix(i, 1);
        }
        // Cross product of the two tangents
        array_1d<double, 3> det_vector = MathUtils<double>::CrossProduct(t1, t2);
        // Norm gives the surface integration factor
        DetJ0 = norm_2(det_vector);
    }
    mIntegrationWeight = r_integration_points[0].Weight() * std::abs(DetJ0);
}

void SupportPressureCondition::CalculateB(
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


void SupportPressureCondition::EquationIdVector(EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType number_of_control_points = GetGeometry().size();
    const unsigned int LocalSize = (mDim + 1) * number_of_control_points;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < number_of_control_points; i++)
    {
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
        rResult[Index++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
        if (mDim > 2) rResult[Index++] = rGeom[i].GetDof(VELOCITY_Z).EquationId();
        rResult[Index++] = rGeom[i].GetDof(PRESSURE).EquationId();
    }
}


void SupportPressureCondition::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const SizeType number_of_control_points = GetGeometry().size();

    rElementalDofList.resize(0);
    rElementalDofList.reserve((mDim + 1) * number_of_control_points);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
        if (mDim > 2) rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(PRESSURE));
    }

    KRATOS_CATCH("")
};


void SupportPressureCondition::GetSolutionCoefficientVector(
        Vector& rValues) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * mDim;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        const array_1d<double, 3>& velocity = GetGeometry()[i].GetSolutionStepValue(VELOCITY);
        IndexType index = i * mDim;

        for (IndexType d = 0; d < mDim; ++d) {
            rValues[index + d] = velocity[d];
        }
    }
}


void SupportPressureCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
#pragma omp critical
{   
    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLaw->InitializeMaterial( r_properties, r_geometry, row(N_values , 0 ));

    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;


    const auto& r_geom = GetGeometry();
    const SizeType n = r_geom.size();

    // 1) Dati di quadratura (GI_GAUSS_1) e funzioni di forma
    const auto& r_ip   = r_geom.IntegrationPoints(this->GetIntegrationMethod());
    const auto& r_dnde = r_geom.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const Matrix& H    = r_geom.ShapeFunctionsValues(this->GetIntegrationMethod());
    const Matrix& DN_DX = r_dnde[0];

    // 2) Normal and |J| for boundary measure
    array_1d<double,3> n_par = -r_geom.Normal(0, this->GetIntegrationMethod());
    if (mDim == 3) {
        r_geom.Calculate(NORMAL, n_par);
    }
    n_par /= MathUtils<double>::Norm(n_par);

    GeometryType::JacobiansType J0;
    r_geom.Jacobian(J0, this->GetIntegrationMethod());
    Matrix J = ZeroMatrix(3,3);
    J(0,0)=J0[0](0,0); J(0,1)=J0[0](0,1);
    J(1,0)=J0[0](1,0); J(1,1)=J0[0](1,1);
    J(2,2)=1.0;

    array_1d<double,3> t_par;
    r_geom.Calculate(LOCAL_TANGENT, t_par);
    Vector det_vec = prod(J, t_par);
    det_vec[2]=0.0;
    double detJ = norm_2(det_vec);

    if (mDim == 3) {
        Matrix tangent_matrix;
        r_geom.Calculate(LOCAL_TANGENT_MATRIX, tangent_matrix);  // 3x2
        array_1d<double,3> t1, t2;
        for (std::size_t i = 0; i < 3; ++i) {
            t1[i] = tangent_matrix(i, 0);
            t2[i] = tangent_matrix(i, 1);
        }
        array_1d<double, 3> det_vector = MathUtils<double>::CrossProduct(t1, t2);
        detJ = norm_2(det_vector);
    }
    const double w = r_ip[0].Weight() * std::abs(detJ);

    // 3) Costitutivo al GP (B*u -> strain-rate), stress (solo viscoso)
    const SizeType strain_size = (mDim == 3) ? 6 : 3;
    Matrix B(strain_size, n*mDim, 0.0);
    CalculateB(B, DN_DX);

    ConstitutiveLaw::Parameters Values(r_geom, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables cv(strain_size);
    Vector u_coeff(n*mDim);
    GetSolutionCoefficientVector(u_coeff);
    Vector strain = prod(B, u_coeff);

    Values.SetStrainVector(strain);
    Values.SetStressVector(cv.StressVector);
    Values.SetConstitutiveMatrix(cv.D);
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);

    const Vector& sigma_voigt = Values.GetStressVector(); // τ in Voigt
    // Build symmetric stress matrix from Voigt.
    Matrix tau = ZeroMatrix(mDim, mDim);
    if (mDim == 2) {
        tau(0,0)=sigma_voigt[0]; tau(1,1)=sigma_voigt[1];
        tau(0,1)=sigma_voigt[2]; tau(1,0)=sigma_voigt[2];
    } else {
        // 3D Voigt order: xx, yy, zz, xy, yz, xz.
        tau(0,0)=sigma_voigt[0]; tau(1,1)=sigma_voigt[1]; tau(2,2)=sigma_voigt[2];
        tau(0,1)=sigma_voigt[3]; tau(1,0)=sigma_voigt[3];
        tau(1,2)=sigma_voigt[4]; tau(2,1)=sigma_voigt[4];
        tau(0,2)=sigma_voigt[5]; tau(2,0)=sigma_voigt[5];
    }

    // 4) Pressione al GP (interpolazione)
    double p_gp = 0.0;
    for (SizeType a=0; a<n; ++a)
        p_gp += H(0,a) * r_geom[a].GetSolutionStepValue(PRESSURE);

    // 5) Normal (size mDim)
    Vector n_vec = ZeroVector(mDim);
    for (IndexType d = 0; d < mDim; ++d) {
        n_vec[d] = n_par[d];
    }

    // 6) Tractions
    Vector t_visc = ZeroVector(mDim);
    Vector t_pres = ZeroVector(mDim);
    for (IndexType i = 0; i < mDim; ++i) {
        for (IndexType j = 0; j < mDim; ++j) {
            t_visc[i] += tau(i,j) * n_vec[j];
        }
        t_pres[i] = -p_gp * n_vec[i];
    }
    Vector t_tot = t_visc + t_pres;

    // 7) Coordinate del GP (qui uso il centro geometrico della condizione)
    const auto center = r_geom.Center();

    // 8) Store results for post-process (1 GP -> 1 row).
    // 2D: [0] w, [1..2] f_tot, [3..4] f_visc, [5..6] f_pres, [7..8] n, [9..10] x,y
    // 3D: [0] w, [1..3] f_tot, [4..6] f_visc, [7..9] f_pres, [10..12] n, [13..15] x,y,z
    const std::size_t results_cols = (mDim == 3) ? 16 : 11;
    Matrix results(1, results_cols, 0.0);
    results(0,0) = w;
    for (IndexType d = 0; d < mDim; ++d) {
        results(0, 1 + d) = t_tot[d];
        results(0, 1 + mDim + d) = t_visc[d];
        results(0, 1 + 2 * mDim + d) = t_pres[d];
        results(0, 1 + 3 * mDim + d) = n_vec[d];
    }
    results(0, 1 + 4 * mDim + 0) = center.X();
    results(0, 1 + 4 * mDim + 1) = center.Y();
    if (mDim == 3) {
        results(0, 1 + 4 * mDim + 2) = center.Z();
    }

    // Usa una chiave diversa da quella del cilindro per non confonderli
    this->SetValue(RESULTS_ON_TRUE_BOUNDARY, results);
}
}

} // Namespace Kratos
