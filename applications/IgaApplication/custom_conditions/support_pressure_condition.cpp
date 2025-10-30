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
    array_1d<double, 3> normal_parameter_space = - r_geometry.Normal(0, GetIntegrationMethod());
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
        
        Vector t_N = ZeroVector(mDim); // Initialize t_N with zero vector of size 2
        double x = r_geometry.Center().X(); 
        double y = r_geometry.Center().Y();

        Matrix grad_u(mDim, mDim); // grad_u(i,j) = du_i/dx_j

        grad_u(0, 0) = sinh(x) * sinh(y);     // ∂u_x / ∂x
        grad_u(0, 1) = cosh(x) * cosh(y);     // ∂u_x / ∂y
        grad_u(1, 0) = -cosh(x) * cosh(y);    // ∂u_y / ∂x
        grad_u(1, 1) = -sinh(x) * sinh(y);     // ∂u_y / ∂y

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

    KRATOS_ERROR_IF(mDim == 3) << "SupportPressureCondition is not implemented in 3D. Current dimension: " << mDim << std::endl;
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
    Jacobian(2,2) = 1.0;
    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);
    array_1d<double, 3> tangent_parameter_space;
    r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space !!
    Vector add_factor = prod(Jacobian, tangent_parameter_space);
    add_factor[2] = 0.0;
    DetJ0 = norm_2(add_factor);
    mIntegrationWeight = r_integration_points[0].Weight() * std::abs(DetJ0);
}

void SupportPressureCondition::CalculateB(
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
    rElementalDofList.reserve(3 * number_of_control_points);

    for (IndexType i = 0; i < number_of_control_points; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(PRESSURE));
    }

    KRATOS_CATCH("")
};


void SupportPressureCondition::GetSolutionCoefficientVector(
        Vector& rValues) const
{
    const SizeType number_of_control_points = GetGeometry().size();
    const SizeType mat_size = number_of_control_points * 2;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);

    for (IndexType i = 0; i < number_of_control_points; ++i)
    {
        const array_1d<double, 3>& velocity = GetGeometry()[i].GetSolutionStepValue(VELOCITY);
        IndexType index = i * 2;

        rValues[index] = velocity[0];
        rValues[index + 1] = velocity[1];
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

    // 2) Normale e |J| per misura di linea (2D)
    array_1d<double,3> n_par = - r_geom.Normal(0, this->GetIntegrationMethod());
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
    const double detJ = norm_2(det_vec);
    const double w = r_ip[0].Weight() * std::abs(detJ);

    // 3) Costitutivo al GP (B*u -> strain-rate), stress (solo viscoso)
    Matrix B(3, n*mDim, 0.0);
    CalculateB(B, DN_DX);

    ConstitutiveLaw::Parameters Values(r_geom, GetProperties(), rCurrentProcessInfo);
    ConstitutiveVariables cv(3);
    Vector u_coeff(n*mDim);
    GetSolutionCoefficientVector(u_coeff);
    Vector strain = prod(B, u_coeff);

    Values.SetStrainVector(strain);
    Values.SetStressVector(cv.StressVector);
    Values.SetConstitutiveMatrix(cv.D);
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);

    const Vector& sigma_voigt = Values.GetStressVector(); // τ in Voigt (xx, yy, xy)
    // 2x2 da Voigt
    Matrix tau = ZeroMatrix(2,2);
    tau(0,0)=sigma_voigt[0]; tau(1,1)=sigma_voigt[1];
    tau(0,1)=sigma_voigt[2]; tau(1,0)=sigma_voigt[2];

    // 4) Pressione al GP (interpolazione)
    double p_gp = 0.0;
    for (SizeType a=0; a<n; ++a)
        p_gp += H(0,a) * r_geom[a].GetSolutionStepValue(PRESSURE);

    // 5) Normale (2D)
    array_1d<double,2> n2; n2[0]=n_par[0]; n2[1]=n_par[1];

    // 6) Tractions
    array_1d<double,2> t_visc = ZeroVector(2);
    t_visc[0] = tau(0,0)*n2[0] + tau(0,1)*n2[1];
    t_visc[1] = tau(1,0)*n2[0] + tau(1,1)*n2[1];

    array_1d<double,2> t_pres = ZeroVector(2);
    t_pres[0] = -p_gp * n2[0];
    t_pres[1] = -p_gp * n2[1];

    array_1d<double,2> t_tot = t_visc + t_pres;

    // 7) Coordinate del GP (qui uso il centro geometrico della condizione)
    const auto center = r_geom.Center();

    // 8) Salva risultati per il post-process (1 GP -> 1 riga)
    //    [0] w, [1] fx_tot, [2] fy_tot, [3] fx_visc, [4] fy_visc, [5] fx_pres, [6] fy_pres,
    //    [7] nx, [8] ny, [9] x_gp, [10] y_gp
    Matrix results(1, 11, 0.0);
    results(0,0)  = w;
    results(0,1)  = t_tot[0];
    results(0,2)  = t_tot[1];
    results(0,3)  = t_visc[0];
    results(0,4)  = t_visc[1];
    results(0,5)  = t_pres[0];
    results(0,6)  = t_pres[1];
    results(0,7)  = n2[0];
    results(0,8)  = n2[1];
    results(0,9)  = center.X();
    results(0,10) = center.Y();

    // Usa una chiave diversa da quella del cilindro per non confonderli
    this->SetValue(RESULTS_ON_TRUE_BOUNDARY, results);
}
}

} // Namespace Kratos
