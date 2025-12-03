
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andea Gorgi
//                  
//

// System includes

// External includes

// Project includes
#include "custom_conditions/load_solid_condition.h"

namespace Kratos
{

void LoadSolidCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
    InitializeMemberVariables();
}


void LoadSolidCondition::InitializeMaterial()
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

void LoadSolidCondition::InitializeMemberVariables()
{
    // calculate the integration weight
    const auto& r_geometry = GetGeometry();
    // Initialize the dimension
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    mDim = r_DN_De[0].size2();
    KRATOS_ERROR_IF(mDim != 2) << "LoadSOlidCOndiyion momentarily only supports 2D elements, but the current condition has dimension " << mDim << std::endl;
    // reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    // Initialize Jacobian
    Matrix InvJ0(3,3);
    GeometryType::JacobiansType J0;
    Matrix delta_position;
    CalculateDeltaPositionMatrix(r_geometry, delta_position);
    r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);

    // Compute the normals
    array_1d<double, 3> tangent_parameter_space;
    array_1d<double, 3> normal_physical_space;
    array_1d<double, 3> normal_parameter_space;

    r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space
    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);

    const double abs_tangent_u = std::abs(tangent_parameter_space[0]);
    const double abs_tangent_v = std::abs(tangent_parameter_space[1]);
    mCharacteristicGeometryLength = mesh_size_uv[0];
    if (abs_tangent_v > abs_tangent_u) {
        mCharacteristicGeometryLength = mesh_size_uv[1];
    }
    mCharacteristicGeometryLength *= 0.5;

    double magnitude = std::sqrt(tangent_parameter_space[0] * tangent_parameter_space[0] + tangent_parameter_space[1] * tangent_parameter_space[1]);
    
    normal_parameter_space[0] = + tangent_parameter_space[1] / magnitude;
    normal_parameter_space[1] = - tangent_parameter_space[0] / magnitude;  // By observations on the result of .Calculate(LOCAL_TANGENT)
    normal_parameter_space[2] = 0.0;

    // compute complete jacobian transformation including parameter->physical space transformation
    double DetJ0;
    Matrix Jacobian = ZeroMatrix(3,3);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);
    Jacobian(2,2) = 1.0;

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

    Vector add_factor = prod(Jacobian, tangent_parameter_space); //additional factor to the determinant of the jacobian for the parameter->physical space transformation
    add_factor[2] = 0.0; 
    DetJ0 = norm_2(add_factor);

    // compute the normal of the physical space
    normal_physical_space = prod(trans(InvJ0),normal_parameter_space);
    normal_physical_space[2] = 0.0;
    normal_physical_space /= norm_2(normal_physical_space);
    SetValue(NORMAL, normal_physical_space);

    const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;

    const double integration_weight = r_integration_points[0].Weight() * std::abs(DetJ0) * thickness;

    SetValue(INTEGRATION_WEIGHT, integration_weight);
}

void LoadSolidCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType mat_size = GetGeometry().size() * 2;

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

void LoadSolidCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_geometry.size();

    const SizeType mat_size = number_of_control_points * mDim;

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size,mat_size); //resetting LHS

    KRATOS_CATCH("")
}


void LoadSolidCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_geometry.size();

     // reading integration points and local gradients
    const Matrix& N = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    const SizeType mat_size = number_of_control_points * mDim;
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);

    // resizing as needed the RHS
    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);
    noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS

    //-------------------------------------------------------------------------
    // Initialize DN_DX
    Matrix DN_DX(number_of_control_points,2);
    Matrix InvJ0(3,3);

    // Initialize Jacobian
    GeometryType::JacobiansType J0;
    Matrix delta_position;
    CalculateDeltaPositionMatrix(r_geometry, delta_position);
    r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);

    // compute complete jacobian transformation including parameter->physical space transformation
    double DetJ0;
    Matrix Jacobian = ZeroMatrix(3,3);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);
    Jacobian(2,2) = 1.0;

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

    // retrieve the normal of the physical space
    array_1d<double, 3> normal_physical_space = GetValue(NORMAL);
    
    // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    Matrix sub_inv_jacobian = ZeroMatrix(2,2);
    sub_inv_jacobian(0,0) = InvJ0(0,0);
    sub_inv_jacobian(1,0) = InvJ0(1,0);
    sub_inv_jacobian(0,1) = InvJ0(0,1);
    sub_inv_jacobian(1,1) = InvJ0(1,1);

    Vector g_N = this->GetValue(FORCE); 

    double t = rCurrentProcessInfo[TIME];
    const double x = r_geometry.Center().X();
    const double y = r_geometry.Center().Y();

    Vector rStressVector = ZeroVector(3);
    // rStressVector[0] = ((std::sqrt(3*std::pow(x, 2) + 4*std::pow(2*x + y, 2))*std::fabs(t) <= 3.0/500.0) ? (
    //     67000*t*(2*x + y)
    //  )
    //  : (
    //     (268.0/3.0)*t*(2*x + y)*(250*std::sqrt(3*std::pow(x, 2) + 4*std::pow(2*x + y, 2))*std::fabs(t) + 3)/(std::sqrt(3*std::pow(x, 2) + 4*std::pow(2*x + y, 2))*std::fabs(t))
    //  ));
    //  rStressVector[1] = ((std::sqrt(3*std::pow(x, 2) + 4*std::pow(2*x + y, 2))*std::fabs(t) <= 3.0/500.0) ? (
    //     0
    //  )
    //  : (
    //     (134.0/3.0)*t*(2*x + y)*(500*std::sqrt(3*std::pow(x, 2) + 4*std::pow(2*x + y, 2))*std::fabs(t) - 3)/(std::sqrt(3*std::pow(x, 2) + 4*std::pow(2*x + y, 2))*std::fabs(t))
    //  ));
    //  rStressVector[2] = ((std::sqrt(3*std::pow(x, 2) + 4*std::pow(2*x + y, 2))*std::fabs(t) <= 3.0/500.0) ? (
    //     33500*t*x
    //  )
    //  : (
    //     201*t*x/(std::sqrt(3*std::pow(x, 2) + 4*std::pow(2*x + y, 2))*std::fabs(t))
    //  ));

    // PLASTICITY IN THE MIDDLE
    // rStressVector[0] = ((std::sqrt(21)*M_PI*std::fabs(t*std::sin((1.0/4.0)*M_PI*x)*std::sin((1.0/2.0)*M_PI*y)) <= 3.0/500.0) ? (
    //     33500*M_PI*t*std::sin((1.0/4.0)*M_PI*x)*std::sin((1.0/2.0)*M_PI*y)
    //  )
    //  : (
    //     (134.0/7.0)*t*(-1750*M_PI + std::sqrt(21)/std::fabs(t*std::sin((1.0/4.0)*M_PI*x)*std::sin((1.0/2.0)*M_PI*y)))*std::sin((1.0/4.0)*M_PI*x)*std::sin((1.0/2.0)*M_PI*y)
    //  ));
    //  rStressVector[1] = ((std::sqrt(21)*M_PI*std::fabs(t*std::sin((1.0/4.0)*M_PI*x)*std::sin((1.0/2.0)*M_PI*y)) <= 3.0/500.0) ? (
    //     -134000*M_PI*t*std::sin((1.0/4.0)*M_PI*x)*std::sin((1.0/2.0)*M_PI*y)
    //  )
    //  : (
    //     -67.0/7.0*t*(3500*M_PI + 3*std::sqrt(21)/std::fabs(t*std::sin((1.0/4.0)*M_PI*x)*std::sin((1.0/2.0)*M_PI*y)))*std::sin((1.0/4.0)*M_PI*x)*std::sin((1.0/2.0)*M_PI*y)
    //  ));
    //  rStressVector[2] = 0;

    //PLASTICITY AT THE BOUNDARY
    // rStressVector[0] = ((std::sqrt(21)*M_PI*std::fabs(t*std::sin((1.0/4.0)*M_PI*x)*std::cos((1.0/2.0)*M_PI*y)) <= 3.0/500.0) ? (
    //     33500*M_PI*t*std::sin((1.0/4.0)*M_PI*x)*std::cos((1.0/2.0)*M_PI*y)
    //  )
    //  : (
    //     (134.0/7.0)*t*(-1750*M_PI + std::sqrt(21)/std::fabs(t*std::sin((1.0/4.0)*M_PI*x)*std::cos((1.0/2.0)*M_PI*y)))*std::sin((1.0/4.0)*M_PI*x)*std::cos((1.0/2.0)*M_PI*y)
    //  ));
    //  rStressVector[1] = ((std::sqrt(21)*M_PI*std::fabs(t*std::sin((1.0/4.0)*M_PI*x)*std::cos((1.0/2.0)*M_PI*y)) <= 3.0/500.0) ? (
    //     -134000*M_PI*t*std::sin((1.0/4.0)*M_PI*x)*std::cos((1.0/2.0)*M_PI*y)
    //  )
    //  : (
    //     -67.0/7.0*t*(3500*M_PI + 3*std::sqrt(21)/std::fabs(t*std::sin((1.0/4.0)*M_PI*x)*std::cos((1.0/2.0)*M_PI*y)))*std::sin((1.0/4.0)*M_PI*x)*std::cos((1.0/2.0)*M_PI*y)
    //  ));
    //  rStressVector[2] = 0;

    rStressVector[0] = ((std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))*std::fabs(t) <= 3.0/1000.0) ? (
        134000*t*y*(x - 2)*(y - 2)
     )
     : (
        (268.0/3.0)*t*y*(x - 2)*(y - 2)*(500*std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))*std::fabs(t) + 3)/(std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))*std::fabs(t))
     ));
     rStressVector[1] = ((std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))*std::fabs(t) <= 3.0/1000.0) ? (
        0
     )
     : (
        (134.0/3.0)*t*y*(y - 2)*(-3*x + 1000*(x - 2)*std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))*std::fabs(t) + 6)/(std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))*std::fabs(t))
     ));
     rStressVector[2] = ((std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))*std::fabs(t) <= 3.0/1000.0) ? (
        67000*t*x*(x - 4)*(y - 1)
     )
     : (
        201*t*x*(x - 4)*(y - 1)/(std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))*std::fabs(t))
     ));

    g_N[0] = rStressVector[0]*normal_physical_space[0] + rStressVector[2]*normal_physical_space[1];
    g_N[1] = rStressVector[2]*normal_physical_space[0] + rStressVector[1]*normal_physical_space[1];

    for (IndexType i = 0; i < number_of_control_points; i++) {
        for (IndexType zdim = 0; zdim < 2; zdim++) {
            
            rRightHandSideVector[2*i+zdim] += N(0,i)*g_N[zdim] * integration_weight;

        }
    }
    KRATOS_CATCH("")
}

    int LoadSolidCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
        return 0;
    }

    void LoadSolidCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        if (rResult.size() != 2 * number_of_control_points)
            rResult.resize(2 * number_of_control_points, false);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * 2;
            const auto& r_node = r_geometry[i];
            rResult[index] = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
        }
    }

    void LoadSolidCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(2 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        }
    };


    void LoadSolidCondition::GetSolutionCoefficientVector(
        Vector& rValues) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 2;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (IndexType i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
            IndexType index = i * 2;

            rValues[index] = displacement[0];
            rValues[index + 1] = displacement[1];
        }
    }

    void LoadSolidCondition::CalculateB(
        Matrix& rB, 
        Matrix& r_DN_DX) const
    {
        const SizeType number_of_control_points = GetGeometry().size();
        const SizeType mat_size = number_of_control_points * 2;

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

    void LoadSolidCondition::ApplyConstitutiveLaw(SizeType matSize, Vector& rStrain, ConstitutiveLaw::Parameters& rValues,
                                        ConstitutiveVariables& rConstitutiVariables)
    {
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=rValues.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        KRATOS_ERROR_IF(mCharacteristicGeometryLength <= 0.0)
            << "Characteristic geometry length must be initialized before applying the constitutive law." << std::endl;
        rValues.SetCharacteristicGeometryLength(mCharacteristicGeometryLength);
        
        rValues.SetStrainVector(rStrain);
        rValues.SetStressVector(rConstitutiVariables.StressVector);
        rValues.SetConstitutiveMatrix(rConstitutiVariables.D);

        mpConstitutiveLaw->CalculateMaterialResponse(rValues, ConstitutiveLaw::StressMeasure_Cauchy); 
    }

    void LoadSolidCondition::CalculateDeltaPositionMatrix(
        const GeometryType& rGeometry,
        Matrix& rDeltaPosition) const
    {
        const SizeType number_of_points = rGeometry.PointsNumber();
        if (rDeltaPosition.size1() != number_of_points || rDeltaPosition.size2() != 3) {
            rDeltaPosition.resize(number_of_points, 3, false);
        }

        for (IndexType i = 0; i < number_of_points; ++i) {
            const auto& r_node = rGeometry[i];
            rDeltaPosition(i, 0) = r_node.X() - r_node.X0();
            rDeltaPosition(i, 1) = r_node.Y() - r_node.Y0();
            rDeltaPosition(i, 2) = r_node.Z() - r_node.Z0();
        }
    }


    void LoadSolidCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        //TODO: build a CalculateOnIntegrationPoints method
        //--------------------------------------------------------------------------------------------
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 2;

        // Initialize Jacobian
        GeometryType::JacobiansType J0;
        // Initialize DN_DX
        Matrix DN_DX(number_of_control_points,2);
        Matrix InvJ0(mDim,mDim);

        const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        Matrix delta_position;
        CalculateDeltaPositionMatrix(r_geometry, delta_position);
        r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);
        double DetJ0;
        Vector old_displacement(mat_size);
        GetSolutionCoefficientVector(old_displacement);
        
        Matrix Jacobian = ZeroMatrix(2,2);
        Jacobian(0,0) = J0[0](0,0);
        Jacobian(0,1) = J0[0](0,1);
        Jacobian(1,0) = J0[0](1,0);
        Jacobian(1,1) = J0[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian,InvJ0,DetJ0);

        // // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[0],InvJ0);

        Matrix B = ZeroMatrix(3,mat_size);

        CalculateB(B, DN_DX);

        array_1d<double, 3> normal_physical_space = GetValue(NORMAL);

        // GET STRESS VECTOR
        ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);
        Vector old_strain = prod(B,old_displacement);

        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        KRATOS_ERROR_IF(mCharacteristicGeometryLength <= 0.0)
            << "Characteristic geometry length must be initialized before finalizing the solution step." << std::endl;
        Values.SetCharacteristicGeometryLength(mCharacteristicGeometryLength);

        ConstitutiveVariables this_constitutive_variables(strain_size);
        Values.SetStrainVector(old_strain);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);
        // update internal variables before requesting stresses
        mpConstitutiveLaw->FinalizeMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);
        mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);

        const Vector sigma = Values.GetStressVector();
        Vector sigma_n(2);

        sigma_n[0] = sigma[0]*normal_physical_space[0] + sigma[2]*normal_physical_space[1];
        sigma_n[1] = sigma[2]*normal_physical_space[0] + sigma[1]*normal_physical_space[1];

        SetValue(NORMAL_STRESS, sigma_n);

        SetValue(CAUCHY_STRESS_XX, sigma[0]);
        SetValue(CAUCHY_STRESS_YY, sigma[1]);
        SetValue(CAUCHY_STRESS_XY, sigma[2]);
        // //---------------------
    }

    void LoadSolidCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
        //--------------------------------------------------------------------------------------------
        // calculate the constitutive law response
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();
        const SizeType mat_size = number_of_control_points * 2;

        // Initialize Jacobian
        GeometryType::JacobiansType J0;
        Matrix delta_position;
        CalculateDeltaPositionMatrix(r_geometry, delta_position);
        r_geometry.Jacobian(J0, this->GetIntegrationMethod(), delta_position);

        // Initialize DN_DX
        const unsigned int dimension = 2;
        Matrix DN_DX(number_of_control_points, 2);
        Matrix InvJ0(dimension, dimension);
        const GeometryType::ShapeFunctionsGradientsType& DN_De =
            r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        double detJ0;

        Matrix Jacobian = ZeroMatrix(2, 2);
        Jacobian(0,0) = J0[0](0,0);
        Jacobian(0,1) = J0[0](0,1);
        Jacobian(1,0) = J0[0](1,0);
        Jacobian(1,1) = J0[0](1,1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian, InvJ0, detJ0);

        // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[0], InvJ0);

        Matrix B = ZeroMatrix(3, mat_size);
        CalculateB(B, DN_DX);

        // Prepare constitutive law parameters
        ConstitutiveLaw::Parameters Values(r_geometry, GetProperties(), rCurrentProcessInfo);

        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        Flags& r_constitutive_law_options = Values.GetOptions();
        r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        ConstitutiveVariables this_constitutive_variables(strain_size);

        Vector displacement_coefficient_vector(mat_size);
        GetSolutionCoefficientVector(displacement_coefficient_vector);
        Vector old_strain = prod(B, displacement_coefficient_vector);

        Values.SetStrainVector(old_strain);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);

        KRATOS_ERROR_IF(mCharacteristicGeometryLength <= 0.0)
            << "Characteristic geometry length must be initialized before initializing the solution step." << std::endl;
        Values.SetCharacteristicGeometryLength(mCharacteristicGeometryLength);

        // update internal variables before requesting stresses
        mpConstitutiveLaw->InitializeMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);
    }

} // Namespace Kratos
