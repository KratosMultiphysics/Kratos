
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
#include "custom_conditions/sbm_load_solid_condition.h"
#include "includes/global_pointer_variables.h"

namespace Kratos
{

void SbmLoadSolidCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    InitializeMaterial();
    InitializeMemberVariables();
    InitializeSbmMemberVariables();
}


void SbmLoadSolidCondition::InitializeMaterial()
{
    KRATOS_TRY
    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLaw->InitializeMaterial( r_properties, r_geometry, row(N_values , 0 ));

        mpConstitutiveLawOnTrue = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLawOnTrue->InitializeMaterial(r_properties, r_geometry, row(N_values , 0 ));

    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH( "" );

}

void SbmLoadSolidCondition::InitializeMemberVariables()
{
    // Compute class memeber variables
    const auto& r_geometry = this->GetGeometry();
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    
    // Initialize DN_DX
    mDim = r_DN_De[0].size2();
    
    Vector mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);

    // Compute the local tangents
    array_1d<double, 3> tangent_parameter_space;
    r_geometry.Calculate(LOCAL_TANGENT, tangent_parameter_space); // Gives the result in the parameter space

    const double abs_tangent_u = std::abs(tangent_parameter_space[0]);
    const double abs_tangent_v = std::abs(tangent_parameter_space[1]);
    mCharacteristicGeometryLength = mesh_size_uv[0];
    if (abs_tangent_v > abs_tangent_u) {
        mCharacteristicGeometryLength = mesh_size_uv[1];
    }
    mCharacteristicGeometryLength *= 0.5;
    
    // Compute basis function order (Note: it is not allow to use different orders in different directions)
    if (mDim == 3) {
        mBasisFunctionsOrder = std::cbrt(r_DN_De[0].size1()) - 1;
    } else {
        mBasisFunctionsOrder = std::sqrt(r_DN_De[0].size1()) - 1;
    }

    // mBasisFunctionsOrder *= 2;

    // Compute the normals
    mNormalParameterSpace = - r_geometry.Normal(0, GetIntegrationMethod());
    mNormalParameterSpace = mNormalParameterSpace / MathUtils<double>::Norm(mNormalParameterSpace);
    mNormalPhysicalSpace = mNormalParameterSpace;

    SetValue(NORMAL, mNormalPhysicalSpace);

    // calculate the integration weight
    // reading integration point
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    // Initialize Jacobian
    Matrix InvJ0(3,3);
    GeometryType::JacobiansType J0;
    Matrix delta_position;
    CalculateDeltaPositionMatrix(r_geometry, delta_position);
    r_geometry.Jacobian(J0,this->GetIntegrationMethod(), delta_position);

    // compute complete jacobian transformation including parameter->physical space transformation
    double detJ0;
    Matrix Jacobian = ZeroMatrix(3,3);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);
    Jacobian(2,2) = 1.0;

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,detJ0);

    Vector add_factor = prod(Jacobian, tangent_parameter_space); //additional factor to the determinant of the jacobian for the parameter->physical space transformation
    add_factor[2] = 0.0; 
    detJ0 = norm_2(add_factor);

    const double thickness = GetProperties().Has(THICKNESS) ? GetProperties()[THICKNESS] : 1.0;

    const double int_to_reference_weight = r_integration_points[0].Weight() * std::abs(detJ0) * thickness;

    SetValue(INTEGRATION_WEIGHT, int_to_reference_weight);
}

void SbmLoadSolidCondition::InitializeSbmMemberVariables()
{
    auto& r_geometry = this->GetGeometry();
    std::string loopIdentifier = this->GetValue(IDENTIFIER);

    // NURBS case
    if (this->GetValue(NEIGHBOUR_NODES).size() != 0) 
    {
        mpProjectionNode = &r_geometry.GetValue(NEIGHBOUR_NODES)[0];

        mTrueNormal = mpProjectionNode->GetValue(NORMAL);

        if (loopIdentifier == "inner")
            mTrueNormal = -mTrueNormal;
            
        mDistanceVector.resize(3);
        noalias(mDistanceVector) = mpProjectionNode->Coordinates() - r_geometry.Center().Coordinates();
        this->SetValue(PROJECTION_NODE_COORDINATES, mpProjectionNode->Coordinates());
        // dot product n dot n_tilde
        mTrueDotSurrogateNormal = inner_prod(mNormalParameterSpace, mTrueNormal);

        return;
    }

    // Retrieve projection
    Condition candidate_closest_skin_segment_1 = this->GetValue(NEIGHBOUR_CONDITIONS)[0] ;
    // Find the closest node in condition
    int closestNodeId;
    if (mDim > 2) {
        double incumbent_dist = 1e16;
        // Loop over the three nodes of the closest skin element
        for (unsigned int i = 0; i < 3; i++) {
            if (norm_2(candidate_closest_skin_segment_1.GetGeometry()[i]-r_geometry.Center()) < incumbent_dist) {
                incumbent_dist = norm_2(candidate_closest_skin_segment_1.GetGeometry()[i]-r_geometry.Center());
                closestNodeId = i;
            }
        }
    } else {
        closestNodeId = 0;
    }
    mpProjectionNode = &candidate_closest_skin_segment_1.GetGeometry()[closestNodeId] ;
    mDistanceVector.resize(3);
    noalias(mDistanceVector) = mpProjectionNode->Coordinates() - r_geometry.Center().Coordinates();

    // calculate the integration weight
    // loopIdentifier is inner or outer
    if (mDim == 2) {
        // Need also the second closest condition in 2D
        Condition candidate_closest_skin_segment_2 = this->GetValue(NEIGHBOUR_CONDITIONS)[1] ;
        array_1d<double,3> vector_skin_segment_1 = candidate_closest_skin_segment_1.GetGeometry()[1] - candidate_closest_skin_segment_1.GetGeometry()[0];
        array_1d<double,3> vector_skin_segment_2 = candidate_closest_skin_segment_2.GetGeometry()[1] - candidate_closest_skin_segment_2.GetGeometry()[0];
        array_1d<double,3> vector_out_of_plane = ZeroVector(3);
        vector_out_of_plane[2] = 1.0;
        
        array_1d<double,3> crossProductSkinSegment1;
        array_1d<double,3> crossProductSkinSegment2; 
        MathUtils<double>::CrossProduct(crossProductSkinSegment1, vector_out_of_plane, vector_skin_segment_1);
        MathUtils<double>::CrossProduct(crossProductSkinSegment2, vector_out_of_plane, vector_skin_segment_2);
        
        mTrueNormal = crossProductSkinSegment1 / MathUtils<double>::Norm(crossProductSkinSegment1) + crossProductSkinSegment2 / MathUtils<double>::Norm(crossProductSkinSegment2);
        if (loopIdentifier == "inner") {
            mTrueNormal = mTrueNormal / MathUtils<double>::Norm(mTrueNormal) ;
        } else { // outer
            mTrueNormal = - mTrueNormal / MathUtils<double>::Norm(mTrueNormal) ;
        }
    } else {
        // 3D CASE
        array_1d<double,3> vector_skin_segment_1 = candidate_closest_skin_segment_1.GetGeometry()[1] - candidate_closest_skin_segment_1.GetGeometry()[0];
        array_1d<double,3> vector_skin_segment_2 = candidate_closest_skin_segment_1.GetGeometry()[2] - candidate_closest_skin_segment_1.GetGeometry()[1];
        MathUtils<double>::CrossProduct(mTrueNormal, vector_skin_segment_1, vector_skin_segment_2);

        if (loopIdentifier == "inner") {
            mTrueNormal = mTrueNormal / MathUtils<double>::Norm(mTrueNormal) ;
        } else { // outer
            mTrueNormal = - mTrueNormal / MathUtils<double>::Norm(mTrueNormal) ;
        }
    }

    this->SetValue(PROJECTION_NODE_COORDINATES, mpProjectionNode->Coordinates());

    // dot product n dot n_tilde
    mTrueDotSurrogateNormal = inner_prod(mNormalParameterSpace, mTrueNormal);

}

void SbmLoadSolidCondition::CalculateLocalSystem(
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

    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_control_points = r_geometry.size();

     // reading integration points and local gradients
    const Matrix& r_N = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    const double int_to_reference_weight = GetValue(INTEGRATION_WEIGHT);

    KRATOS_ERROR_IF(mDim != 2) << "SbmLoadSolidCondition momentarily only supports 2D conditions, but the current dimension is" << mDim << std::endl;

    //-------------------------------------------------------------------------
    // Initialize DN_DX
    Matrix DN_DX(number_of_control_points,2);
    Matrix InvJ0(3,3);

    // Initialize Jacobian
    GeometryType::JacobiansType J0;
    Matrix delta_position;
    CalculateDeltaPositionMatrix(r_geometry, delta_position);
    r_geometry.Jacobian(J0,this->GetIntegrationMethod(), delta_position);

    // compute complete jacobian transformation including parameter->physical space transformation
    double detJ0;
    Matrix Jacobian = ZeroMatrix(3,3);
    Jacobian(0,0) = J0[0](0,0);
    Jacobian(0,1) = J0[0](0,1);
    Jacobian(1,0) = J0[0](1,0);
    Jacobian(1,1) = J0[0](1,1);
    Jacobian(2,2) = 1.0;

    // Calculating inverse jacobian and jacobian determinant
    MathUtils<double>::InvertMatrix(Jacobian,InvJ0,detJ0);
    
    // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
    Matrix sub_inv_jacobian = ZeroMatrix(2,2);
    sub_inv_jacobian(0,0) = InvJ0(0,0);
    sub_inv_jacobian(1,0) = InvJ0(1,0);
    sub_inv_jacobian(0,1) = InvJ0(0,1);
    sub_inv_jacobian(1,1) = InvJ0(1,1);
    noalias(DN_DX) = prod(r_DN_De[0],sub_inv_jacobian);

    Matrix B = ZeroMatrix(mDim,mat_size);
    CalculateB(B, DN_DX);

    // Obtain the tangent costitutive law matrix
    ConstitutiveLaw::Parameters values_surrogate(r_geometry, GetProperties(), rCurrentProcessInfo);

    Vector old_displacement_coefficient_vector(mat_size);
    GetSolutionCoefficientVector(old_displacement_coefficient_vector);
    Vector old_strain = prod(B,old_displacement_coefficient_vector);

    const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables(strain_size);
    ApplyConstitutiveLaw(mat_size, old_strain, values_surrogate, this_constitutive_variables);

    // const Matrix& r_D = values_surrogate.GetConstitutiveMatrix();

    Matrix r_D = ZeroMatrix(3,3);
    AnalyticalConstitutiveMatrix(r_D, old_strain);

    const Matrix DB = prod(r_D,B);

    // compute Taylor expansion contribution: H_sum_vec
    Matrix grad_H_sum_transposed = ZeroMatrix(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(mDistanceVector, grad_H_sum_transposed);

    Matrix grad_H_sum = trans(grad_H_sum_transposed);

    Matrix B_sum = ZeroMatrix(mDim,mat_size);
    CalculateB(B_sum, grad_H_sum);

    // obtain the old stress vector on the true boundary (on the projection node)
    ConstitutiveLaw::Parameters values_true(r_geometry, GetProperties(), rCurrentProcessInfo);

    Vector old_strain_on_true = prod(B_sum,old_displacement_coefficient_vector);

    const SizeType strain_size_true = mpConstitutiveLawOnTrue->GetStrainSize();
    ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
    ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true, mpConstitutiveLawOnTrue);
    // const Matrix& r_D_on_true = values_true.GetConstitutiveMatrix();

    Matrix r_D_on_true = ZeroMatrix(3,3);
    AnalyticalConstitutiveMatrix(r_D_on_true, old_strain_on_true);
    
    Matrix DB_sum = prod(r_D_on_true, B_sum); //

    for (IndexType i = 0; i < number_of_control_points; i++) {
        for (IndexType j = 0; j < number_of_control_points; j++) {
            
            for (IndexType idim = 0; idim < 2; idim++) {
                const int id1 = 2*idim;
                const int iglob = 2*i+idim;

                for (IndexType jdim = 0; jdim < 2; jdim++) {
                    const int id2 = (id1+2)%3;
                    const int jglob = 2*j+jdim;

                    // // FLUX STANDARD TERM
                    rLeftHandSideMatrix(iglob, jglob) -= r_N(0,i)*(DB(id1, jglob)* mNormalPhysicalSpace[0] + DB(id2, jglob)* mNormalPhysicalSpace[1]) * int_to_reference_weight;
                
                    // SBM TERM
                    rLeftHandSideMatrix(iglob, jglob) += r_N(0,i)*(DB_sum(id1, jglob)* mTrueNormal[0] + DB_sum(id2, jglob)* mTrueNormal[1]) 
                                                        * int_to_reference_weight * mTrueDotSurrogateNormal;
                }
            }
        }
    }

    // const Vector& r_stress_vector = values_surrogate.GetStressVector();
    Vector r_stress_vector = ZeroVector(3);
    AnalyticalStress(r_stress_vector, old_strain);

    // const Vector& r_stress_vector_on_true = values_true.GetStressVector();
    Vector r_stress_vector_on_true = ZeroVector(3);
    AnalyticalStress(r_stress_vector_on_true, old_strain_on_true);

    // Assembly
    double nu = this->GetProperties().GetValue(POISSON_RATIO);
    double E = this->GetProperties().GetValue(YOUNG_MODULUS);
    Vector g_N = ZeroVector(3);
    const double x = mpProjectionNode->X();
    const double y = mpProjectionNode->Y();

    // g_N[0] = -E*(1+4*nu)/(1-nu*nu)*(sin(x/80)*sin(y/20))/80 * mTrueNormal[0] + E/32/(1+nu) *cos(x/80)*cos(y/20) * mTrueNormal[1]; 
    // g_N[1] = -E*(4+nu)/(1-nu*nu)*(sin(x/80)*sin(y/20))/80 * mTrueNormal[1]   + E/32/(1+nu) *cos(x/80)*cos(y/20) * mTrueNormal[0];


    double t = rCurrentProcessInfo[TIME];

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

    // rStressVector[0] = ((std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))*std::fabs(t*x*y*(x - 4)*(y - 2)) <= 3.0/2000.0) ? (
    //     268000*t*x*std::pow(y, 2)*(2 - x)*(x - 4)*std::pow(y - 2, 2)
    //  )
    //  : (
    //     (268.0/3.0)*t*x*std::pow(y, 2)*(2 - x)*(x - 4)*std::pow(y - 2, 2)*(1000*std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))*std::fabs(t*x*y*(x - 4)*(y - 2)) + 3)*std::fabs(1/(t*x*y*(x - 4)*(y - 2)))/std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))
    //  ));
    //  rStressVector[1] = ((std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))*std::fabs(t*x*y*(x - 4)*(y - 2)) <= 3.0/2000.0) ? (
    //     0
    //  )
    //  : (
    //     (134.0/3.0)*t*x*std::pow(y, 2)*(x - 4)*std::pow(y - 2, 2)*(3*x + 2000*(2 - x)*std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))*std::fabs(t*x*y*(x - 4)*(y - 2)) - 6)*std::fabs(1/(t*x*y*(x - 4)*(y - 2)))/std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))
    //  ));
    //  rStressVector[2] = ((std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))*std::fabs(t*x*y*(x - 4)*(y - 2)) <= 3.0/2000.0) ? (
    //     134000*t*std::pow(x, 2)*y*(1 - y)*std::pow(x - 4, 2)*(y - 2)
    //  )
    //  : (
    //     201*t*std::pow(x, 2)*y*(1 - y)*std::pow(x - 4, 2)*(y - 2)*std::fabs(1/(t*x*y*(x - 4)*(y - 2)))/std::sqrt(3*std::pow(x, 2)*std::pow(x - 4, 2)*std::pow(y - 1, 2) + 4*std::pow(y, 2)*std::pow(x - 2, 2)*std::pow(y - 2, 2))
    //  ));

    // // MIDDLE PLASTICITY    
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
    //  rStressVector[0] = ((std::sqrt(21)*M_PI*std::fabs(t*std::sin((1.0/4.0)*M_PI*x)*std::sin(M_PI*((1.0/2.0)*y + 0.75))) <= 3.0/500.0) ? (
    //     33500*M_PI*t*std::sin((1.0/4.0)*M_PI*x)*std::sin(M_PI*((1.0/2.0)*y + 0.75))
    //  )
    //  : (
    //     (134.0/7.0)*t*(-1750*M_PI + std::sqrt(21)/std::fabs(t*std::sin((1.0/4.0)*M_PI*x)*std::sin(M_PI*((1.0/2.0)*y + 0.75))))*std::sin((1.0/4.0)*M_PI*x)*std::sin(M_PI*((1.0/2.0)*y + 0.75))
    //  ));
    //  rStressVector[1] = ((std::sqrt(21)*M_PI*std::fabs(t*std::sin((1.0/4.0)*M_PI*x)*std::sin(M_PI*((1.0/2.0)*y + 0.75))) <= 3.0/500.0) ? (
    //     -134000*M_PI*t*std::sin((1.0/4.0)*M_PI*x)*std::sin(M_PI*((1.0/2.0)*y + 0.75))
    //  )
    //  : (
    //     -67.0/7.0*t*(3500*M_PI + 3*std::sqrt(21)/std::fabs(t*std::sin((1.0/4.0)*M_PI*x)*std::sin(M_PI*((1.0/2.0)*y + 0.75))))*std::sin((1.0/4.0)*M_PI*x)*std::sin(M_PI*((1.0/2.0)*y + 0.75))
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
     
    g_N[0] = rStressVector[0]*mTrueNormal[0] + rStressVector[2]*mTrueNormal[1];
    g_N[1] = rStressVector[2]*mTrueNormal[0] + rStressVector[1]*mTrueNormal[1];

    // // // cosinusoidal
    // g_N[0] = E/(1-nu)*(sin(x)*sinh(y)) * mTrueNormal[0]; 
    // g_N[1] = E/(1-nu)*(sin(x)*sinh(y)) * mTrueNormal[1]; 


    // g_N = mpProjectionNode->GetValue(FORCE);
    Vector normal_stress_old = ZeroVector(3);
    normal_stress_old[0] = (r_stress_vector[0] * mNormalPhysicalSpace[0] + r_stress_vector[2] * mNormalPhysicalSpace[1]);
    normal_stress_old[1] = (r_stress_vector[2] * mNormalPhysicalSpace[0] + r_stress_vector[1] * mNormalPhysicalSpace[1]);
    Vector normal_stress_true_old = ZeroVector(3);
    normal_stress_true_old[0] = (r_stress_vector_on_true[0] * mTrueNormal[0] + r_stress_vector_on_true[2] * mTrueNormal[1]);
    normal_stress_true_old[1] = (r_stress_vector_on_true[2] * mTrueNormal[0] + r_stress_vector_on_true[1] * mTrueNormal[1]);

    for (IndexType i = 0; i < number_of_control_points; i++) {
        
        for (IndexType idim = 0; idim < 2; idim++) {
            const int iglob = 2*i+idim;

            // // External load term
            rRightHandSideVector[2*i+idim] += r_N(0,i) * g_N[idim] * mTrueDotSurrogateNormal * int_to_reference_weight;

            // // // Residual terms
            // // // FLUX STANDARD TERM
            rRightHandSideVector(iglob) += r_N(0,i) * normal_stress_old[idim] * int_to_reference_weight;
        
            // // SBM TERM
            rRightHandSideVector(iglob) -= r_N(0,i)*normal_stress_true_old[idim] * int_to_reference_weight * mTrueDotSurrogateNormal;
        }
    }
    KRATOS_CATCH("")
}

void SbmLoadSolidCondition::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{   
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

void SbmLoadSolidCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

    int SbmLoadSolidCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(PENALTY_FACTOR))
            << "No penalty factor (PENALTY_FACTOR) defined in property of SupportPenaltyLaplacianCondition" << std::endl;
        return 0;
    }

    void SbmLoadSolidCondition::EquationIdVector(
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

    void SbmLoadSolidCondition::GetDofList(
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


    void SbmLoadSolidCondition::GetSolutionCoefficientVector(
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

    void SbmLoadSolidCondition::CalculateB(
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

    void SbmLoadSolidCondition::ApplyConstitutiveLaw(
        SizeType matSize,
        Vector& rStrain,
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rConstitutiVariables,
        ConstitutiveLaw::Pointer pConstitutiveLaw)
    {
        ConstitutiveLaw::Pointer p_constitutive_law = pConstitutiveLaw ? pConstitutiveLaw : mpConstitutiveLaw;
        KRATOS_ERROR_IF(p_constitutive_law == nullptr)
            << "Constitutive law pointer must be initialized before applying the constitutive law." << std::endl;

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=rValues.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        
        rValues.SetStrainVector(rStrain);
        rValues.SetStressVector(rConstitutiVariables.StressVector);
        rValues.SetConstitutiveMatrix(rConstitutiVariables.D);

        KRATOS_ERROR_IF(mCharacteristicGeometryLength <= 0.0)
            << "Characteristic geometry length must be initialized before applying the constitutive law." << std::endl;
        rValues.SetCharacteristicGeometryLength(mCharacteristicGeometryLength);

        p_constitutive_law->CalculateMaterialResponse(rValues, ConstitutiveLaw::StressMeasure_Cauchy); 
    }

void SbmLoadSolidCondition::CalculateDeltaPositionMatrix(
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


    void SbmLoadSolidCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
        //---------- SET STRESS VECTOR VALUE ----------------------------------------------------------------
        //TODO: build a CalculateOnIntegrationPoints method
        //--------------------------------------------------------------------------------------------
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

        Vector displacement_coefficient_vector(mat_size);
        GetSolutionCoefficientVector(displacement_coefficient_vector);

        Matrix Jacobian = ZeroMatrix(2, 2);
        Jacobian(0, 0) = J0[0](0, 0);
        Jacobian(0, 1) = J0[0](0, 1);
        Jacobian(1, 0) = J0[0](1, 0);
        Jacobian(1, 1) = J0[0](1, 1);

        // Calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(Jacobian, InvJ0, detJ0);

        // Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[0], InvJ0);

        Matrix B = ZeroMatrix(3, mat_size);
        CalculateB(B, DN_DX);

        // Prepare constitutive law parameters
        ConstitutiveLaw::Parameters values_surrogate(r_geometry, GetProperties(), rCurrentProcessInfo);

        const SizeType strain_size = mpConstitutiveLaw->GetStrainSize();
        Flags& r_constitutive_law_options = values_surrogate.GetOptions();
        r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        ConstitutiveVariables this_constitutive_variables(strain_size);

        Vector old_strain = prod(B, displacement_coefficient_vector);

        values_surrogate.SetStrainVector(old_strain);
        values_surrogate.SetStressVector(this_constitutive_variables.StressVector);
        values_surrogate.SetConstitutiveMatrix(this_constitutive_variables.D);

        KRATOS_ERROR_IF(mCharacteristicGeometryLength <= 0.0)
            << "Characteristic geometry length must be initialized before finalizing the solution step." << std::endl;
        values_surrogate.SetCharacteristicGeometryLength(mCharacteristicGeometryLength);

        // update internal variables before requesting stresses
        mpConstitutiveLaw->FinalizeMaterialResponse(values_surrogate, ConstitutiveLaw::StressMeasure_Cauchy);
        mpConstitutiveLaw->CalculateMaterialResponse(values_surrogate, ConstitutiveLaw::StressMeasure_Cauchy);

        {
            Matrix grad_H_sum_transposed = ZeroMatrix(3, number_of_control_points);
            ComputeGradientTaylorExpansionContribution(mDistanceVector, grad_H_sum_transposed);
            Matrix grad_H_sum = trans(grad_H_sum_transposed);

            Matrix B_sum_true = ZeroMatrix(mDim, mat_size);
            CalculateB(B_sum_true, grad_H_sum);

            ConstitutiveLaw::Parameters values_true(r_geometry, GetProperties(), rCurrentProcessInfo);

            Flags& r_constitutive_law_options_true = values_true.GetOptions();
            r_constitutive_law_options_true.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
            r_constitutive_law_options_true.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
            r_constitutive_law_options_true.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

            const SizeType strain_size_true = mpConstitutiveLawOnTrue->GetStrainSize();
            ConstitutiveVariables constitutive_variables_true(strain_size_true);

            Vector old_strain_on_true = prod(B_sum_true, displacement_coefficient_vector);
            values_true.SetStrainVector(old_strain_on_true);
            values_true.SetStressVector(constitutive_variables_true.StressVector);
            values_true.SetConstitutiveMatrix(constitutive_variables_true.D);

            KRATOS_ERROR_IF(mCharacteristicGeometryLength <= 0.0)
                << "Characteristic geometry length must be initialized before finalizing the solution step on the true boundary." << std::endl;
            values_true.SetCharacteristicGeometryLength(mCharacteristicGeometryLength);

            mpConstitutiveLawOnTrue->FinalizeMaterialResponse(values_true, ConstitutiveLaw::StressMeasure_Cauchy);
            mpConstitutiveLawOnTrue->CalculateMaterialResponse(values_true, ConstitutiveLaw::StressMeasure_Cauchy);
        }

        // const Vector sigma = values_surrogate.GetStressVector();

        Vector sigma = ZeroVector(3);
        AnalyticalStress(sigma, old_strain);
        // const Matrix& r_D = values_surrogate.GetConstitutiveMatrix();

        Matrix r_D = ZeroMatrix(3,3);
        AnalyticalConstitutiveMatrix(r_D, old_strain);
        Vector sigma_n(2);

        sigma_n[0] = sigma[0]*mNormalPhysicalSpace[0] + sigma[2]*mNormalPhysicalSpace[1];
        sigma_n[1] = sigma[2]*mNormalPhysicalSpace[0] + sigma[1]*mNormalPhysicalSpace[1];

        SetValue(NORMAL_STRESS, sigma_n);

        SetValue(CAUCHY_STRESS_XX, sigma[0]);
        SetValue(CAUCHY_STRESS_YY, sigma[1]);
        SetValue(CAUCHY_STRESS_XY, sigma[2]);
        // //---------------------
        // // Set the stress vector on the true boundary
        // //---------------------
        std::vector<double> integration_weight_list = GetValue(INTEGRATION_WEIGHTS);
        std::vector<Vector> integration_point_list = GetValue(INTEGRATION_POINTS);
        std::vector<Vector> integration_normal_list = GetValue(INTEGRATION_POINTS_NORMAL);
        const std::string loopIdentifier = this->GetValue(IDENTIFIER);
        const SizeType number_of_integration_points_on_true = integration_weight_list.size();
        Matrix values_on_true_boundary(number_of_integration_points_on_true, 7);
        for (IndexType i = 0; i < number_of_integration_points_on_true; ++i) {
            // Get the integration point
            const Vector& r_integration_point = integration_point_list[i];
            const double weight = integration_weight_list[i];

            Vector normal_on_true = mTrueNormal;
            if (integration_normal_list.size() > i) {
                Vector candidate_normal = integration_normal_list[i];
                if (loopIdentifier == "inner") {
                    candidate_normal *= -1.0;
                }
                if (norm_2(candidate_normal) > 1e-12) {
                    normal_on_true = candidate_normal;
                }
            }
            const double normal_on_true_norm = norm_2(normal_on_true);
            if (normal_on_true_norm > 1e-12) {
                normal_on_true /= normal_on_true_norm;
            }

            Vector distance_vector = ZeroVector(3);
            distance_vector = r_integration_point - r_geometry.Center().Coordinates();
            // compute Taylor expansion contribution: H_sum_vec
            Matrix grad_H_sum_transposed = ZeroMatrix(3, number_of_control_points);
            ComputeGradientTaylorExpansionContribution(distance_vector, grad_H_sum_transposed);

            Matrix grad_H_sum = trans(grad_H_sum_transposed);

            Matrix B_sum = ZeroMatrix(mDim,mat_size);
            CalculateB(B_sum, grad_H_sum);

            // obtain the stress vector on the true boundary 
            ConstitutiveLaw::Parameters values_true(r_geometry, GetProperties(), rCurrentProcessInfo);

            Vector old_strain_on_true = prod(B_sum,displacement_coefficient_vector);

            const SizeType strain_size_true = mpConstitutiveLawOnTrue->GetStrainSize();
            ConstitutiveVariables this_constitutive_variables_true(strain_size_true);
            ApplyConstitutiveLaw(mat_size, old_strain_on_true, values_true, this_constitutive_variables_true, mpConstitutiveLawOnTrue);

            // const Vector& r_stress_vector_on_true = values_true.GetStressVector();

            Vector r_stress_vector_on_true = ZeroVector(3);
            AnalyticalStress(r_stress_vector_on_true, old_strain_on_true);

            Vector normal_stress_true = ZeroVector(3); //FIXME:  correct the normal and use the ones at the new projection points
            normal_stress_true[0] = (r_stress_vector_on_true[0] * normal_on_true[0] + r_stress_vector_on_true[2] * normal_on_true[1]);
            normal_stress_true[1] = (r_stress_vector_on_true[2] * normal_on_true[0] + r_stress_vector_on_true[1] * normal_on_true[1]);

            const double normal_stress = (normal_stress_true[0] * normal_on_true[0] + normal_stress_true[1] * normal_on_true[1]);
            const double shear_stress = (-normal_stress_true[0] * normal_on_true[1] + normal_stress_true[1] * normal_on_true[0]);
            
            values_on_true_boundary(i, 0) = weight;
            values_on_true_boundary(i, 1) = r_integration_point[0];
            values_on_true_boundary(i, 2) = r_integration_point[1];
            values_on_true_boundary(i, 3) = r_integration_point[2];
            values_on_true_boundary(i, 4) = normal_stress;    
            values_on_true_boundary(i, 5) = shear_stress;
        }
        this->SetValue(RESULTS_ON_TRUE_BOUNDARY, values_on_true_boundary);

    }

void SbmLoadSolidCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo){
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
    Jacobian(0, 0) = J0[0](0, 0);
    Jacobian(0, 1) = J0[0](0, 1);
    Jacobian(1, 0) = J0[0](1, 0);
    Jacobian(1, 1) = J0[0](1, 1);

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

    mpConstitutiveLaw->InitializeMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);
    mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_Cauchy);

    Matrix grad_H_sum_transposed = ZeroMatrix(3, number_of_control_points);
    ComputeGradientTaylorExpansionContribution(mDistanceVector, grad_H_sum_transposed);
    Matrix grad_H_sum = trans(grad_H_sum_transposed);

    Matrix B_sum = ZeroMatrix(mDim,mat_size);
    CalculateB(B_sum, grad_H_sum);

    ConstitutiveLaw::Parameters values_true(r_geometry, GetProperties(), rCurrentProcessInfo);

    Flags& r_constitutive_law_options_true = values_true.GetOptions();
    r_constitutive_law_options_true.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_constitutive_law_options_true.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_constitutive_law_options_true.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    const SizeType strain_size_true = mpConstitutiveLawOnTrue->GetStrainSize();
    ConstitutiveVariables constitutive_variables_true(strain_size_true);

    Vector old_strain_on_true = prod(B_sum, displacement_coefficient_vector);
    values_true.SetStrainVector(old_strain_on_true);
    values_true.SetStressVector(constitutive_variables_true.StressVector);
    values_true.SetConstitutiveMatrix(constitutive_variables_true.D);

    KRATOS_ERROR_IF(mCharacteristicGeometryLength <= 0.0)
        << "Characteristic geometry length must be initialized before initializing the solution step on the true boundary." << std::endl;
    values_true.SetCharacteristicGeometryLength(mCharacteristicGeometryLength);

    mpConstitutiveLawOnTrue->InitializeMaterialResponse(values_true, ConstitutiveLaw::StressMeasure_Cauchy);
    mpConstitutiveLawOnTrue->CalculateMaterialResponse(values_true, ConstitutiveLaw::StressMeasure_Cauchy);
}

void SbmLoadSolidCondition::ComputeGradientTaylorExpansionContribution(const Vector& rDistanceVector, Matrix& grad_H_sum)
{
    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_control_points = r_geometry.PointsNumber();
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
                    H_taylor_term_X += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
                }
                for (IndexType k = 0; k <= n-1; k++) {
                    IndexType n_k = n - 1 - k;
                    double derivative = shapeFunctionDerivatives(i,k+1); 
                    // Compute the Taylor term for this derivative
                    H_taylor_term_Y += ComputeTaylorTerm(derivative, rDistanceVector[0], n_k, rDistanceVector[1], k);
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
                            H_taylor_term_X += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x-1, rDistanceVector[1], k_y, rDistanceVector[2], k_z);
                        }
                        if (k_y >= 1) {
                            H_taylor_term_Y += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y-1, rDistanceVector[2], k_z);
                        }
                        if (k_z >= 1) {
                            H_taylor_term_Z += ComputeTaylorTerm3D(derivative, rDistanceVector[0], k_x, rDistanceVector[1], k_y, rDistanceVector[2], k_z-1);
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
double SbmLoadSolidCondition::ComputeTaylorTerm(
    const double derivative, 
    const double dx, 
    const IndexType n_k, 
    const double dy, 
    const IndexType k)
{
    return derivative * std::pow(dx, n_k) * std::pow(dy, k) / (MathUtils<double>::Factorial(k) * MathUtils<double>::Factorial(n_k));    
}

double SbmLoadSolidCondition::ComputeTaylorTerm3D(
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

} // Namespace Kratos
