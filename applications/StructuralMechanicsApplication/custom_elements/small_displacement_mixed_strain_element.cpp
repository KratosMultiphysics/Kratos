// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/small_displacement_mixed_strain_element.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

Element::Pointer SmallDisplacementMixedStrainElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedStrainElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedStrainElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedStrainElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedStrainElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    KRATOS_TRY

    SmallDisplacementMixedStrainElement::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementMixedStrainElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto &r_geometry = GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    const unsigned int dim = r_geometry.WorkingSpaceDimension();
    const unsigned int dof_size  = n_nodes*(dim+1);

    if (rResult.size() != dof_size){
        rResult.resize(dof_size);
    }

    if (dim == 2) {
        for(unsigned int i = 0; i < n_nodes; ++i) {
            rResult[i * (dim + 1)] = this->GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[i * (dim + 1) + 1] = this->GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[i * (dim + 1) + 2] = this->GetGeometry()[i].GetDof(VOLUMETRIC_STRAIN).EquationId();
        }
    } else if (dim == 3) {
        for(unsigned int i = 0; i < n_nodes; ++i){
            rResult[i * (dim + 1)] = this->GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[i * (dim + 1) + 1] = this->GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[i * (dim + 1) + 2] = this->GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[i * (dim + 1) + 3] = this->GetGeometry()[i].GetDof(VOLUMETRIC_STRAIN).EquationId();
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto &r_geometry = GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    const unsigned int dim = r_geometry.WorkingSpaceDimension();
    const unsigned int dof_size  = n_nodes*(dim+1);

    if (rElementalDofList.size() != dof_size){
        rElementalDofList.resize(dof_size);
    }

    if (dim == 2) {
        for(unsigned int i = 0; i < n_nodes; ++i) {
            rElementalDofList[i * (dim + 1)] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * (dim + 1) + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * (dim + 1) + 2] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN);
        }
    } else if (dim == 3) {
        for(unsigned int i = 0; i < n_nodes; ++i){
            rElementalDofList[i * (dim + 1)] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * (dim + 1) + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * (dim + 1) + 2] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
            rElementalDofList[i * (dim + 1) + 3] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{

}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)
{

}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{

}

/***********************************************************************************/
/***********************************************************************************/

bool SmallDisplacementMixedStrainElement::UseElementProvidedStrain() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY;

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false );

        rRightHandSideVector = ZeroVector( mat_size ); //resetting RHS
    }

    // Reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    if ( CalculateStiffnessMatrixFlag ) {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    } else {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    }

    // If strain has to be computed inside of the constitutive law with PK2
    Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

    // Some declarations
    array_1d<double, 3> body_force;
    double int_to_reference_weight;

    // Computing in all integrations points
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        // Contribution to external forces
        noalias(body_force) = this->GetBodyForce(integration_points, point_number);

        // Compute element kinematics B, F, DN_DX ...
        // CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material reponse
        // CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

        // Calculating weights for integration on the reference configuration
        // int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if ( dimension == 2 && GetProperties().Has( THICKNESS ))
            int_to_reference_weight *= GetProperties()[THICKNESS];

        // if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
        //     // Contributions to stiffness matrix calculated on the reference config
        //     this->CalculateAndAddKm( rLeftHandSideMatrix, this_kinematic_variables.B, this_constitutive_variables.D, int_to_reference_weight );
        // }

        // if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
        //     this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, int_to_reference_weight);
        // }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;

    // Displacements vector
    Vector displacements(mat_size);
    GetValuesVector(displacements);

    // Compute strain
    noalias(rThisConstitutiveVariables.StrainVector) = prod(rThisKinematicVariables.B, displacements);

    // Here we essentially set the input parameters
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF); //assuming the determinant is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); //assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure)
{
    // Set the constitutive variables
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> SmallDisplacementMixedStrainElement::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const IndexType PointNumber) const
{

}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::CalculateB(
    Matrix& rB,
    const Matrix& rDN_DX,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const IndexType PointNumber) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    rB.clear();

    if(dimension == 2) {
        for ( SizeType i = 0; i < number_of_nodes; ++i ) {
            rB( 0, i*2     ) = rDN_DX( i, 0 );
            rB( 1, i*2 + 1 ) = rDN_DX( i, 1 );
            rB( 2, i*2     ) = rDN_DX( i, 1 );
            rB( 2, i*2 + 1 ) = rDN_DX( i, 0 );
        }
    } else if(dimension == 3) {
        for ( SizeType i = 0; i < number_of_nodes; ++i ) {
            rB( 0, i*3     ) = rDN_DX( i, 0 );
            rB( 1, i*3 + 1 ) = rDN_DX( i, 1 );
            rB( 2, i*3 + 2 ) = rDN_DX( i, 2 );
            rB( 3, i*3     ) = rDN_DX( i, 1 );
            rB( 3, i*3 + 1 ) = rDN_DX( i, 0 );
            rB( 4, i*3 + 1 ) = rDN_DX( i, 2 );
            rB( 4, i*3 + 2 ) = rDN_DX( i, 1 );
            rB( 5, i*3     ) = rDN_DX( i, 2 );
            rB( 5, i*3 + 2 ) = rDN_DX( i, 0 );
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::CalculateDeviatoricStrainOperator(Matrix& rDevStrainOp) const
{
    KRATOS_TRY;

    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    rDevStrainOp = ZeroMatrix(strain_size, strain_size);

    const double two_thirds = 2.0/3.0;
    const double minus_one_third = -1.0/3.0;

    if (strain_size == 3) {
        rDevStrainOp(0,0) = two_thirds;
        rDevStrainOp(0,1) = minus_one_third;
        rDevStrainOp(1,0) = minus_one_third;
        rDevStrainOp(1,1) = two_thirds;
        rDevStrainOp(2,2) = 1.0;
    } else if (strain_size == 6) {
        rDevStrainOp(0,0) = two_thirds;
        rDevStrainOp(0,1) = minus_one_third;
        rDevStrainOp(0,2) = minus_one_third;
        rDevStrainOp(1,0) = minus_one_third;
        rDevStrainOp(1,1) = two_thirds;
        rDevStrainOp(1,2) = minus_one_third;
        rDevStrainOp(2,0) = minus_one_third;
        rDevStrainOp(2,1) = minus_one_third;
        rDevStrainOp(2,2) = two_thirds;
        rDevStrainOp(3,3) = 1.0;
        rDevStrainOp(4,4) = 1.0;
        rDevStrainOp(5,5) = 1.0;
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

int  SmallDisplacementMixedStrainElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int base_element_check = SmallDisplacementMixedStrainElement::BaseType::Check(rCurrentProcessInfo);

    return base_element_check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacementMixedStrainElement::BaseType);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacementMixedStrainElement::BaseType);
}

} // Namespace Kratos
