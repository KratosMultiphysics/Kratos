// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrándiz
//

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/small_displacement.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
SmallDisplacement::SmallDisplacement( IndexType NewId, GeometryType::Pointer pGeometry )
        : BaseSolidElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacement::SmallDisplacement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : BaseSolidElement( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_shared<SmallDisplacement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacement::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_shared<SmallDisplacement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacement::~SmallDisplacement()
{
}

/***********************************************************************************/
/***********************************************************************************/

bool SmallDisplacement::UseElementProvidedStrain()
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacement::CalculateAll( 
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag 
    )
{
    KRATOS_TRY;
    
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const unsigned int mat_size = number_of_nodes * dimension;

    if ( CalculateStiffnessMatrixFlag == true ) { // Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false );

        rRightHandSideVector = ZeroVector( mat_size ); //resetting RHS
    }

    // Reading integration points and local gradients
    const GeometryType::IntegrationMethod integration_method =
        GetGeometry().GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(integration_method);
    
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    
    // If strain has to be computed inside of the constitutive law with PK2
    Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter
    
    for ( unsigned int point_number = 0; point_number < integration_points.size(); point_number++ ) {
        // Contribution to external forces
        const Vector body_force = this->GetBodyForce(integration_points, point_number);
        
        // Compute element kinematics B, F, DN_DX ...
        CalculateKinematicVariables(this_kinematic_variables, point_number, integration_method);
        
        // Compute material reponse
        CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());
        
        // Calculating weights for integration on the reference configuration
        double int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0); 

        if ( dimension == 2 && GetProperties().Has( THICKNESS )) 
            int_to_reference_weight *= GetProperties()[THICKNESS];

        if ( CalculateStiffnessMatrixFlag == true ) { // Calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config
            this->CalculateAndAddKm( rLeftHandSideMatrix, this_kinematic_variables.B, this_constitutive_variables.D, int_to_reference_weight );
        }

        if ( CalculateResidualVectorFlag == true ) { // Calculation of the matrix is required
            this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, int_to_reference_weight);
        }
    }
    
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacement::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables, 
    const unsigned int PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{        
    const GeometryType::IntegrationPointsArrayType& r_integration_points =
        GetGeometry().IntegrationPoints(rIntegrationMethod);
    // Shape functions
    rThisKinematicVariables.N = GetGeometry().ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());
    
    rThisKinematicVariables.detJ0 = CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod); 
    
    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;
    
    // Compute B
    //CalculateB( rThisKinematicVariables.B,/* rThisKinematicVariables.J0,*/ rThisKinematicVariables.DN_DX, r_integration_points, PointNumber );
    CalculateB_DSG( rThisKinematicVariables.B, rThisKinematicVariables.J0, rThisKinematicVariables.DN_DX, r_integration_points, PointNumber );

    // Compute equivalent F
    Vector displacements;
    GetValuesVector(displacements);
    Vector strain_vector = prod(rThisKinematicVariables.B, displacements);
    rThisKinematicVariables.F = ComputeEquivalentF(strain_vector);
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacement::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables, 
    ConstitutiveVariables& rThisConstitutiveVariables, 
    ConstitutiveLaw::Parameters& rValues,
    const unsigned int PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
{        
    // Displacements vector
    Vector displacements;
    GetValuesVector(displacements);
    
    // Compute strain
    noalias(rThisConstitutiveVariables.StrainVector) = prod(rThisKinematicVariables.B, displacements);

    // Here we essentially set the input parameters
    rValues.SetDeterminantF(rThisKinematicVariables.detF); //assuming the determinant is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else
    
    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); //assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
    
    // Actually do the computations in the ConstitutiveLaw    
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done 
}

/***********************************************************************************/
/***********************************************************************************/
Vector getTLine( const Matrix& J_con, int a, int b )
{
    Vector Line;
    Line = ZeroVector( 6 );
    Line(0) = J_con(0,a)*J_con(0,b);
    Line(1) = J_con(1,a)*J_con(1,b);
    Line(2) = J_con(2,a) * J_con(2,b);
    Line(3) = J_con(1,b) * J_con(2,a) + J_con(1,a) * J_con(2,b);
    Line(4) = J_con(0,b) * J_con(2,a) + J_con(0,a) * J_con(2,b);
    Line(5) = J_con(0,b) * J_con(1,a) + J_con(0,a) * J_con(1,b);

    return Line;

}

void SmallDisplacement::CalculateB_DSG(
    Matrix& rB,
    Matrix J,
    const Matrix& rDN_DX,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const unsigned int PointNumber
    )
{
    KRATOS_TRY;
    
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rB.clear();
    
    //Gcov
    Matrix G_cov;
    Vector a;
    Vector b;
    a = ZeroVector( dimension );
    b = ZeroVector( dimension );
    G_cov = ZeroMatrix( dimension, dimension );

    for( int i = 0; i < 3; i++ )
    {
        for( int j = 0; j < 3; j++ )
        {
            a(0) = J(i,0);
            a(1) = J(i,1);
            a(2) = J(i,2);

            b(0) = J(j,0);
            b(1) = J(j,1);
            b(2) = J(j,2);

            G_cov(i,j) = a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
        }
    }
    //Gcon
    Matrix G_con;
    G_con = ZeroMatrix( dimension, dimension );
    double detG_cov = MathUtils<double>::Det( G_cov );
    MathUtils<double>::InvertMatrix( G_cov, G_con, detG_cov );

    //Jcon
    Matrix J_con;
    J_con = ZeroMatrix( dimension, dimension );
    for( int i = 0; i < 3; i++ )
    {
        for( int j = 0; j < 3; j++ )
        {
            J_con( i, 0 ) += G_con( i, j ) * J( j, 0 );
            J_con( i, 1 ) += G_con( i, j ) * J( j, 1 );
            J_con( i, 2 ) += G_con( i, j ) * J( j, 2 );
        }
    }

    //getT
    Matrix T;
    T = ZeroMatrix( 6 , 6 );
    // 1. Line
    Vector TLine = getTLine( J_con, 0, 0 );
    for( int i = 0; i < 6; i++ )
    {
        T(0,i) = TLine(i);
    }
    // 2.Line
    TLine = getTLine( J_con, 1, 1 );
    for( int i = 0; i < 6; i++ )
    {
        T(1,i) = TLine(i);
    }
    //3. Line
    TLine = getTLine( J_con, 2, 2 );
    for( int i = 0; i < 6; i++ )
    {
        T(2,i) = TLine(i);
    }
    // 4. Line
    TLine = getTLine( J_con, 1, 2 );
    for( int i = 0; i < 6; i++ )
    {
        T(3,i) = TLine(i);
    }
    // 5. Line
    TLine = getTLine( J_con, 2, 0 );
    for( int i = 0; i < 6; i++ )
    {
        T(4,i) = TLine(i);
    }
    // 6. Line
    TLine = getTLine( J_con, 0, 1 );
    for( int i = 0; i < 6; i++ )
    {
        T(5,i) = TLine(i);
    }

    if(dimension == 2) {
        for ( unsigned int i = 0; i < number_of_nodes; ++i ) {
            rB( 0, i*2     ) = rDN_DX( i, 0 );
            rB( 1, i*2 + 1 ) = rDN_DX( i, 1 );
            rB( 2, i*2     ) = rDN_DX( i, 1 );
            rB( 2, i*2 + 1 ) = rDN_DX( i, 0 );
        }
    } else if(dimension == 3) {

        rB( 0, 0 ) = -1.0;
        rB( 0, 3 ) =  1.0;

        rB( 1, 1 ) = -1.0;
        rB( 1, 7 ) =  1.0;

        rB( 2, 2 ) = -1.0;
        rB( 2, 11 ) = 1.0;

        rB( 3, 1 ) = -0.5;
        rB( 3, 2 ) = -0.5;
        rB( 3, 8 ) = 0.5;
        rB( 3, 10 ) = 0.5;

        rB( 4, 0 ) = -0.5;
        rB( 4, 2 ) = -0.5;
        rB( 4, 5 ) = 0.5;
        rB( 4, 9 ) = 0.5;

        rB( 5, 0 ) = -0.5;
        rB( 5, 1 ) = -0.5;
        rB( 5, 4 ) = 0.5;
        rB( 5, 6 ) = 0.5;
            
        rB = prod( T, rB );

        /* All others are zero
        for ( unsigned int i = 0; i < number_of_nodes; ++i ) {
            rB( 0, i*3     ) = rDN_DX( i, 0 );
            rB( 1, i*3 + 1 ) = rDN_DX( i, 1 );
            rB( 2, i*3 + 2 ) = rDN_DX( i, 2 );
            rB( 3, i*3     ) = rDN_DX( i, 1 );
            rB( 3, i*3 + 1 ) = rDN_DX( i, 0 );
            rB( 4, i*3 + 1 ) = rDN_DX( i, 2 );
            rB( 4, i*3 + 2 ) = rDN_DX( i, 1 );
            rB( 5, i*3     ) = rDN_DX( i, 2 );
            rB( 5, i*3 + 2 ) = rDN_DX( i, 0 );
        }*/
    }

    KRATOS_CATCH( "" )
}

void SmallDisplacement::CalculateB(
    Matrix& rB,
    const Matrix& rDN_DX,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const unsigned int PointNumber
    )
{
    KRATOS_TRY;
    
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rB.clear();
    
    if(dimension == 2) {
        for ( unsigned int i = 0; i < number_of_nodes; ++i ) {
            rB( 0, i*2     ) = rDN_DX( i, 0 );
            rB( 1, i*2 + 1 ) = rDN_DX( i, 1 );
            rB( 2, i*2     ) = rDN_DX( i, 1 );
            rB( 2, i*2 + 1 ) = rDN_DX( i, 0 );
        }
    } else if(dimension == 3) {

         for ( unsigned int i = 0; i < number_of_nodes; ++i ) {
            rB( 0, i*3     ) = rDN_DX( i, 0 )*10.0;
            rB( 1, i*3 + 1 ) = rDN_DX( i, 1 )*10.0;
            rB( 2, i*3 + 2 ) = rDN_DX( i, 2 )*10.0;
            rB( 3, i*3     ) = rDN_DX( i, 1 )*10.0;
            rB( 3, i*3 + 1 ) = rDN_DX( i, 0 )*10.0;
            rB( 4, i*3 + 1 ) = rDN_DX( i, 2 )*10.0;
            rB( 4, i*3 + 2 ) = rDN_DX( i, 1 )*10.0;
            rB( 5, i*3     ) = rDN_DX( i, 2 )*10.0;
            rB( 5, i*3 + 2 ) = rDN_DX( i, 0 )*10.0;
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

Matrix SmallDisplacement::ComputeEquivalentF(const Vector& rStrainTensor)
{
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    Matrix F(dim,dim);
    
    if(dim == 2) {
        F(0,0) = 1.0+rStrainTensor(0);  
        F(0,1) = 0.5*rStrainTensor(2); 
        F(1,0) = 0.5*rStrainTensor(2);   
        F(1,1) = 1.0+rStrainTensor(1);
    } else {
        F(0,0) = 1.0+rStrainTensor(0);     
        F(0,1) = 0.5*rStrainTensor(3); 
        F(0,2) = 0.5*rStrainTensor(5);
        F(1,0) = 0.5*rStrainTensor(3);   
        F(1,1) = 1.0+rStrainTensor(1);   
        F(1,2) = 0.5*rStrainTensor(4);
        F(2,0) = 0.5*rStrainTensor(5);   
        F(2,1) = 0.5*rStrainTensor(4); 
        F(2,2) = 1.0+rStrainTensor(2);
    }
    
    return F;
}

/***********************************************************************************/
/***********************************************************************************/

int  SmallDisplacement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    int ier = BaseSolidElement::Check(rCurrentProcessInfo);

    return ier;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseSolidElement );
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseSolidElement );
}

} // Namespace Kratos


