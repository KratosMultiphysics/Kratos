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
#include "custom_elements/total_lagrangian.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"

namespace Kratos
{
TotalLagrangian::TotalLagrangian( IndexType NewId, GeometryType::Pointer pGeometry )
        : BaseSolidElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

TotalLagrangian::TotalLagrangian( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : BaseSolidElement( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer TotalLagrangian::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_shared<TotalLagrangian>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

TotalLagrangian::~TotalLagrangian()
{
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::CalculateAll( 
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag 
    )
{
    KRATOS_TRY;

    const unsigned int number_of_nodes = this->GetGeometry().size();
    const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
    const unsigned int strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

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

    // Reading integration points
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

    for ( unsigned int point_number = 0; point_number < integration_points.size(); ++point_number ){
        // Contribution to external forces
        const Vector body_force = this->GetBodyForce(integration_points, point_number);
        
        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematicVariables(this_kinematic_variables, point_number, integration_method);
        
        // Compute material reponse
        this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, this->GetStressMeasure());

        // Calculating weights for integration on the reference configuration
        double int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0); 

        if ( dimension == 2 && this->GetProperties().Has( THICKNESS )) 
            int_to_reference_weight *= this->GetProperties()[THICKNESS];

        if ( CalculateStiffnessMatrixFlag == true ) { // Calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config
            /* Material stiffness matrix */
            this->CalculateAndAddKm( rLeftHandSideMatrix, this_kinematic_variables.B, this_constitutive_variables.D, int_to_reference_weight );

            /* Geometric stiffness matrix */
            this->CalculateAndAddKg( rLeftHandSideMatrix, this_kinematic_variables.DN_DX, this_constitutive_variables.StressVector, int_to_reference_weight );
        }

        if ( CalculateResidualVectorFlag == true ) { // Calculation of the matrix is required
            this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, int_to_reference_weight);
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const unsigned int PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{
    // Shape functions
    rThisKinematicVariables.N = row(GetGeometry().ShapeFunctionsValues(rIntegrationMethod), PointNumber);

    // Calculating jacobian
    Matrix J;
    J = this->GetGeometry().Jacobian( J, PointNumber, rIntegrationMethod );
    
    rThisKinematicVariables.detJ0 = this->CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod);
    
    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;
    
    // Deformation gradient
    GeometryUtils::DeformationGradient(J, rThisKinematicVariables.InvJ0,
                                       rThisKinematicVariables.F);

    // Axisymmetric case
    const unsigned int strain_size = (rThisKinematicVariables.B).size1();
    if (strain_size == 4)
    {
        BoundedMatrix<double, 2, 2> F2x2 = rThisKinematicVariables.F;
        rThisKinematicVariables.F.resize(3, 3, false);
        for (unsigned i = 0; i < 2; ++i)
        {
            for (unsigned j = 0; j < 2; ++j)
                rThisKinematicVariables.F(i, j) = F2x2(i, j);
            rThisKinematicVariables.F(i, 2) = rThisKinematicVariables.F(2, i) = 0.0;
        }
        const double current_radius = StructuralMechanicsMathUtilities::CalculateRadius(
            rThisKinematicVariables.N, this->GetGeometry(), Current);
        const double initial_radius = StructuralMechanicsMathUtilities::CalculateRadius(
            rThisKinematicVariables.N, this->GetGeometry(), Initial);
        rThisKinematicVariables.F(2, 2) = current_radius / initial_radius;

        CalculateAxisymmetricB(
            rThisKinematicVariables.B, rThisKinematicVariables.F,
            rThisKinematicVariables.DN_DX, rThisKinematicVariables.N);
    }
    else
    {
        CalculateB(rThisKinematicVariables.B, rThisKinematicVariables.F,
                   rThisKinematicVariables.DN_DX);
    }

    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    KRATOS_TRY
    
    const unsigned int number_of_nodes = this->GetGeometry().PointsNumber();
    const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
    const unsigned strain_size = GetProperties()[CONSTITUTIVE_LAW]->GetStrainSize();
    
    for ( unsigned int i = 0; i < number_of_nodes; ++i ){
        const unsigned int index = dimension * i;

        if ( strain_size == 3 ) {
            rB( 0, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 0 );
            rB( 0, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 0 );
            rB( 1, index + 0 ) = rF( 0, 1 ) * rDN_DX( i, 1 );
            rB( 1, index + 1 ) = rF( 1, 1 ) * rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 1 ) + rF( 0, 1 ) * rDN_DX( i, 0 );
            rB( 2, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 1 ) + rF( 1, 1 ) * rDN_DX( i, 0 );
        } else {
            rB( 0, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 0 );
            rB( 0, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 0 );
            rB( 0, index + 2 ) = rF( 2, 0 ) * rDN_DX( i, 0 );
            rB( 1, index + 0 ) = rF( 0, 1 ) * rDN_DX( i, 1 );
            rB( 1, index + 1 ) = rF( 1, 1 ) * rDN_DX( i, 1 );
            rB( 1, index + 2 ) = rF( 2, 1 ) * rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rF( 0, 2 ) * rDN_DX( i, 2 );
            rB( 2, index + 1 ) = rF( 1, 2 ) * rDN_DX( i, 2 );
            rB( 2, index + 2 ) = rF( 2, 2 ) * rDN_DX( i, 2 );
            rB( 3, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 1 ) + rF( 0, 1 ) * rDN_DX( i, 0 );
            rB( 3, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 1 ) + rF( 1, 1 ) * rDN_DX( i, 0 );
            rB( 3, index + 2 ) = rF( 2, 0 ) * rDN_DX( i, 1 ) + rF( 2, 1 ) * rDN_DX( i, 0 );
            rB( 4, index + 0 ) = rF( 0, 1 ) * rDN_DX( i, 2 ) + rF( 0, 2 ) * rDN_DX( i, 1 );
            rB( 4, index + 1 ) = rF( 1, 1 ) * rDN_DX( i, 2 ) + rF( 1, 2 ) * rDN_DX( i, 1 );
            rB( 4, index + 2 ) = rF( 2, 1 ) * rDN_DX( i, 2 ) + rF( 2, 2 ) * rDN_DX( i, 1 );
            rB( 5, index + 0 ) = rF( 0, 2 ) * rDN_DX( i, 0 ) + rF( 0, 0 ) * rDN_DX( i, 2 );
            rB( 5, index + 1 ) = rF( 1, 2 ) * rDN_DX( i, 0 ) + rF( 1, 0 ) * rDN_DX( i, 2 );
            rB( 5, index + 2 ) = rF( 2, 2 ) * rDN_DX( i, 0 ) + rF( 2, 0 ) * rDN_DX( i, 2 );
        }
    }

    KRATOS_CATCH( "" )
}

void TotalLagrangian::CalculateAxisymmetricB(Matrix& rB,
                                             const Matrix& rF,
                                             const Matrix& rDN_DX,
                                             const Vector& rN)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = this->GetGeometry().PointsNumber();
    const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
    double radius = 0.0;
    radius = StructuralMechanicsMathUtilities::CalculateRadius(rN, this->GetGeometry());
    for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
        const unsigned int index = dimension * i;

            rB(0, index + 0) = rF(0, 0) * rDN_DX(i, 0);
            rB(0, index + 1) = rF(1, 0) * rDN_DX(i, 0);
            rB(1, index + 1) = rF(0, 1) * rDN_DX(i, 1);
            rB(1, index + 1) = rF(1, 1) * rDN_DX(i, 1);
            rB(2, index + 0) = rN[i] / radius;
            rB(3, index + 0) = rF(0, 0) * rDN_DX(i, 1) + rF(0, 1) * rDN_DX(i, 0);
            rB(3, index + 1) = rF(1, 0) * rDN_DX(i, 1) + rF(1, 1) * rDN_DX(i, 0);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

int  TotalLagrangian::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    int ier = BaseSolidElement::Check(rCurrentProcessInfo);

    return ier;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseSolidElement );
}

/***********************************************************************************/
/***********************************************************************************/

void TotalLagrangian::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseSolidElement );
}

} // Namespace Kratos


