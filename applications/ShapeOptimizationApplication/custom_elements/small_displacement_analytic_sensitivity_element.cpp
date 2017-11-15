// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

// Includes
#include "small_displacement_analytic_sensitivity_element.hpp"
#include "shape_optimization_application.h"

namespace Kratos
{

// =============================================================================================================================================
// CONSTRUCTOR
// =============================================================================================================================================

SmallDisplacementAnalyticSensitivityElement::SmallDisplacementAnalyticSensitivityElement( IndexType NewId, GeometryType::Pointer pGeometry )
: SmallDisplacement( NewId, pGeometry )
{
	//DO NOT ADD DOFS HERE!!!
}


// =============================================================================================================================================
// CONSTRUCTOR
// =============================================================================================================================================

SmallDisplacementAnalyticSensitivityElement::SmallDisplacementAnalyticSensitivityElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
: SmallDisplacement( NewId, pGeometry, pProperties )
{
}


// =============================================================================================================================================
// COPY CONSTRUCTOR
// =============================================================================================================================================

SmallDisplacementAnalyticSensitivityElement::SmallDisplacementAnalyticSensitivityElement( SmallDisplacementAnalyticSensitivityElement const& rOther)
:SmallDisplacement(rOther)
{
}


// =============================================================================================================================================
// OPERATIONS
// =============================================================================================================================================

Element::Pointer SmallDisplacementAnalyticSensitivityElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
	return Element::Pointer( new SmallDisplacementAnalyticSensitivityElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


// =============================================================================================================================================
// CLONE
// =============================================================================================================================================

Element::Pointer SmallDisplacementAnalyticSensitivityElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

	SmallDisplacementAnalyticSensitivityElement NewElement ( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

	//-----------//

	if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
	{
		NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

		if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
			KRATOS_ERROR << "Constitutive law not has the correct size. Size is: " << NewElement.mConstitutiveLawVector.size() << std::endl;
	}


	return Element::Pointer( new SmallDisplacementAnalyticSensitivityElement(NewElement) );
}

// =============================================================================================================================================
// DESTRUCTOR
// =============================================================================================================================================

SmallDisplacementAnalyticSensitivityElement::~SmallDisplacementAnalyticSensitivityElement()
{
}


// =============================================================================================================================================
// STARTING / ENDING METHODS
// =============================================================================================================================================

void SmallDisplacementAnalyticSensitivityElement::Calculate(const Variable<Vector> &rVariable, Vector& rU, const ProcessInfo &rCurrentProcessInfo)
{
	KRATOS_TRY

	// Calculation of derivative of stiffness matrix w.r.t. x,y,z multiplied by a given vector U
	if (rVariable == DKDXU)
	{
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
	    const unsigned int number_of_dofs = dimension*GetGeometry().size();
		const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
		Vector result_x = ZeroVector(number_of_dofs);
		Vector result_y = ZeroVector(number_of_dofs);
		Vector result_z = ZeroVector(number_of_dofs);

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);		
        
        // Reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(  );

        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);		
        
        // Displacements vector
        Vector displacements;
        GetValuesVector(displacements);
        
		//ask for node for which DKDXU shall be computed
		int active_node_index = this->GetValue(ACTIVE_NODE_INDEX);

        for ( unsigned int point_number = 0; point_number < integration_points.size(); point_number++ )
        {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, integration_points);

            // Compute material reponse
            CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure(), displacements);

            // Calculating weights for integration on the reference configuration
            double IntegrationWeight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0); 

            if ( dimension == 2 && GetProperties().Has( THICKNESS )) 
            {
                IntegrationWeight *= GetProperties()[THICKNESS];
            }

			// Computation of X-component of DKDXU
			Matrix DB_DX = ZeroMatrix(this_kinematic_variables.B.size1(),this_kinematic_variables.B.size2());
			this->CalculateDerivedDeformationMatrix(DB_DX,this_kinematic_variables.DN_DX,active_node_index,0);
			result_x += IntegrationWeight*prod( trans( DB_DX ), Vector(prod(this_constitutive_variables.D ,Vector(prod(this_kinematic_variables.B, rU)))) );
			result_x += IntegrationWeight*prod( trans( this_kinematic_variables.B ), Vector(prod(this_constitutive_variables.D ,Vector(prod(DB_DX, rU)))) );
			result_x += IntegrationWeight*prod( trans( this_kinematic_variables.B ), Vector(prod(this_constitutive_variables.D ,Vector(prod(this_kinematic_variables.B, rU)))) )*this_kinematic_variables.DN_DX(active_node_index,0);
			this->SetValue(DKDXU_X,result_x);

			// Computation of Y-component of DKDXU
			DB_DX.clear();
			this->CalculateDerivedDeformationMatrix(DB_DX,this_kinematic_variables.DN_DX,active_node_index,1);
			result_y += IntegrationWeight*prod( trans( DB_DX ), Vector(prod(this_constitutive_variables.D ,Vector(prod(this_kinematic_variables.B, rU)))) );
			result_y += IntegrationWeight*prod( trans( this_kinematic_variables.B ), Vector(prod(this_constitutive_variables.D ,Vector(prod(DB_DX, rU)))) );
			result_y += IntegrationWeight*prod( trans( this_kinematic_variables.B ), Vector(prod(this_constitutive_variables.D ,Vector(prod(this_kinematic_variables.B, rU)))) )*this_kinematic_variables.DN_DX(active_node_index,1);
			this->SetValue(DKDXU_Y,result_y);

			// Computation of Z-component of DKDXU
			DB_DX.clear();
			this->CalculateDerivedDeformationMatrix(DB_DX,this_kinematic_variables.DN_DX,active_node_index,2);
			result_z += IntegrationWeight*prod( trans( DB_DX ), Vector(prod(this_constitutive_variables.D ,Vector(prod(this_kinematic_variables.B, rU)))) );
			result_z += IntegrationWeight*prod( trans( this_kinematic_variables.B ), Vector(prod(this_constitutive_variables.D ,Vector(prod(DB_DX, rU)))) );
			result_z += IntegrationWeight*prod( trans( this_kinematic_variables.B ), Vector(prod(this_constitutive_variables.D ,Vector(prod(this_kinematic_variables.B, rU)))) )*this_kinematic_variables.DN_DX(active_node_index,2);
			this->SetValue(DKDXU_Z,result_z);
		}
	}
	
	KRATOS_CATCH( "" )
}

// =============================================================================================================================================
// =============================================================================================================================================


void SmallDisplacementAnalyticSensitivityElement::save( Serializer& rSerializer ) const
{
	KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallDisplacement )
}

void SmallDisplacementAnalyticSensitivityElement::load( Serializer& rSerializer )
{
	KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallDisplacement )
}

// ----------------------------------------------------------------------------------------------------
void SmallDisplacementAnalyticSensitivityElement::CalculateDerivedDeformationMatrix( Matrix& rDB_DX,
        																			 const Matrix& rDN_DX, 
																					 const int node_index, 
																					 const int direction )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int voigt_size = dimension * (dimension +1) * 0.5;

    if ( rDB_DX.size1() != voigt_size || rDB_DX.size2() != dimension*number_of_nodes )
    	rDB_DX.resize(voigt_size, dimension*number_of_nodes, false );

    if( dimension == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 2 * i;

            rDB_DX( 0, index + 0 ) = - rDN_DX( i, direction ) * rDN_DX( node_index, 0 );
            rDB_DX( 1, index + 1 ) = - rDN_DX( i, direction ) * rDN_DX( node_index, 1 );
            rDB_DX( 2, index + 0 ) = - rDN_DX( i, direction ) * rDN_DX( node_index, 1 );
            rDB_DX( 2, index + 1 ) = - rDN_DX( i, direction ) * rDN_DX( node_index, 0 );
        }
    }
    else if( dimension == 3 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 3 * i;

            rDB_DX( 0, index + 0 ) = - rDN_DX( i, direction ) * rDN_DX( node_index, 0 );
            rDB_DX( 1, index + 1 ) = - rDN_DX( i, direction ) * rDN_DX( node_index, 1 );
            rDB_DX( 2, index + 2 ) = - rDN_DX( i, direction ) * rDN_DX( node_index, 2 );
            rDB_DX( 3, index + 0 ) = - rDN_DX( i, direction ) * rDN_DX( node_index, 1 );
            rDB_DX( 3, index + 1 ) = - rDN_DX( i, direction ) * rDN_DX( node_index, 0 );
            rDB_DX( 4, index + 1 ) = - rDN_DX( i, direction ) * rDN_DX( node_index, 2 );
            rDB_DX( 4, index + 2 ) = - rDN_DX( i, direction ) * rDN_DX( node_index, 1 );
            rDB_DX( 5, index + 0 ) = - rDN_DX( i, direction ) * rDN_DX( node_index, 2 );
            rDB_DX( 5, index + 2 ) = - rDN_DX( i, direction ) * rDN_DX( node_index, 0 );
        }
    }
    else
        KRATOS_ERROR << "Wrong dimension specified." << std::endl;
		
    KRATOS_CATCH( "" )
}

} // Namespace Kratos
