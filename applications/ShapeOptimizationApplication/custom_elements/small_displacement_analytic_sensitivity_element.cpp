// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:		 BSD License
//					 license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//					 Geiser Armin, https://github.com/armingeiser
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
: SmallDisplacementElement( NewId, pGeometry )
{
	//DO NOT ADD DOFS HERE!!!
}


// =============================================================================================================================================
// CONSTRUCTOR
// =============================================================================================================================================

SmallDisplacementAnalyticSensitivityElement::SmallDisplacementAnalyticSensitivityElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
: SmallDisplacementElement( NewId, pGeometry, pProperties )
{
	mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}


// =============================================================================================================================================
// COPY CONSTRUCTOR
// =============================================================================================================================================

SmallDisplacementAnalyticSensitivityElement::SmallDisplacementAnalyticSensitivityElement( SmallDisplacementAnalyticSensitivityElement const& rOther)
:SmallDisplacementElement(rOther)
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

	NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

	if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
	{
		NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

		if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
			KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() );
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
		Vector result_x = ZeroVector(number_of_dofs);
		Vector result_y = ZeroVector(number_of_dofs);
		Vector result_z = ZeroVector(number_of_dofs);

		//create and initialize element variables:
		GeneralVariables Variables;
		this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

		//create constitutive law parameters:
		ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

		//set constitutive law flags:
		Flags &ConstitutiveLawOptions=Values.GetOptions();

		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

		//reading integration points
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

		//auxiliary terms
		Vector VolumeForce = ZeroVector(dimension);

		//ask for node for which DKDXU shall be computed
		int active_node_index = this->GetValue(ACTIVE_NODE_INDEX);

		for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
		{
			//compute element kinematics B, F, DN_DX ...
			this->CalculateKinematics(Variables,PointNumber);

			//set general variables to constitutive law parameters
			this->SetGeneralVariables(Variables,Values,PointNumber);

			//compute stresses and constitutive parameters
			mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

			//calculating weights for integration on the "reference configuration"
			double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
			IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );

			// Computation of X-component of DKDXU
			Matrix DB_DX = ZeroMatrix(Variables.B.size1(),Variables.B.size2());
			this->CalculateDerivedDeformationMatrix(DB_DX,Variables.DN_DX,active_node_index,0);
			result_x += IntegrationWeight*prod( trans( DB_DX ), Vector(prod(Variables.ConstitutiveMatrix ,Vector(prod(Variables.B, rU)))) );
			result_x += IntegrationWeight*prod( trans( Variables.B ), Vector(prod(Variables.ConstitutiveMatrix ,Vector(prod(DB_DX, rU)))) );
			result_x += IntegrationWeight*prod( trans( Variables.B ), Vector(prod(Variables.ConstitutiveMatrix ,Vector(prod(Variables.B, rU)))) )*Variables.DN_DX(active_node_index,0);
			this->SetValue(DKDXU_X,result_x);

			// Computation of Y-component of DKDXU
			DB_DX.clear();
			this->CalculateDerivedDeformationMatrix(DB_DX,Variables.DN_DX,active_node_index,1);
			result_y += IntegrationWeight*prod( trans( DB_DX ), Vector(prod(Variables.ConstitutiveMatrix ,Vector(prod(Variables.B, rU)))) );
			result_y += IntegrationWeight*prod( trans( Variables.B ), Vector(prod(Variables.ConstitutiveMatrix ,Vector(prod(DB_DX, rU)))) );
			result_y += IntegrationWeight*prod( trans( Variables.B ), Vector(prod(Variables.ConstitutiveMatrix ,Vector(prod(Variables.B, rU)))) )*Variables.DN_DX(active_node_index,1);
			this->SetValue(DKDXU_Y,result_y);

			// Computation of Z-component of DKDXU
			DB_DX.clear();
			this->CalculateDerivedDeformationMatrix(DB_DX,Variables.DN_DX,active_node_index,2);
			result_z += IntegrationWeight*prod( trans( DB_DX ), Vector(prod(Variables.ConstitutiveMatrix ,Vector(prod(Variables.B, rU)))) );
			result_z += IntegrationWeight*prod( trans( Variables.B ), Vector(prod(Variables.ConstitutiveMatrix ,Vector(prod(DB_DX, rU)))) );
			result_z += IntegrationWeight*prod( trans( Variables.B ), Vector(prod(Variables.ConstitutiveMatrix ,Vector(prod(Variables.B, rU)))) )*Variables.DN_DX(active_node_index,2);
			this->SetValue(DKDXU_Z,result_z);
		}
	}
	
	KRATOS_CATCH( "" )
}

// =============================================================================================================================================
// =============================================================================================================================================


void SmallDisplacementAnalyticSensitivityElement::save( Serializer& rSerializer ) const
{
	KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallDisplacementElement )
}

void SmallDisplacementAnalyticSensitivityElement::load( Serializer& rSerializer )
{
	KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallDisplacementElement )
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
        KRATOS_THROW_ERROR( std::invalid_argument, "Wrong dimension specified.", "" )
		
    KRATOS_CATCH( "" )
}

} // Namespace Kratos
