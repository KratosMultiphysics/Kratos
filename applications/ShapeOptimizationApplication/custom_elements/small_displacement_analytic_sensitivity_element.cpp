// ==============================================================================
/*
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//==============================================================================
//
//   Project Name:        KratosShape                            $
//   Created by:          $Author:    daniel.baumgaertner@tum.de $
//                        $Author:           armin.geiser@tum.de $
//   Date:                $Date:                   December 2016 $
//   Revision:            $Revision:                         0.0 $
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
