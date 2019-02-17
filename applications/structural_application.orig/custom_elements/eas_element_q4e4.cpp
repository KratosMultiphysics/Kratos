/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
giang.bui@rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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

==============================================================================
 */
/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 2013-04-09 14:47:00 $
 *   Revision:            $Revision: 1.1 $
 *
 * ***********************************************************/

// System includes


// External includes

// Project includes
// #include "includes/define.h"
#include "custom_elements/eas_element_q4e4.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application.h"


//#define ENABLE_Q4E4_DEBUG_LEVEL1


namespace Kratos
{

	EASElementQ4E4::EASElementQ4E4( IndexType NewId, GeometryType::Pointer pGeometry ) : Element( NewId, pGeometry )
	{
		mIsInitialized = false;
		//DO NOT ADD DOFS HERE!!!
		//THIS IS THE DEFAULT CONSTRUCTOR
	}

	/**
	 * A simple kinematic linear 3D element for the solution
	 * of the momentum balance in structural mechanics.
	 * This element is used for students training at the Ruhr University Bochum.
	 * Therefore it may includes comments that are obvious for the
	 * experienced user.
	 * hbui modified to work in 2D
	 */
	EASElementQ4E4::EASElementQ4E4( IndexType NewId,
			GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
		: Element( NewId, pGeometry, pProperties )
	{
		mIsInitialized = false;
		
		mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();//default method

        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        
        if(dim != 2)
            KRATOS_THROW_ERROR(std::logic_error, "This element only works in 2D", "");
        
        if ( GetGeometry().size() != 4 && GetGeometry().size() != 8 && GetGeometry().size() != 9 )
            KRATOS_THROW_ERROR(std::logic_error, "This element only works with quadrilateral geometry", "");
        
	}

	Element::Pointer EASElementQ4E4::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
	{
		return Element::Pointer( new EASElementQ4E4( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
	}

	EASElementQ4E4::~EASElementQ4E4()
	{
	}

	/**
	 * Initialization of the element, called at the begin of each simulation.
	 * Member variables and the Material law are initialized here
	 */
	void EASElementQ4E4::Initialize()
	{
		KRATOS_TRY

		//dimension of the problem
	    unsigned int dim = GetGeometry().WorkingSpaceDimension();
	    
		if ( mIsInitialized )
		{
			//Set Up Initial displacement for StressFreeActivation of Elements
			mInitialDisp.resize( GetGeometry().size(), dim, false );

			for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
				for ( unsigned int i = 0; i < dim; i++ )
					mInitialDisp( node, i ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT )[i];

			return;
		}

        //hbui remark: the code below will be called if mIsInitialized == false
        
		//number of integration points used, mThisIntegrationMethod refers to the
		//integration method defined in the constructor
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

		//initializing the Jacobian, the inverse Jacobian and Jacobians determinant in the reference
		// configuration
		GeometryType::JacobiansType J0( integration_points.size() );

		mInvJ0.resize( integration_points.size() );

		mTotalDomainInitialSize = 0.00;

		for ( unsigned int i = 0; i < integration_points.size(); i++ )
		{
			mInvJ0[i].resize( dim, dim, false );
			noalias( mInvJ0[i] ) = ZeroMatrix( dim, dim );
		}

		mDetJ0.resize( integration_points.size(), false );

		noalias( mDetJ0 ) = ZeroVector( integration_points.size() );

		//calculating the Jacobian
		J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod );

		//calculating the inverse Jacobian
		for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
		{
			//getting informations for integration
			double IntegrationWeight = integration_points[PointNumber].Weight();
			//calculating and storing inverse of the jacobian and the parameters needed
			MathUtils<double>::InvertMatrix( J0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber] );
			//calculating the total area
			mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
		}

		//Set Up Initial displacement for StressFreeActivation of Elements
		mInitialDisp.resize( GetGeometry().size(), dim, false );

		for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
			for ( unsigned int i = 0; i < dim; i++ )
				mInitialDisp( node, i ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT )[i];

		//                         mInitialDisp(node, i)= 0.0;
		//Initialization of the constitutive law vector and
		// declaration, definition and initialization of the material
		// lwas at each integration point
		if ( mConstitutiveLawVector.size() != integration_points.size() )
		{
			mConstitutiveLawVector.resize( integration_points.size() );
		}

		InitializeMaterial();

        //initialize incompatible mode
        mIncompatibleMode.resize( 4 );
        noalias( mIncompatibleMode ) = ZeroVector( 4 );
        
        //calculate Jacobian at centre
        Matrix CentreJ( dim, dim );
        noalias( CentreJ ) = ZeroMatrix( dim, dim );
        Matrix F0 (3, 3);
        noalias( F0 ) = ZeroMatrix( 3, 3 );
        double temp;
	    
	    GeometryType::PointType p0(0.00 , 0.00 , 0.00);
        CentreJ = GetGeometry().Jacobian( CentreJ, p0 );
        
        CalculateF0operator(F0 , CentreJ);
        
        mDetJcentre = MathUtils<double>::Det(CentreJ);
        MathUtils<double>::InvertMatrix( F0, mInverseF0operator, temp );
        
		mIsInitialized = true;


//        KRATOS_WATCH(CentreJ);
//        KRATOS_WATCH(F0);
//        KRATOS_WATCH(mDetJcentre);
//        KRATOS_WATCH(mInverseF0operator);


		KRATOS_CATCH( "" )
	}

	/**
	 * Calculate double Variables at each integration point, used for postprocessing etc.
	 * @param rVariable Global name of the variable to be calculated
	 * @param output Vector to store the values on the qudrature points, output of the method
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& Output, const ProcessInfo& rCurrentProcessInfo )
	{
		if ( Output.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
			Output.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(), false );

		for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
			Output[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, Output[ii] );
	}

	/**
	 * Calculate Vector Variables at each integration point, used for postprocessing etc.
	 * @param rVariable Global name of the variable to be calculated
	 * @param output Vector to store the values on the qudrature points, output of the method
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable,
			std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo )
	{
		GetValueOnIntegrationPoints( rVariable, Output, rCurrentProcessInfo );
	}

	/**
	 * Calculate Matrix Variables at each integration point, used for postprocessing etc.
	 * @param rVariable Global name of the variable to be calculated
	 * @param output Vector to store the values on the qudrature points, output of the method
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable,
			std::vector<Matrix>& Output, const ProcessInfo& rCurrentProcessInfo )
	{
		KRATOS_TRY

        //TODO add incompatible mode
//        
//		unsigned int number_of_nodes = GetGeometry().size();
//		unsigned int dim = GetGeometry().WorkingSpaceDimension();;
//		unsigned int StrainSize = dim * ( dim + 1 ) / 2;

//		//Initialize local variables
//		Matrix B( StrainSize, number_of_nodes*dim );
//		Matrix TanC( StrainSize, StrainSize );
//		Vector StrainVector( StrainSize );
//		Vector StressVector( StrainSize );
//		Matrix DN_DX( number_of_nodes, dim );
//		Matrix CurrentDisp( number_of_nodes, dim );

//		//reading integration points and local gradients
//		const GeometryType::IntegrationPointsArrayType& integration_points =
//			GetGeometry().IntegrationPoints( mThisIntegrationMethod );
//		const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

//		if ( Output.size() != integration_points.size() )
//			Output.resize( integration_points.size() );

//		const GeometryType::ShapeFunctionsGradientsType& DN_De =
//			GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

//		//Current displacements
//		for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
//			noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

//		//Declaration of the integration weight
//		//    double Weight;

//		//loop over all integration points
//		for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
//		{
//			noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );
//			//Initializing B_Operator at the current integration point
//			CalculateBoperator( B, DN_DX );

//			//calculate strain
//			CalculateStrain( B, CurrentDisp, StrainVector );
//			//assign the integration weight at the current integration point
//			//        Weight = integration_points[PointNumber].Weight();

//			//calculate material response
//			mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
//					StrainVector,
//					ZeroMatrix( 1 ),
//					StressVector,
//					TanC,
//					rCurrentProcessInfo,
//					GetProperties(),
//					GetGeometry(),
//					row( Ncontainer, PointNumber ),
//					true,
//					0,
//					true );

//			if ( Output[PointNumber].size2() != StrainVector.size() )
//				Output[PointNumber].resize( 1, StrainVector.size(), false );

//            // Add in variable if needed here. Refer to KinematicLinear2 element.

//		}

		//         std::cout << std::endl;
		KRATOS_CATCH( "" )
	}

	/**
	 * Initialization of the Material law at each integration point
	 */
	void EASElementQ4E4::InitializeMaterial()
	{
		KRATOS_TRY

			for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
			{
				mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
				mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(),
						row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
			}

		KRATOS_CATCH( "" )
	}

	void EASElementQ4E4::ResetConstitutiveLaw()
	{
		KRATOS_TRY

			for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
			{
				mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
			}

		KRATOS_CATCH( "" )
	}



	/**
	 * THIS is the main method here the integration in space (loop over the integration points) is done,
	 * the algorithmic tangent and the (inner and outer) load vector is computed
	 * @param rLeftHandSideMatrix algorithmic tangent, size (number_of_nodes*dim)*(number_of_nodes*dim)
	 * @param rRightHandSideVector (inner and outer) load vector, size (number_of_nodes*dim)
	 * @param rCurrentProcessInfo
	 * @param CalculateStiffnessMatrixFlag true: algorithmic tangent has to be computed
	 * @param CalculateResidualVectorFlag true: load vector has to be computed
	 */
	void EASElementQ4E4::CalculateAll( MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,
			bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag )
	{
		KRATOS_TRY

		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();;
		unsigned int StrainSize = dim * ( dim + 1 ) / 2;

		//this is the size of the elements stiffness matrix/force vector
		unsigned int MatSize = GetGeometry().size() * dim;

        
//        for(int i = 0; i < number_of_nodes; i++)
//        {
//            KRATOS_WATCH(GetGeometry().GetPoint(i).X());
//            KRATOS_WATCH(GetGeometry().GetPoint(i).Y());
            //Note that X() and Y() are absolute coordinate of a point (which includes displacement)
//        }


		//Initialize local variables
		Matrix B( StrainSize, MatSize );
		Matrix TanC( StrainSize, StrainSize );
		Vector StrainVector( StrainSize );
		Vector StressVector( StrainSize );
		Matrix DN_DX( number_of_nodes, dim );
		Matrix CurrentDisp( number_of_nodes, dim );
		Matrix G( 3, 4 );
		Matrix Kbg( MatSize, 4 );
		Matrix Kgg( 4, 4 );
		Matrix Kgb( 4, MatSize );


		if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
		{
			//resize the LHS=StiffnessMatrix if its size is not correct
			if ( rLeftHandSideMatrix.size1() != MatSize )
				rLeftHandSideMatrix.resize( MatSize, MatSize, false );

			noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
			
			noalias( Kbg ) = ZeroMatrix( MatSize, 4 );
			noalias( Kgg ) = ZeroMatrix( 4, 4 );
			noalias( Kgb ) = ZeroMatrix( 4 , MatSize );
		}

		//resizing as needed the RHS
		if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
		{
			//resize the RHS=force vector if its size is not correct
			if ( rRightHandSideVector.size() != MatSize )
				rRightHandSideVector.resize( MatSize, false );

			noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
		}

		//reading integration points and local gradients
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

		const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

		const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

		//Current displacements
		for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
			noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );
        
		//auxiliary terms
		Vector BodyForce;

        //compute incompatible mode
        CalculateIncompatibleMode( mIncompatibleMode , rCurrentProcessInfo );
        
        #ifdef ENABLE_Q4E4_DEBUG_LEVEL1
        KRATOS_WATCH(mIncompatibleMode);
        #endif

		/////////////////////////////////////////////////////////////////////////
		//// Integration in space over quadrature points
		/////////////////////////////////////////////////////////////////////////
		for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
		{
			noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );
			//Initializing B_Operator at the current integration point
			CalculateBoperator( B, DN_DX );
            
            //Compute enhanced strain operator
			CalculateGoperator( G, PointNumber);
            
			//calculate strain
			CalculateStrain( B, CurrentDisp, G, mIncompatibleMode, StrainVector );

//            KRATOS_WATCH(StrainVector);
//            KRATOS_WATCH(StressVector);

			//calculate stress and tangent operator
			mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
					StrainVector,
					ZeroMatrix( 1 ),
					StressVector,
					TanC,
					rCurrentProcessInfo,
					GetProperties(),
					GetGeometry(),
					row( Ncontainer, PointNumber ),
					true,
//					2, //retains old tangent and update stress incrementally for DruckerPragerExtended
					1,
					true );

//            KRATOS_WATCH(StressVector);
//            KRATOS_WATCH(TanC);

			//calculating weights for integration on the reference configuration
			double IntToReferenceWeight = integration_points[PointNumber].Weight();

            //multiply with thickness in case of 2D
			if ( dim == 2 ) IntToReferenceWeight *= GetProperties()[THICKNESS];

			if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
			{
			    double Temp = ( IntToReferenceWeight * mDetJ0[PointNumber] );
			    
				//calculate stiffness matrix
				noalias( rLeftHandSideMatrix ) += prod( trans( B ), Temp * Matrix( prod( TanC, B ) ) );
					
				noalias( Kbg ) += prod( trans( B ), Temp * Matrix( prod( TanC, G ) ) );
				noalias( Kgg ) += prod( trans( G ), Temp * Matrix( prod( TanC, G ) ) );
				noalias( Kgb ) += prod( trans( G ), Temp * Matrix( prod( TanC, B ) ) );
			}

			if ( CalculateResidualVectorFlag == true )
			{
				//contribution to external forces
				BodyForce = GetProperties()[BODY_FORCE];
				
				CalculateAndAdd_ExtForceContribution( row( Ncontainer, PointNumber ), rCurrentProcessInfo, BodyForce, rRightHandSideVector, IntToReferenceWeight, mDetJ0[PointNumber]);
				
				//contribution of gravity (if there is)
				AddBodyForcesToRHS( rRightHandSideVector, row( Ncontainer, PointNumber ), IntToReferenceWeight, mDetJ0[PointNumber] );

				AddInternalForcesToRHS( rRightHandSideVector, B, StressVector, IntToReferenceWeight, mDetJ0[PointNumber] );
			}
			
//			KRATOS_WATCH(StrainVector);
//			KRATOS_WATCH(mIncompatibleMode);
//			KRATOS_WATCH(CurrentDisp);
			
//			//update new tangent for the material
//			mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
//					StrainVector,
//					ZeroMatrix( 1 ),
//					StressVector,
//					TanC,
//					rCurrentProcessInfo,
//					GetProperties(),
//					GetGeometry(),
//					row( Ncontainer, PointNumber ),
//					true,
//					1,
//					true );
					
			
		}//loop over integration points
	    
	    
	    if ( CalculateStiffnessMatrixFlag == true ) //modify stiffness matrix is required
		{
		    Matrix InverseKgg( 4, 4 );
            SD_MathUtils<double>::InvertMatrix(Kgg, InverseKgg);
            #ifdef ENABLE_Q4E4_DEBUG_LEVEL1
            KRATOS_WATCH(prod( InverseKgg, Kgb ));
            KRATOS_WATCH(prod( Kbg, Matrix( prod( InverseKgg, Kgb ) ) ));
            #endif
		    noalias( rLeftHandSideMatrix ) += ( -1 ) * prod( Kbg, Matrix( prod( InverseKgg, Kgb ) ) );
		}
        
//        KRATOS_WATCH("------");        
        
		KRATOS_CATCH( "" )
	}

	/**
	 * THIS method is called from the scheme during the iteration loop, it calls the CalculateAll()
	 * method with CalculateStiffnessMatrixFlag = false and CalculateResidualVectorFlag = true
	 * @param rRightHandSideVector (inner and outer) load vector, size (number_of_nodes*dim)
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::CalculateRightHandSide( VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo )
	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = false;
		bool CalculateResidualVectorFlag = true;
		MatrixType temp = Matrix();

		CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );
	}

	/**
	 * THIS method is called from the scheme during the iteration loop, it calls the CalculateAll()
	 * method with CalculateStiffnessMatrixFlag = true and CalculateResidualVectorFlag = true
	 * @param rLeftHandSideMatrix algorithmic tangent, size (number_of_nodes*dim)*(number_of_nodes*dim)
	 * @param rRightHandSideVector (inner and outer) load vector, size (number_of_nodes*dim)
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = true;
		bool CalculateResidualVectorFlag = true;
		CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
				CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
	}

	/**
	 * THIS method is called from the scheme after each solution step, here the time step
	 * start and end point variables can be transferred n --> n+1
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
	{
		//reading integration points and local gradients
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

		for ( unsigned int Point = 0; Point < integration_points.size(); Point++ )
		{
			mConstitutiveLawVector[Point]->FinalizeSolutionStep( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), Point ), CurrentProcessInfo );
		}

	}

	/**
	 * THIS method is called from the scheme at the start of each solution step
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
	{
		//reading integration points and local gradients
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

		for ( unsigned int Point = 0; Point < integration_points.size(); Point++ )
		{
			mConstitutiveLawVector[Point]->InitializeSolutionStep( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), Point ), CurrentProcessInfo );
		}
	}
	
	void EASElementQ4E4::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
        //reset all resistant forces at node
        for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
        {
            GetGeometry()[i].GetSolutionStepValue( REACTION_X ) = 0.0;
            GetGeometry()[i].GetSolutionStepValue( REACTION_Y ) = 0.0;
        }
    }

	void EASElementQ4E4::MassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
	{
	    // this function need to be reimplemented or not support
	}

	//************************************************************************************
	//************************************************************************************
	void EASElementQ4E4::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
	{
	    // this function need to be reimplemented or not support
	}

	//************************************************************************************
	//************************************************************************************
	void EASElementQ4E4::GetValuesVector( Vector& values, int Step )
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes * dim;

		if ( values.size() != MatSize )    values.resize( MatSize, false );

		for ( unsigned int i = 0; i < number_of_nodes; i++ )
		{
			unsigned int index = i * dim;
			values[index] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

			if ( dim == 3 )
				values[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
		}
	}


	//************************************************************************************
	//************************************************************************************
	void EASElementQ4E4::GetFirstDerivativesVector( Vector& values, int Step )
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes * dim;

		if ( values.size() != MatSize )   values.resize( MatSize, false );

		for ( unsigned int i = 0; i < number_of_nodes; i++ )
		{
			unsigned int index = i * dim;
			values[index] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

			if ( dim == 3 )
				values[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
		}
	}

	//************************************************************************************
	//************************************************************************************
	void EASElementQ4E4::GetSecondDerivativesVector( Vector& values, int Step )
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes * dim;

		if ( values.size() != MatSize ) values.resize( MatSize, false );

		for ( unsigned int i = 0; i < number_of_nodes; i++ )
		{
			unsigned int index = i * dim;
			values[index] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

			if ( dim == 3 )
				values[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
		}
	}

	/**
	 * returns the used integration method
	 */
	EASElementQ4E4::IntegrationMethod EASElementQ4E4::GetIntegrationMethod() const
	{
		return mThisIntegrationMethod;
	}

	/**
	 * not used
	 */
	void EASElementQ4E4::CalculateAndAddExtForceContribution(
	        const Vector& N,
			const ProcessInfo& CurrentProcessInfo,
			Vector& BodyForce,
			VectorType& rRightHandSideVector,
			double weight )
	{
		KRATOS_TRY
        
        //TODO ?
        
		KRATOS_CATCH( "" )
	}

	/**
	 * Informations for the builder and solver to assemble the global vectors and matrices.
	 * Here a Vector containing the EquationIds of the differnt Dofs is created
	 * @param rResult Vector of the EquationIds
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
	{
		unsigned int dim = ( GetGeometry().WorkingSpaceDimension() );
		unsigned int MatSize = GetGeometry().size() * dim;

		if ( rResult.size() != MatSize )
			rResult.resize( MatSize, false );

		for ( unsigned int i = 0 ; i < GetGeometry().size() ; i++ )
		{
			int index = i * dim;
			rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
			rResult[index+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
			if(dim == 3)
				rResult[index+2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
		}
	}

	/**
	 * Informations for the builder and solver to assemble the global vectors and matrices.
	 * Here a Container containing the pointers rto the DOFs of this element is created
	 * @param ElementalDofList Container with of the DOFs associated with the nodes
	 *                           of this element
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo )
	{
		unsigned int dim = GetGeometry().WorkingSpaceDimension();

		ElementalDofList.resize( 0 );

		for ( unsigned int i = 0 ; i < GetGeometry().size() ; i++ )
		{
			ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
			ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
			if(dim == 3)
				ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
		}
	}

	/**
	 * Adds the Body Forces to the load vector
	 * @param R RHS Vector
	 * @param N_DISP shape function values at the current integration points
	 * @param Weight current integration weight
	 * @param detJ current Determinant of the Jacobian
	 */
	inline void EASElementQ4E4::AddBodyForcesToRHS( Vector& R, const Vector& N_DISP, double Weight, double detJ )
	{
		KRATOS_TRY

		unsigned int dim = GetGeometry().WorkingSpaceDimension();

		Vector gravity( dim );

		double density = 0.0;
		if( GetValue( USE_DISTRIBUTED_PROPERTIES ) )
		{
			noalias( gravity ) = GetValue(GRAVITY);
			density = GetValue(DENSITY);
		}
		else
		{
			noalias( gravity ) = GetProperties()[GRAVITY];
			density = GetProperties()[DENSITY];
		}

		for ( unsigned int prim = 0; prim < GetGeometry().size(); prim++ )
		{
			for ( unsigned int i = 0; i < dim; i++ )
			{
				R( prim*dim + i ) += N_DISP( prim ) * density * gravity( i ) * detJ * Weight;
			}
		}

		KRATOS_CATCH( "" )
	}

	inline void EASElementQ4E4::CalculateAndAdd_ExtForceContribution(
			const Vector& N,
			const ProcessInfo& CurrentProcessInfo,
			Vector& BodyForce,
			VectorType& rRightHandSideVector,
			double weight,
			double detJ
			)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		for ( unsigned int i = 0; i < number_of_nodes; i++ )
		{
			int index = dimension * i;

			for ( unsigned int j = 0; j < dimension; j++ ) rRightHandSideVector[index + j] += weight * detJ * N[i] * BodyForce[j];
		}

		KRATOS_CATCH( "" )
	}

	/**
	 * Adds the Internal Forces to the load vector
	 * @param R RHS Vector
	 * @param B_Operator B-Operator at the current integration point
	 * @param StressVector current stress vector
	 * @param Weight current integration weight
	 * @param detJ current Determinant of the Jacobian
	 */
//	void EASElementQ4E4::AddInternalForcesToRHS( Vector& R, const Matrix& B_Operator, Vector& StressVector, double Weight, double detJ )
//	{
//		KRATOS_TRY

//		unsigned int dim = GetGeometry().WorkingSpaceDimension();
//		unsigned int StrainSize = dim * (dim + 1) / 2;

//		for ( unsigned int prim = 0; prim < GetGeometry().size(); prim++ )
//		{
//			for ( unsigned int i = 0; i < dim; i++ )
//			{
//				for ( unsigned int gamma = 0; gamma < StrainSize; gamma++ )
//				{
//					R( prim * dim + i ) += ( -1 ) * ( B_Operator( gamma, prim * dim + i ) * StressVector( gamma ) * detJ * Weight );
//				}
//			}
//		}
//		//         noalias(R) -= detJ*Weight* prod(trans(B_Operator),StressVector);

//		KRATOS_CATCH( "" )
//	}

    void EASElementQ4E4::AddInternalForcesToRHS( Vector& R, const Matrix& B_Operator, Vector& StressVector, double Weight, double detJ )
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int StrainSize = dim * (dim + 1) / 2;
        double aux;
        
        for ( unsigned int prim = 0; prim < GetGeometry().size(); prim++ )
        {
            for ( unsigned int i = 0; i < dim; i++ )
            {
                aux = 0.0;
                for ( unsigned int gamma = 0; gamma < StrainSize; gamma++ )
                {
                    aux += B_Operator( gamma, prim * dim + i ) * StressVector( gamma ) * detJ * Weight;
                }
                
                R( prim * dim + i ) -= aux;
                
                if( i == 0 )
                    GetGeometry()[prim].GetSolutionStepValue( REACTION_X ) += aux;
                
                if( i == 1 )
                    GetGeometry()[prim].GetSolutionStepValue( REACTION_Y ) += aux;
                    
                if( i == 2 )
                    GetGeometry()[prim].GetSolutionStepValue( REACTION_Z ) += aux;
            }
        }
        
        KRATOS_CATCH( "" )
    }

	/**
	 * Computes the strain vector
	 */
	void EASElementQ4E4::CalculateStrain(
	    const Matrix& B,
	    const Matrix& Displacements,
	    const Matrix& G,
	    const Vector& IncompatibleMode,
	    Vector& StrainVector
	    )
	{
		KRATOS_TRY
		unsigned int Dim = GetGeometry().WorkingSpaceDimension();
		unsigned int StrainSize = Dim * (Dim + 1) / 2;
		noalias( StrainVector ) = ZeroVector( StrainSize );

		for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
		{
			for ( unsigned int item = 0; item < StrainSize; item++ )
				for ( unsigned int dim = 0; dim < Dim; dim++ )
					StrainVector[item] += B( item, Dim * node + dim ) * ( Displacements( node, dim ) - mInitialDisp( node, dim ) );
		}

        noalias( StrainVector ) += prod( G , IncompatibleMode );

		KRATOS_CATCH( "" )
	}

	/**
	 * Computes the standard B-Operator at the current quadrature point
	 * @param B_Operator current B-operator
	 * @param DN_DX shape function values at the current integration point
	 */
	void EASElementQ4E4::CalculateBoperator ( Matrix& B_Operator, const Matrix& DN_DX )
    {
		KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().PointsNumber();

		unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int StrainSize = dim * (dim + 1) / 2;
		
		//         if(B_Operator.size() != number_of_nodes)
		//             B_Operator.resize(number_of_nodes);
		noalias( B_Operator ) = ZeroMatrix( StrainSize, number_of_nodes * dim );

		if(dim == 2) //hbui added
		{
			for ( unsigned int i = 0; i < number_of_nodes; i++ )
			{
				B_Operator( 0, i*2 ) = DN_DX( i, 0 );
				B_Operator( 1, i*2 + 1 ) = DN_DX( i, 1 );
				B_Operator( 2, i*2 ) = DN_DX( i, 1 );
				B_Operator( 2, i*2 + 1 ) = DN_DX( i, 0 );
			}
		}
		else if(dim == 3)
		{
			for ( unsigned int i = 0; i < number_of_nodes; i++ )
			{
				//             if(B_Operator[i].size1() != 6 || B_Operator[i].size2() != 3)
				//                 B_Operator[i].resize(6,3);
				//             noalias(B_Operator[i])= ZeroMatrix(6,3);

				B_Operator( 0, i*3 ) = DN_DX( i, 0 );
				B_Operator( 1, i*3 + 1 ) = DN_DX( i, 1 );
				B_Operator( 2, i*3 + 2 ) = DN_DX( i, 2 );
				B_Operator( 3, i*3 ) = DN_DX( i, 1 );
				B_Operator( 3, i*3 + 1 ) = DN_DX( i, 0 );
				B_Operator( 4, i*3 + 1 ) = DN_DX( i, 2 );
				B_Operator( 4, i*3 + 2 ) = DN_DX( i, 1 );
				B_Operator( 5, i*3 ) = DN_DX( i, 2 );
				B_Operator( 5, i*3 + 2 ) = DN_DX( i, 0 );
			}
		}

		KRATOS_CATCH( "" )
	}

    void EASElementQ4E4::CalculateF0operator ( Matrix& F0_Operator, Matrix& J )
    {
        F0_Operator(0,0) = J(0,0) * J(0,0);
        F0_Operator(0,1) = J(1,0) * J(0,1); /*J(0,1) * J(0,1);*/
        F0_Operator(0,2) = J(0,0) * J(0,1) * 2.00;
        
        F0_Operator(1,0) = J(0,1) * J(1,0); /*J(1,0) * J(1,0);*/
        F0_Operator(1,1) = J(1,1) * J(1,1);
        F0_Operator(1,2) = J(1,0) * J(1,1) * 2.00;
        
        F0_Operator(2,0) = J(0,0) * J(1,0);
        F0_Operator(2,1) = J(0,1) * J(1,1);
        F0_Operator(2,2) = J(0,0) * J(1,1) + J(0,1) * J(1,0);
    }
    
    void EASElementQ4E4::CalculateGoperator ( Matrix& G_Operator, IndexType IntegrationPointIndex)
    {
//         unsigned int dim = GetGeometry().WorkingSpaceDimension();
        noalias( G_Operator ) = ZeroMatrix( 3, 4 );
        
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
        const GeometryType::IntegrationPointType IntegrationPoint = integration_points[IntegrationPointIndex];
        
        Matrix E_Operator( 3, 4 );
        noalias( E_Operator ) = ZeroMatrix( 3, 4 );
        double xi = IntegrationPoint.X();
        double eta = IntegrationPoint.Y();
        
        E_Operator( 0, 0 ) = xi;
        E_Operator( 1, 1 ) = eta;
        E_Operator( 2, 2 ) = xi;
        E_Operator( 2, 3 ) = eta;
        
        noalias( G_Operator ) = (mDetJcentre / mDetJ0[IntegrationPointIndex]) * prod( trans( mInverseF0operator ) , E_Operator );
    }
    
    /**
	 * Computes the incompatible mode
	 */
	void EASElementQ4E4::CalculateIncompatibleMode (
	    Vector& rIncompatibleMode,
	    const ProcessInfo& rCurrentProcessInfo
	    )
	{
		KRATOS_TRY
        unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();;
		unsigned int StrainSize = dim * (dim + 1) / 2;
		
		//this is the size of the elements stiffness matrix/force vector
		unsigned int MatSize = GetGeometry().size() * dim;
        
        //reset incompatible mode
        
		//Initialize local variables
		Matrix B( StrainSize, MatSize );
		Matrix TanC( StrainSize, StrainSize );
		Vector StrainVector( StrainSize );
		Vector StressVector( StrainSize );
		Matrix DN_DX( number_of_nodes, dim );
		Matrix CurrentDisp( number_of_nodes, dim );
		Matrix G( 3, 4 );
		Matrix Kgg( 4, 4 );
		Matrix InverseKgg( 4, 4 );
		Vector h( 4 );
		#ifdef ENABLE_Q4E4_DEBUG_LEVEL1
		Vector RHS_test( 4 ); //a test value for rIncompatibleMode
		#endif
        
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
	    
	    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

		const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

		//Current displacements
		for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
			noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

        double eps = 1.00;
        double tol = 1e-9;
        int max_cnt = 15;
        int cnt = 0;
//         bool CalculateStiffnessMatrixFlag = true;

        #ifdef ENABLE_Q4E4_DEBUG_LEVEL1
        KRATOS_WATCH(CurrentDisp);
        #endif

        while(cnt++ < max_cnt)
        {
            noalias( Kgg ) = ZeroMatrix( 4, 4 );
            noalias( h ) = ZeroVector( 4 );
            #ifdef ENABLE_Q4E4_DEBUG_LEVEL1
            noalias( RHS_test ) = ZeroVector( 4 );
            #endif
        
		    /////////////////////////////////////////////////////////////////////////
		    //// Integration in space over quadrature points
		    /////////////////////////////////////////////////////////////////////////
		    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
		    {
		        noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );
			    //Initializing B_Operator at the current integration point
			    CalculateBoperator( B, DN_DX );
                
                //Compute enhanced strain operator
			    CalculateGoperator( G, PointNumber);
			    
//			    KRATOS_WATCH(integration_points[PointNumber]);
//                KRATOS_WATCH(G);
                
			    //calculate strain
			    CalculateStrain( B, CurrentDisp, G, rIncompatibleMode, StrainVector );

//                KRATOS_WATCH(StrainVector);
//                KRATOS_WATCH(StressVector);

			    //calculate stress and tangent operator
			    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
					    StrainVector,
					    ZeroMatrix( 1 ),
					    StressVector,
					    TanC,
					    rCurrentProcessInfo,
					    GetProperties(),
					    GetGeometry(),
					    row( Ncontainer, PointNumber ),
					    true,
//					    2, //retain old tangent, do stress update incrementally
					    1,
					    false );
                
//                KRATOS_WATCH(StressVector);
//                KRATOS_WATCH(TanC);
                
//                KRATOS_WATCH(StrainVector);
//                KRATOS_WATCH(TanC);
//                KRATOS_WATCH(StressVector);
//                Vector tmpStress( 3 );
//                noalias( tmpStress ) = prod( TanC , StrainVector );
//                KRATOS_WATCH(tmpStress);
                
                
			    //calculating weights for integration on the reference configuration
			    double IntToReferenceWeight = integration_points[PointNumber].Weight();

                //multiply with thickness in case of 2D
			    if ( dim == 2 ) IntToReferenceWeight *= GetProperties()[THICKNESS];

                double Temp = ( IntToReferenceWeight * mDetJ0[PointNumber] );
                
                //articulate to orthogonal condition
                noalias( h ) += Temp * prod( trans( G ), StressVector );
			    noalias( Kgg ) += prod( trans( G ), Temp * Matrix( prod( TanC, G ) ) );
    
//                KRATOS_WATCH(h);
                
		        #ifdef ENABLE_Q4E4_DEBUG_LEVEL1
                //compute test value
		        CalculateStrain( B, CurrentDisp, G, ZeroVector( 4 ), StrainVector ); //StrainVector = B * CurrentDisp
		        
		        KRATOS_WATCH(StrainVector);
		        KRATOS_WATCH(Temp);
		        KRATOS_WATCH(G);
		        KRATOS_WATCH(TanC);
		        Vector product = prod( trans( G ), Temp * Vector( prod( TanC, StrainVector ) ) );
		        KRATOS_WATCH(product);
		        noalias( RHS_test ) += product;
		        #endif
		    }

            #ifdef ENABLE_Q4E4_DEBUG_LEVEL1
//		    KRATOS_WATCH(Kgg);
//		    KRATOS_WATCH(InverseKgg);
            #endif

		    //incrementally update incomaptible mode
		    SD_MathUtils<double>::InvertMatrix(Kgg , InverseKgg);
            noalias( rIncompatibleMode ) += ( -1.00 ) * prod( InverseKgg , h );

		    //compute stopping criteria
		    double NormalizedCoefficient = fabs( (StressVector[0] + StressVector[1]) / 2 );
		    if(NormalizedCoefficient > tol)
		        eps = sqrt(inner_prod(h, h)) / NormalizedCoefficient;
		    else
		        eps = sqrt(inner_prod(h, h));
		    if(eps < tol)
		        break;
		    
            
            #ifdef ENABLE_Q4E4_DEBUG_LEVEL1
		    KRATOS_WATCH(eps);
//		    KRATOS_WATCH(h);
//		    KRATOS_WATCH(InverseKgg);
//		    KRATOS_WATCH(rIncompatibleMode);
            KRATOS_WATCH("--------------");
            #endif
            
		}
		
		if(cnt >= max_cnt)
		{
		    KRATOS_WATCH(eps);
		    KRATOS_THROW_ERROR(std::logic_error, "Iteration to update incompatible mode does not converge in 15 loops", "");
		}
		
		#ifdef ENABLE_Q4E4_DEBUG_LEVEL1
		KRATOS_WATCH(rIncompatibleMode);
		KRATOS_WATCH(RHS_test);
		Vector test = ( -1.00 ) * prod( InverseKgg , RHS_test );
		KRATOS_WATCH(test);
		#endif
		
		
		KRATOS_CATCH( "" )
	}
	
	/**
	 * Calculate Matrix Variables at each integration point, used for postprocessing etc.
	 * @param rVariable Global name of the variable to be calculated
	 * @param rValues Vector to store the values on the qudrature points, output of the method
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
	{
	    //add matrix variables if needed
		return;
	}

	/**
	 * Calculate Vector Variables at each integration point, used for postprocessing etc.
	 * @param rVariable Global name of the variable to be calculated
	 * @param rValues Vector to store the values on the qudrature points, output of the method
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
	{
		// std::cout<<"GetValue On Integration Points"<<std::endl;
		if ( rValues.size() != mConstitutiveLawVector.size() )
			rValues.resize( mConstitutiveLawVector.size() );

		unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int StrainSize = dim * (dim + 1) / 2;

		if ( rVariable == MATERIAL_PARAMETERS )
		{
			if ( rValues.size() !=
					GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
				rValues.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() );

			for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
				rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
		}

		if ( rVariable == PLASTIC_STRAIN_VECTOR )
		{
			for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
			{
				if ( rValues[i].size() != StrainSize )
					rValues[i].resize( StrainSize );

				noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( PLASTIC_STRAIN_VECTOR, rValues[i] );
			}
		}

		if ( rVariable == STRESSES )
		{
			for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
			{
				if ( rValues[i].size() != StrainSize )
					rValues[i].resize( StrainSize );

				noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( STRESSES, rValues[i] );
			}
		}

		if ( rVariable == INTERNAL_VARIABLES )
		{
			//        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
			//        {
			//            if ( rValues[i].size() != 8 ) //hbui: why 8?
			//                rValues[i].resize( 8 );

			//            noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( INTERNAL_VARIABLES, rValues[i] );
			//        }

			KRATOS_WATCH("INTERNAL_VARIABLE are disabled for this element");

		}

		// std::cout<<"END::GetValue On Integration Points"<<std::endl;
	}

	/**
	 * Calculate double Variables at each integration point, used for postprocessing etc.
	 * @param rVariable Global name of the variable to be calculated
	 * @param rValues Vector to store the values on the qudrature points, output of the method
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
	{
		if ( rValues.size() != GetGeometry().IntegrationPoints().size() )
			rValues.resize( GetGeometry().IntegrationPoints().size() );

		//reading integration points and local gradients
		for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
		{
			rValues[Point] = mConstitutiveLawVector[Point]->GetValue( rVariable, rValues[Point] );
		}
	}

	/**
	 * Set a Matrix Variable from outside
	 * @param rVariable Global name of the variable to be calculated
	 * @param rValues Vector of the values on the quadrature points
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
	{
	    // This function is not necessary. Check for KinematicLinear2 for way to add in variables.
	}

	/**
	 * Set a Vector Variable from outside
	 * @param rVariable Global name of the variable to be calculated
	 * @param rValues Vector of the values on the quadrature points
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
	{
	    // This function is not necessary. Check for KinematicLinear2 for way to add in variables.
	}

	/**
	 * Set a Double Variable from outside
	 * @param rVariable Global name of the variable to be calculated
	 * @param rValue value on the quadrature points
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
			std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
	{
	    // This function is not necessary. Check for KinematicLinear2 for way to add in variables.
	}

    /**
	 * Set ConstitutiveLaw value from outside
	 * @param rVariable Global name of the variable to be calculated
	 * @param rValue value on the quadrature points
	 * @param rCurrentProcessInfo
	 */
	void EASElementQ4E4::SetValueOnIntegrationPoints( const Variable< ConstitutiveLaw::Pointer >& rVariable, std::vector< ConstitutiveLaw::Pointer >& rValues, const Kratos::ProcessInfo& rCurrentProcessInfo )
	{
	    // This function is not necessary. Check for KinematicLinear2 for way to add in variables.
	}

	int EASElementQ4E4::Check( const Kratos::ProcessInfo& rCurrentProcessInfo )
	{
		KRATOS_TRY

	    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

		if ( this->Id() < 1 )
		{
			KRATOS_THROW_ERROR( std::logic_error, "Element found with Id 0 or negative", "" );
		}

		if ( this->GetGeometry().Area() < 0 )
		{
			std::cout << "error on element -> " << this->Id() << std::endl;
			KRATOS_THROW_ERROR( std::logic_error, "Area can not be less than 0", "" );
		}

		if ( this->GetGeometry().Volume() < 0 ) //hbui added
		{
			std::cout << "error on element -> " << this->Id() << std::endl;
			KRATOS_THROW_ERROR( std::logic_error, "Volume can not be less than 0", "" );
		}


		//verify that the constitutive law exists
		if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
		{
			KRATOS_THROW_ERROR( std::logic_error, "constitutive law not provided for property ", this->GetProperties().Id() );
		}

		//Verify that the body force is defined
		if ( this->GetProperties().Has( BODY_FORCE ) == false )
		{
			KRATOS_THROW_ERROR( std::logic_error, "BODY_FORCE not provided for property ", this->GetProperties().Id() )
		}

		//verify that the constitutive law has the correct dimension
		if ( dimension == 2 )
		{
			if ( this->GetProperties().Has( THICKNESS ) == false )
				KRATOS_THROW_ERROR( std::logic_error, "THICKNESS not provided for element ", this->Id() );

			if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 3 )
				KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law used. This is a 2D element! expected strain size is 3 (el id = ) ", this->Id() );
		}
		else if(dimension == 3)
		{
			if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 6 )
				KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) ", this->Id() );
		}

		//check constitutive law
		for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
		{
			return mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );
		}



		int ok = 0;
		for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
		{
			ok = mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );
			if( ok != 0 ) break;
//			if( mConstitutiveLawVector[i]->IsIncremental() )
//				KRATOS_THROW_ERROR( std::logic_error, "This element does not provide incremental strains!", "" );
//			if( mConstitutiveLawVector[i]->GetStrainMeasure() != ConstitutiveLaw::StrainMeasure_Linear )
//				KRATOS_THROW_ERROR( std::logic_error, "This element formulated in linear strain measure", "" );
//			if( mConstitutiveLawVector[i]->GetStressMeasure() != ConstitutiveLaw::StressMeasure_PK1 )
//				KRATOS_THROW_ERROR( std::logic_error, "This element is formulated in PK1 stresses", "" );
		}

		return ok;

		KRATOS_CATCH( "" );

	}



} // Namespace Kratos


