//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "custom_elements/linear_solid_element.hpp"

#include "solid_mechanics_application_variables.h"


namespace Kratos
{

//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************

LinearSolidElement::LinearSolidElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearSolidElement::LinearSolidElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{
    //BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LinearSolidElement::LinearSolidElement( LinearSolidElement const& rOther)
    :Element(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    ,mConstitutiveLawVector(rOther.mConstitutiveLawVector)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT AFTER COPYING AN ELEMENT HAVE TO BE DEFINED HERE
    //IF NO ASSIGMENT OPERATOR IS DEFINED THE COPY CONSTRUCTOR WILL DEFINE IT BY DEFFAULT
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

LinearSolidElement&  LinearSolidElement::operator=(LinearSolidElement const& rOther)
{

    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Element::operator=(rOther);

    mThisIntegrationMethod = rOther.mThisIntegrationMethod;

    mConstitutiveLawVector.clear();
    mConstitutiveLawVector.resize( rOther.mConstitutiveLawVector.size());

    for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
    {
        mConstitutiveLawVector[i] = rOther.mConstitutiveLawVector[i];
    }


    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer LinearSolidElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    //NEEDED TO CREATE AN ELEMENT   
    return Element::Pointer( new LinearSolidElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer LinearSolidElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    //YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
    //ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

    LinearSolidElement NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    //-----------//

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());
	
	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() );
      }
    
       
    for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
      {
	NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
      }

    //-----------//

    return Element::Pointer( new LinearSolidElement(NewElement) );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LinearSolidElement::~LinearSolidElement()
{
}


//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

LinearSolidElement::IntegrationMethod LinearSolidElement::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}

//************************************************************************************
//************************************************************************************

void LinearSolidElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    //NEEDED TO DEFINE THE DOFS OF THE ELEMENT 
    rElementalDofList.resize( 0 );
    const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}


//************************************************************************************
//************************************************************************************

void LinearSolidElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    //NEEDED TO DEFINE GLOBAL IDS FOR THE CORRECT ASSEMBLY
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int element_size          = number_of_nodes * dimension;

    if ( rResult.size() != element_size )
        rResult.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension;
        rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        if( dimension == 3)
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void LinearSolidElement::GetValuesVector( Vector& rValues, int Step )
{
    //GIVES THE VECTOR WITH THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT DISPLACEMENTS)
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension;

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );

    }
}


//************************************VELOCITY****************************************
//************************************************************************************

void LinearSolidElement::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    //GIVES THE VECTOR WITH THE TIME DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT VELOCITIES)
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension;

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

        if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
    }
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void LinearSolidElement::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    //GIVES THE VECTOR WITH THE TIME SECOND DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT ACCELERATIONS)
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension;

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
    }

}


//*********************************SET VALUE METHODS**********************************
//************************************************************************************
//(see element.h)

//*********************************GET VALUE METHODS**********************************
//************************************************************************************
//(see element.h)

void LinearSolidElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
    {
      CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
  
  return;

  KRATOS_CATCH( "" )

}


//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************


void LinearSolidElement::Initialize()
{
    KRATOS_TRY

    //INITIALIZE THE CONSTITUTIVE LAW DEFINED IN PROPERTIES
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //Constitutive Law initialisation
    if ( mConstitutiveLawVector.size() != integration_points.size() )
    {
        mConstitutiveLawVector.resize( integration_points.size() );
    }

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(),
                    row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }
    }
    else{
        KRATOS_THROW_ERROR( std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id() )
    }


    KRATOS_CATCH( "" )
}



////************************************************************************************
////************************************************************************************

void LinearSolidElement::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //call the constitutive law to initialize the solution step
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->InitializeSolutionStep( GetProperties(),
                GetGeometry(),
                row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),
                rCurrentProcessInfo );

    KRATOS_CATCH( "" )
}


////************************************************************************************
////************************************************************************************
void LinearSolidElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{

}

////************************************************************************************
////************************************************************************************

void LinearSolidElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{

}

////************************************************************************************
////************************************************************************************

void LinearSolidElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //call the constitutive law to finalize the solution step
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
      mConstitutiveLawVector[i]->FinalizeSolutionStep( GetProperties(),
	      GetGeometry(),
              row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),
	      rCurrentProcessInfo );

    //IF INTERNAL VARIABLES EXIST; THE MATERIAL RESPONSE MUST BE FINALIZED ALSO
    //see constitutive_law.h and other solid elements


    //explicit case:
    this->ClearNodalForces();

    KRATOS_CATCH( "" )
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************



//METHOD TO BUILD ALL TERMA OF THE LOCAL SYSTEM: FOR SINGLE TERMS IT MUST BE SIMPLIFIED 
//IS NEEDED IN IMPLICIT AND IN EXPLICIT CALCULATIONS

//************************************************************************************
//************************************************************************************

void LinearSolidElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY

    //1.-Initialize sizes for the system components:

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int system_size = number_of_nodes * dimension;

    if ( rLeftHandSideMatrix.size1() != system_size )
      rLeftHandSideMatrix.resize( system_size, system_size, false );

    noalias( rLeftHandSideMatrix ) = ZeroMatrix( system_size, system_size ); //resetting LHS

    //resizing as needed the RHS
    if ( rRightHandSideVector.size() != system_size )
      rRightHandSideVector.resize( system_size, false );

    rRightHandSideVector = ZeroVector( system_size ); //resetting RHS

    //2.- Initialize local variables
    unsigned int voigt_size = dimension * (dimension +1) * 0.5;
    Vector StrainVector(voigt_size);
    Vector StressVector(voigt_size);
    Matrix ConstitutiveMatrix(voigt_size, voigt_size);
    Matrix B(voigt_size,voigt_size);
    Matrix DN_DX(number_of_nodes, voigt_size);
    
    //deffault values for the infinitessimal theory
    double detF = 1;
    Matrix F = identity_matrix<double>(dimension);

    //3.-Calculate elemental system:

    //reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //get the shape functions [N] (for the order of the default integration method)
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

    //get the shape functions parent coodinates derivative [dN/d£] (for the order of the default integration method)
    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    //calculate delta position (here coincides with the current displacement)
    Matrix DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    GeometryType::JacobiansType J = GetGeometry().Jacobian( J, mThisIntegrationMethod, DeltaPosition );

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {
	//a.-compute element kinematics

	//calculating the inverse of the jacobian for this integration point[d£/dx_n]
	Matrix InvJ;
	double detJ;
	MathUtils<double>::InvertMatrix( J[PointNumber], InvJ, detJ);

	if(detJ<0)
	  KRATOS_THROW_ERROR( std::invalid_argument," SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0 ) detJ = ", detJ )
	    
	//compute cartesian derivatives for this integration point  [dN/dx_n]
	noalias( DN_DX ) = prod( DN_De[PointNumber] , InvJ );
	
	//set shape functions for this integration point
        Vector N=row( Ncontainer, PointNumber);


 	//b.-compute infinitessimal strain

	this->CalculateInfinitesimalStrain(StrainVector, DN_DX);

	//c.-set general variables to constitutivelaw structure

	//create constitutive law parameters structure:
	ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

	//set constitutive law variables: (it passes only references to this local variables)
	Values.SetStrainVector(StrainVector);
	Values.SetStressVector(StressVector);
	Values.SetConstitutiveMatrix(ConstitutiveMatrix);
	Values.SetShapeFunctionsDerivatives(DN_DX);
	Values.SetShapeFunctionsValues(N);
	//values to be set:
	Values.SetDeterminantF(detF);
	Values.SetDeformationGradientF(F);

	//set constitutive law flags:
	Flags &ConstitutiveLawOptions=Values.GetOptions();

	//compute stress and constitutive matrix
	ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
	ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

	//CALL THE CONSTITUTIVE LAW (for this integration point)
	//(after calling the constitutive law StressVector and ConstitutiveMatrix are set and can be used)
	// see consitutive_law.h
        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

	//d.- compute system component integrals

        //calculating weights for integration on the "reference configuration"
        double IntegrationWeight = integration_points[PointNumber].Weight() * detJ;

	if( dimension == 2 ) IntegrationWeight *= GetProperties()[THICKNESS];

	//compute the deformation matrix B
	this->CalculateDeformationMatrix(B,DN_DX);

	
	//compute and add stiffness matrix (LHS = rLeftHandSideMatrix = K)
	noalias( rLeftHandSideMatrix ) += prod( trans( B ),  IntegrationWeight * Matrix( prod( ConstitutiveMatrix, B ) ) );


	//compute and add external forces 
	Vector VolumeForce  = this->CalculateVolumeForce( VolumeForce, N );

	for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
	    int index = dimension * i;
	    for ( unsigned int j = 0; j < dimension; j++ )
	      {
		rRightHandSideVector[index + j]  += IntegrationWeight * N[i] * VolumeForce[j];
		
	      }
	  }

	//compute and add internal forces (RHS = rRightHandSideVector = Fext - Fint)
	noalias( rRightHandSideVector ) -= IntegrationWeight * prod( trans( B ), StressVector );

      }


    KRATOS_CATCH( "" )
}


//TWO METHODS ACCOUNTING FOR EXPLICIT CASES:
//AFTER COMPUTING RHS VECTORS IT WILL BE SET EXPLICITLY AS NODAL VECTORS

//***********************************************************************************
//***********************************************************************************

void LinearSolidElement::AddExplicitContribution(const VectorType& rRHSVector, 
						    const Variable<VectorType>& rRHSVariable, 
						    Variable<array_1d<double,3> >& rDestinationVariable, 
						    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( rRHSVariable == EXTERNAL_FORCES_VECTOR && rDestinationVariable == EXTERNAL_FORCE )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		ExternalForce[j] += rRHSVector[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }

    if( rRHSVariable == INTERNAL_FORCES_VECTOR && rDestinationVariable == INTERNAL_FORCE )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &InternalForce = GetGeometry()[i].FastGetSolutionStepValue(INTERNAL_FORCE);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		InternalForce[j] += rRHSVector[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }


    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		ForceResidual[j] += rRHSVector[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }

    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************

void LinearSolidElement::ClearNodalForces()
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      if( GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_FORCE) && GetGeometry()[i].SolutionStepsDataHas(INTERNAL_FORCE) ){

        array_1d<double, 3 > & ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
        array_1d<double, 3 > & InternalForce = GetGeometry()[i].FastGetSolutionStepValue(INTERNAL_FORCE);
  
    	GetGeometry()[i].SetLock();
        ExternalForce.clear();
        InternalForce.clear();
    	GetGeometry()[i].UnSetLock();

      }

    }

    KRATOS_CATCH( "" )
}


//*************************COMPUTE DELTA POSITION*************************************
//************************************************************************************


Matrix& LinearSolidElement::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY

    //KRATOS NODAL CURRENT POSITION (X = X0 + DISPLACEMENT_X) IS ALWAYS COMPUTED
    GeometryType& geom = GetGeometry();
    const unsigned int number_of_nodes = geom.PointsNumber();
    unsigned int dimension = geom.WorkingSpaceDimension();

    rDeltaPosition = zero_matrix<double>( number_of_nodes , dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        rDeltaPosition(i, 0) = GetGeometry()[i].X() - GetGeometry()[i].X0();
	rDeltaPosition(i, 1) = GetGeometry()[i].Y() - GetGeometry()[i].Y0();
	if(dimension == 3)
	  rDeltaPosition(i, 2) = GetGeometry()[i].Z() - GetGeometry()[i].Z0();
    }

    return rDeltaPosition;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LinearSolidElement::CalculateInfinitesimalStrain( Vector& rStrainVector, const Matrix& rDN_DX )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    Matrix H = zero_matrix<double> ( dimension ); //[dU/dx_n]

    if( dimension == 2 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {

            array_1d<double, 3 > & Displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

            H ( 0 , 0 ) += Displacement[0]*rDN_DX ( i , 0 );
            H ( 0 , 1 ) += Displacement[0]*rDN_DX ( i , 1 );
            H ( 1 , 0 ) += Displacement[1]*rDN_DX ( i , 0 );
            H ( 1 , 1 ) += Displacement[1]*rDN_DX ( i , 1 );
        }


        //Infinitesimal Strain Calculation
        if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

        rStrainVector[0] = H( 0, 0 );
        rStrainVector[1] = H( 1, 1 );
        rStrainVector[2] = ( H( 0, 1 ) + H( 1, 0 ) ); // xy

    }
    else if( dimension == 3 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            array_1d<double, 3 > & Displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

            H ( 0 , 0 ) += Displacement[0]*rDN_DX ( i , 0 );
            H ( 0 , 1 ) += Displacement[0]*rDN_DX ( i , 1 );
            H ( 0 , 2 ) += Displacement[0]*rDN_DX ( i , 2 );

            H ( 1 , 0 ) += Displacement[1]*rDN_DX ( i , 0 );
            H ( 1 , 1 ) += Displacement[1]*rDN_DX ( i , 1 );
            H ( 1 , 2 ) += Displacement[1]*rDN_DX ( i , 2 );

            H ( 2 , 0 ) += Displacement[2]*rDN_DX ( i , 0 );
            H ( 2 , 1 ) += Displacement[2]*rDN_DX ( i , 1 );
            H ( 2 , 2 ) += Displacement[2]*rDN_DX ( i , 2 );
        }


        //Infinitesimal Strain Calculation
        if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

        rStrainVector[0] = H( 0, 0 );
        rStrainVector[1] = H( 1, 1 );
        rStrainVector[2] = H( 2, 2 );
        rStrainVector[3] = ( H( 0, 1 ) + H( 1, 0 ) ); // xy
        rStrainVector[4] = ( H( 1, 2 ) + H( 2, 1 ) ); // yz
        rStrainVector[5] = ( H( 0, 2 ) + H( 2, 0 ) ); // zx
    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

    }


    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************

void LinearSolidElement::CalculateDeformationMatrix(Matrix& rB, const Matrix& rDN_DX)
{
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rB.clear();
    
    if( dimension == 2 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 2 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 0, index + 1 ) = 0.0;
            rB( 1, index + 0 ) = 0.0;
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rDN_DX( i, 1 );
            rB( 2, index + 1 ) = rDN_DX( i, 0 );

        }

    }
    else if( dimension == 3 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 3 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 2 ) = rDN_DX( i, 2 );

            rB( 3, index + 0 ) = rDN_DX( i, 1 );
            rB( 3, index + 1 ) = rDN_DX( i, 0 );

            rB( 4, index + 1 ) = rDN_DX( i, 2 );
            rB( 4, index + 2 ) = rDN_DX( i, 1 );

            rB( 5, index + 0 ) = rDN_DX( i, 2 );
            rB( 5, index + 2 ) = rDN_DX( i, 0 );

        }

    }
    else
    {

        KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

    }

    KRATOS_CATCH( "" )
}


//************************************CALCULATE VOLUME ACCELERATION*******************
//************************************************************************************

Vector& LinearSolidElement::CalculateVolumeForce( Vector& rVolumeForce, const Vector& rN )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rVolumeForce = ZeroVector(dimension);
    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
      if( GetGeometry()[j].SolutionStepsDataHas(VOLUME_ACCELERATION) ) // it must be checked once at the begining only
	rVolumeForce += rN[j] * GetGeometry()[j].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    rVolumeForce *= GetProperties()[DENSITY];

    return rVolumeForce;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LinearSolidElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //lumped
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int system_size = dimension * number_of_nodes;

    if ( rMassMatrix.size1() != system_size )
        rMassMatrix.resize( system_size, system_size, false );

    rMassMatrix = ZeroMatrix( system_size, system_size );

    double TotalMass = GetGeometry().DomainSize() * GetProperties()[DENSITY];

    if( dimension == 2 ) TotalMass *= GetProperties()[THICKNESS];

    Vector LumpFact = ZeroVector(number_of_nodes);

    LumpFact = GetGeometry().LumpingFactors( LumpFact );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        double temp = LumpFact[i] * TotalMass;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            unsigned int index = i * dimension + j;
            rMassMatrix( index, index ) = temp;
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LinearSolidElement::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //0.-Initialize the DampingMatrix:
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int system_size = number_of_nodes * dimension;

    if ( rDampingMatrix.size1() != system_size )
        rDampingMatrix.resize( system_size, system_size, false );

    noalias( rDampingMatrix ) = ZeroMatrix( system_size, system_size );


    //1.-Calculate StiffnessMatrix:

    MatrixType StiffnessMatrix  = Matrix();

    this->CalculateLeftHandSide( StiffnessMatrix, rCurrentProcessInfo );

    //2.-Calculate MassMatrix:

    MatrixType MassMatrix  = Matrix();

    this->CalculateMassMatrix ( MassMatrix, rCurrentProcessInfo );
    
    
    //3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) ){ //if is stored in material properties
      alpha = GetProperties()[RAYLEIGH_ALPHA];
    }
    else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) ){ //if is stored in the process info parameters
      alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
    }

    double beta  = 0;
    if( GetProperties().Has(RAYLEIGH_BETA) ){ //if is stored in material properties
      beta = GetProperties()[RAYLEIGH_BETA];
    }
    else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) ){ //if is stored in the process info parameters
      beta = rCurrentProcessInfo[RAYLEIGH_BETA];
    }

    //4.-Compose the Damping Matrix:
   
    //Rayleigh Damping Matrix: alpha*M + beta*K
    rDampingMatrix  = alpha * MassMatrix;
    rDampingMatrix += beta  * StiffnessMatrix;


    KRATOS_CATCH( "" )
}


//*******************************CALCULATE VALUE METHODS******************************
//************************************************************************************
//(see element.h)

void LinearSolidElement::CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //TO WRITE SOME VARIABLE:
    //a.- It can be stored for all integration points as a member variables an supplied here (avoid it: memory consuming)
    //b.- It can be calculated again (recomendable: is it preferable to recompute than store it)

    const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();

    if ( rOutput.size() != mConstitutiveLawVector.size() )
        rOutput.resize( mConstitutiveLawVector.size() );


    if ( rVariable == CAUCHY_STRESS_TENSOR )
    {
      const unsigned int number_of_nodes = GetGeometry().size();
      
      //2.- Initialize local variables
      unsigned int voigt_size = dimension * (dimension +1) * 0.5;
      Vector StrainVector(voigt_size);
      Vector StressVector(voigt_size);
      Matrix DN_DX(number_of_nodes, voigt_size);
    
      //deffault values for the infinitessimal theory
      double detF = 1;
      Matrix F = identity_matrix<double>(dimension);

      //get the shape functions [N] (for the order of the default integration method)
      const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

      //get the shape functions parent coodinates derivative [dN/d£] (for the order of the default integration method)
      const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

      //calculate delta position (here coincides with the current displacement)
      Matrix DeltaPosition = zero_matrix<double>( number_of_nodes , dimension);
      DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

      //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
      //std::cout<<" DeltaPosition "<<DeltaPosition<<std::endl;
      //std::cout<<" mThisIntegrationMethod "<<mThisIntegrationMethod<<std::endl;
 
      GeometryType::JacobiansType J = GetGeometry().Jacobian( J, mThisIntegrationMethod, DeltaPosition );

      //loop integration points
      for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
	  //a.-compute element kinematics

	  //calculating the inverse of the jacobian for this integration point[d£/dx_n]
	  Matrix InvJ;
	  double detJ;
	  //std::cout<<" Jacobian "<<J[PointNumber]<<std::endl;
	  MathUtils<double>::InvertMatrix( J[PointNumber], InvJ, detJ);
	  
	  if(detJ<0)
	    KRATOS_THROW_ERROR( std::invalid_argument," SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0 ) detJ = ", detJ )
	    
          //compute cartesian derivatives for this integration point  [dN/dx_n]
	  noalias( DN_DX ) = prod( DN_De[PointNumber], InvJ );
	
	  //set shape functions for this integration point
	  Vector N=row( Ncontainer, PointNumber);


	  //b.-compute infinitessimal strain

	  this->CalculateInfinitesimalStrain(StrainVector, DN_DX);

	  //c.-set general variables to constitutivelaw structure

	  //create constitutive law parameters structure:
	  ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

	  //set constitutive law variables: (it passes only references to this local variables)
	  Values.SetStrainVector(StrainVector);
	  Values.SetStressVector(StressVector);
	  Values.SetShapeFunctionsDerivatives(DN_DX);
	  Values.SetShapeFunctionsValues(N);
	  //values to be set:
	  Values.SetDeterminantF(detF);
	  Values.SetDeformationGradientF(F);

	  //set constitutive law flags:
	  Flags &ConstitutiveLawOptions=Values.GetOptions();

	  //compute stress and constitutive matrix
	  ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

	  //CALL THE CONSTITUTIVE LAW (for this integration point)
	  //(after calling the constitutive law StressVector and ConstitutiveMatrix are set and can be used)
	  // see consitutive_law.h
	  mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

	  //b.-compute infinitessimal strain
	  this->CalculateInfinitesimalStrain(StrainVector, DN_DX);

	  if ( rOutput[PointNumber].size2() != dimension )
	    rOutput[PointNumber].resize( dimension, dimension, false );
	  
	  rOutput[PointNumber] =  MathUtils<double>::StressVectorToTensor(StressVector);
        }

    }
    else if( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR)
    {

      const unsigned int number_of_nodes = GetGeometry().size();
      
      //2.- Initialize local variables
      unsigned int voigt_size = dimension * (dimension +1) * 0.5;
      Matrix DN_DX(number_of_nodes, voigt_size);
   
      //Calculate Strain (in this case the infinitessimal strain vector)

      //get the shape functions parent coodinates derivative [dN/d£] (for the order of the default integration method)
      const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
      
      //calculate delta position (here coincides with the current displacement)
      Matrix DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);
      
      //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
      GeometryType::JacobiansType J = GetGeometry().Jacobian( J, mThisIntegrationMethod, DeltaPosition );

      //loop integration points
      for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {

	  //a.-compute element kinematics

	  //calculating the inverse of the jacobian for this integration point[d£/dx_n]
	  Matrix InvJ;
	  double detJ;
	  MathUtils<double>::InvertMatrix( J[ PointNumber], InvJ, detJ);

	  if(detJ<0)
	    KRATOS_THROW_ERROR( std::invalid_argument," SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0 ) detJ = ", detJ )
	  

	  //compute cartesian derivatives for this integration point  [dN/dx_n]
	  noalias( DN_DX ) = prod( DN_De[PointNumber] , InvJ );
	
	  //b.-compute infinitessimal strain
	  Vector StrainVector;
	  this->CalculateInfinitesimalStrain(StrainVector, DN_DX);

	  if ( rOutput[PointNumber].size2() != dimension )
	    rOutput[PointNumber].resize( dimension, dimension, false );
	  
	  rOutput[PointNumber] = MathUtils<double>::StrainVectorToTensor(StrainVector);
        }
    }
    else
    {
        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        {
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable , rOutput[ii] );
        }
    }


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

int LinearSolidElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    //verify that the variables are correctly initialized

    if ( VELOCITY.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "" )

    if ( DISPLACEMENT.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )

    if ( ACCELERATION.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "" )

    if ( DENSITY.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DENSITY has Key zero! (check if the application is correctly registered", "" )


    if ( VOLUME_ACCELERATION.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "VOLUME_ACCELERATION has Key zero! (check if the application is correctly registered", "" )

    for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
    {
        if ( this->GetGeometry()[i].SolutionStepsDataHas( VOLUME_ACCELERATION ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable VOLUME_ACCELERATION on node ", this->GetGeometry()[i].Id() )
    }


    //verify that the dofs exist
    for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
    {
        if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id() )

        if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id() )
    }

    //verify that the constitutive law exists
    if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
    {
        KRATOS_THROW_ERROR( std::logic_error, "constitutive law not provided for property ", this->GetProperties().Id() )
    }

    //verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

    bool correct_strain_measure = false;
    for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    {
	    if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Infinitesimal)
		    correct_strain_measure = true;
    }

    if( correct_strain_measure == false )
	    KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the element type ", " Small Displacements " );


    //verify that the constitutive law has the correct dimension
    if ( dimension == 2 )
    {
        if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 3 )
	    KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law used. This is a 2D element! expected strain size is 3 (el id = ) ", this->Id() )

        if ( THICKNESS.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "" )

        if ( this->GetProperties().Has( THICKNESS ) == false )
            KRATOS_THROW_ERROR( std::logic_error, "THICKNESS not provided for element ", this->Id() )


    }
    else
    {
        if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 6 )
            KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) ", this->Id() )
    }

    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->Check( this->GetProperties(), this->GetGeometry(), rCurrentProcessInfo );

    return 0;

    KRATOS_CATCH( "" );
}


//************************************************************************************
//************************************************************************************

//SERIALIZER, TO STORE AND LOAD VARIABLES IN CASE OF A CALL TO THE PROBLEM RESTART

void LinearSolidElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    int IntMethod = int(mThisIntegrationMethod);
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
}

void LinearSolidElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
}


} // Namespace Kratos


