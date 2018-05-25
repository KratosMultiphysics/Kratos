//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               April 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/solid_elements/large_displacement_V_element.hpp"
#include "solid_mechanics_application_variables.h"


namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementVElement::LargeDisplacementVElement()
    : LargeDisplacementElement()
{
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementVElement::LargeDisplacementVElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : LargeDisplacementElement( NewId, pGeometry )
{
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementVElement::LargeDisplacementVElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacementElement( NewId, pGeometry, pProperties )
{
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LargeDisplacementVElement::LargeDisplacementVElement( LargeDisplacementVElement const& rOther)
    :LargeDisplacementElement(rOther)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

LargeDisplacementVElement&  LargeDisplacementVElement::operator=(LargeDisplacementVElement const& rOther)
{
    LargeDisplacementElement::operator=(rOther);

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer LargeDisplacementVElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new LargeDisplacementVElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer LargeDisplacementVElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    KRATOS_THROW_ERROR( std::logic_error, "calling the default constructor for a large displacement 3D element ... illegal operation!!", "" )

    LargeDisplacementVElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    //-----------//

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;


    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());
	
	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() )
      }
    
    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());
       
    return Element::Pointer( new LargeDisplacementVElement(NewElement) );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementVElement::~LargeDisplacementVElement()
{
}


//************* GETTING METHODS
//************************************************************************************
//************************************************************************************



void LargeDisplacementVElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    rElementalDofList.resize( 0 );

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_X ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_Y ) );

        if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_Z ) );
    }
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementVElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = GetDofsSize();

    if ( rResult.size() != dofs_size )
        rResult.resize( dofs_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension;
        rResult[index]     = GetGeometry()[i].GetDof( VELOCITY_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( VELOCITY_Y ).EquationId();

        if( dimension == 3)
          rResult[index + 2] = GetGeometry()[i].GetDof( VELOCITY_Z ).EquationId();

    }

}


//************************************************************************************
//************************************************************************************

void LargeDisplacementVElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementVariables& rVariables, double& rIntegrationWeight)
{
    KRATOS_TRY
       
    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix(); 

    // operation performed: add Km to the rLefsHandSideMatrix
    this->CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
   
    // operation performed: add Kg to the rLefsHandSideMatrix
    this->CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    rLeftHandSideMatrix *= rVariables.GetProcessInfo()[DELTA_TIME]; // backward Euler Approach (BDF order 1)  
    
    //KRATOS_WATCH( rLeftHandSideMatrix )
  
    KRATOS_CATCH( "" )  
}




//************************************************************************************
//************************************************************************************

unsigned int LargeDisplacementVElement::GetDofsSize()
{
  KRATOS_TRY
     
  const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
  const unsigned int number_of_nodes = GetGeometry().PointsNumber();    
  
  unsigned int size = number_of_nodes * dimension; //usual size for velocity based elements
  
  return size;   
  
  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementVElement::SetElementVariables(ElementVariables& rVariables,
                                                    ConstitutiveLaw::Parameters& rValues,
                                                    const int & rPointNumber)
{

    //to take in account previous step for output print purposes
    unsigned int step = 0;
    if( mFinalizedStep ){
      step = 1;
      this->GetHistoricalVariables(rVariables,rPointNumber);
    }
  
    if(rVariables.detF<0){
        
	std::cout<<" Element: "<<this->Id()<<std::endl;

	unsigned int number_of_nodes = GetGeometry().PointsNumber();

	for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
	    array_1d<double, 3> & CurrentPosition  = GetGeometry()[i].Coordinates();
	    array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
	    array_1d<double, 3> & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
	    array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
	    std::cout<<" NODE ["<<GetGeometry()[i].Id()<<"]: "<<PreviousPosition<<" (Cur: "<<CurrentPosition<<") "<<std::endl;
	    std::cout<<" ---Disp: "<<CurrentDisplacement<<" (Pre: "<<PreviousDisplacement<<")"<<std::endl;
	  }

	for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
	    if( GetGeometry()[i].SolutionStepsDataHas(CONTACT_FORCE) ){
	      array_1d<double, 3 > & PreContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE,1);
	      array_1d<double, 3 > & ContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
	      std::cout<<" ---Contact_Force: (Pre:"<<PreContactForce<<", Cur:"<<ContactForce<<") "<<std::endl;
	    }
	    else{
	      std::cout<<" ---Contact_Force: NULL "<<std::endl;
	    }
	  }
	
        KRATOS_THROW_ERROR( std::invalid_argument," LARGE DISPLACEMENT ELEMENT INVERTED: |F|<0  detF = ", rVariables.detF )
    }



    //Compute strain rate measures if they are required by the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    mConstitutiveLawVector[rPointNumber]->GetLawFeatures(LawFeatures);
    
    bool strain_rate_measure = false;
    for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    {
      if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Velocity_Gradient)
	strain_rate_measure = true;
    }

    if( strain_rate_measure ){     
      //Compute symmetric spatial velocity gradient [DN_DX = dN/dx_n*1] stored in a vector
      this->CalculateVelocityGradientVector( rVariables.StrainVector, rVariables.DN_DX, step );
      Flags &ConstitutiveLawOptions=rValues.GetOptions();
      ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    }

    //Compute F and detF (from 0 to n+1) : store it in H variable and detH
    rVariables.detH = rVariables.detF * rVariables.detF0;
    noalias(rVariables.H) = prod( rVariables.F, rVariables.F0 );
    
    rValues.SetDeterminantF(rVariables.detH);
    rValues.SetDeformationGradientF(rVariables.H);
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
    rValues.SetShapeFunctionsValues(rVariables.N);

}

//************************************************************************************
//************************************************************************************

void LargeDisplacementVElement::CalculateVelocityGradient(Matrix& rH,
                                                          const Matrix& rDN_DX,
                                                          unsigned int step)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rH = zero_matrix<double> ( dimension );

    if( dimension == 2 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            array_1d<double,3>& rCurrentVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,step);
          
            rH ( 0 , 0 ) += rCurrentVelocity[0]*rDN_DX ( i , 0 );
            rH ( 0 , 1 ) += rCurrentVelocity[0]*rDN_DX ( i , 1 );
            rH ( 1 , 0 ) += rCurrentVelocity[1]*rDN_DX ( i , 0 );
            rH ( 1 , 1 ) += rCurrentVelocity[1]*rDN_DX ( i , 1 );
        }

    }
    else if( dimension == 3)
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
          array_1d<double,3>& rCurrentVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,step);
            
            rH ( 0 , 0 ) += rCurrentVelocity[0]*rDN_DX ( i , 0 );
            rH ( 0 , 1 ) += rCurrentVelocity[0]*rDN_DX ( i , 1 );
            rH ( 0 , 2 ) += rCurrentVelocity[0]*rDN_DX ( i , 2 );
            rH ( 1 , 0 ) += rCurrentVelocity[1]*rDN_DX ( i , 0 );
            rH ( 1 , 1 ) += rCurrentVelocity[1]*rDN_DX ( i , 1 );
            rH ( 1 , 2 ) += rCurrentVelocity[1]*rDN_DX ( i , 2 );
            rH ( 2 , 0 ) += rCurrentVelocity[2]*rDN_DX ( i , 0 );
            rH ( 2 , 1 ) += rCurrentVelocity[2]*rDN_DX ( i , 1 );
            rH ( 2 , 2 ) += rCurrentVelocity[2]*rDN_DX ( i , 2 );
        }

    }
    else
    {
      KRATOS_ERROR << " something is wrong with the dimension when computing velocity gradient " << std::endl;
    }
    
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementVElement::CalculateVelocityGradientVector(Vector& rH,
                                                                const Matrix& rDN_DX,
                                                                unsigned int step)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rH = ZeroVector( dimension * dimension );

    if( dimension == 2 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            array_1d<double,3>& rCurrentVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,step);
          
            rH[0] += rCurrentVelocity[0]*rDN_DX ( i , 0 );
            rH[1] += rCurrentVelocity[1]*rDN_DX ( i , 1 );
            rH[2] += rCurrentVelocity[0]*rDN_DX ( i , 1 );
            rH[3] += rCurrentVelocity[1]*rDN_DX ( i , 0 );

        }

    }
    else if( dimension == 3)
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
          array_1d<double,3>& rCurrentVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,step);
            
            rH[0] += rCurrentVelocity[0]*rDN_DX ( i , 0 );
            rH[1] += rCurrentVelocity[1]*rDN_DX ( i , 1 );
            rH[2] += rCurrentVelocity[2]*rDN_DX ( i , 2 );
            
            rH[3] += rCurrentVelocity[0]*rDN_DX ( i , 1 );
            rH[4] += rCurrentVelocity[1]*rDN_DX ( i , 2 );
            rH[5] += rCurrentVelocity[2]*rDN_DX ( i , 0 );

            rH[6] += rCurrentVelocity[1]*rDN_DX ( i , 0 );
            rH[7] += rCurrentVelocity[2]*rDN_DX ( i , 1 );
            rH[8] += rCurrentVelocity[0]*rDN_DX ( i , 2 );
        }

    }
    else
    {
      KRATOS_ERROR << " something is wrong with the dimension when computing velocity gradient " << std::endl;
    }
    
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementVElement::CalculateSymmetricVelocityGradient(const Matrix& rH,
                                                                  Vector& rStrainVector)
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
    {
        if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

        rStrainVector[0] = rH( 0, 0 );
        rStrainVector[1] = rH( 1, 1 );
        rStrainVector[2] = (rH( 0, 1 ) + rH( 1, 0 )); // xy

    }
    else if( dimension == 3 )
    {
        if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

        rStrainVector[0] = rH( 0, 0 );
        rStrainVector[1] = rH( 1, 1 );
        rStrainVector[2] = rH( 2, 2 );
        rStrainVector[3] = ( rH( 0, 1 ) + rH( 1, 0 ) ); // xy
        rStrainVector[4] = ( rH( 1, 2 ) + rH( 2, 1 ) ); // yz
        rStrainVector[5] = ( rH( 0, 2 ) + rH( 2, 0 ) ); // xz

    }
    else
    {
        KRATOS_ERROR << " something is wrong with the dimension symmetric velocity gradient " << std::endl;
    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementVElement::CalculateSkewSymmetricVelocityGradient(const Matrix& rH,
                                                                      Vector& rStrainVector)
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
    {
        if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

        rStrainVector[0] = 0.0;
        rStrainVector[1] = 0.0;
        rStrainVector[2] = (rH( 0, 1 ) - rH( 1, 0 )); // xy

    }
    else if( dimension == 3 )
    {
        if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

        rStrainVector[0] = 0.0;
        rStrainVector[1] = 0.0;
        rStrainVector[2] = 0.0;
        rStrainVector[3] = ( rH( 0, 1 ) - rH( 1, 0 ) ); // xy
        rStrainVector[4] = ( rH( 1, 2 ) - rH( 2, 1 ) ); // yz
        rStrainVector[5] = ( rH( 0, 2 ) - rH( 2, 0 ) ); // xz

    }
    else
    {
        KRATOS_ERROR << " something is wrong with the dimension symmetric velocity gradient " << std::endl;
    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

int  LargeDisplacementVElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = LargeDisplacementElement::Check(rCurrentProcessInfo);

    // Check that the element nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
      {
	// Nodal data
	Node<3> &rNode = this->GetGeometry()[i];
	KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rNode);
	//KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUME_ACCELERATION,rNode);
	
	// Nodal dofs
	KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X,rNode);
	KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y,rNode);
	if( rCurrentProcessInfo[SPACE_DIMENSION] == 3)
	  KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Z,rNode);
      }
    // Check compatibility with the constitutive law

    // Check that all required variables have been registered

    return ErrorCode;

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************


void LargeDisplacementVElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
}

void LargeDisplacementVElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
}


} // Namespace Kratos


