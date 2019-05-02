//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/beam_elements/beam_element.hpp"

#include "solid_mechanics_application_variables.h"


namespace Kratos
{

  /**
   * Flags related to the element computation
   */
  KRATOS_CREATE_LOCAL_FLAG( BeamElement, COMPUTE_RHS_VECTOR,                 0 );
  KRATOS_CREATE_LOCAL_FLAG( BeamElement, COMPUTE_LHS_MATRIX,                 1 );
  KRATOS_CREATE_LOCAL_FLAG( BeamElement, FINALIZED_STEP,                     2 );

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  BeamElement::BeamElement(IndexType NewId,GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
  {

  }

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************


  BeamElement::BeamElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
  {
    this->Set(STRUCTURE);
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  BeamElement::BeamElement(BeamElement const& rOther)
    :Element(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    ,mInitialLocalQuaternion(rOther.mInitialLocalQuaternion)
  {
  }

  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Element::Pointer BeamElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
  {
    return Kratos::make_shared<BeamElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  BeamElement::~BeamElement()
  {
  }


  //************************************************************************************
  //************************************************************************************

  BeamElement::IntegrationMethod  BeamElement::GetIntegrationMethod() const
  {
    return mThisIntegrationMethod;
  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::IncreaseIntegrationMethod(IntegrationMethod& rThisIntegrationMethod, unsigned int increment) const
  {
    int IntMethod = int(rThisIntegrationMethod);
    IntMethod += increment;
    rThisIntegrationMethod = IntegrationMethod(IntMethod);
  }

  //************************************************************************************
  //************************************************************************************

  void BeamElement::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
  {
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    rElementalDofList.resize(0);

    for ( SizeType i = 0; i < GetGeometry().size(); i++ )
      {

	rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
	rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));

	if( dimension == 2 ){
	  rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
	}
	else{
          rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
	  rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
	  rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
	  rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
	}

      }
  }

  //************************************************************************************
  //************************************************************************************

  void BeamElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
  {

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = number_of_nodes * ( (dimension-1) * 3 );

    if ( rResult.size() != dofs_size )
      rResult.resize( dofs_size, false );

    unsigned int index = 0;

    if( dimension == 2 ){
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  index = i * ( (dimension-1) * 3 );
	  rResult[index]   = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
	  rResult[index+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
	  rResult[index+2] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
	}
    }
    else{
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  index = i * ( (dimension-1) * 3 );
	  rResult[index]   = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
	  rResult[index+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
	  rResult[index+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();

	  rResult[index+3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
	  rResult[index+4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
	  rResult[index+5] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
	}
    }

  }


  //*********************************DISPLACEMENT***************************************
  //************************************************************************************

  void BeamElement::GetValuesVector(Vector& rValues, int Step)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = number_of_nodes * ( (dimension-1) * 3 );

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    unsigned int index = 0;
    if( dimension == 2 ){
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  index = i * ( (dimension-1) * 3 );
	  rValues[index]   = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
	  rValues[index+1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
	  rValues[index+2] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Z, Step );
	}
    }
    else{
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  index = i * ( (dimension-1) * 3 );
	  rValues[index]   = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
	  rValues[index+1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
	  rValues[index+2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );

	  rValues[index+3] = GetGeometry()[i].GetSolutionStepValue( ROTATION_X, Step );
	  rValues[index+4] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Y, Step );
	  rValues[index+5] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Z, Step );
	}
    }

    KRATOS_CATCH( "" )
  }

  //************************************VELOCITY****************************************
  //************************************************************************************

  void BeamElement::GetFirstDerivativesVector(Vector& rValues, int Step)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = number_of_nodes * ( (dimension-1) * 3 );

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    unsigned int index = 0;
    if( dimension == 2 ){
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  index = i * ( (dimension-1) * 3 );
	  rValues[index]   = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
	  rValues[index+1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
	  rValues[index+2] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_Z, Step );
	}
    }
    else{
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  index = i * ( (dimension-1) * 3 );
	  rValues[index]   = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
	  rValues[index+1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
	  rValues[index+2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );

	  rValues[index+3] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_X, Step );
	  rValues[index+4] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_Y, Step );
	  rValues[index+5] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_Z, Step );
	}
    }

    KRATOS_CATCH( "" )
  }

  //*********************************ACCELERATION***************************************
  //************************************************************************************

  void BeamElement::GetSecondDerivativesVector(Vector& rValues, int Step)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = number_of_nodes * ( (dimension-1) * 3 );

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    unsigned int index = 0;
    if( dimension == 2 ){

      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  index = i * ( (dimension-1) * 3 );
	  rValues[index]   = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
	  rValues[index+1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );
	  rValues[index+2] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Z, Step );
	}
    }
    else{
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  index = i * ( (dimension-1) * 3 );
	  rValues[index]   = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
	  rValues[index+1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );
	  rValues[index+2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );

	  rValues[index+3] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_X, Step );
	  rValues[index+4] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Y, Step );
	  rValues[index+5] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Z, Step );
	}
    }

    KRATOS_CATCH( "" )
  }


  //*********************************GET DOUBLE VALUE***********************************
  //************************************************************************************

  void  BeamElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
						  std::vector<double>& rValues,
						  const ProcessInfo& rCurrentProcessInfo )
  {
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
  }

  //**********************************GET VECTOR VALUE**********************************
  //************************************************************************************

  void BeamElement::GetValueOnIntegrationPoints( const Variable<array_1d<double, 3 > >& rVariable,
						 std::vector< array_1d<double, 3 > >& rValues,
						 const ProcessInfo& rCurrentProcessInfo )
  {
    this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
  }

  //************************************************************************************
  //************************************************************************************

  void BeamElement::Initialize()
  {
    KRATOS_TRY

    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    Matrix LocalTransformationMatrix(dimension,dimension);
    noalias(LocalTransformationMatrix) = ZeroMatrix(dimension,dimension);
    this->CalculateLocalAxesMatrix( LocalTransformationMatrix );

    //Quaternions initialization
    mInitialLocalQuaternion = QuaternionType::FromRotationMatrix( LocalTransformationMatrix );

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    this->Set(BeamElement::FINALIZED_STEP,false);

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void BeamElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    this->Set(BeamElement::FINALIZED_STEP,true);

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  unsigned int BeamElement::GetDofsSize()
  {
    KRATOS_TRY

    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    const SizeType number_of_nodes  = GetGeometry().PointsNumber();

    unsigned int size = number_of_nodes * (dimension-1) * 3;

    return size;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
					     VectorType& rRightHandSideVector,
					     Flags& rCalculationFlags)

  {
    KRATOS_TRY

    const unsigned int MatSize = this->GetDofsSize();

    if ( rCalculationFlags.Is(BeamElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
      {
        if ( rLeftHandSideMatrix.size1() != MatSize )
	  rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
      }

    if ( rCalculationFlags.Is(BeamElement::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
      {
        if ( rRightHandSideVector.size() != MatSize )
	  rRightHandSideVector.resize( MatSize, false );

	noalias(rRightHandSideVector) = ZeroVector( MatSize ); //resetting RHS

      }

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  Vector& BeamElement::MapToInitialLocalFrame(Vector& rVariable, unsigned int PointNumber)
  {
    KRATOS_TRY

    BeamMathUtilsType::MapToCurrentLocalFrame(mInitialLocalQuaternion, rVariable);

    return rVariable;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void BeamElement::MapLocalToGlobal(ElementDataType& rVariables, MatrixType& rMatrix)
  {
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
      BeamMathUtilsType::MapLocalToGlobal2D(mInitialLocalQuaternion, rMatrix);
    else
      BeamMathUtilsType::MapLocalToGlobal3D(mInitialLocalQuaternion, rMatrix);

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::MapLocalToGlobal(ElementDataType& rVariables, VectorType& rVector)
  {
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
      BeamMathUtilsType::MapLocalToGlobal2D(mInitialLocalQuaternion, rVector);
    else
      BeamMathUtilsType::MapLocalToGlobal3D(mInitialLocalQuaternion, rVector);

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  Vector& BeamElement::GetLocalCurrentValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const Vector& rN)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);

    //strains due to displacements and rotations
    for ( SizeType i = 0; i <number_of_nodes ; i++ )
      {
	CurrentValueVector = GetNodalCurrentValue( rVariable, CurrentValueVector, i );
	for( SizeType j = 0; j < dimension; j++ )
	  rValue[j] += rN[i] * CurrentValueVector[j];
      }

    //Current Frame is the Local Frame
    rValue = this->MapToInitialLocalFrame( rValue );

    return rValue;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  Vector& BeamElement::GetLocalPreviousValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const Vector& rN)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    Vector PreviousValueVector(3);
    noalias(PreviousValueVector) = ZeroVector(3);

    //strains due to displacements and rotations
    for( SizeType i = 0; i < number_of_nodes ; i++ )
      {
	PreviousValueVector = GetNodalPreviousValue( rVariable, PreviousValueVector, i );

	for( SizeType j = 0; j < dimension; j++ )
	  rValue[j] += rN[i] * PreviousValueVector[j];
      }

    //Current Frame is the Local Frame
    rValue = this->MapToInitialLocalFrame( rValue );

    return rValue;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  Vector& BeamElement::GetNodalCurrentValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode)
  {
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    if( rValue.size() != dimension )
      rValue.resize(dimension, false);

    rValue = GetGeometry()[rNode].FastGetSolutionStepValue( rVariable );

    return rValue;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  Vector& BeamElement::GetNodalPreviousValue(const Variable<array_1d<double,3> >&rVariable, Vector& rValue, const unsigned int& rNode)
  {
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    if( rValue.size() != dimension )
      rValue.resize(dimension, false);

    rValue = GetGeometry()[rNode].FastGetSolutionStepValue( rVariable, 1 );

    return rValue;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::InitializeElementData(ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    const unsigned int voigt_size      = dimension * (dimension +1 ) * 0.5;

    rVariables.Initialize(voigt_size,dimension,number_of_nodes);

    //Compute Section Properties:
    this->CalculateSectionProperties(rVariables.Section);

    rVariables.Length = GetGeometry().Length();

    if(rVariables.Length == 0.00)
      KRATOS_ERROR << "Zero length found in element #" << this->Id() << std::endl;


    //set variables including all integration points values

    //reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

    //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );

    //Calculate Delta Position
    ElementUtilities::CalculateDeltaPosition(rVariables.DeltaPosition,this->GetGeometry());

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

    KRATOS_CATCH( "" )
  }


  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void BeamElement::CalculateKinematics(ElementDataType& rVariables, const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    KRATOS_ERROR << " calling the default method CalculateKinematics for a beam element " << std::endl;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************


  void BeamElement::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
					      ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    //create and initialize element variables:
    ElementDataType Variables;
    this->InitializeElementData(Variables,rCurrentProcessInfo);

    //reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //auxiliary terms
    Vector VolumeForce(3);
    noalias(VolumeForce) = ZeroVector(3);

    double IntegrationWeight = 1.0;

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {
        //compute element kinematics  ...
        this->CalculateKinematics(Variables,PointNumber);

	//compute element ConstitutiveTensor
	this->CalculateConstitutiveMatrix(Variables);

	//compute element Strain and Stress Resultants and Couples
	this->CalculateStressResultants(Variables, PointNumber);

	IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
	IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );


	if ( rLocalSystem.CalculationFlags.Is(BeamElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
	  {
	    this->CalculateAndAddLHS( rLocalSystem, Variables, IntegrationWeight );


	    //std::cout<<"["<<this->Id()<<"] Beam Rotated rLeftHandSideMatrix "<<rLocalSystem.GetLeftHandSideMatrix()<<std::endl;
	  }

	if ( rLocalSystem.CalculationFlags.Is(BeamElement::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
	  {
	    //contribution to external forces
	    VolumeForce  = this->CalculateVolumeForce( VolumeForce, Variables.N );

	    this->CalculateAndAddRHS( rLocalSystem , Variables, VolumeForce, IntegrationWeight );

	    //std::cout<<"["<<this->Id()<<"] Beam Rotated rRightHandSideVector "<<rLocalSystem.GetRightHandSideVector()<<std::endl;
	  }

      }


    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateDynamicSystem( LocalSystemComponents& rLocalSystem,
					    ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY


    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateSectionProperties(SectionProperties & rSection)
  {
    KRATOS_TRY

      if( GetProperties().Has(CROSS_SECTION_AREA) ){
        rSection.Area = GetProperties()[CROSS_SECTION_AREA];
      }
      else{
        rSection.Area = GetValue(CROSS_SECTION_AREA);
      }


    if( GetProperties().Has(LOCAL_INERTIA_TENSOR) )
      {
        Matrix& inertia = GetProperties()[LOCAL_INERTIA_TENSOR];
        rSection.Inertia_z = inertia(0,0);
        rSection.Inertia_y = inertia(1,1);
        rSection.Polar_Inertia = inertia(0,1);
      }
    else
      {
        Matrix& inertia = GetValue(LOCAL_INERTIA_TENSOR);
        rSection.Inertia_z = inertia(0,0);
        rSection.Inertia_y = inertia(1,1);
        rSection.Polar_Inertia = inertia(0,1);
      }

    rSection.Rotational_Inertia = rSection.Polar_Inertia;

    KRATOS_CATCH( "" )

  }

  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateInertiaDyadic(SectionProperties& rSection, Matrix& rInertiaDyadic)
  {
    KRATOS_TRY

    if( rInertiaDyadic.size1() != 3 )
      rInertiaDyadic.resize(3, 3, false);

    noalias(rInertiaDyadic) = ZeroMatrix(3,3);

    //if the local axes are the principal axes of the cross section

    //Axis Local E1
    rInertiaDyadic(0,0) = rSection.Rotational_Inertia; //beam axis
    rInertiaDyadic(1,1) = rSection.Inertia_y; //vertial axis
    rInertiaDyadic(2,2) = rSection.Inertia_z; //horizontal axis


    rInertiaDyadic *= GetProperties()[DENSITY];

    KRATOS_CATCH( "" )
  }


  //************************************CALCULATE TOTAL MASS****************************
  //************************************************************************************

  double& BeamElement::CalculateTotalMass( SectionProperties& Section, double& rTotalMass )
  {
    KRATOS_TRY

    rTotalMass = ( Section.Area ) * GetProperties()[DENSITY];

    return rTotalMass;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateMaterialConstitutiveMatrix(Matrix& rConstitutiveMatrix, ElementDataType& rVariables)
  {
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 ){

     if( rConstitutiveMatrix.size1() != 3 || rConstitutiveMatrix.size2() != 3)
	rConstitutiveMatrix.resize(3,3, false);

     noalias(rConstitutiveMatrix) = ZeroMatrix(3,3);

     if(GetProperties().Has(LOCAL_CONSTITUTIVE_MATRIX)){

       //Axis local E1
       rConstitutiveMatrix = GetProperties()[LOCAL_CONSTITUTIVE_MATRIX];

     }
     else{

       const double PoissonCoefficient = GetProperties()[POISSON_RATIO];
       const double YoungModulus       = GetProperties()[YOUNG_MODULUS];
       const double ShearModulus       = YoungModulus*0.5/(1.0 + PoissonCoefficient);

       //Axis local E1
       rConstitutiveMatrix( 0, 0 ) = YoungModulus * rVariables.Section.Area;         //local beam axis
       rConstitutiveMatrix( 1, 1 ) = ShearModulus * rVariables.Section.Area;         //local vertical axis
       rConstitutiveMatrix( 2, 2 ) = YoungModulus * rVariables.Section.Inertia_z;    //local horizontal axis

     }

    }
    else{

      if( rConstitutiveMatrix.size1() != 6 || rConstitutiveMatrix.size2() != 6)
	rConstitutiveMatrix.resize(6,6, false);

      noalias(rConstitutiveMatrix) = ZeroMatrix(6,6);

      if(GetProperties().Has(LOCAL_CONSTITUTIVE_MATRIX)){

	//Axis local E1
	rConstitutiveMatrix = GetProperties()[LOCAL_CONSTITUTIVE_MATRIX];

      }
      else{

	const double PoissonCoefficient = GetProperties()[POISSON_RATIO];
	const double YoungModulus       = GetProperties()[YOUNG_MODULUS];
	const double ShearModulus       = YoungModulus*0.5/(1.0 + PoissonCoefficient);


	//Axis local E1
	rConstitutiveMatrix( 0, 0 ) = YoungModulus * rVariables.Section.Area;          //local beam axis
	rConstitutiveMatrix( 1, 1 ) = ShearModulus * rVariables.Section.Area;          //local vertial axis
	rConstitutiveMatrix( 2, 2 ) = ShearModulus * rVariables.Section.Area;          //local horizontal axis

	rConstitutiveMatrix( 3, 3 ) = ShearModulus * rVariables.Section.Polar_Inertia; //local torsion
	rConstitutiveMatrix( 4, 4 ) = YoungModulus * rVariables.Section.Inertia_y;     //local vertical axis
	rConstitutiveMatrix( 5, 5 ) = YoungModulus * rVariables.Section.Inertia_z;     //local horizontal axis

      }

    }
    //std::cout<<" ConstitutiveMatrix "<<rConstitutiveMatrix<<std::endl;

    KRATOS_CATCH( "" )
  }



  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateConstitutiveMatrix(ElementDataType& rVariables)
  {
    KRATOS_TRY

    //Material Elastic constitutive matrix
    this->CalculateMaterialConstitutiveMatrix(rVariables.ConstitutiveMatrix, rVariables);

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateStressResultants(ElementDataType& rVariables, const unsigned int& rPointNumber)
  {
    KRATOS_TRY

    KRATOS_ERROR << " calling the default method CalculateStressResultants for a beam element " << std::endl;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  double& BeamElement::CalculateIntegrationWeight(double& rIntegrationWeight)
  {
    KRATOS_TRY

    return rIntegrationWeight;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
  {
    KRATOS_TRY

    Flags       LocalFlags;
    MatrixType  LocalLeftHandSideMatrix;
    VectorType  LocalRightHandSideVector;

    //Initialize sizes for the local system components:
    LocalFlags.Set(BeamElement::COMPUTE_LHS_MATRIX);
    this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, LocalFlags );

    // Local material stiffness
    this->CalculateAndAddKuum( LocalLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Local geometrical stiffness
    this->CalculateAndAddKuug( LocalLeftHandSideMatrix, rVariables, rIntegrationWeight );

    //// Local geometrical stiffness (not used in quasi-static additive rotations it degradates convergence)
    //if( mIterationCounter < 2 || mIterationCounter > 5 ) //TEST
    //  this->CalculateAndAddKuug( LocalLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // LocalToGlobalSystem for the correct assembly
    this->MapLocalToGlobal(rVariables, LocalLeftHandSideMatrix);

    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();
    rLeftHandSideMatrix += LocalLeftHandSideMatrix;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
  {
    KRATOS_TRY

    Flags       LocalFlags;
    MatrixType  LocalLeftHandSideMatrix;
    VectorType  LocalRightHandSideVector;

    //Initialize sizes for the local system components:
    LocalFlags.Set(BeamElement::COMPUTE_RHS_VECTOR);
    this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, LocalFlags );

    // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
    this->CalculateAndAddExternalForces( LocalRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

    // operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
    this->CalculateAndAddInternalForces( LocalRightHandSideVector, rVariables, rIntegrationWeight );

    // LocalToGlobalSystem for the correct assembly
    this->MapLocalToGlobal(rVariables, LocalRightHandSideVector);

    VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
    rRightHandSideVector += LocalRightHandSideVector;

    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  void  BeamElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
					    ProcessInfo& rCurrentProcessInfo)
  {

    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(BeamElement::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate elemental system
    this->CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(BeamElement::COMPUTE_LHS_MATRIX);

    VectorType RightHandSideVector = Vector();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector,  LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(RightHandSideVector);

    //Calculate elemental system
    this->CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(BeamElement::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(BeamElement::COMPUTE_LHS_MATRIX);

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate elemental system
    this->CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
						  ElementDataType& rVariables,
						  Vector& rVolumeForce,
						  double& rIntegrationWeight)
  {
    KRATOS_TRY

    SizeType number_of_nodes  = GetGeometry().PointsNumber();
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    double DomainSize = rVariables.Section.Area;

    //gravity load
    Vector GravityForce(3);
    noalias(GravityForce) = ZeroVector(3);

    Vector GravityCouple(3);
    noalias(GravityCouple) = ZeroVector(3);

    Matrix SkewSymMatrix(3,3);
    noalias(SkewSymMatrix) = ZeroMatrix(3,3);

    Vector IntegrationPointPosition(3);
    noalias(IntegrationPointPosition) = ZeroVector(3);

    Vector CurrentValueVector(3);
    noalias(CurrentValueVector) = ZeroVector(3);

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	CurrentValueVector = GetGeometry()[i].Coordinates();
	CurrentValueVector = this->MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );
	IntegrationPointPosition += rVariables.N[i] * CurrentValueVector;
      }

    unsigned int RowIndex = 0;
    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
      	RowIndex = i * ( (dimension-1) * 3 );

        GravityForce  = rIntegrationWeight * rVariables.N[i] * rVolumeForce * DomainSize;

	//integration for the moment force vector
	BeamMathUtilsType::VectorToSkewSymmetricTensor(GravityForce,SkewSymMatrix); // m = f x r = skewF · r
	CurrentValueVector = GetGeometry()[i].Coordinates();
	CurrentValueVector = this->MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );
	CurrentValueVector -= IntegrationPointPosition;
        GravityCouple = prod(SkewSymMatrix,CurrentValueVector);

	if( dimension == 2 ){
	  GravityForce[2] = GravityCouple[2];
	  BeamMathUtilsType::AddVector(GravityForce,  rRightHandSideVector, RowIndex);
	}
	else{
	  BeamMathUtilsType::AddVector(GravityForce,  rRightHandSideVector, RowIndex);
	  BeamMathUtilsType::AddVector(GravityCouple, rRightHandSideVector, RowIndex+dimension);
	}
      }


    KRATOS_CATCH( "" )
  }



  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
						  ElementDataType & rVariables,
						  double& rIntegrationWeight)
  {
    KRATOS_TRY

    KRATOS_ERROR << " calling the default method Fint for a beam element " << std::endl;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void BeamElement::AddExplicitContribution(const VectorType& rRHSVector,
					    const Variable<VectorType>& rRHSVariable,
					    Variable<array_1d<double,3> >& rDestinationVariable,
					    const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    if( (rRHSVariable == RESIDUAL_VECTOR) ){

      if ( rDestinationVariable == FORCE_RESIDUAL )
	{

	  for(SizeType i=0; i< number_of_nodes; i++)
	    {
	      int index = ((dimension-1) * 3) * i;

	      GetGeometry()[i].SetLock();

	      array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

	      for(SizeType j=0; j<dimension; j++)
		{
		  ForceResidual[j] += rRHSVector[index + j];
		}

	      GetGeometry()[i].UnSetLock();
	    }
	}
      else if( rDestinationVariable == MOMENT_RESIDUAL )
	{

	  for(SizeType i=0; i< number_of_nodes; i++)
	    {
	      int index = dimension + ( (dimension-1) * 3 ) * i;

	      GetGeometry()[i].SetLock();

	      array_1d<double, 3 > &MomentResidual = GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);

	      for(SizeType j=0; j<dimension; j++)
		{
		  MomentResidual[j] += rRHSVector[index + j];
		}

	      GetGeometry()[i].UnSetLock();
	    }

	}


    }


    KRATOS_CATCH( "" )
  }



  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
					ElementDataType& rVariables,
					double& rIntegrationWeight)
  {
    KRATOS_TRY

    KRATOS_ERROR << " calling the default method Kuum for a beam element " << std::endl;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
					ElementDataType& rVariables,
					double& rIntegrationWeight)
  {
    KRATOS_TRY


    KRATOS_CATCH( "" )
  }


  //****************GID DEFINITION OF THE AUTOMATIC LOCAL AXES******************
  //*****************************************************************************

  //Local E1 beam axis
  void BeamElement::CalculateLocalAxesMatrix(Matrix& rRotationMatrix)
  {

    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    Vector LocalX(3);
    noalias(LocalX) = ZeroVector(3);
    Vector ReferenceCoordinates(6);
    noalias(ReferenceCoordinates) = ZeroVector(6);

    ReferenceCoordinates[0] = GetGeometry()[0].X();
    ReferenceCoordinates[1] = GetGeometry()[0].Y();
    ReferenceCoordinates[2] = GetGeometry()[0].Z();

    int k = number_of_nodes - 1 ;

    ReferenceCoordinates[3] = GetGeometry()[k].X();
    ReferenceCoordinates[4] = GetGeometry()[k].Y();
    ReferenceCoordinates[5] = GetGeometry()[k].Z();

    for( SizeType i = 0; i < dimension; i++ )
      {
	LocalX[i]  = (ReferenceCoordinates[i+3] - ReferenceCoordinates[i]);
      }


    if(this->Has(LOCAL_AXIS_2)){
      Vector LocalY(3);
      noalias(LocalY) = this->GetValue(LOCAL_AXIS_2);
      BeamMathUtilsType::CalculateLocalAxesMatrix(LocalX,LocalY,rRotationMatrix);
    }
    else{
      BeamMathUtilsType::CalculateLocalAxesMatrix(LocalX,rRotationMatrix);
    }

    KRATOS_CATCH( "" )

  }


  //************************************CALCULATE VOLUME ACCELERATION*******************
  //************************************************************************************

  Vector&  BeamElement::CalculateVolumeForce(Vector& rVolumeForce, const Vector &rN)
  {
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();

    if( rVolumeForce.size() != 3 )
      rVolumeForce.resize(3, false);

    noalias(rVolumeForce) = ZeroVector(3);

    for ( SizeType j = 0; j < number_of_nodes; j++ )
      {
	//temporary, will be checked once at the beginning only
	if( GetGeometry()[j].SolutionStepsDataHas(VOLUME_ACCELERATION) ){
	  rVolumeForce += rN[j] * GetGeometry()[j].FastGetSolutionStepValue(VOLUME_ACCELERATION);
	}
      }

    rVolumeForce *= GetProperties()[DENSITY];

    //Current Frame is the local frame
    rVolumeForce = this->MapToInitialLocalFrame( rVolumeForce );

    return rVolumeForce;

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    bool ComputeDynamicTangent = false;
    if( rCurrentProcessInfo.Has(COMPUTE_DYNAMIC_TANGENT) ){
      if(rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true){
	ComputeDynamicTangent = true;
      }
    }

    if( ComputeDynamicTangent == true ){

      //create local system components
      LocalSystemComponents LocalSystem;

      //calculation flags
      LocalSystem.CalculationFlags.Set(BeamElement::COMPUTE_RHS_VECTOR);
      LocalSystem.CalculationFlags.Set(BeamElement::COMPUTE_LHS_MATRIX);

      //Initialize sizes for the system components:
      this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

      //Set Variables to Local system components
      LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
      LocalSystem.SetRightHandSideVector(rRightHandSideVector);

      //Calculate elemental system
      this->CalculateDynamicSystem( LocalSystem, rCurrentProcessInfo );

    }
    else{

      //1.-Calculate Tangent Inertia Matrix:
      this->CalculateMassMatrix( rLeftHandSideMatrix, rCurrentProcessInfo );

      double MatSize = rLeftHandSideMatrix.size1();

      //2.-Calculate Inertial Forces:
      if ( rRightHandSideVector.size() != MatSize )
	rRightHandSideVector.resize( MatSize, false );

      noalias(rRightHandSideVector) = ZeroVector( MatSize ); //resetting RHS
    }


    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    bool ComputeDynamicTangent = false;
    if( rCurrentProcessInfo.Has(COMPUTE_DYNAMIC_TANGENT) ){
      if(rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true){
	ComputeDynamicTangent = true;
      }
    }

    if( ComputeDynamicTangent == true ){

      //create local system components
      LocalSystemComponents LocalSystem;

      //calculation flags
      LocalSystem.CalculationFlags.Set(BeamElement::COMPUTE_LHS_MATRIX);

      VectorType RightHandSideVector = Vector();

      //Initialize sizes for the system components:
      this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector,  LocalSystem.CalculationFlags );

      //Set Variables to Local system components
      LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
      LocalSystem.SetRightHandSideVector(RightHandSideVector);

      //Calculate elemental system
      this->CalculateDynamicSystem( LocalSystem, rCurrentProcessInfo );

    }
    else{

      //1.-Calculate Tangent Inertia Matrix:
      this->CalculateMassMatrix( rLeftHandSideMatrix, rCurrentProcessInfo );

    }

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    bool ComputeDynamicTangent = false;
    if( rCurrentProcessInfo.Has(COMPUTE_DYNAMIC_TANGENT) )
      if(rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true)
	ComputeDynamicTangent = true;

    if( ComputeDynamicTangent == true ){

      //create local system components
      LocalSystemComponents LocalSystem;

      //calculation flags
      LocalSystem.CalculationFlags.Set(BeamElement::COMPUTE_RHS_VECTOR);

      MatrixType LeftHandSideMatrix = Matrix();

      //Initialize sizes for the system components:
      this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

      //Set Variables to Local system components
      LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
      LocalSystem.SetRightHandSideVector(rRightHandSideVector);

      //Calculate elemental system
      this->CalculateDynamicSystem( LocalSystem, rCurrentProcessInfo );

    }
    else{

      unsigned int MatSize = this->GetDofsSize();

      //2.-Calculate Inertial Forces:
      if ( rRightHandSideVector.size() != MatSize )
	rRightHandSideVector.resize( MatSize, false );

      noalias(rRightHandSideVector) = ZeroVector( MatSize ); //resetting RHS
    }


    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix, ElementDataType& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight )
  {

    KRATOS_TRY


    KRATOS_CATCH( "" )

  }

  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector, ElementDataType& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
  {
    KRATOS_TRY


    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
  {

    KRATOS_TRY


    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void BeamElement::CalculateOnIntegrationPoints(  const Variable<array_1d<double, 3 > >& rVariable,
						   std::vector< array_1d<double, 3 > >& rOutput,
						   const ProcessInfo& rCurrentProcessInfo )
  {

    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
      rOutput.resize( integration_points_number );

    if(rVariable==MOMENT || rVariable==FORCE){

      //create and initialize element variables:
      ElementDataType Variables;
      this->InitializeElementData(Variables,rCurrentProcessInfo);

      //reading integration points (in fact is the two nodes beam element, only one integration point)
      const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

      const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

      for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
	{
	  //compute element kinematics  ...
	  this->CalculateKinematics(Variables,PointNumber);

	  //compute element Strain and Stress Resultants and Couples
	  this->CalculateStressResultants(Variables, PointNumber);

	  // LocalToGlobalSystem for the correct assembly
	  this->MapLocalToGlobal(Variables, Variables.StressVector);

	  if(rVariable==MOMENT)
	    {
	      if( dimension == 2 ){
		rOutput[PointNumber][dimension] = Variables.StressVector[dimension];
	      }
	      else{
		for( SizeType i=0; i<dimension; i++ )
		  {
		    rOutput[PointNumber][i] = Variables.StressVector[i+dimension];
		  }
	      }
	    }

	  if(rVariable==FORCE)
	    {
	      for( SizeType i=0; i<dimension; i++ )
		{
		  rOutput[PointNumber][i] = Variables.StressVector[i];
		}
	    }
	}
    }



    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************


  /**
   * This function provides the place to perform checks on the completeness of the input.
   * It is designed to be called only once (or anyway, not often) typically at the beginning
   * of the calculations, so to verify that nothing is missing from the input
   * or that no common error is found.
   * @param rCurrentProcessInfo
   */
  int  BeamElement::Check(const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = Element::Check(rCurrentProcessInfo);

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(ROTATION);
    KRATOS_CHECK_VARIABLE_KEY(ANGULAR_VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ANGULAR_ACCELERATION);

    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(CROSS_SECTION_AREA);
    KRATOS_CHECK_VARIABLE_KEY(LOCAL_INERTIA_TENSOR);
    //KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION);

    // Check that the element nodes contain all required SolutionStepData and Degrees of freedom
    for(SizeType i=0; i<this->GetGeometry().size(); ++i)
      {
	// Nodal data
	Node<3> &rNode = this->GetGeometry()[i];
	KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rNode);
	KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION,rNode);
	//KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUME_ACCELERATION,rNode);

	// Nodal dofs
	KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X,rNode);
	KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y,rNode);
	if( rCurrentProcessInfo[SPACE_DIMENSION] == 3)
	  KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z,rNode);

	KRATOS_CHECK_DOF_IN_NODE(ROTATION_Z,rNode);
	if( rCurrentProcessInfo[SPACE_DIMENSION] == 3){
	  KRATOS_CHECK_DOF_IN_NODE(ROTATION_X,rNode);
	  KRATOS_CHECK_DOF_IN_NODE(ROTATION_Y,rNode);
	}
      }

    //verify that the area is given by properties
    if ( this->GetProperties().Has(CROSS_SECTION_AREA) == false )
      {
        if( this->GetProperties()[CROSS_SECTION_AREA] == 0.0 )
	  KRATOS_ERROR << "CROSS_SECTION_AREA not provided for this element " << this->Id() << std::endl;
      }

    if ( this->GetProperties().Has(LOCAL_INERTIA_TENSOR) == false )
      {
        if( this->GetProperties()[LOCAL_INERTIA_TENSOR](1,1) == 0.0 )
	  KRATOS_ERROR << "LOCAL_INERTIA_TENSOR not provided for this element " << this->Id() << std::endl;
      }

    return ErrorCode;

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void BeamElement::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    int IntMethod = int(mThisIntegrationMethod);
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("InitialLocalQuaternion",mInitialLocalQuaternion);
  }

  void BeamElement::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("InitialLocalQuaternion",mInitialLocalQuaternion);
  }


} // Namespace Kratos
