//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2015 $
//   Revision:            $Revision:                  0.0 $
//
// 

// System includes

// External includes

// Project includes
#include "custom_elements/beam_elements/large_displacement_beam_semc_element.hpp"

#include "solid_mechanics_application_variables.h"


namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementBeamSEMCElement::LargeDisplacementBeamSEMCElement(IndexType NewId,GeometryType::Pointer pGeometry)
    : LargeDisplacementBeamElement(NewId, pGeometry)
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************


LargeDisplacementBeamSEMCElement::LargeDisplacementBeamSEMCElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : LargeDisplacementBeamElement(NewId, pGeometry, pProperties)
{
    KRATOS_TRY

    mThisIntegrationMethod = GeometryData::GI_GAUSS_1;

    KRATOS_CATCH( "" )
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LargeDisplacementBeamSEMCElement::LargeDisplacementBeamSEMCElement(LargeDisplacementBeamSEMCElement const& rOther)
    :LargeDisplacementBeamElement(rOther)
    ,mPreviousCurvatureVectors(rOther.mPreviousCurvatureVectors)
{

}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer LargeDisplacementBeamSEMCElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new LargeDisplacementBeamSEMCElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementBeamSEMCElement::~LargeDisplacementBeamSEMCElement()
{
}

//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamSEMCElement::Initialize()
{
    KRATOS_TRY
      
    LargeDisplacementBeamElement::Initialize();

    //------------- REDUCED QUADRATURE INTEGRATION 
      
    IntegrationMethod ReducedIntegrationMethod = GeometryData::GI_GAUSS_1;

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( ReducedIntegrationMethod );

    //Curvature Initialization
    if ( mPreviousCurvatureVectors.size() != integration_points.size() )
      {
        mPreviousCurvatureVectors.resize( integration_points.size() );
      }
    
    for ( unsigned int i = 0; i < mPreviousCurvatureVectors.size(); i++ )
      {
	mPreviousCurvatureVectors[i] = ZeroVector(3);
      }
  

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamSEMCElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    LargeDisplacementBeamElement::InitializeSolutionStep(rCurrentProcessInfo);

    for ( unsigned int i = 0; i < mCurrentCurvatureVectors.size(); i++ )
      {
	mPreviousCurvatureVectors[i] = mCurrentCurvatureVectors[i] ;
      }

    mEnergy.clear();

    KRATOS_CATCH( "" )
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


void LargeDisplacementBeamSEMCElement::CalculateStressResultants(ElementVariables& rVariables, const unsigned int& rPointNumber, double alpha)
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //set current CURVATURES
    rVariables.CurrentCurvatureVector  = mPreviousCurvatureVectors[rPointNumber];
    rVariables.PreviousCurvatureVector = mPreviousCurvatureVectors[rPointNumber];


    //compute Strain Resultants and Couples
    Vector StrainResultants = ZeroVector(3); 
    Vector StrainCouples = ZeroVector(3); 

    CalculateStrainResultants(StrainResultants, rVariables, alpha);
    CalculateStrainCouples(StrainCouples, rVariables, alpha);
    
    for ( unsigned int i = 0; i < dimension; i++ )
      {
	rVariables.StrainVector[i]   = StrainResultants[i];
	rVariables.StrainVector[i+3] = StrainCouples[i];
      }
    
    //----------------

    Matrix ConstitutiveMatrix  = ZeroMatrix(6,6);
    this->GetMaterialConstitutiveMatrix(ConstitutiveMatrix, rVariables);

    //Reference Stress Vector
    rVariables.StressVector = prod( ConstitutiveMatrix, rVariables.StrainVector );


    Vector StressResultants = ZeroVector(dimension);
    Vector StressCouples    = ZeroVector(dimension);

    for ( unsigned int i = 0; i < dimension; i++ )
      {
    	StressResultants[i] = rVariables.StressVector[i];
    	StressCouples[i]    = rVariables.StressVector[i+3];
      }
	
    //----------------

    //Current frame given by the Frame Rotation
    StressResultants = prod(rVariables.CurrentRotationMatrix, StressResultants);
    StressCouples    = prod(rVariables.CurrentRotationMatrix, StressCouples);

    for ( unsigned int i = 0; i < dimension; i++ )
      {
    	rVariables.StressVector[i]   = StressResultants[i];
    	rVariables.StressVector[i+3] = StressCouples[i];
      }


    //set variables for the initialization update
    this->UpdateStrainVariables(rVariables,rPointNumber);


    // std::cout<<" Element ["<<this->Id()<<"]"<<std::endl;
    // std::cout<<" StressResultants "<<StressResultants<<std::endl;
    // std::cout<<" StressCouples "<<StressCouples<<std::endl;

    // Vector CurrentStepRotation = ZeroVector(3);
    // this->GetLocalCurrentValue(STEP_ROTATION, CurrentStepRotation, rVariables.N);
    // std::cout<<" StepRotation "<<CurrentStepRotation<<std::endl;

    KRATOS_CATCH( "" )
}


//*********************************SET STRAIN VARIABLES*******************************
//************************************************************************************

void LargeDisplacementBeamSEMCElement::UpdateStrainVariables(ElementVariables& rVariables, const unsigned int& rPointNumber)
{
    KRATOS_TRY

    LargeDisplacementBeamElement::UpdateStrainVariables(rVariables, rPointNumber);

    mCurrentCurvatureVectors[rPointNumber] = rVariables.CurrentCurvatureVector;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamSEMCElement::CalculateStrainResultants(Vector& rStrainResultants, ElementVariables& rVariables, double alpha)
{
    KRATOS_TRY

    // current Strain Resultants
    if( rStrainResultants.size() != 3 )
      rStrainResultants.resize(3, false);

    rStrainResultants = ZeroVector(3); 

    Vector E3 = ZeroVector(3);
    E3[2] = 1.0;
    
    rStrainResultants = prod( trans(rVariables.CurrentRotationMatrix), rVariables.CurrentAxisPositionDerivatives ) - E3;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamSEMCElement::CalculateStrainCouples(Vector& rStrainCouples, ElementVariables& rVariables, double alpha)
{
    KRATOS_TRY

    // current Strain Couples
    if( rStrainCouples.size() != 3 )
      rStrainCouples.resize(3, false);
    
    rStrainCouples = ZeroVector(3); 

    this->CalculateCurrentCurvature(rVariables,STEP_ROTATION);

    rStrainCouples = rVariables.CurrentCurvatureVector;

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamSEMCElement::CalculateCurrentCurvature(ElementVariables& rVariables, const Variable<array_1d<double, 3 > >& rVariable)
{
    KRATOS_TRY

    //------------------[option 0]:start
    // LargeDisplacementBeamElement::CalculateCurrentCurvature(rVariables,rVariable);
    //------------------[option 0]:end
      
    //------------------[option 1]:start  
    const unsigned int number_of_nodes = GetGeometry().size();

    Vector CurrentStepRotationDerivativesVector = ZeroVector(3);
    Vector CurrentValueVector                   = ZeroVector(3);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
    	//Current Rotation Derivatives	
    	CurrentValueVector = GetNodalCurrentValue( rVariable, CurrentValueVector, i );
    	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector, rVariables.PointNumber );

    	CurrentStepRotationDerivativesVector += rVariables.DN_DX(i,0) * ( CurrentValueVector );
      }

    rVariables.CurrentCurvatureVector  = rVariables.PreviousCurvatureVector;
    rVariables.CurrentCurvatureVector += prod( rVariables.CurrentRotationMatrix, CurrentStepRotationDerivativesVector );
    //------------------[option 1]:end


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

//Strain Energy Calculation 
void LargeDisplacementBeamSEMCElement::CalculateStrainEnergy(ElementVariables& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
{
    KRATOS_TRY

    //Internal Energy Calculation: alpha = 1     
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //compute Strain Resultants and Couples
    Vector StrainResultants = ZeroVector(3); 
    Vector StrainCouples    = ZeroVector(3); 

    double Alpha = 1.0;
    CalculateStrainResultants(StrainResultants, rVariables, Alpha);
    CalculateStrainCouples(StrainCouples, rVariables, Alpha);

    //----------------

    Matrix ConstitutiveMatrix  = ZeroMatrix(6,6);
    this->GetMaterialConstitutiveMatrix(ConstitutiveMatrix, rVariables);

    //----------------

    Vector StrainVector = ZeroVector(6);
    for ( unsigned int i = 0; i < dimension; i++ )
      {
    	StrainVector[i]   = StrainResultants[i];
    	StrainVector[i+3] = StrainCouples[i];
      }

    //Reference Stress Vector
    Vector StressVector = prod( ConstitutiveMatrix, StrainVector );

    //std::cout<<" Strain "<<StrainVector<<" Stress "<<StressVector<<std::endl;

    mEnergy.Deformation += 0.5 * (inner_prod(StressVector, StrainVector)) * rIntegrationWeight ;
    
    //std::cout<<" mEnergy.Deformation "<<mEnergy.Deformation<<" rIntegrationWeight "<<rIntegrationWeight<<std::endl;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamSEMCElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
      
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize               = number_of_nodes * ( dimension * 2 );

    if(rMassMatrix.size1() != MatSize)
      rMassMatrix.resize (MatSize, MatSize, false);
    
    rMassMatrix = ZeroMatrix( MatSize, MatSize );
      
    bool exact_integration = false;
    if( exact_integration ){

      //create and initialize element variables:
      ElementVariables Variables;

      IntegrationMethod ThisIntegrationMethod = mThisIntegrationMethod;
      //full quadrature integration:
      mThisIntegrationMethod = GeometryData::GI_GAUSS_2;

      const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

      // initialize variables short version:

      //Compute Section Properties:
      this->CalculateSectionProperties(Variables.Section);
      Variables.Length = GetGeometry().Length();

      //Set equilibrium point initial:0/mid:0.5/final:1
      if( rCurrentProcessInfo.Has(EQUILIBRIUM_POINT) )
	Variables.Alpha = rCurrentProcessInfo[EQUILIBRIUM_POINT];
      else
	Variables.Alpha = 1;

      //reading shape functions
      Variables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

      //get the shape functions for the order of the integration method [N]
      const Matrix& Ncontainer = Variables.GetShapeFunctions();

      // initialize variables short version;
 
      MatrixType  LocalLeftHandSideMatrix = ZeroMatrix( MatSize, MatSize );

      //(in fact is the two nodes beam element, full quadrature 2 integration points)
      for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
	{
	  //integration point number
	  Variables.PointNumber = PointNumber;

	  //set shape functions values for this integration point
	  Variables.N=row( Ncontainer, PointNumber);

	  //compute local to global frame
	  this->CalculateFrameMapping( Variables, PointNumber );
	
	  Variables.detJ = mDetJ0Full[PointNumber];
 
	  double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
	  IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );
	
	  LocalLeftHandSideMatrix.clear();
	    
	  LargeDisplacementBeamElement::CalculateAndAddInertiaLHS( LocalLeftHandSideMatrix, Variables, rCurrentProcessInfo, IntegrationWeight ); // (R_N+1, R_N)

	  BeamMathUtilsType::MapLocalToGlobal3D(Variables.LocalTransformationMatrix,LocalLeftHandSideMatrix);
	    
	  rMassMatrix += LocalLeftHandSideMatrix;
		
	}

      mThisIntegrationMethod = ThisIntegrationMethod;

    }
    else{


      SectionProperties Section;
      this->CalculateSectionProperties(Section);

      //block m(1,1) of the mass matrix

      MatrixType m11 = ZeroMatrix(3,3);

      double TotalMass = 0;
      TotalMass  = this->CalculateTotalMass( Section, TotalMass );
      TotalMass *= GetGeometry().Length();

      Vector LumpFact = ZeroVector(number_of_nodes); 

      LumpFact = GetGeometry().LumpingFactors(LumpFact);

      //block m(2,2) of the mass matrix
    
      MatrixType m22 = ZeroMatrix(3,3);

      Matrix InertiaDyadic = ZeroMatrix(3,3);
      this->CalculateInertiaDyadic( Section, InertiaDyadic );


      for( unsigned int i=0; i < number_of_nodes; i++ )
	{

	  double temp = LumpFact[i] * TotalMass;

	  int RowIndex = i * (dimension * 2);

	  for( unsigned int k=0; k < dimension; k++ )
	    {
	      m11(k,k) = temp;
	    }

	  m22 = InertiaDyadic * temp;
	
	  //Building the Local Tangent Inertia Matrix
	  MathUtils<double>::AddMatrix( rMassMatrix, m11, RowIndex, RowIndex );
	  MathUtils<double>::AddMatrix( rMassMatrix, m22, RowIndex+3, RowIndex+3 );
	
	}
 
      // Note: 
      // Variables.LocalTransformationMatrix is a rotation matrix with new base in columns
      // That means that the standard rotation K = Q·K'·QT and F = Q·F' is the correct transformation

      // initialize local transformation/rotation matrix
      Matrix LocalTransformationMatrix = ZeroMatrix(3,3);

      mInitialLocalQuaternionsReduced[0].ToRotationMatrix(LocalTransformationMatrix);

      //std::cout<<"(Element:"<<this->Id()<<") Mass "<<rMassMatrix<<std::endl;

      // Transform Local to Global LHSMatrix:
      BeamMathUtilsType::MapLocalToGlobal3D(LocalTransformationMatrix,rMassMatrix);
    }


    std::cout<<" DYNAMIC_TANGENT   ("<<this->Id()<<"): "<<rMassMatrix<<std::endl;

    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

//Inertia in the SPATIAL configuration 
void LargeDisplacementBeamSEMCElement::CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight )
{

    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize               = number_of_nodes * ( dimension * 2 );

    if(rLeftHandSideMatrix.size1() != MatSize)
        rLeftHandSideMatrix.resize (MatSize, MatSize, false);

    rLeftHandSideMatrix = ZeroMatrix( MatSize, MatSize );

    SectionProperties Section;
    this->CalculateSectionProperties(Section);


    //rCurrentProcessInfo must give it:
    //double DeltaTime = rCurrentProcessInfo[DELTA_TIME];
 
    //double Newmark = (1.0/ ( DeltaTime * DeltaTime * 0.5));

    //block m(1,1) of the mass matrix

    MatrixType m11 = ZeroMatrix(3,3);

    double TotalMass = 0;
    TotalMass  = this->CalculateTotalMass( Section, TotalMass );
    
    //block m(2,2) of the mass matrix
     
    MatrixType m22 = ZeroMatrix(3,3);

    Vector CurrentStepRotationVector = ZeroVector(3);
    Vector CurrentValueVector        = ZeroVector(3);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	//Current Compound Rotation Vector 
	CurrentValueVector = GetNodalCurrentValue( STEP_ROTATION, CurrentValueVector, i );
	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector );
	CurrentStepRotationVector += rVariables.N[i] * CurrentValueVector;
      }

    Matrix DiagonalMatrix = IdentityMatrix(3);   

    Matrix Tangent   = IdentityMatrix(3);
    double tolerance = 1e-12; 
    double rotation  = norm_2(CurrentStepRotationVector);

    if( rotation >= tolerance ){

      double sinus   = std::sin(rotation);
      double cosinus = std::cos(rotation);
      
      double c1 = ( sinus/rotation );
      double c2 = ( 1.0 - cosinus )/( rotation * rotation );
      double c3 = ( cosinus * rotation - sinus )/( rotation * rotation * rotation );
      double c4 = ( sinus * rotation - 2.0 * ( 1.0 - cosinus ) )/( rotation * rotation * rotation * rotation );
    
      // (1) get inertia dyadic
      Matrix InertiaDyadic = ZeroMatrix(3,3);
      this->CalculateInertiaDyadic( Section, InertiaDyadic );
      // InertiaDyadic = prod(rVariables.CurrentRotationMatrix,InertiaDyadic);
      // InertiaDyadic = prod(InertiaDyadic,trans(rVariables.CurrentRotationMatrix));
      InertiaDyadic = prod(rVariables.PreviousRotationMatrix,InertiaDyadic);
      InertiaDyadic = prod(InertiaDyadic,trans(rVariables.PreviousRotationMatrix));


      // (2) rotation and inertia terms:
      Matrix RotationTensor = ZeroMatrix(3,3);
      BeamMathUtilsType::VectorToSkewSymmetricTensor( CurrentStepRotationVector, RotationTensor );

      Vector InertiaxRotationVector = prod(InertiaDyadic, CurrentStepRotationVector);

      Matrix InertiaxRotationTensor = ZeroMatrix(3,3);
      BeamMathUtilsType::VectorToSkewSymmetricTensor( InertiaxRotationVector, InertiaxRotationTensor );

      Vector RotationxInertiaRotationVector = MathUtils<double>::CrossProduct(CurrentStepRotationVector, InertiaxRotationVector);
      //Vector RotationxInertiaRotationVector = prod(RotationTensor,InertiaxRotationVector);

      // (3) tangent:
      Tangent  = c1 * InertiaDyadic;
      Tangent += c2 * ( (-1) * InertiaxRotationTensor+ prod( RotationTensor, InertiaDyadic ) );
      Tangent += c3 * outer_prod( InertiaxRotationVector, CurrentStepRotationVector );
      Tangent += c4 * outer_prod( RotationxInertiaRotationVector, CurrentStepRotationVector );

    }
    
    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	m11 = ZeroMatrix(3,3);
	m22 = ZeroMatrix(3,3);

	RowIndex = i * (dimension * 2);


	for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {

	    ColIndex = j * (dimension * 2);

	    m11 = TotalMass * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight * DiagonalMatrix;

	    m22 = rVariables.N[i] * rVariables.N[j] * rIntegrationWeight * Tangent;
	    

	    //Building the Local Tangent Inertia Matrix
	    MathUtils<double>::AddMatrix( rLeftHandSideMatrix, m11, RowIndex, ColIndex );
	    MathUtils<double>::AddMatrix( rLeftHandSideMatrix, m22, RowIndex+3, ColIndex+3 );
	    
	  }
      }
         

    //std::cout<<"(Element:"<<this->Id()<<") Tangent "<<Tangent<<" LHS "<< rLeftHandSideMatrix<<std::endl;

    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************

//Inertia in the SPATIAL configuration 
void LargeDisplacementBeamSEMCElement::CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize               = number_of_nodes * ( dimension * 2 );

    if(rRightHandSideVector.size() != MatSize)
        rRightHandSideVector.resize(MatSize, false);

    rRightHandSideVector = ZeroVector( MatSize );

    SectionProperties Section;
    this->CalculateSectionProperties(Section);
  
    double TotalMass = 0;
    TotalMass = this->CalculateTotalMass( Section, TotalMass );     
 

    Vector CurrentLinearAccelerationVector = ZeroVector(3);
    Vector CurrentStepRotationVector       = ZeroVector(3);
    Vector CurrentValueVector              = ZeroVector(3);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	//Current Compound Rotation Vector 
	CurrentValueVector = GetNodalCurrentValue( STEP_ROTATION, CurrentValueVector, i );
	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector );
	CurrentStepRotationVector += rVariables.N[i] * CurrentValueVector;
	
	//Current Linear Acceleration Vector
	CurrentValueVector = GetNodalCurrentValue( ACCELERATION, CurrentValueVector, i );
	CurrentValueVector = MapToInitialLocalFrame( CurrentValueVector );
	CurrentLinearAccelerationVector += rVariables.N[i] * CurrentValueVector;
      }


    //-----------------
    //block m(1) of the inertial force vector

    //Compute Linear Term
    Vector LinearInertialForceVector = ZeroVector(3);
    LinearInertialForceVector = TotalMass * CurrentLinearAccelerationVector;


    //-----------------
    //block m(2,2) of the inertial force vector (rotations part::to be defined)
    

    Vector Residual  = ZeroVector(3);
    double tolerance = 1e-12; 
    double rotation  = norm_2(CurrentStepRotationVector);

    if( rotation >= tolerance ){

      double sinus   = std::sin(rotation);
      double cosinus = std::cos(rotation);
      
      double c1 = ( sinus/rotation );
      double c2 = ( 1.0 - cosinus )/( rotation * rotation );
    
      // (1) get inertia dyadic
      Matrix InertiaDyadic = ZeroMatrix(3,3);
      this->CalculateInertiaDyadic( Section, InertiaDyadic );
      // InertiaDyadic = prod(rVariables.CurrentRotationMatrix,InertiaDyadic);
      // InertiaDyadic = prod(InertiaDyadic,trans(rVariables.CurrentRotationMatrix));
      InertiaDyadic = prod(rVariables.PreviousRotationMatrix,InertiaDyadic);
      InertiaDyadic = prod(InertiaDyadic,trans(rVariables.PreviousRotationMatrix));

      // (2) rotation and inertia terms:
      Vector InertiaxRotationVector = prod(InertiaDyadic, CurrentStepRotationVector);
      Vector RotationxInertiaRotationVector = MathUtils<double>::CrossProduct(CurrentStepRotationVector, InertiaxRotationVector);
 
      // (3) residual:
      Residual  = c1 * InertiaxRotationVector + c2 * RotationxInertiaRotationVector;

    }

    //compose total acceleration integral function:
    Vector TotalInertialForceVector = ZeroVector(6);
    
    BeamMathUtilsType::AddVector(LinearInertialForceVector,  TotalInertialForceVector, 0);
    BeamMathUtilsType::AddVector(Residual, TotalInertialForceVector, 3);

    //Initialize Local Matrices
    VectorType Fi         = ZeroVector(6);
    MatrixType OperatorI  = ZeroMatrix(6,6);
    unsigned int RowIndex = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	RowIndex = i * (dimension * 2);

	Fi = ZeroVector(6);
	
 	this->CalculateOperator( OperatorI, rVariables.N, i );

    	//nodal force vector
    	Fi  = prod( OperatorI, TotalInertialForceVector );

	Fi *= rIntegrationWeight;
	
	BeamMathUtilsType::AddVector(Fi, rRightHandSideVector, RowIndex);
      }

    //std::cout<<"(Element:"<<this->Id()<<") Residual "<<Residual<<" RHS "<< rRightHandSideVector<<" intweight "<<rIntegrationWeight<<std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamSEMCElement::AddExplicitContribution(const MatrixType& rLHSMatrix, 
							       const Variable<MatrixType>& rLHSVariable, 
							       Variable<Matrix >& rDestinationVariable, 
							       const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
      
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
 

    if( (rLHSVariable == TANGENT_MATRIX) ){

      if ( rDestinationVariable == TANGENT_LYAPUNOV )
	{

	  for(unsigned int i=0; i< number_of_nodes; i++)
	    {
	      int index = dimension + (dimension * 2) * i;

	      GetGeometry()[i].SetLock();

	      Matrix& Tangent= GetGeometry()[i].FastGetSolutionStepValue(TANGENT_LYAPUNOV);

	      for(unsigned int j=0; j<dimension; j++)
		{
		  for(unsigned int k=0; k<dimension; k++)
		    {
		      Tangent(j,k) += rLHSMatrix(index+j,index+k);
		    }
		}

	      GetGeometry()[i].UnSetLock();
	    }
	}
      
    }


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementBeamSEMCElement::AddExplicitContribution(const VectorType& rRHSVector, 
							      const Variable<VectorType>& rRHSVariable, 
							      Variable<array_1d<double,3> >& rDestinationVariable, 
							      const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY
      
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( (rRHSVariable == RESIDUAL_VECTOR) ){

      if ( rDestinationVariable == RESIDUAL_LYAPUNOV )
	{

	  for(unsigned int i=0; i< number_of_nodes; i++)
	    {
	      int index = dimension + (dimension * 2) * i;

	      GetGeometry()[i].SetLock();

	      array_1d<double, 3 > &Residual = GetGeometry()[i].FastGetSolutionStepValue(RESIDUAL_LYAPUNOV);

	      for(unsigned int j=0; j<dimension; j++)
		{
		  Residual[j] += rRHSVector[index + j];
		}

	      GetGeometry()[i].UnSetLock();
	    }
	}
      else{

	LargeDisplacementBeamElement::AddExplicitContribution(rRHSVector, rRHSVariable, rDestinationVariable, rCurrentProcessInfo);
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
int LargeDisplacementBeamSEMCElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (GetGeometry().WorkingSpaceDimension() != 3 || GetGeometry().size()!=2 )
    {
      KRATOS_THROW_ERROR( std::invalid_argument, "This element works only in 3D and with 2 noded linear elements", "")
    }

    //verify that the variables are correctly initialized
    if(VELOCITY.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered", "" )
    if(DISPLACEMENT.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )
    if(ACCELERATION.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered", "" )
    if(DENSITY.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"DENSITY has Key zero! (check if the application is correctly registered", "" )
     if(VOLUME_ACCELERATION.Key() == 0)
        // KRATOS_THROW_ERROR( std::invalid_argument,"VOLUME_ACCELERATION has Key zero! (check if the application is correctly registered", "" )
    if(CROSS_SECTION_AREA.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"CROSS_SECTION_AREA has Key zero! (check if the application is correctly registered", "" )
    if(LOCAL_INERTIA_TENSOR.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"LOCAL_INERTIA_TENSOR has Key zero! (check if the application is correctly registered", "" )
    if(LOCAL_CONSTITUTIVE_MATRIX.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"LOCAL_CONSTITUTIVE_MATRIX has Key zero! (check if the application is correctly registered", "" )
    if(ROTATION.Key() == 0)
        KRATOS_THROW_ERROR( std::invalid_argument,"ROTATION has Key zero! (check if the application is correctly registered", "" )

    //verify that the dofs exist
    for(unsigned int i=0; i<this->GetGeometry().size(); i++)
    {
        if(this->GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false)
            KRATOS_THROW_ERROR( std::invalid_argument,"missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id() )
        if(this->GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false)
            KRATOS_THROW_ERROR( std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id() )
    }

    //verify that the area is given by properties
    if (this->GetProperties().Has(CROSS_SECTION_AREA)==false)
    {
        if( GetValue(CROSS_SECTION_AREA) == 0.0 )
            KRATOS_THROW_ERROR( std::logic_error,"CROSS_SECTION_AREA not provided for this element", this->Id() )
    }

    //verify that the inertia is given by properties
    if (this->GetProperties().Has(LOCAL_INERTIA_TENSOR)==false)
    {
        if( GetValue(LOCAL_INERTIA_TENSOR)(0,0) == 0.0 )
	  KRATOS_THROW_ERROR( std::logic_error,"LOCAL_INERTIA_TENSOR not provided for this element ", this->Id() )
    }


    return 0;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************


void LargeDisplacementBeamSEMCElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementBeamElement )
    rSerializer.save("PreviousCurvatureVecors",mPreviousCurvatureVectors);

}

void LargeDisplacementBeamSEMCElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementBeamElement )
    rSerializer.load("PreviousCurvatureVecors",mPreviousCurvatureVectors);

}


} // Namespace Kratos


