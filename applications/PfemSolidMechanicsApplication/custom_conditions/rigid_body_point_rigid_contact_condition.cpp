//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:           September 2014 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes


// External includes


// Project includes
#include "custom_conditions/rigid_body_point_rigid_contact_condition.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{


/**
 * Flags related to the condition computation
 */
KRATOS_CREATE_LOCAL_FLAG( RigidBodyPointRigidContactCondition, COMPUTE_RHS_VECTOR,                 0 );
KRATOS_CREATE_LOCAL_FLAG( RigidBodyPointRigidContactCondition, COMPUTE_LHS_MATRIX,                 1 );
KRATOS_CREATE_LOCAL_FLAG( RigidBodyPointRigidContactCondition, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
KRATOS_CREATE_LOCAL_FLAG( RigidBodyPointRigidContactCondition, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );


//***********************************************************************************
//***********************************************************************************
RigidBodyPointRigidContactCondition::RigidBodyPointRigidContactCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//***********************************************************************************
//***********************************************************************************
RigidBodyPointRigidContactCondition::RigidBodyPointRigidContactCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{

    //DO NOT ADD DOFS HERE!!!
}


//************************************************************************************
//************************************************************************************
RigidBodyPointRigidContactCondition::RigidBodyPointRigidContactCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
  : Condition(NewId, pGeometry, pProperties)
{
    mpRigidWall = pRigidWall;
    
    const unsigned int inode = GetGeometry().PointsNumber()-1;
    
    mMasterElements = GetGeometry()[inode].GetValue(MASTER_ELEMENTS);
    
    //DO NOT ADD DOFS HERE!!!
}


//************************************************************************************
//************************************************************************************
RigidBodyPointRigidContactCondition::RigidBodyPointRigidContactCondition( RigidBodyPointRigidContactCondition const& rOther )
    : Condition(rOther)
    , mMasterElements(rOther.mMasterElements)
      
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer RigidBodyPointRigidContactCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new RigidBodyPointRigidContactCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}


//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer RigidBodyPointRigidContactCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  return this->Create( NewId, rThisNodes, pGetProperties() );
}


//***********************************************************************************
//***********************************************************************************
RigidBodyPointRigidContactCondition::~RigidBodyPointRigidContactCondition()
{
}

//************* GETTING METHODS

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::GetDofList(DofsVectorType& rConditionDofList,
				    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rConditionDofList.resize(0);

    Element& MasterElement = mMasterElements.back();
    MasterElement.GetDofList(rConditionDofList, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::EquationIdVector(EquationIdVectorType& rResult,
					  ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rResult.resize( 0, false );

    Element& MasterElement = mMasterElements.back();
    MasterElement.EquationIdVector(rResult, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::GetValuesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    Element& MasterElement = mMasterElements.back();
    MasterElement.GetValuesVector(rValues, Step);


    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    KRATOS_TRY

    Element& MasterElement = mMasterElements.back();
    MasterElement.GetFirstDerivativesVector(rValues, Step);

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    KRATOS_TRY

    Element& MasterElement = mMasterElements.back();
    MasterElement.GetSecondDerivativesVector(rValues, Step);

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
void RigidBodyPointRigidContactCondition::ClearNodalForces()
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	GetGeometry()[i].SetLock();
	array_1d<double, 3> & ContactForce  = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
	ContactForce.clear();
	GetGeometry()[i].UnSetLock();
      }

    //KRATOS_WATCH( " CLEAR NODAL FORCE " )

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::AddExplicitContribution(const VectorType& rRHSVector, 
						 const Variable<VectorType>& rRHSVariable, 
						 Variable<array_1d<double,3> >& rDestinationVariable, 
						 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( rRHSVariable == CONTACT_FORCES_VECTOR && rDestinationVariable == CONTACT_FORCE )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = i * (dimension * (dimension-1));

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		ContactForce[j] += rRHSVector[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }


    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = i * (dimension * (dimension-1));

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

//************* STARTING - ENDING  METHODS
//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    GeneralVariables ContactVariables;
    int ContactFace = 0;
       
    if ( this->mpRigidWall->IsInside( GetGeometry()[0], ContactVariables.Gap.Normal, ContactVariables.Gap.Tangent, ContactVariables.Surface.Normal, ContactVariables.Surface.Tangent, ContactFace ) ) {

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      array_1d<double, 3> &ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);

      mTangentialVariables.PreviousTangentForceModulus = 0.0;
      for (unsigned int i = 0; i < dimension; ++i) {
	mTangentialVariables.PreviousTangentForceModulus += ContactForce[i] * ContactVariables.Surface.Tangent[i];
      }


    }
    else {
      mTangentialVariables.PreviousTangentForceModulus = 0.0;
    } 
    
    mTangentialVariables.DeltaTime = rCurrentProcessInfo[DELTA_TIME];

    mTangentialVariables.Sign = 1;
    
    mTangentialVariables.FrictionCoefficient = 0.3;
    mTangentialVariables.DynamicFrictionCoefficient = 0.3;
    mTangentialVariables.StaticFrictionCoefficient  = 0.4;
    
    ClearNodalForces();
    
    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
void RigidBodyPointRigidContactCondition::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
  CurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] = 0;
  CurrentProcessInfo[NUMBER_OF_STICK_CONTACTS]  = 0;
  CurrentProcessInfo[NUMBER_OF_SLIP_CONTACTS]   = 0;
  
  ClearNodalForces();
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointRigidContactCondition::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointRigidContactCondition::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
  KRATOS_TRY
    
  CurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] = 0;
  CurrentProcessInfo[NUMBER_OF_STICK_CONTACTS]  = 0;
  CurrentProcessInfo[NUMBER_OF_SLIP_CONTACTS]   = 0;
  
  KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
								   VectorType& rRightHandSideVector,
								   Flags& rCalculationFlags)

{

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * (dimension * 2);

    if ( rCalculationFlags.Is(RigidBodyPointRigidContactCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }


    //resizing as needed the RHS
    if ( rCalculationFlags.Is(RigidBodyPointRigidContactCondition::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
	    rRightHandSideVector.resize( MatSize, false );
      
	rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
	  
    }
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::InitializeGeneralVariables(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

 
    KRATOS_CATCH( "" )

}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void RigidBodyPointRigidContactCondition::CalculateKinematics(GeneralVariables& rVariables,
					     const double& rPointNumber)
{
    KRATOS_TRY

    int ContactFace = 0; //free surface
      
    if( this->mpRigidWall->IsInside( GetGeometry()[0], rVariables.Gap.Normal, rVariables.Gap.Tangent, rVariables.Surface.Normal, rVariables.Surface.Tangent, ContactFace ) ){

      rVariables.Options.Set(ACTIVE,true);

      //get contact properties and parameters
      this->CalculateContactFactors( rVariables );

      //get correct tangent vector
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      const array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
      const array_1d<double, 3> & PreviousDisplacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, 1);
      array_1d<double, 3 > DeltaDisplacement           = CurrentDisplacement-PreviousDisplacement;
 
      for (unsigned int i = 0; i < dimension ; ++i) {
	rVariables.Surface.Tangent[i] = DeltaDisplacement[i] - mTangentialVariables.DeltaTime * this->mpRigidWall->Velocity()[i];
      }

      if( norm_2(rVariables.Surface.Normal) )
	rVariables.Surface.Normal /=  norm_2(rVariables.Surface.Normal);

      rVariables.Surface.Tangent -= inner_prod(rVariables.Surface.Tangent, rVariables.Surface.Normal) * rVariables.Surface.Normal;

      if( norm_2(rVariables.Surface.Tangent) )
	rVariables.Surface.Tangent /=  norm_2(rVariables.Surface.Tangent);

      rVariables.Surface.Tangent *= (-1); // for consistency in the tangent contact force direction

    }
    else{
      
      rVariables.Options.Set(ACTIVE,false);
      
    }

    KRATOS_CATCH( "" )
}

  //************************************************************************************
  //************************************************************************************


  void RigidBodyPointRigidContactCondition::CalculateContactFactors(GeneralVariables &rVariables)
  {

    KRATOS_TRY
      
    WeakPointerVector<Node<3> >& rN = GetGeometry()[0].GetValue(NEIGHBOUR_NODES);

    array_1d<double,3> Contact_Point = GetGeometry()[0].Coordinates();
    array_1d<double,3> Neighb_Point;

    double distance = 0;
    double counter = 0;

    for(unsigned int i = 0; i < rN.size(); i++)
      {
	if(rN[i].Is(BOUNDARY)){
	    
	  Neighb_Point[0] = rN[i].X();
	  Neighb_Point[1] = rN[i].Y();
	  Neighb_Point[2] = rN[i].Z();
	    
	  distance += norm_2(Contact_Point-Neighb_Point);

	  counter ++;
	}
      }

    if( counter != 0 )
      distance /= counter;

    if( distance == 0 )
      distance = 1;
    
    //get contact properties and parameters
    double PenaltyParameter = GetProperties()[PENALTY_PARAMETER];
    double ElasticModulus   = mMasterElements.front().GetProperties()[YOUNG_MODULUS];

    //reduction of the penalty parameter:
    PenaltyParameter *=1e-6;
    
    double factor = 4;
    if( distance < 1.0 ){ //take a number bigger than 1.0 (length units)
      int order = (int)((-1) * std::log10(distance) + 1) ;
      distance *= factor * pow(10,order);
    }

    rVariables.Penalty.Normal  = distance * PenaltyParameter * ElasticModulus;
    rVariables.Penalty.Tangent = rVariables.Penalty.Normal;  
    

    rVariables.CentroidPosition =  GetGeometry()[0].Coordinates() - mMasterElements.front().GetGeometry()[0].Coordinates();
    rVariables.CentroidDistance = norm_2(rVariables.CentroidPosition);


    rVariables.SkewSymDistance = ZeroMatrix(3);

    //compute the skewsymmmetric tensor of the distance
    this->VectorToSkewSymmetricTensor(rVariables.CentroidPosition, rVariables.SkewSymDistance);


    //std::cout<<" Node "<<GetGeometry()[0].Id()<<" Contact Factors "<<rVariables.Penalty.Normal<<" Gap Normal "<<rVariables.Gap.Normal<<" Gap Tangent "<<rVariables.Gap.Tangent<<" Surface.Normal "<<rVariables.Surface.Normal<<" Surface.Tangent "<<rVariables.Surface.Tangent<<" distance "<<distance<<" ElasticModulus "<<ElasticModulus<<" PenaltyParameter "<<PenaltyParameter<<std::endl;
    
    // std::cout<<" Penalty.Normal "<<rVariables.Penalty.Normal<<" Penalty.Tangent "<<rVariables.Penalty.Tangent<<std::endl;

    KRATOS_CATCH( "" )
      }



//************************************************************************************
//************************************************************************************

void RigidBodyPointRigidContactCondition::CalculateConditionSystem(LocalSystemComponents& rLocalSystem,
							  ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //create and initialize condition variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    //reading integration points
    for ( unsigned int PointNumber = 0; PointNumber < 1; PointNumber++ )
    {
        //compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables,PointNumber);

        //calculating weights for integration on the "reference configuration"
        double IntegrationWeight = 1;
        IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );

	if( Variables.Options.Is(ACTIVE) )
	  rCurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] += 1; 

        if ( rLocalSystem.CalculationFlags.Is(RigidBodyPointRigidContactCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
        {
            //contributions to stiffness matrix calculated on the reference config
	    this->CalculateAndAddLHS ( rLocalSystem, Variables, IntegrationWeight );
        }

        if ( rLocalSystem.CalculationFlags.Is(RigidBodyPointRigidContactCondition::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
        {
            //contribution to external forces 
	    this->CalculateAndAddRHS ( rLocalSystem, Variables, IntegrationWeight );
        }

    }

    KRATOS_CATCH( "" )
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void RigidBodyPointRigidContactCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
{

  //contributions of the stiffness matrix calculated on the reference configuration
  if( rLocalSystem.CalculationFlags.Is( RigidBodyPointRigidContactCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
    {
      std::vector<MatrixType>& rLeftHandSideMatrices = rLocalSystem.GetLeftHandSideMatrices();
      const std::vector< Variable< MatrixType > >& rLeftHandSideVariables = rLocalSystem.GetLeftHandSideVariables();

      for( unsigned int i=0; i<rLeftHandSideVariables.size(); i++ )
	{
	  bool calculated = false;
	  
	  if( rLeftHandSideVariables[i] == GEOMETRIC_STIFFNESS_MATRIX ){
	    // operation performed: add Kg to the rLefsHandSideMatrix
	    this->CalculateAndAddKuug( rLeftHandSideMatrices[i], rVariables, rIntegrationWeight );
	    calculated = true;
	  }

	  if(calculated == false)
	    {
	      KRATOS_THROW_ERROR(std::logic_error, " ELEMENT can not supply the required local system variable: ",rLeftHandSideVariables[i])
	    }

	}
    } 
  else{
    
    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix(); 

    // operation performed: add Kg to the rLefsHandSideMatrix
    this->CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    //KRATOS_WATCH( rLeftHandSideMatrix )
  }

}


//************************************************************************************
//************************************************************************************

void RigidBodyPointRigidContactCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
{
    //contribution of the internal and external forces
    if( rLocalSystem.CalculationFlags.Is( RigidBodyPointRigidContactCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
    {

      std::vector<VectorType>& rRightHandSideVectors = rLocalSystem.GetRightHandSideVectors();
      const std::vector< Variable< VectorType > >& rRightHandSideVariables = rLocalSystem.GetRightHandSideVariables();
      for( unsigned int i=0; i<rRightHandSideVariables.size(); i++ )
	{
	  bool calculated = false;
	  if( rRightHandSideVariables[i] == CONTACT_FORCES_VECTOR ){
	    // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
	    this->CalculateAndAddContactForces( rRightHandSideVectors[i], rVariables, rIntegrationWeight );
	    calculated = true;
	  }
	  
	  if(calculated == false)
	    {
	      KRATOS_THROW_ERROR(std::logic_error, " ELEMENT can not supply the required local system variable: ",rRightHandSideVariables[i])
	    }

	}
    }
    else{
      
      VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector(); 

      // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
      this->CalculateAndAddContactForces( rRightHandSideVector, rVariables, rIntegrationWeight );

      //KRATOS_WATCH( rRightHandSideVector )

    }
    

}

//***********************************************************************************
//************************************************************************************

double& RigidBodyPointRigidContactCondition::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    return rIntegrationWeight;
}


//************************************************************************************
//************************************************************************************

void RigidBodyPointRigidContactCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyPointRigidContactCondition::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );

}

//************************************************************************************
//************************************************************************************

void RigidBodyPointRigidContactCondition::CalculateRightHandSide( std::vector< VectorType >& rRightHandSideVectors, const std::vector< Variable< VectorType > >& rRHSVariables, ProcessInfo& rCurrentProcessInfo )
{
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyPointRigidContactCondition::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(RigidBodyPointRigidContactCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    if( rRHSVariables.size() != rRightHandSideVectors.size() )
      rRightHandSideVectors.resize(rRHSVariables.size());
    
    for( unsigned int i=0; i<rRightHandSideVectors.size(); i++ )
      {
	this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVectors[i], LocalSystem.CalculationFlags );
      }

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

    LocalSystem.SetRightHandSideVariables(rRHSVariables);

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );


}



//************************************************************************************
//************************************************************************************

void RigidBodyPointRigidContactCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyPointRigidContactCondition::COMPUTE_LHS_MATRIX);
    LocalSystem.CalculationFlags.Set(RigidBodyPointRigidContactCondition::COMPUTE_RHS_VECTOR);

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );

    //KRATOS_WATCH( rLeftHandSideMatrix )
    //KRATOS_WATCH( rRightHandSideVector )

}


//************************************************************************************
//************************************************************************************

void RigidBodyPointRigidContactCondition::CalculateLocalSystem( std::vector< MatrixType >& rLeftHandSideMatrices,
					       const std::vector< Variable< MatrixType > >& rLHSVariables,
					       std::vector< VectorType >& rRightHandSideVectors,
					       const std::vector< Variable< VectorType > >& rRHSVariables,
					       ProcessInfo& rCurrentProcessInfo )
{
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyPointRigidContactCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS);
    LocalSystem.CalculationFlags.Set(RigidBodyPointRigidContactCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);


    //Initialize sizes for the system components:
    if( rLHSVariables.size() != rLeftHandSideMatrices.size() )
      rLeftHandSideMatrices.resize(rLHSVariables.size());

    if( rRHSVariables.size() != rRightHandSideVectors.size() )
      rRightHandSideVectors.resize(rRHSVariables.size());
    
    LocalSystem.CalculationFlags.Set(RigidBodyPointRigidContactCondition::COMPUTE_LHS_MATRIX);
    for( unsigned int i=0; i<rLeftHandSideMatrices.size(); i++ )
      {
	//Note: rRightHandSideVectors.size() > 0
	this->InitializeSystemMatrices( rLeftHandSideMatrices[i], rRightHandSideVectors[0], LocalSystem.CalculationFlags );
      }

    LocalSystem.CalculationFlags.Set(RigidBodyPointRigidContactCondition::COMPUTE_RHS_VECTOR);
    LocalSystem.CalculationFlags.Set(RigidBodyPointRigidContactCondition::COMPUTE_LHS_MATRIX,false);

    for( unsigned int i=0; i<rRightHandSideVectors.size(); i++ )
      {
	//Note: rLeftHandSideMatrices.size() > 0
    	this->InitializeSystemMatrices( rLeftHandSideMatrices[0], rRightHandSideVectors[i], LocalSystem.CalculationFlags );
      }

    LocalSystem.CalculationFlags.Set(RigidBodyPointRigidContactCondition::COMPUTE_LHS_MATRIX,true);


    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrices(rLeftHandSideMatrices);
    LocalSystem.SetRightHandSideVectors(rRightHandSideVectors);

    LocalSystem.SetLeftHandSideVariables(rLHSVariables);
    LocalSystem.SetRightHandSideVariables(rRHSVariables);

    //Calculate condition system
    this->CalculateConditionSystem( LocalSystem, rCurrentProcessInfo );

}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rMassMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rDampingMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
							      GeneralVariables& rVariables,
							      double& rIntegrationWeight)

{
    KRATOS_TRY


    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();


    if( rVariables.Options.Is(ACTIVE)){

      //Force
      Matrix ForceMatrix  = ZeroMatrix(3);

      noalias(ForceMatrix) = rVariables.Penalty.Normal * rIntegrationWeight  * custom_outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal);
      

      for(unsigned int i=0; i<dimension; i++)
	{
	  for(unsigned int j=0; j<dimension; j++)
	    {
	      rLeftHandSideMatrix(i,j) += ForceMatrix(i,j);
	    }
	}
	      

      //Moment
      Matrix MomentMatrix = ZeroMatrix(3);
      
      MomentMatrix = prod(ForceMatrix,rVariables.SkewSymDistance);
      
      for(unsigned int i=0; i<dimension; i++)
	{
	  for(unsigned int j=0; j<dimension; j++)
	    {
	      rLeftHandSideMatrix(i+dimension,j) -= MomentMatrix(i,j);
	    }
	}

      

      // std::cout<<std::endl;
      // std::cout<<" Penalty.Normal "<<rVariables.Penalty.Normal<<" rVariables.Gap.Normal "<<rVariables.Gap.Normal<<" rVariables.Surface.Normal "<<rVariables.Surface.Normal<<" rVariables.Surface.Tangent "<<rVariables.Surface.Tangent<<" rIntegrationWeight "<<rIntegrationWeight<<" nxn : "<<custom_outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal)<<std::endl;

      this->CalculateAndAddKuugTangent( rLeftHandSideMatrix,  rVariables, rIntegrationWeight);
      // std::cout<<std::endl;
      // std::cout<<" Kcont "<<rLeftHandSideMatrix<<std::endl;

    }
    else{

      rLeftHandSideMatrix= ZeroMatrix(dimension*2,dimension*2);

    }
 
    //KRATOS_WATCH( rLeftHandSideMatrix )

    KRATOS_CATCH( "" )
}


//************* Tangent Contact Force constitutive matrix      **********************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables, double& rIntegrationWeight)
{

  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();


  double NormalForceModulus = 0;
  NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

  double TangentRelativeMovement = 0;
  TangentRelativeMovement = this->CalculateTangentRelativeMovement( TangentRelativeMovement, rVariables );
       
  double TangentForceModulus = this->CalculateCoulombsFrictionLaw( TangentRelativeMovement, NormalForceModulus, rVariables );


  //Force
  Matrix ForceMatrix  = ZeroMatrix(3);


  if( fabs(TangentForceModulus) >= 1e-25 ){
       
    if ( mTangentialVariables.Slip ) {
      //simpler expression:
      noalias(ForceMatrix) =  mTangentialVariables.Sign * mTangentialVariables.FrictionCoefficient * rVariables.Penalty.Normal * rIntegrationWeight * ( custom_outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Normal) );

      //added extra term
      //ForceMatrix +=  mTangentialVariables.FrictionCoefficient * rVariables.Penalty.Normal * rIntegrationWeight * fabs(rVariables.Gap.Normal) * custom_outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Tangent);

      //added extra term, maybe not necessary
      //ForceMatrix -=  mTangentialVariables.Sign * mTangentialVariables.FrictionCoefficient * rVariables.Penalty.Normal * rIntegrationWeight * ( custom_outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Normal) + rVariables.Gap.Normal * custom_outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal) );

    }
    else {
      noalias(ForceMatrix) =  rVariables.Penalty.Tangent * rIntegrationWeight * custom_outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Tangent);

    }

  }
       

  //noalias(ForceMatrix) = rVariables.Penalty.Normal * rIntegrationWeight  * custom_outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal);
      

  for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
	{
	  rLeftHandSideMatrix(i,j) += ForceMatrix(i,j);
	}
    }
	      
  // std::cout<<" KuuT "<<ForceMatrix<<std::endl;

  //Moment
  Matrix MomentMatrix = ZeroMatrix(3);
      
  MomentMatrix = prod(ForceMatrix,rVariables.SkewSymDistance);
      
  for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
	{
	  rLeftHandSideMatrix(i+dimension,j) -= MomentMatrix(i,j);
	}
    }


}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::CalculateAndAddContactForces(VectorType& rRightHandSideVector,
							      GeneralVariables& rVariables,
							      double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( rVariables.Options.Is(ACTIVE)){

       this->CalculateAndAddNormalContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );
       this->CalculateAndAddTangentContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );

    }
    else{

      rRightHandSideVector = ZeroVector(dimension*2);
    
    }

    KRATOS_CATCH( "" )
}


//**************************** Calculate Normal Contact Force ***********************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::CalculateAndAddNormalContactForce(VectorType& rRightHandSideVector,
									    GeneralVariables& rVariables,
									    double& rIntegrationWeight)
{
      
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();


  double NormalForceModulus = 0;
  NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

  NormalForceModulus *= (-1) * rIntegrationWeight;
   
  VectorType ContactForceVector = ZeroVector(3);

  for (unsigned int i = 0; i < dimension ; ++i) {
    ContactForceVector[i]    = NormalForceModulus * rVariables.Surface.Normal[i];
    rRightHandSideVector[i] += ContactForceVector[i];
  }

  GetGeometry()[0].SetLock();

  array_1d<double, 3 >& ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);


  for(unsigned int j = 0; j < dimension; j++)
    {
      ContactForce[j] += NormalForceModulus * rVariables.Surface.Normal[j];
    }

  GetGeometry()[0].UnSetLock();

  VectorType ContactTorque = ZeroVector(3);

  //ContactTorque = MathUtils<double>::CrossProduct( rVariables.CentroidPosition, ContactForceVector);      
  //std::cout<<" [ContactTorqueA]: "<<ContactTorque;

  ContactTorque = prod(rVariables.SkewSymDistance,ContactForceVector); //  = (D x F)
       
  // std::cout<<" [ContactTorqueB]: "<<ContactTorque;
  // std::cout<<" [ContactForce]: "<<ContactForceVector;
  // std::cout<<" [Normal]: "<<rVariables.Surface.Normal;
  // std::cout<<" [Distance]: "<<rVariables.SkewSymDistance;
  // std::cout<<std::endl;

  //Contact torque due to contact force on beam surface
  for (unsigned int i =0; i < dimension; ++i) {
    rRightHandSideVector[i+dimension] += ContactTorque[i];
  }
       

  //std::cout<<" Fcont "<<rRightHandSideVector<<std::endl;

}

//**************************** Calculate Tangent Contact Force **********************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::CalculateAndAddTangentContactForce(VectorType& rRightHandSideVector,
									     GeneralVariables& rVariables,
									     double& rIntegrationWeight)
{

  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

  double NormalForceModulus = 0;
  NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

  double TangentRelativeMovement = 0;
  TangentRelativeMovement = this->CalculateTangentRelativeMovement( TangentRelativeMovement, rVariables );
       
  double TangentForceModulus =  this->CalculateCoulombsFrictionLaw( TangentRelativeMovement, NormalForceModulus, rVariables );

  TangentForceModulus *= (-1) * rIntegrationWeight;


  // std::cout<< "["<<mTangentialVariables.Sign<<"] Tangent Force Node ["<<GetGeometry()[0].Id()<<" ]:"<<TangentForceModulus<<" RelativeMovement: "<<TangentRelativeMovement<<" Tangent Gap: "<<rVariables.Gap.Tangent<<" SLIP ["<<mTangentialVariables.Slip<<"]"<<std::endl; 


  VectorType ContactForceVector = ZeroVector(3);

  for (unsigned int i = 0; i < dimension ; ++i) {
    ContactForceVector[i]    = TangentForceModulus * rVariables.Surface.Tangent[i];
    rRightHandSideVector[i] += ContactForceVector[i];
  }


  GetGeometry()[0].SetLock();

  array_1d<double, 3 >& ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);


  for(unsigned int j = 0; j < dimension; j++)
    {
      ContactForce[j] += TangentForceModulus * rVariables.Surface.Tangent[j];
    }

  GetGeometry()[0].UnSetLock();

  // if( fabs(inner_prod(rVariables.Surface.Tangent, rVariables.Surface.Normal)) > 0.001 ){
  //   std::cout<<" Tangent "<<rVariables.Surface.Tangent<<" Normal "<<rVariables.Surface.Normal<<" prod "<<inner_prod(rVariables.Surface.Tangent, rVariables.Surface.Normal)<<std::endl;
  // }

  VectorType ContactTorque = ZeroVector(3);

  //ContactTorque = MathUtils<double>::CrossProduct( rVariables.CentroidPosition, ContactForceVector);      
  ContactTorque = prod(rVariables.SkewSymDistance,ContactForceVector); //  = (D x F)

  // std::cout<<" [ContactTorque]: "<<ContactTorque;
  // std::cout<<" [ContactForce]:  "<<ContactForceVector;
  // std::cout<<" [Normal]:  "<<rVariables.Surface.Normal;
  // std::cout<<std::endl;
  
  //Contact torque due to contact tangent force on beam surface
  for (unsigned int i =0; i < dimension; ++i) {
    rRightHandSideVector[i+dimension] += ContactTorque[i];
  }
     
}

//**************************** Calculate Normal Force Modulus ***********************
//***********************************************************************************

double& RigidBodyPointRigidContactCondition::CalculateNormalForceModulus ( double& rNormalForceModulus, GeneralVariables& rVariables )
{

  rNormalForceModulus = (rVariables.Penalty.Normal * rVariables.Gap.Normal); 

  return rNormalForceModulus;

}

//**************************** Calculate Tangent Force Modulus **********************
//***********************************************************************************

double& RigidBodyPointRigidContactCondition::CalculateTangentRelativeMovement( double& rTangentRelativeMovement, GeneralVariables& rVariables )
{
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

  const array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
  const array_1d<double, 3> & PreviousDisplacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, 1);
  array_1d<double, 3 > DeltaDisplacement           = CurrentDisplacement-PreviousDisplacement;

  // bool Regularization = false;

  // if( Regularization == true ){

  // 	 //regularization taking the contigous boundary segments
  // 	 WeakPointerVector<Node<3> >& rN = GetGeometry()[0].GetValue(NEIGHBOUR_NODES);

  // 	 array_1d<double, 3 > NeighbDeltaDisplacement;

  // 	 double counter = 0;
       
  // 	 for(unsigned int i = 0; i < rN.size(); i++)
  // 	   {
  // 	     if(rN[i].Is(BOUNDARY)){
	          
  // 	       const array_1d<double, 3> & NeighbCurrentDisplacement  =  GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
  // 	       const array_1d<double, 3> & NeighbPreviousDisplacement =  GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, 1);

  // 	       array_1d<double, 3 > NeighbDeltaDisplacement           =  NeighbCurrentDisplacement-NeighbPreviousDisplacement;

  // 	       DeltaDisplacement += NeighbDeltaDisplacement;
     
  // 	       counter ++;
  // 	     }
  // 	   }
  //        if( counter!= 0)
  // 	    DeltaDisplacement /= counter;

  // }

  VectorType WallDisplacement = mTangentialVariables.DeltaTime * this->mpRigidWall->Velocity();
       
  rTangentRelativeMovement = 0.0;
  double WallTangentRelativeMovement    = 0.0;

  for (unsigned int i = 0; i < dimension; ++i)
    {
      rTangentRelativeMovement    += DeltaDisplacement[i] * rVariables.Surface.Tangent[i];
      WallTangentRelativeMovement += WallDisplacement[i]  * rVariables.Surface.Tangent[i];  
    }


  rTangentRelativeMovement -= WallTangentRelativeMovement;      

  rVariables.Gap.Tangent = rTangentRelativeMovement;


  return rTangentRelativeMovement;

}

//**************************** Check Coulomb law for Tangent Contact Force **********
//***********************************************************************************

double RigidBodyPointRigidContactCondition::CalculateCoulombsFrictionLaw(double & rTangentRelativeMovement, double & rNormalForceModulus , GeneralVariables& rVariables)
{
  mTangentialVariables.FrictionCoefficient = this->CalculateFrictionCoefficient(rTangentRelativeMovement);

 
  double TangentForceModulus = rVariables.Penalty.Tangent * rVariables.Gap.Tangent; //+ mTangentialVariables.PreviousTangentForceModulus; 
     
  //std::cout<<" Gap.Tangent "<<rVariables.Gap.Tangent<<std::endl;

 
  if ( fabs(TangentForceModulus) >  mTangentialVariables.FrictionCoefficient * fabs(rNormalForceModulus) && fabs(rVariables.Gap.Tangent) > 1e-200) {

    mTangentialVariables.Sign =  rVariables.Gap.Tangent/ fabs( rVariables.Gap.Tangent ) ; 

    TangentForceModulus =  mTangentialVariables.Sign * mTangentialVariables.FrictionCoefficient * fabs(rNormalForceModulus) ;
    mTangentialVariables.Slip = true;

  }
  else {
    mTangentialVariables.Slip = false;
  }


  return TangentForceModulus;
}

  

//**************************** Check friction coefficient ***************************
//***********************************************************************************

double RigidBodyPointRigidContactCondition::CalculateFrictionCoefficient(double & rTangentRelativeMovement)
{

  //---FRICTION LAW in function of the relative sliding velocity ---//

  double Velocity = 0;
  Velocity = rTangentRelativeMovement / mTangentialVariables.DeltaTime;

 
  //Addicional constitutive parameter  C
  //which describes how fast the static coefficient approaches the dynamic:
  double C = 2;//0.1;

  //Addicional constitutive parameter  E
  //regularization parameter (->0, classical Coulomb law)
  double E = 0;//0.01;


  double FrictionCoefficient = mTangentialVariables.DynamicFrictionCoefficient + ( mTangentialVariables.StaticFrictionCoefficient-  mTangentialVariables.DynamicFrictionCoefficient ) * exp( (-1) * C * fabs(Velocity) );


  //Square root regularization
  FrictionCoefficient *= fabs(Velocity)/sqrt( ( Velocity * Velocity ) + ( E * E ) );
       
  //Hyperbolic regularization
  //FrictionCoefficient *= tanh( fabs(Velocity)/E );

  return FrictionCoefficient;

}

//************************************************************************************
//************************************************************************************

void RigidBodyPointRigidContactCondition::VectorToSkewSymmetricTensor( const Vector& rVector, 
							       Matrix& rSkewSymmetricTensor )
{
    KRATOS_TRY

    //Initialize Local Matrices
    if( rSkewSymmetricTensor.size1() != 3 )
      rSkewSymmetricTensor.resize(3, 3, false);
    
    rSkewSymmetricTensor = ZeroMatrix(3);

    rSkewSymmetricTensor( 0, 1 ) = -rVector[2];
    rSkewSymmetricTensor( 0, 2 ) =  rVector[1];
    rSkewSymmetricTensor( 1, 2 ) = -rVector[0];

    rSkewSymmetricTensor( 1, 0 ) =  rVector[2];
    rSkewSymmetricTensor( 2, 0 ) = -rVector[1];
    rSkewSymmetricTensor( 2, 1 ) =  rVector[0];


    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

inline Condition::MatrixType RigidBodyPointRigidContactCondition::custom_outer_prod(const array_1d<double, 3>& a, const array_1d<double, 3>& b)
{
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

  Condition::MatrixType A(dimension,dimension);
    
  A(0,0)=a[0]*b[0];
  A(0,1)=a[0]*b[1];
  A(1,0)=a[1]*b[0];
  A(1,1)=a[1]*b[1];
  if( dimension == 3 ){
    A(0,2)=a[0]*b[2];
    A(1,2)=a[1]*b[2];
    A(2,0)=a[2]*b[0];
    A(2,1)=a[2]*b[1];
    A(2,2)=a[2]*b[2];
  }

  return A;
}

//***********************************************************************************
//***********************************************************************************


int RigidBodyPointRigidContactCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
  return 0;
}

//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
