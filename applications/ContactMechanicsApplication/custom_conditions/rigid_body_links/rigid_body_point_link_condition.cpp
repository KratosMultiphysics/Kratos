//
//   Project Name:        KratosStringDynamicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:             October 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/rigid_body_point_link_condition.hpp"

#include "contact_mechanics_application_variables.h"

namespace Kratos
{

/**
 * Flags related to the condition computation
 */
KRATOS_CREATE_LOCAL_FLAG( RigidBodyPointLinkCondition, COMPUTE_RHS_VECTOR,                 0 );
KRATOS_CREATE_LOCAL_FLAG( RigidBodyPointLinkCondition, COMPUTE_LHS_MATRIX,                 1 );
KRATOS_CREATE_LOCAL_FLAG( RigidBodyPointLinkCondition, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
KRATOS_CREATE_LOCAL_FLAG( RigidBodyPointLinkCondition, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );

//***********************************************************************************
//***********************************************************************************
RigidBodyPointLinkCondition::RigidBodyPointLinkCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{

   for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
	GetGeometry()[i].Set(SLAVE); //Flag to set MASTER_ELEMENTS in that nodes (if is SLAVE, a MASTER is required)
      }
    //DO NOT ADD DOFS HERE!!!
}

//***********************************************************************************
//***********************************************************************************
RigidBodyPointLinkCondition::RigidBodyPointLinkCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
	GetGeometry()[i].Set(SLAVE); //Flag to set MASTER_ELEMENTS in that nodes (if is SLAVE, a MASTER is required)
      }

    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
RigidBodyPointLinkCondition::RigidBodyPointLinkCondition( RigidBodyPointLinkCondition const& rOther )
    : Condition(rOther)
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer RigidBodyPointLinkCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new RigidBodyPointLinkCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}


//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer RigidBodyPointLinkCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  return this->Create( NewId, rThisNodes, pGetProperties() );
}


//***********************************************************************************
//***********************************************************************************
RigidBodyPointLinkCondition::~RigidBodyPointLinkCondition()
{
}



//************* GETTING METHODS

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::GetDofList(DofsVectorType& rConditionDofList,
					     ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rConditionDofList.clear();
    rConditionDofList.resize(0);

    const unsigned int inode = GetGeometry().PointsNumber()-1;

    WeakPointerVector<Element>& SlaveElements  = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

    Element& MasterElement = (GetGeometry()[inode].GetValue(MASTER_ELEMENTS)).back();

    for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ie++)
      {

	DofsVectorType SlaveConditionDofList;
	ie->GetDofList(SlaveConditionDofList, rCurrentProcessInfo);

	for(unsigned int i=0; i<SlaveConditionDofList.size(); i++)
	  rConditionDofList.push_back(SlaveConditionDofList[i]);

	DofsVectorType MasterConditionDofList;
	MasterElement.GetDofList(MasterConditionDofList, rCurrentProcessInfo);

	for(unsigned int i=0; i<MasterConditionDofList.size(); i++)
	  rConditionDofList.push_back(MasterConditionDofList[i]);
      }



    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::EquationIdVector(EquationIdVectorType& rResult,
						   ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rResult.clear();
    rResult.resize( 0, false );

    const unsigned int inode = GetGeometry().PointsNumber()-1;

    WeakPointerVector<Element>& SlaveElements  = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));
    Element& MasterElement = (GetGeometry()[inode].GetValue(MASTER_ELEMENTS)).back();


    for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ie++)
      {
	EquationIdVectorType SlaveResult;
	ie->EquationIdVector(SlaveResult, rCurrentProcessInfo);

	for(unsigned int i=0; i<SlaveResult.size(); i++)
	  rResult.push_back(SlaveResult[i]);

	EquationIdVectorType MasterResult;
	MasterElement.EquationIdVector(MasterResult, rCurrentProcessInfo);

	for(unsigned int i=0; i<MasterResult.size(); i++)
	  rResult.push_back(MasterResult[i]);
      }



    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::GetValuesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    rValues.clear();
    rValues.resize(0);

    const unsigned int inode = GetGeometry().PointsNumber()-1;

    WeakPointerVector<Element>& SlaveElements  = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));
    Element& MasterElement = (GetGeometry()[inode].GetValue(MASTER_ELEMENTS)).back();

    unsigned int indexi = 0;
    unsigned int sizei  = 0;
    for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ie++)
      {
	Vector SlaveValues;
	ie->GetValuesVector(SlaveValues,Step);

	sizei += SlaveValues.size();
	rValues.resize(sizei,true);

	for(unsigned int i=0; i<SlaveValues.size(); i++)
	  rValues[indexi+i] = SlaveValues[i];

	indexi += SlaveValues.size();

	Vector MasterValues;
	MasterElement.GetValuesVector(MasterValues, Step);

	sizei += MasterValues.size();

	for(unsigned int i=0; i<MasterValues.size(); i++)
	  rValues[indexi+i] = MasterValues[i];

	indexi += MasterValues.size();
      }


    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    KRATOS_TRY

    rValues.clear();
    rValues.resize(0);

    const unsigned int inode = GetGeometry().PointsNumber()-1;

    WeakPointerVector<Element>& SlaveElements  = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

    Element& MasterElement = (GetGeometry()[inode].GetValue(MASTER_ELEMENTS)).back();

    unsigned int indexi = 0;
    unsigned int sizei  = 0;
    for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ie++)
      {
	Vector SlaveValues;
	ie->GetFirstDerivativesVector(SlaveValues,Step);

	sizei += SlaveValues.size();
	rValues.resize(sizei,true);

	for(unsigned int i=0; i<SlaveValues.size(); i++)
	  rValues[indexi+i] = SlaveValues[i];

	indexi += SlaveValues.size();

	Vector MasterValues;
	MasterElement.GetFirstDerivativesVector(MasterValues, Step);

	sizei += MasterValues.size();

	for(unsigned int i=0; i<MasterValues.size(); i++)
	  rValues[indexi+i] = MasterValues[i];

	indexi += MasterValues.size();
      }

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    KRATOS_TRY

    rValues.clear();
    rValues.resize(0);

    const unsigned int inode = GetGeometry().PointsNumber()-1;

    WeakPointerVector<Element>& SlaveElements  = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));
    Element& MasterElement = (GetGeometry()[inode].GetValue(MASTER_ELEMENTS)).back();

    unsigned int indexi = 0;
    unsigned int sizei  = 0;
    for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ie++)
      {
	Vector SlaveValues;
	ie->GetSecondDerivativesVector(SlaveValues,Step);

	sizei += SlaveValues.size();
	rValues.resize(sizei,true);

	for(unsigned int i=0; i<SlaveValues.size(); i++)
	  rValues[indexi+i] = SlaveValues[i];

	indexi += SlaveValues.size();

	Vector MasterValues;
	MasterElement.GetSecondDerivativesVector(MasterValues, Step);

	sizei += MasterValues.size();

	for(unsigned int i=0; i<MasterValues.size(); i++)
	  rValues[indexi+i] = MasterValues[i];

	indexi += MasterValues.size();
      }




    KRATOS_CATCH( "" )
}



//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::AddExplicitContribution(const VectorType& rRHSVector,
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

void RigidBodyPointLinkCondition::Initialize()
{
    KRATOS_TRY

      //std::cout<<" INITIALIZE LINK "<<std::endl;

    //Fix linked deformable point dofs
    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
      {
	(GetGeometry()[i].pGetDof(DISPLACEMENT_X))->FixDof();
	(GetGeometry()[i].pGetDof(DISPLACEMENT_Y))->FixDof();
	(GetGeometry()[i].pGetDof(DISPLACEMENT_Z))->FixDof();

	(GetGeometry()[i].pGetDof(ROTATION_X))->FixDof();
	(GetGeometry()[i].pGetDof(ROTATION_Y))->FixDof();
	(GetGeometry()[i].pGetDof(ROTATION_Z))->FixDof();

	//std::cout<<" NODE:"<<GetGeometry()[i].Id()<<" FIXED: "<<(GetGeometry()[i].pGetDof(DISPLACEMENT_X))->IsFixed()<<std::endl;
      }


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    mIteration = 0;

    const unsigned int inode = GetGeometry().PointsNumber()-1;

    PointPointerType  mpSlaveNode = GetGeometry()(inode);


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
void RigidBodyPointLinkCondition::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

      mIteration +=1;

    // const unsigned int inode = GetGeometry().PointsNumber()-1;

    // PointPointerType  mpSlaveNode = GetGeometry()(inode);

    // array_1d<double, 3 >&  Displacement = mpSlaveNode->FastGetSolutionStepValue(DISPLACEMENT);
    // array_1d<double, 3 >&  Rotation     = mpSlaveNode->FastGetSolutionStepValue(ROTATION);
    // array_1d<double, 3 >&  StepRotation = mpSlaveNode->FastGetSolutionStepValue(STEP_ROTATION);


    // std::cout<<" Link Initialize RB iteration [Rotation:"<<Rotation<<"; StepRota:"<<StepRotation<<"; Disp:"<<Displacement<<"]"<<std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY



    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
  KRATOS_TRY


  KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
							  VectorType& rRightHandSideVector,
							  const unsigned int& rSlaveElementSize,
							  Flags& rCalculationFlags)

{

    //resizing as needed the LHS
    unsigned int MatSize = rSlaveElementSize + 6; //slave body element + rigid body element

    if ( rCalculationFlags.Is(RigidBodyPointLinkCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }


    //resizing as needed the RHS
    if ( rCalculationFlags.Is(RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
	    rRightHandSideVector.resize( MatSize, false );

	rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS

    }
}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::InitializeGeneralVariables(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int inode = GetGeometry().PointsNumber()-1;

    // set the slave node
    GeometryType& SlaveGeometry = rVariables.pSlaveElement->GetGeometry();

    for ( unsigned int i = 0; i < SlaveGeometry.size(); i++ )
      {
	if( SlaveGeometry[i].Id() == GetGeometry()[inode].Id() )
	  rVariables.SlaveNode = i;
      }

    rVariables.Distance = ZeroVector(3);
    rVariables.SkewSymDistance = ZeroMatrix(3,3);
    rVariables.LagrangeMultipliers   = ZeroVector(3);


    //compute distance from the slave node to the master rigid body element node (center of gravity)
    Element& MasterElement = (GetGeometry()[inode].GetValue(MASTER_ELEMENTS)).back();

    //rVariables.Distance = MasterElement.GetGeometry()[0].Coordinates() - GetGeometry()[inode].Coordinates();

    rVariables.Distance = GetGeometry()[inode].Coordinates() - MasterElement.GetGeometry()[0].Coordinates();
    //std::cout<<" Distance "<<norm_2(rVariables.Distance)<<std::endl;

    //compute the skewsymmmetric tensor of the distance
    this->VectorToSkewSymmetricTensor(rVariables.Distance, rVariables.SkewSymDistance);

    KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateConditionSystem(LocalSystemComponents& rLocalSystem,
							   LocalSystemComponents& rLinkedSystem,
							   ElementType::Pointer& pSlaveElement,
							   ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


      //std::cout<<" [ID:"<<this->Id()<<"]"<<std::endl;

    //create and initialize condition variables:
    GeneralVariables Variables;
    Variables.pSlaveElement = pSlaveElement;

    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    if ( rLocalSystem.CalculationFlags.Is(RigidBodyPointLinkCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
      {
	//contributions to stiffness matrix calculated on the reference config
	this->CalculateAndAddLHS ( rLocalSystem, rLinkedSystem, Variables );
      }

    if ( rLocalSystem.CalculationFlags.Is(RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
      {
	//contribution to external forces
	this->CalculateAndAddRHS ( rLocalSystem, rLinkedSystem, Variables );
      }


    KRATOS_CATCH( "" )
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, LocalSystemComponents& rLinkedSystem, GeneralVariables& rVariables)
{

  //contributions of the stiffness matrix calculated on the reference configuration
  if( rLocalSystem.CalculationFlags.Is( RigidBodyPointLinkCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS ) )
    {
      std::vector<MatrixType>& rLeftHandSideMatrices = rLocalSystem.GetLeftHandSideMatrices();
      const std::vector< Variable< MatrixType > >& rLeftHandSideVariables = rLocalSystem.GetLeftHandSideVariables();

      std::vector<MatrixType>& rLinkedLeftHandSideMatrices = rLinkedSystem.GetLeftHandSideMatrices();

      for( unsigned int i=0; i<rLeftHandSideVariables.size(); i++ )
	{
	  bool calculated = false;

	  if( rLeftHandSideVariables[i] == GEOMETRIC_STIFFNESS_MATRIX ){
	    // operation performed: add Kg to the rLefsHandSideMatrix
	    this->CalculateAndAddStiffness( rLeftHandSideMatrices[i], rLinkedLeftHandSideMatrices[i], rVariables );
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
    MatrixType& rLinkedLeftHandSideMatrix = rLinkedSystem.GetLeftHandSideMatrix();

    // operation performed: add Kg to the rLefsHandSideMatrix
    this->CalculateAndAddStiffness( rLeftHandSideMatrix, rLinkedLeftHandSideMatrix, rVariables );

    //KRATOS_WATCH( rLeftHandSideMatrix )
  }

}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, LocalSystemComponents& rLinkedSystem, GeneralVariables& rVariables)
{
    //contribution of the internal and external forces
    if( rLocalSystem.CalculationFlags.Is( RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS ) )
    {

      std::vector<VectorType>& rRightHandSideVectors = rLocalSystem.GetRightHandSideVectors();
      const std::vector< Variable< VectorType > >& rRightHandSideVariables = rLocalSystem.GetRightHandSideVariables();

      std::vector<VectorType>& rLinkedRightHandSideVectors = rLinkedSystem.GetRightHandSideVectors();

      for( unsigned int i=0; i<rRightHandSideVariables.size(); i++ )
	{
	  bool calculated = false;
	  if( rRightHandSideVariables[i] == INTERNAL_FORCES_VECTOR ){
	    // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
	    this->CalculateAndAddForces( rRightHandSideVectors[i], rLinkedRightHandSideVectors[i], rVariables);
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
      VectorType& rLinkedRightHandSideVector = rLinkedSystem.GetRightHandSideVector();

      // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
      this->CalculateAndAddForces( rRightHandSideVector, rLinkedRightHandSideVector, rVariables);

      //KRATOS_WATCH( rRightHandSideVector )

    }


}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
							VectorType& rRightHandSideVector,
							ProcessInfo& rCurrentProcessInfo )
{
    //set sizes to zero
    rLeftHandSideMatrix.clear();
    rLeftHandSideMatrix.resize(0,0);

    rRightHandSideVector.clear();
    rRightHandSideVector.resize(0);


    //Ask to the linked deformable element the LocalRightHandSide (only one element link)
    const unsigned int inode = GetGeometry().PointsNumber()-1;
    WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

    //std::cout<<" [ID:"<<this->Id()<<"] [iter:"<<mIteration<<"] [SlaveElements: "<<SlaveElements.size()<<"] [NodeId "<<GetGeometry()[0].Id()<<"]"<<std::endl;

    for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ie++)
      {
	//create local system components
	LocalSystemComponents LocalSystem;

	//calculation flags
	LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_LHS_MATRIX);
	LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR);

	MatrixType LocalLeftHandSideMatrix;
	VectorType LocalRightHandSideVector;

	//get dofs size of the slave body element
	EquationIdVectorType EquationIds;
	ie->EquationIdVector(EquationIds,rCurrentProcessInfo);
	unsigned int SlaveElementSize = EquationIds.size();

	//Initialize sizes for the system components:
	this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, SlaveElementSize, LocalSystem.CalculationFlags );

	//Set Variables to Local system components
	LocalSystem.SetLeftHandSideMatrix(LocalLeftHandSideMatrix);
	LocalSystem.SetRightHandSideVector(LocalRightHandSideVector);

	//create linked system components
	LocalSystemComponents LinkedSystem;

	MatrixType LeftHandSideMatrix;
	VectorType RightHandSideVector;
	//std::cout<<"[LINK]"<<std::endl;
	ie->CalculateLocalSystem(LeftHandSideMatrix,RightHandSideVector,rCurrentProcessInfo);

	LinkedSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
	LinkedSystem.SetRightHandSideVector(RightHandSideVector);

	//Calculate condition system
	ElementType::Pointer SlaveElement = ElementType::Pointer(*(ie.base()));
	this->CalculateConditionSystem( LocalSystem, LinkedSystem, SlaveElement, rCurrentProcessInfo );

	//assemble the local system into the global system

	//resizing as needed the LHS
	unsigned int GlobalSize1 = rLeftHandSideMatrix.size1();
	unsigned int LocalSize1  = LocalLeftHandSideMatrix.size1();
	unsigned int MatSize1    = GlobalSize1+LocalSize1;

	unsigned int GlobalSize2 = rLeftHandSideMatrix.size2();
	unsigned int LocalSize2  = LocalLeftHandSideMatrix.size2();
	unsigned int MatSize2    = GlobalSize2+LocalSize2;

	rLeftHandSideMatrix.resize( MatSize1, MatSize2, true );

	unsigned int indexi = 0;
	for(unsigned int i=GlobalSize1; i<MatSize1; i++)
	  {
	    unsigned int indexj = 0;
	    for(unsigned int j=GlobalSize2; j<MatSize2; j++)
	      {
		rLeftHandSideMatrix(i,j) = LocalLeftHandSideMatrix(indexi,indexj);
		indexj++;
	      }
	    indexi++;
	  }

	//resizing as needed the RHS
	GlobalSize1 = rRightHandSideVector.size();
	LocalSize1  = LocalRightHandSideVector.size();
	MatSize1    = GlobalSize1+LocalSize1;

	rRightHandSideVector.resize( MatSize1, true );

	indexi = 0;
	for(unsigned int i=GlobalSize1; i<MatSize1; i++)
	  {
	    rRightHandSideVector[i]=LocalRightHandSideVector[indexi];
	    indexi++;
	  }

      }

    //std::cout<<" LINK RighHandSide "<<rRightHandSideVector<<std::endl;

    // KRATOS_WATCH( rLeftHandSideMatrix )
    // KRATOS_WATCH( rRightHandSideVector )

}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
									  VectorType& rRightHandSideVector,
									  ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //set sizes to zero
    rLeftHandSideMatrix.clear();
    rLeftHandSideMatrix.resize(0,0);

    rRightHandSideVector.clear();
    rRightHandSideVector.resize(0);


    //Ask to the linked deformable element the LocalRightHandSide (only one element link)
    const unsigned int inode = GetGeometry().PointsNumber()-1;
    WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

    for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ie++)
      {

	//create local system components
	LocalSystemComponents LocalSystem;

	//calculation flags
	LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_LHS_MATRIX);
	LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR);

	MatrixType LocalLeftHandSideMatrix;
	VectorType LocalRightHandSideVector;

	//get dofs size of the slave body element
	EquationIdVectorType EquationIds;
	//std::cout<<"  [ Node :"<<GetGeometry()[inode].Id() <<" Linked to Elements :"<<ie->Id()<<" ] "<<std::endl;
	ie->EquationIdVector(EquationIds,rCurrentProcessInfo);
	unsigned int SlaveElementSize = EquationIds.size();

	//Initialize sizes for the system components:
	this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, SlaveElementSize, LocalSystem.CalculationFlags );

	//Set Variables to Local system components
	LocalSystem.SetLeftHandSideMatrix(LocalLeftHandSideMatrix);
	LocalSystem.SetRightHandSideVector(LocalRightHandSideVector);


	//create linked system components
	LocalSystemComponents LinkedSystem;

	MatrixType LeftHandSideMatrix;
	VectorType RightHandSideVector;
	//std::cout<<"[LINK]"<<std::endl;
	ie->CalculateSecondDerivativesContributions(LeftHandSideMatrix,RightHandSideVector,rCurrentProcessInfo);

	LinkedSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
	LinkedSystem.SetRightHandSideVector(RightHandSideVector);

	//Calculate condition system
	ElementType::Pointer SlaveElement = ElementType::Pointer(*(ie.base()));
	this->CalculateConditionSystem( LocalSystem, LinkedSystem, SlaveElement, rCurrentProcessInfo );

	//assemble the local system into the global system

	//resizing as needed the LHS
	unsigned int GlobalSize1 = rLeftHandSideMatrix.size1();
	unsigned int LocalSize1  = LocalLeftHandSideMatrix.size1();
	unsigned int MatSize1    = GlobalSize1+LocalSize1;

	unsigned int GlobalSize2 = rLeftHandSideMatrix.size2();
	unsigned int LocalSize2  = LocalLeftHandSideMatrix.size2();
	unsigned int MatSize2    = GlobalSize2+LocalSize2;

	rLeftHandSideMatrix.resize( MatSize1, MatSize2, true );

	unsigned int indexi = 0;
	for(unsigned int i=GlobalSize1; i<MatSize1; i++)
	  {
	    unsigned int indexj = 0;
	    for(unsigned int j=GlobalSize2; j<MatSize2; j++)
	      {
		rLeftHandSideMatrix(i,j) = LocalLeftHandSideMatrix(indexi,indexj);
		indexj++;
	      }
	    indexi++;
	  }

	//resizing as needed the RHS
	GlobalSize1 = rRightHandSideVector.size();
	LocalSize1  = LocalRightHandSideVector.size();
	MatSize1    = GlobalSize1+LocalSize1;

	rRightHandSideVector.resize( MatSize1, true );

	indexi = 0;
	for(unsigned int i=GlobalSize1; i<MatSize1; i++)
	  {
	    rRightHandSideVector[i]=LocalRightHandSideVector[indexi];
	    indexi++;
	  }

      }

    //std::cout<<" LINK Dynamic RighHandSide "<<rRightHandSideVector<<std::endl;

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
							  ProcessInfo& rCurrentProcessInfo )
{
    //set sizes to zero
    rRightHandSideVector.clear();
    rRightHandSideVector.resize(0);

    MatrixType LocalLeftHandSideMatrix = Matrix();

    //Ask to the linked deformable element the LocalRightHandSide (only one element link)
    const unsigned int inode = GetGeometry().PointsNumber()-1;
    WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

    for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ie++)
      {

	//create local system components
	LocalSystemComponents LocalSystem;

	//calculation flags
	LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR);

	VectorType LocalRightHandSideVector;

	//get dofs size of the slave body element
	EquationIdVectorType EquationIds;
	ie->EquationIdVector(EquationIds,rCurrentProcessInfo);
	unsigned int SlaveElementSize = EquationIds.size();

	//Initialize sizes for the system components:
	this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, SlaveElementSize, LocalSystem.CalculationFlags );
	//Set Variables to Local system components
	LocalSystem.SetRightHandSideVector(LocalRightHandSideVector);

	//create linked system components
	LocalSystemComponents LinkedSystem;

	VectorType RightHandSideVector;
	ie->CalculateRightHandSide(RightHandSideVector,rCurrentProcessInfo);

	LinkedSystem.SetRightHandSideVector(RightHandSideVector);

	//Calculate condition system
	ElementType::Pointer SlaveElement = ElementType::Pointer(*(ie.base()));
	this->CalculateConditionSystem( LocalSystem, LinkedSystem, SlaveElement, rCurrentProcessInfo );


	//assemble the local system into the global system

	//resizing as needed the RHS
	unsigned int GlobalSize1 = rRightHandSideVector.size();
	unsigned int LocalSize1  = LocalRightHandSideVector.size();
	unsigned int MatSize1    = GlobalSize1+LocalSize1;

	rRightHandSideVector.resize( MatSize1, true );

	unsigned int indexi = 0;
	for(unsigned int i=GlobalSize1; i<MatSize1; i++)
	  {
	    rRightHandSideVector[i]=LocalRightHandSideVector[indexi];
	    indexi++;
	  }

      }
    //KRATOS_WATCH( rLeftHandSideMatrix )
    //KRATOS_WATCH( rRightHandSideVector )

}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateRightHandSide( std::vector< VectorType >& rRightHandSideVectors,
							  const std::vector< Variable< VectorType > >& rRHSVariables,
							  ProcessInfo& rCurrentProcessInfo )
{

    for( unsigned int i=0; i< rRightHandSideVectors.size(); i++)
      {
	rRightHandSideVectors[i].clear();
	rRightHandSideVectors[i].resize(0);
      }

    MatrixType LocalLeftHandSideMatrix = Matrix();


    //Ask to the linked deformable element the LocalRightHandSide (only one element link)
    const unsigned int inode = GetGeometry().PointsNumber()-1;
    WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

    for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ie++)
      {

	//create local system components
	LocalSystemComponents LocalSystem;

	//calculation flags
	LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);


	std::vector< VectorType > LocalRightHandSideVectors;

	//Initialize sizes for the system components:
	if( rRHSVariables.size() != LocalRightHandSideVectors.size() )
	  LocalRightHandSideVectors.resize(rRHSVariables.size());

	LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR);

	//get dofs size of the slave body element
	EquationIdVectorType EquationIds;
	ie->EquationIdVector(EquationIds,rCurrentProcessInfo);
	unsigned int SlaveElementSize = EquationIds.size();

	for( unsigned int i=0; i<LocalRightHandSideVectors.size(); i++ )
	  {
	    //Note: LocalLeftHandSideMatrices.size() > 0
	    this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVectors[i], SlaveElementSize, LocalSystem.CalculationFlags );
	  }

	//Set Variables to Local system components
	LocalSystem.SetRightHandSideVectors(LocalRightHandSideVectors);
	LocalSystem.SetRightHandSideVariables(rRHSVariables);

	//create linked system components
	LocalSystemComponents LinkedSystem;

	std::vector< VectorType > RightHandSideVectors;
	ie->CalculateRightHandSide(RightHandSideVectors,rRHSVariables,rCurrentProcessInfo);

	LinkedSystem.SetRightHandSideVectors(RightHandSideVectors);

	LinkedSystem.SetRightHandSideVariables(rRHSVariables);


	//Calculate condition system
	ElementType::Pointer SlaveElement = ElementType::Pointer(*(ie.base()));
	this->CalculateConditionSystem( LocalSystem, LinkedSystem, SlaveElement, rCurrentProcessInfo );

	//assemble the local system into the global system

	//resizing as needed the RHSs
	for( unsigned int i=0; i< rRightHandSideVectors.size(); i++)
	  {

	    unsigned int GlobalSize1 = rRightHandSideVectors[i].size();
	    unsigned int LocalSize1  = LocalRightHandSideVectors[i].size();
	    unsigned int MatSize1    = GlobalSize1+LocalSize1;

	    rRightHandSideVectors[i].resize( MatSize1, true );

	    unsigned int indexi = 0;
	    for(unsigned int j=GlobalSize1; j<MatSize1; j++)
	      {
		rRightHandSideVectors[i][j]=LocalRightHandSideVectors[i][indexi];
		indexi++;
	      }

	  }
      }

}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateLocalSystem( std::vector< MatrixType >& rLeftHandSideMatrices,
							const std::vector< Variable< MatrixType > >& rLHSVariables,
							std::vector< VectorType >& rRightHandSideVectors,
							const std::vector< Variable< VectorType > >& rRHSVariables,
							ProcessInfo& rCurrentProcessInfo )
{

    for( unsigned int i=0; i< rLeftHandSideMatrices.size(); i++)
      {
	rLeftHandSideMatrices[i].clear();
	rLeftHandSideMatrices[i].resize(0,0);
      }

    for( unsigned int i=0; i< rRightHandSideVectors.size(); i++)
      {
	rRightHandSideVectors[i].clear();
	rRightHandSideVectors[i].resize(0);
      }

   //Ask to the linked deformable element the LocalRightHandSide (only one element link)
    const unsigned int inode = GetGeometry().PointsNumber()-1;
    WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

    for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ie++)
      {

	//create local system components
	LocalSystemComponents LocalSystem;

	//calculation flags
	LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_LHS_MATRIX_WITH_COMPONENTS);
	LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR_WITH_COMPONENTS);

	std::vector< MatrixType > LocalLeftHandSideMatrices;
	std::vector< VectorType > LocalRightHandSideVectors;

	//Initialize sizes for the system components:
	if( rLHSVariables.size() != LocalLeftHandSideMatrices.size() )
	  LocalLeftHandSideMatrices.resize(rLHSVariables.size());

	if( rRHSVariables.size() != LocalRightHandSideVectors.size() )
	  LocalRightHandSideVectors.resize(rRHSVariables.size());

	LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_LHS_MATRIX);

	//get dofs size of the slave body element
	EquationIdVectorType EquationIds;
	ie->EquationIdVector(EquationIds,rCurrentProcessInfo);
	unsigned int SlaveElementSize = EquationIds.size();

	for( unsigned int i=0; i<LocalLeftHandSideMatrices.size(); i++ )
	  {
	    //Note: LocalRightHandSideVectors.size() > 0
	    this->InitializeSystemMatrices( LocalLeftHandSideMatrices[i], LocalRightHandSideVectors[0], SlaveElementSize, LocalSystem.CalculationFlags );
	  }

	LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR);
	LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_LHS_MATRIX,false);

	for( unsigned int i=0; i<LocalRightHandSideVectors.size(); i++ )
	  {
	    //Note: LocalLeftHandSideMatrices.size() > 0
	    this->InitializeSystemMatrices( LocalLeftHandSideMatrices[0], LocalRightHandSideVectors[i], SlaveElementSize, LocalSystem.CalculationFlags );
	  }
	LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_LHS_MATRIX,true);


	//Set Variables to Local system components
	LocalSystem.SetLeftHandSideMatrices(LocalLeftHandSideMatrices);
	LocalSystem.SetRightHandSideVectors(LocalRightHandSideVectors);

	LocalSystem.SetLeftHandSideVariables(rLHSVariables);
	LocalSystem.SetRightHandSideVariables(rRHSVariables);


	//create linked system components
	LocalSystemComponents LinkedSystem;

	std::vector< MatrixType > LeftHandSideMatrices;
	std::vector< VectorType > RightHandSideVectors;
	ie->CalculateLocalSystem(LeftHandSideMatrices,rLHSVariables,RightHandSideVectors,rRHSVariables,rCurrentProcessInfo);

	LinkedSystem.SetLeftHandSideMatrices(LeftHandSideMatrices);
	LinkedSystem.SetRightHandSideVectors(RightHandSideVectors);

	LinkedSystem.SetLeftHandSideVariables(rLHSVariables);
	LinkedSystem.SetRightHandSideVariables(rRHSVariables);


	//Calculate condition system
	ElementType::Pointer SlaveElement = ElementType::Pointer(*(ie.base()));
	this->CalculateConditionSystem( LocalSystem, LinkedSystem, SlaveElement, rCurrentProcessInfo );

	//assemble the local system into the global system

	//resizing as needed the LHSs
	for( unsigned int i=0; i< rLeftHandSideMatrices.size(); i++)
	  {
	    unsigned int GlobalSize1 = rLeftHandSideMatrices[i].size1();
	    unsigned int LocalSize1  = LocalLeftHandSideMatrices[i].size1();
	    unsigned int MatSize1    = GlobalSize1+LocalSize1;

	    unsigned int GlobalSize2 = rLeftHandSideMatrices[i].size2();
	    unsigned int LocalSize2  = LocalLeftHandSideMatrices[i].size2();
	    unsigned int MatSize2    = GlobalSize2+LocalSize2;

	    rLeftHandSideMatrices[i].resize( MatSize1, MatSize2, true );

	    unsigned int indexi = 0;
	    for(unsigned int i=GlobalSize1; i<MatSize1; i++)
	      {
		unsigned int indexj = 0;
		for(unsigned int j=GlobalSize2; j<MatSize2; j++)
		  {
		    rLeftHandSideMatrices[i](i,j)=LocalLeftHandSideMatrices[i](indexi,indexj);
		    indexj++;
		  }
		indexi++;
	      }

	  }


	//resizing as needed the RHSs
	for( unsigned int i=0; i< rRightHandSideVectors.size(); i++)
	  {

	    unsigned int GlobalSize1 = rRightHandSideVectors[i].size();
	    unsigned int LocalSize1  = LocalRightHandSideVectors[i].size();
	    unsigned int MatSize1    = GlobalSize1+LocalSize1;

	    rRightHandSideVectors[i].resize( MatSize1, true );

	    unsigned int indexi = 0;
	    for(unsigned int j=GlobalSize1; j<MatSize1; j++)
	      {
		rRightHandSideVectors[i][j]=LocalRightHandSideVectors[i][indexi];
		indexi++;
	      }

	  }
      }

}



//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
								ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //set sizes to zero
    rLeftHandSideMatrix.clear();
    rLeftHandSideMatrix.resize(0,0);

    VectorType LocalRightHandSideVector = Vector(0);

    //Ask to the linked deformable element the LocalRightHandSide (only one element link)
    const unsigned int inode = GetGeometry().PointsNumber()-1;
    WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

    for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ie++)
      {

	//create local system components
	LocalSystemComponents LocalSystem;

	//calculation flags
	LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_LHS_MATRIX);

	MatrixType LocalLeftHandSideMatrix;

	//get dofs size of the slave body element
	EquationIdVectorType EquationIds;
	ie->EquationIdVector(EquationIds,rCurrentProcessInfo);
	unsigned int SlaveElementSize = EquationIds.size();

	//Initialize sizes for the system components:
	this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, SlaveElementSize, LocalSystem.CalculationFlags );

	//Set Variables to Local system components
	LocalSystem.SetLeftHandSideMatrix(LocalLeftHandSideMatrix);
	LocalSystem.SetRightHandSideVector(LocalRightHandSideVector);


	//create linked system components
	LocalSystemComponents LinkedSystem;

	MatrixType LeftHandSideMatrix;
	ie->CalculateSecondDerivativesLHS(LeftHandSideMatrix,rCurrentProcessInfo);

	LinkedSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);

	//Calculate condition system
	ElementType::Pointer SlaveElement = ElementType::Pointer(*(ie.base()));
	this->CalculateConditionSystem( LocalSystem, LinkedSystem, SlaveElement, rCurrentProcessInfo );

	//assemble the local system into the global system

	//resizing as needed the LHS
	unsigned int GlobalSize1 = rLeftHandSideMatrix.size1();
	unsigned int LocalSize1  = LocalLeftHandSideMatrix.size1();
	unsigned int MatSize1    = GlobalSize1+LocalSize1;

	unsigned int GlobalSize2 = rLeftHandSideMatrix.size2();
	unsigned int LocalSize2  = LocalLeftHandSideMatrix.size2();
	unsigned int MatSize2    = GlobalSize2+LocalSize2;

	rLeftHandSideMatrix.resize( MatSize1, MatSize2, true );

	unsigned int indexi = 0;
	for(unsigned int i=GlobalSize1; i<MatSize1; i++)
	  {
	    unsigned int indexj = 0;
	    for(unsigned int j=GlobalSize2; j<MatSize2; j++)
	      {
		rLeftHandSideMatrix(i,j) = LocalLeftHandSideMatrix(indexi,indexj);
		indexj++;
	      }
	    indexi++;
	  }


      }


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
								ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //set sizes to zero
    rRightHandSideVector.clear();
    rRightHandSideVector.resize(0);

    MatrixType LocalLeftHandSideMatrix = Matrix();

    //Ask to the linked deformable element the LocalRightHandSide (only one element link)
    const unsigned int inode = GetGeometry().PointsNumber()-1;
    WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

    for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ie++)
      {

	//create local system components
	LocalSystemComponents LocalSystem;

	//calculation flags
	LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR);

	VectorType LocalRightHandSideVector;

	//get dofs size of the slave body element
	EquationIdVectorType EquationIds;
	ie->EquationIdVector(EquationIds,rCurrentProcessInfo);
	unsigned int SlaveElementSize = EquationIds.size();

	//Initialize sizes for the system components:
	this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, SlaveElementSize, LocalSystem.CalculationFlags );
	//Set Variables to Local system components
	LocalSystem.SetRightHandSideVector(LocalRightHandSideVector);

	//create linked system components
	LocalSystemComponents LinkedSystem;

	VectorType RightHandSideVector;
	ie->CalculateSecondDerivativesRHS(RightHandSideVector,rCurrentProcessInfo);

	LinkedSystem.SetRightHandSideVector(RightHandSideVector);

	//Calculate condition system
	ElementType::Pointer SlaveElement = ElementType::Pointer(*(ie.base()));
	this->CalculateConditionSystem( LocalSystem, LinkedSystem, SlaveElement, rCurrentProcessInfo );

	//assemble the local system into the global system

	//resizing as needed the RHS
	unsigned int GlobalSize1 = rRightHandSideVector.size();
	unsigned int LocalSize1  = LocalRightHandSideVector.size();
	unsigned int MatSize1    = GlobalSize1+LocalSize1;

	rRightHandSideVector.resize( MatSize1, true );

	unsigned int indexi = 0;
	for(unsigned int i=GlobalSize1; i<MatSize1; i++)
	  {
	    rRightHandSideVector[i]=LocalRightHandSideVector[indexi];
	    indexi++;
	  }

      }


    KRATOS_CATCH( "" )
}



//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rMassMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rDampingMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::CalculateAndAddStiffness(MatrixType& rLeftHandSideMatrix,
							   MatrixType& rLinkedLeftHandSideMatrix,
							   GeneralVariables& rVariables)

{
    KRATOS_TRY

    //Set variables from the slave linked elements (deformable elements)
    //Get the stiffness tangent matrix from the linearization of the internal forces and the external forces
    //Predict the lagrange multiplier due to the link and add it to the stiffness components

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    const unsigned int size = dimension * 2;

    // std::cout<<"  SLAVE NODE: "<< rVariables.SlaveNode <<std::endl;
    // std::cout<<"  Distance "<< rVariables.Distance <<" Skew "<< rVariables.SkewSymDistance<<std::endl;

    unsigned int start_master = rLeftHandSideMatrix.size1() - (dimension * 2);
    unsigned int start_slave  = rVariables.SlaveNode * (dimension * 2);
    unsigned int end_slave    = start_slave + (dimension * 2);


    if( end_slave == rLinkedLeftHandSideMatrix.size1() )
      end_slave = 0;

    // WriteMatrixInRows( "rLeftHandSideMatrix", rLeftHandSideMatrix );

    // WriteMatrixInRows( "rLinkedLeftHandSideMatrix", rLinkedLeftHandSideMatrix );


    for(unsigned int i=0; i<size; i++)
      {
    	for(unsigned int j=0; j<size; j++)
    	  {
	    //nodal block 11
    	    rLeftHandSideMatrix(start_master+i,start_master+j) += rLinkedLeftHandSideMatrix(start_slave+i,start_slave+j);
    	    //nodal block 21
    	    rLeftHandSideMatrix(end_slave+i,start_master+j)    += rLinkedLeftHandSideMatrix(end_slave+i,start_slave+j);
	    //nodal block 12
    	    rLeftHandSideMatrix(start_master+i,end_slave+j)    += rLinkedLeftHandSideMatrix(start_slave+i,end_slave+j);
   	  }
      }


    //WriteMatrixInRows( "rLeftHandSideMatrix S", rLeftHandSideMatrix );

    Matrix ForceMatrix  = ZeroMatrix(3,3);
    Matrix MomentMatrix = ZeroMatrix(3,3);



    //WriteMatrixInRows( "rLeftHandSideMatrix S", rLeftHandSideMatrix );

    //******************************//

    //add column matrices due to the KINEMATIC LINK

    //column matrices related to the rotation constraint P (U1U2)
    ForceMatrix = ZeroMatrix(3,3);

    for(unsigned int i=0; i<dimension; i++)
      {
    	for(unsigned int j=0; j<dimension; j++)
    	  {
    	    ForceMatrix(i,j)= rLinkedLeftHandSideMatrix(end_slave+i,start_slave+j);
    	  }
      }

    MomentMatrix = prod(ForceMatrix,rVariables.SkewSymDistance);

    for(unsigned int i=0; i<dimension; i++)
      {
    	for(unsigned int j=0; j<dimension; j++)
    	  {
    	    rLeftHandSideMatrix(end_slave+i,start_master+dimension+j) -= MomentMatrix(i,j);
    	  }
      }


    //column matrices related to the rotation constraint H (U1O2)
    ForceMatrix = ZeroMatrix(3,3);

    for(unsigned int i=0; i<dimension; i++)
      {
    	for(unsigned int j=0; j<dimension; j++)
    	  {
    	    ForceMatrix(i,j)= rLinkedLeftHandSideMatrix(end_slave+dimension+i,start_slave+j);
    	  }
      }

    MomentMatrix = prod(ForceMatrix,rVariables.SkewSymDistance);

    for(unsigned int i=0; i<dimension; i++)
      {
    	for(unsigned int j=0; j<dimension; j++)
    	  {
    	    rLeftHandSideMatrix(end_slave+dimension+i,start_master+dimension+j) -= MomentMatrix(i,j);
    	  }
      }



    //******************************//

    //add matrices in the rigid body positions:

    //column matrices related to the rotation constraint P (U1U1)
    ForceMatrix = ZeroMatrix(3,3);

    for(unsigned int i=0; i<dimension; i++)
      {
    	for(unsigned int j=0; j<dimension; j++)
    	  {
    	    ForceMatrix(i,j)= rLinkedLeftHandSideMatrix(start_slave+i,start_slave+j);
    	  }
      }

    MomentMatrix = prod(ForceMatrix,rVariables.SkewSymDistance);

    for(unsigned int i=0; i<dimension; i++)
      {
    	for(unsigned int j=0; j<dimension; j++)
    	  {
    	    rLeftHandSideMatrix(start_master+i,start_master+dimension+j) -= MomentMatrix(i,j);
    	  }
      }


    //column matrices related to the rotation constraint H (U1O1)
    ForceMatrix = ZeroMatrix(3,3);

    for(unsigned int i=0; i<dimension; i++)
      {
    	for(unsigned int j=0; j<dimension; j++)
    	  {
    	    ForceMatrix(i,j)= rLinkedLeftHandSideMatrix(start_slave+dimension+i,start_slave+j);
    	  }
      }

    MomentMatrix = prod(ForceMatrix,rVariables.SkewSymDistance);

    for(unsigned int i=0; i<dimension; i++)
      {
    	for(unsigned int j=0; j<dimension; j++)
    	  {
    	    rLeftHandSideMatrix(start_master+dimension+i,start_master+dimension+j) -= MomentMatrix(i,j);
    	  }
      }


   //******************************//

    // //add row matrices due to the FORCE-MOMENTUM LINK

    // //row matrices related to the rotation constraint P (U1U2)
    // for(unsigned int i=0; i<dimension; i++)
    //   {
    // 	for(unsigned int j=0; j<dimension; j++)
    // 	  {
    // 	    ForceMatrix(i,j)= rLinkedLeftHandSideMatrix(start_slave+i,end_slave+j);
    // 	  }
    //   }

    // MomentMatrix = prod(rVariables.SkewSymDistance,ForceMatrix);

    // for(unsigned int i=0; i<dimension; i++)
    //   {
    // 	for(unsigned int j=0; j<dimension; j++)
    // 	  {
    // 	    rLeftHandSideMatrix(start_master+dimension+i,end_slave+j) += MomentMatrix(i,j);
    // 	  }
    //   }


    // //row matrices related to the rotation constraint H (U1O2)
    // ForceMatrix = ZeroMatrix(3,3);

    // for(unsigned int i=0; i<dimension; i++)
    //   {
    // 	for(unsigned int j=0; j<dimension; j++)
    // 	  {
    // 	    ForceMatrix(i,j)= rLinkedLeftHandSideMatrix(start_slave+i,end_slave+dimension+j);
    // 	  }
    //   }

    // MomentMatrix = prod(rVariables.SkewSymDistance,ForceMatrix);

    // for(unsigned int i=0; i<dimension; i++)
    //   {
    // 	for(unsigned int j=0; j<dimension; j++)
    // 	  {
    // 	    rLeftHandSideMatrix(start_master+dimension+i,end_slave+dimension+j) += MomentMatrix(i,j);
    // 	  }
    //   }



    // //******************************//

    // //add matrices in the rigid body positions:

    // //row matrices related to the rotation constraint P (U1U1)
    // ForceMatrix = ZeroMatrix(3,3);

    // for(unsigned int i=0; i<dimension; i++)
    //   {
    // 	for(unsigned int j=0; j<dimension; j++)
    // 	  {
    // 	    ForceMatrix(i,j)= rLinkedLeftHandSideMatrix(start_slave+i,start_slave+j);
    // 	  }
    //   }

    // MomentMatrix = prod(rVariables.SkewSymDistance,ForceMatrix);

    // for(unsigned int i=0; i<dimension; i++)
    //   {
    // 	for(unsigned int j=0; j<dimension; j++)
    // 	  {
    // 	    rLeftHandSideMatrix(start_master+dimension+i,start_master+j) += MomentMatrix(i,j);
    // 	  }
    //   }


    // //row matrices related to the rotation constraint H (U1O1)
    // ForceMatrix = ZeroMatrix(3,3);

    // for(unsigned int i=0; i<dimension; i++)
    //   {
    // 	for(unsigned int j=0; j<dimension; j++)
    // 	  {
    // 	    ForceMatrix(i,j)= rLinkedLeftHandSideMatrix(start_slave+i,start_slave+dimension+j);
    // 	  }
    //   }

    // MomentMatrix = prod(rVariables.SkewSymDistance,ForceMatrix);

    // for(unsigned int i=0; i<dimension; i++)
    //   {
    // 	for(unsigned int j=0; j<dimension; j++)
    // 	  {
    // 	    rLeftHandSideMatrix(start_master+dimension+i,start_master+dimension+j) += MomentMatrix(i,j);
    // 	  }
    //   }



    // //******************************//


    // //extra row matrices related to the rotation constraint P
    // ForceMatrix = ZeroMatrix(3,3);

    // for(unsigned int i=0; i<dimension; i++)
    //   {
    // 	for(unsigned int j=0; j<dimension; j++)
    // 	  {
    // 	    ForceMatrix(i,j)= rLinkedLeftHandSideMatrix(start_slave+i,start_slave+j);
    // 	  }
    //   }

    // MomentMatrix = prod(rVariables.SkewSymDistance,ForceMatrix);
    // MomentMatrix = prod(MomentMatrix, trans(rVariables.SkewSymDistance));

    // for(unsigned int i=0; i<dimension; i++)
    //   {
    // 	for(unsigned int j=0; j<dimension; j++)
    // 	  {
    // 	    rLeftHandSideMatrix(start_master+dimension+i,start_master+dimension+j) -= MomentMatrix(i,j);
    // 	  }
    //   }


    // WriteMatrixInRows( "rLeftHandSideMatrix S+R", rLeftHandSideMatrix );


    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::CalculateAndAddForces(VectorType& rRightHandSideVector,
							VectorType& rLinkedRightHandSideVector,
							GeneralVariables& rVariables)

{
    KRATOS_TRY

    //Set variables from the slave linked elements (deformable elements)
    //Get the internal forces and the external forces, compute the residual
    //Predict the lagrange multiplier due to the link and add it to the forces

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int start_master = rRightHandSideVector.size() - (dimension * 2);
    unsigned int start_slave  = rVariables.SlaveNode * (dimension * 2);

    const unsigned int size = dimension * 2;

    //rLinkedRightHandSideVector.clear();

    for(unsigned int i=0; i<size; i++)
      {
    	rRightHandSideVector[start_master+i] += rLinkedRightHandSideVector[start_slave+i];
      }


    // Account for the moment of the shear force
    // Vector Force = ZeroVector(3);
    // for(unsigned int i=0; i<dimension; i++)
    //   {
    //  	Force[i] = rLinkedRightHandSideVector[start_slave+i];
    //   }

    // Vector Moment = prod(rVariables.SkewSymDistance,Force);  //this is equal to (D x F)

    // for(unsigned int i=0; i<dimension; i++)
    //   {
    // 	rRightHandSideVector[start_master+dimension+i] += Moment[i];
    //   }


    //std::cout<<" [LINK ID:"<<this->Id()<<"]  RHS Vector "<<rRightHandSideVector<<std::endl;

    //std::cout<<" [ID:"<<this->Id()<<"][iter:"<<mIteration<<"][Link Force "<<Force<<"][Link Moment "<<Moment<<"][NodeId "<<GetGeometry()[0].Id()<<"]"<<std::endl;
    //std::cout<<" RHS Vector "<<rRightHandSideVector<<std::endl;


    //std::cout<<" LINK: "<<this->Id()<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;



    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::VectorToSkewSymmetricTensor( const Vector& rVector,
							       Matrix& rSkewSymmetricTensor )
{
    KRATOS_TRY

    //Initialize Local Matrices
    if( rSkewSymmetricTensor.size1() != 3 )
      rSkewSymmetricTensor.resize(3, 3, false);

    rSkewSymmetricTensor = ZeroMatrix(3,3);

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

void RigidBodyPointLinkCondition::WriteMatrixInRows( std::string MatrixName, const MatrixType& rMatrix )
{

  std::cout<<MatrixName<<" ["<<rMatrix.size1()<<","<<rMatrix.size2()<<"]:"<<std::endl;

  for(unsigned int i=0; i<rMatrix.size1(); i++)
    {
      std::cout<<"(";
      for(unsigned int j=0; j<rMatrix.size2()-1; j++)
	{
	  std::cout<<rMatrix(i,j)<<", ";
	}

      std::cout<<rMatrix(i,rMatrix.size2()-1)<<")"<<std::endl;
    }
}

//***********************************************************************************
//***********************************************************************************


int RigidBodyPointLinkCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
  return 0;
}

//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
