//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:             September 2018 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/rigid_body_links/rigid_body_point_link_segregated_V_condition.hpp"

#include "contact_mechanics_application_variables.h"

namespace Kratos
{

/**
 * Flags related to the condition computation
 */
KRATOS_CREATE_LOCAL_FLAG(RigidBodyPointLinkSegregatedVCondition, COMPUTE_RHS_VECTOR, 0);
KRATOS_CREATE_LOCAL_FLAG(RigidBodyPointLinkSegregatedVCondition, COMPUTE_LHS_MATRIX, 1);


//***********************************************************************************
//***********************************************************************************
RigidBodyPointLinkSegregatedVCondition::RigidBodyPointLinkSegregatedVCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
  for ( SizeType i = 0; i < GetGeometry().size(); i++ )
  {
    GetGeometry()[i].Set(SLAVE); //Flag to set MASTER_ELEMENTS in that nodes (if is SLAVE, a MASTER is required)
  }
}

//***********************************************************************************
//***********************************************************************************
RigidBodyPointLinkSegregatedVCondition::RigidBodyPointLinkSegregatedVCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
  for ( SizeType i = 0; i < GetGeometry().size(); i++ )
  {
    GetGeometry()[i].Set(SLAVE); //Flag to set MASTER_ELEMENTS in that nodes (if is SLAVE, a MASTER is required)
  }
}

//************************************************************************************
//************************************************************************************
RigidBodyPointLinkSegregatedVCondition::RigidBodyPointLinkSegregatedVCondition( RigidBodyPointLinkSegregatedVCondition const& rOther )
    : Condition(rOther)
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer RigidBodyPointLinkSegregatedVCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
  return Kratos::make_shared<RigidBodyPointLinkSegregatedVCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}


//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer RigidBodyPointLinkSegregatedVCondition::Clone(IndexType NewId, NodesArrayType const& rThisNodes) const
{
  return this->Create( NewId, rThisNodes, pGetProperties() );
}

//***********************************************************************************
//***********************************************************************************
RigidBodyPointLinkSegregatedVCondition::~RigidBodyPointLinkSegregatedVCondition()
{
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  rConditionDofList.resize(0);

  switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
  {
    case VELOCITY_STEP:
      {
        RigidBodyPointLinkCondition::GetDofList(rConditionDofList, rCurrentProcessInfo);
        break;
      }
    case PRESSURE_STEP:
      {
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }

  KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  rResult.resize( 0, false );

  switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
  {
    case VELOCITY_STEP:
      {
        RigidBodyPointLinkCondition::EquationIdVector(rResult, rCurrentProcessInfo);
        break;
      }
    case PRESSURE_STEP:
      {
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }

  KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::GetValuesVector(Vector& rValues, int Step)
{
  KRATOS_TRY

  rValues.resize(0);

  const SizeType inode = GetGeometry().PointsNumber()-1;

  WeakPointerVector<Element>& SlaveElements  = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

  SizeType indexi = 0;
  SizeType sizei  = 0;
  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    Vector SlaveValues;
    ie->GetValuesVector(SlaveValues,Step);

    sizei += SlaveValues.size();
    rValues.resize(sizei,true);

    for(SizeType i=0; i<SlaveValues.size(); i++)
      rValues[indexi+i] = SlaveValues[i];

    indexi += SlaveValues.size();
  }

  Element& MasterElement = (GetGeometry()[inode].GetValue(MASTER_ELEMENTS)).back();

  Vector MasterValues;
  MasterElement.GetValuesVector(MasterValues, Step);

  sizei += MasterValues.size();

  for(SizeType i=0; i<MasterValues.size(); i++)
    rValues[indexi+i] = MasterValues[i];

  KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
{
  KRATOS_TRY

  rValues.resize(0);

  const SizeType inode = GetGeometry().PointsNumber()-1;

  WeakPointerVector<Element>& SlaveElements  = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

  SizeType indexi = 0;
  SizeType sizei  = 0;
  for (WeakPointerVector<Element>::iterator ie = SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    Vector SlaveValues;
    ie->GetFirstDerivativesVector(SlaveValues,Step);

    sizei += SlaveValues.size();
    rValues.resize(sizei,true);

    for(SizeType i=0; i<SlaveValues.size(); i++)
      rValues[indexi+i] = SlaveValues[i];

    indexi += SlaveValues.size();
  }

  Element& MasterElement = (GetGeometry()[inode].GetValue(MASTER_ELEMENTS)).back();

  Vector MasterValues;
  MasterElement.GetFirstDerivativesVector(MasterValues, Step);

  sizei += MasterValues.size();

  for(SizeType i=0; i<MasterValues.size(); i++)
    rValues[indexi+i] = MasterValues[i];


  KRATOS_CATCH("")
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
{
  KRATOS_TRY

  rValues.resize(0);

  const SizeType inode = GetGeometry().PointsNumber()-1;

  WeakPointerVector<Element>& SlaveElements  = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

  SizeType indexi = 0;
  SizeType sizei  = 0;
  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    Vector SlaveValues;
    ie->GetSecondDerivativesVector(SlaveValues,Step);

    sizei += SlaveValues.size();
    rValues.resize(sizei,true);

    for(SizeType i=0; i<SlaveValues.size(); i++)
      rValues[indexi+i] = SlaveValues[i];

    indexi += SlaveValues.size();
  }

  Element& MasterElement = (GetGeometry()[inode].GetValue(MASTER_ELEMENTS)).back();

  Vector MasterValues;
  MasterElement.GetSecondDerivativesVector(MasterValues, Step);

  sizei += MasterValues.size();

  for(SizeType i=0; i<MasterValues.size(); i++)
    rValues[indexi+i] = MasterValues[i];

  KRATOS_CATCH("")
}



//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::AddExplicitContribution(const VectorType& rRHSVector,
                                                                     const Variable<VectorType>& rRHSVariable,
                                                                     Variable<array_1d<double,3> >& rDestinationVariable,
                                                                     const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  const SizeType number_of_nodes = GetGeometry().PointsNumber();
  const SizeType dimension       = GetGeometry().WorkingSpaceDimension();

  if( rRHSVariable == CONTACT_FORCES_VECTOR && rDestinationVariable == CONTACT_FORCE )
  {

    for(SizeType i=0; i< number_of_nodes; i++)
    {
      int index = i * (dimension * (dimension-1));

      GetGeometry()[i].SetLock();

      array_1d<double, 3 > &ContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
      for(SizeType j=0; j<dimension; j++)
      {
        ContactForce[j] += rRHSVector[index + j];
      }

      GetGeometry()[i].UnSetLock();
    }
  }


  if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
  {

    for(SizeType i=0; i< number_of_nodes; i++)
    {
      int index = i * (dimension * (dimension-1));

      GetGeometry()[i].SetLock();

      array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
      for(SizeType j=0; j<dimension; j++)
      {
        ForceResidual[j] += rRHSVector[index + j];
      }

      GetGeometry()[i].UnSetLock();
    }
  }

  KRATOS_CATCH("")
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::Initialize()
{
  KRATOS_TRY

  RigidBodyPointLinkCondition::Initialize();

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  this->SetProcessInformation(rCurrentProcessInfo);

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        RigidBodyPointLinkCondition::InitializeSolutionStep(rCurrentProcessInfo);
        break;
      }
    case PRESSURE_STEP:
      {
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************
void RigidBodyPointLinkSegregatedVCondition::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
  KRATOS_TRY

  // const SizeType inode = GetGeometry().PointsNumber()-1;

  // PointPointerType  mpSlaveNode = GetGeometry()(inode);

  // array_1d<double, 3 >&  Displacement = mpSlaveNode->FastGetSolutionStepValue(DISPLACEMENT);
  // array_1d<double, 3 >&  Rotation     = mpSlaveNode->FastGetSolutionStepValue(ROTATION);

  // std::cout<<" Link Initialize RB iteration [Rotation:"<<Rotation<<"; Disp:"<<Displacement<<"]"<<std::endl;

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
  KRATOS_TRY

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
  KRATOS_TRY

  KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
                                                           VectorType& rRightHandSideVector,
                                                           const SizeType& rSlaveElementSize,
                                                           Flags& rCalculationFlags)

{
  //resizing as needed the LHS
  SizeType MatSize = rSlaveElementSize + this->GetDofsSize(); //slave elements + rigid-body element

  if ( rCalculationFlags.Is(RigidBodyPointLinkSegregatedVCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
  {
    if ( rLeftHandSideMatrix.size1() != MatSize )
      rLeftHandSideMatrix.resize( MatSize, MatSize, false );

    noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize); //resetting LHS
  }

  //resizing as needed the RHS
  if ( rCalculationFlags.Is(RigidBodyPointLinkSegregatedVCondition::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
  {
    if ( rRightHandSideVector.size() != MatSize )
      rRightHandSideVector.resize( MatSize, false );

    rRightHandSideVector = ZeroVector(MatSize); //resetting RHS
  }
}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::InitializeGeneralVariables(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  const SizeType inode = GetGeometry().PointsNumber()-1;

  // set the slave node
  GeometryType& SlaveGeometry = rVariables.pSlaveElement->GetGeometry();

  for ( SizeType i = 0; i < SlaveGeometry.size(); i++ )
  {
    if( SlaveGeometry[i].Id() == GetGeometry()[inode].Id() )
      rVariables.SlaveNode = i;
  }

  rVariables.Distance.resize(3);
  noalias(rVariables.Distance) = ZeroVector(3);
  rVariables.SkewSymDistance.resize(3,3,false);
  noalias(rVariables.SkewSymDistance) = ZeroMatrix(3,3);
  rVariables.LagrangeMultipliers.resize(3);
  noalias(rVariables.LagrangeMultipliers) = ZeroVector(3);

  //compute distance from the slave node to the master rigid body element node (center of gravity)
  Element& MasterElement = (GetGeometry()[inode].GetValue(MASTER_ELEMENTS)).back();

  //rVariables.Distance = MasterElement.GetGeometry()[0].Coordinates() - GetGeometry()[inode].Coordinates();

  rVariables.Distance = GetGeometry()[inode].Coordinates() - MasterElement.GetGeometry()[0].Coordinates();
  //std::cout<<" Distance "<<norm_2(rVariables.Distance)<<std::endl;

  //compute the skewsymmmetric tensor of the distance
  this->VectorToSkewSymmetricTensor(rVariables.Distance, rVariables.SkewSymDistance);

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::CalculateConditionSystem(LocalSystemComponents& rLocalSystem,
							   LocalSystemComponents& rLinkedSystem,
							   ElementType::Pointer& pSlaveElement,
							   ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //create and initialize condition variables:
  GeneralVariables Variables;
  Variables.pSlaveElement = pSlaveElement;

  this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

  if ( rLocalSystem.CalculationFlags.Is(RigidBodyPointLinkSegregatedVCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
  {
    //contributions to stiffness matrix calculated on the reference config
    this->CalculateAndAddLHS ( rLocalSystem, rLinkedSystem, Variables );
  }

  if ( rLocalSystem.CalculationFlags.Is(RigidBodyPointLinkSegregatedVCondition::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
  {
    //contribution to external forces
    this->CalculateAndAddRHS ( rLocalSystem, rLinkedSystem, Variables );
  }

  KRATOS_CATCH("")
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, LocalSystemComponents& rLinkedSystem, GeneralVariables& rVariables)
{

  MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();
  MatrixType& rLinkedLeftHandSideMatrix = rLinkedSystem.GetLeftHandSideMatrix();

  // operation performed: add Kg to the rLefsHandSideMatrix
  this->CalculateAndAddStiffness( rLeftHandSideMatrix, rLinkedLeftHandSideMatrix, rVariables );

  //KRATOS_WATCH( rLeftHandSideMatrix )

}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, LocalSystemComponents& rLinkedSystem, GeneralVariables& rVariables)
{
  //contribution of the internal and external forces
  VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
  VectorType& rLinkedRightHandSideVector = rLinkedSystem.GetRightHandSideVector();

  // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
  this->CalculateAndAddForces( rRightHandSideVector, rLinkedRightHandSideVector, rVariables);

  //KRATOS_WATCH( rRightHandSideVector )
}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
							VectorType& rRightHandSideVector,
							ProcessInfo& rCurrentProcessInfo )
{
  //set sizes to zero
  rLeftHandSideMatrix.clear();
  rLeftHandSideMatrix.resize(0,0);

  rRightHandSideVector.clear();
  rRightHandSideVector.resize(0);

  //Ask to the linked deformable element the LocalRightHandSide (only one element link)
  const SizeType inode = GetGeometry().PointsNumber()-1;
  WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

  //std::cout<<" [ID:"<<this->Id()<<"] [SlaveElements: "<<SlaveElements.size()<<"] [NodeId "<<GetGeometry()[0].Id()<<"]"<<std::endl;

  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyPointLinkSegregatedVCondition::COMPUTE_LHS_MATRIX);
    LocalSystem.CalculationFlags.Set(RigidBodyPointLinkSegregatedVCondition::COMPUTE_RHS_VECTOR);

    MatrixType LocalLeftHandSideMatrix;
    VectorType LocalRightHandSideVector;

    //get dofs size of the slave body element
    EquationIdVectorType EquationIds;
    ie->EquationIdVector(EquationIds,rCurrentProcessInfo);
    SizeType SlaveElementSize = EquationIds.size();

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
    SizeType GlobalSize1 = rLeftHandSideMatrix.size1();
    SizeType LocalSize1  = LocalLeftHandSideMatrix.size1();
    SizeType MatSize1    = GlobalSize1+LocalSize1;

    SizeType GlobalSize2 = rLeftHandSideMatrix.size2();
    SizeType LocalSize2  = LocalLeftHandSideMatrix.size2();
    SizeType MatSize2    = GlobalSize2+LocalSize2;

    rLeftHandSideMatrix.resize( MatSize1, MatSize2, true );

    SizeType indexi = 0;
    for(SizeType i=GlobalSize1; i<MatSize1; i++)
    {
      SizeType indexj = 0;
      for(SizeType j=GlobalSize2; j<MatSize2; j++)
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
    for(SizeType i=GlobalSize1; i<MatSize1; i++)
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

void RigidBodyPointLinkSegregatedVCondition::CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
									  VectorType& rRightHandSideVector,
									  ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //set sizes to zero
  rLeftHandSideMatrix.resize(0,0);
  rRightHandSideVector.resize(0);


  //Ask to the linked deformable element the LocalRightHandSide (only one element link)
  const SizeType inode = GetGeometry().PointsNumber()-1;
  WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyPointLinkSegregatedVCondition::COMPUTE_LHS_MATRIX);
    LocalSystem.CalculationFlags.Set(RigidBodyPointLinkSegregatedVCondition::COMPUTE_RHS_VECTOR);

    MatrixType LocalLeftHandSideMatrix;
    VectorType LocalRightHandSideVector;

    //get dofs size of the slave body element
    EquationIdVectorType EquationIds;
    //std::cout<<"  [ Node :"<<GetGeometry()[inode].Id() <<" Linked to Elements :"<<ie->Id()<<" ] "<<std::endl;
    ie->EquationIdVector(EquationIds,rCurrentProcessInfo);
    SizeType SlaveElementSize = EquationIds.size();

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
    SizeType GlobalSize1 = rLeftHandSideMatrix.size1();
    SizeType LocalSize1  = LocalLeftHandSideMatrix.size1();
    SizeType MatSize1    = GlobalSize1+LocalSize1;

    SizeType GlobalSize2 = rLeftHandSideMatrix.size2();
    SizeType LocalSize2  = LocalLeftHandSideMatrix.size2();
    SizeType MatSize2    = GlobalSize2+LocalSize2;

    rLeftHandSideMatrix.resize( MatSize1, MatSize2, true );

    SizeType indexi = 0;
    for(SizeType i=GlobalSize1; i<MatSize1; i++)
    {
      SizeType indexj = 0;
      for(SizeType j=GlobalSize2; j<MatSize2; j++)
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
    for(SizeType i=GlobalSize1; i<MatSize1; i++)
    {
      rRightHandSideVector[i]=LocalRightHandSideVector[indexi];
      indexi++;
    }

  }

  //std::cout<<" LINK Dynamic RighHandSide "<<rRightHandSideVector<<std::endl;

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
							  ProcessInfo& rCurrentProcessInfo )
{
  //set sizes to zero
  rRightHandSideVector.clear();
  rRightHandSideVector.resize(0);

  MatrixType LocalLeftHandSideMatrix = Matrix();

  //Ask to the linked deformable element the LocalRightHandSide (only one element link)
  const SizeType inode = GetGeometry().PointsNumber()-1;
  WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyPointLinkSegregatedVCondition::COMPUTE_RHS_VECTOR);

    VectorType LocalRightHandSideVector;

    //get dofs size of the slave body element
    EquationIdVectorType EquationIds;
    ie->EquationIdVector(EquationIds,rCurrentProcessInfo);
    SizeType SlaveElementSize = EquationIds.size();

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
    SizeType GlobalSize1 = rRightHandSideVector.size();
    SizeType LocalSize1  = LocalRightHandSideVector.size();
    SizeType MatSize1    = GlobalSize1+LocalSize1;

    rRightHandSideVector.resize( MatSize1, true );

    SizeType indexi = 0;
    for(SizeType i=GlobalSize1; i<MatSize1; i++)
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

void RigidBodyPointLinkSegregatedVCondition::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
								ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //set sizes to zero
  rLeftHandSideMatrix.resize(0,0);

  VectorType LocalRightHandSideVector = Vector(0);

  //Ask to the linked deformable element the LocalRightHandSide (only one element link)
  const SizeType inode = GetGeometry().PointsNumber()-1;
  WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyPointLinkSegregatedVCondition::COMPUTE_LHS_MATRIX);

    MatrixType LocalLeftHandSideMatrix;

    //get dofs size of the slave body element
    EquationIdVectorType EquationIds;
    ie->EquationIdVector(EquationIds,rCurrentProcessInfo);
    SizeType SlaveElementSize = EquationIds.size();

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
    SizeType GlobalSize1 = rLeftHandSideMatrix.size1();
    SizeType LocalSize1  = LocalLeftHandSideMatrix.size1();
    SizeType MatSize1    = GlobalSize1+LocalSize1;

    SizeType GlobalSize2 = rLeftHandSideMatrix.size2();
    SizeType LocalSize2  = LocalLeftHandSideMatrix.size2();
    SizeType MatSize2    = GlobalSize2+LocalSize2;

    rLeftHandSideMatrix.resize( MatSize1, MatSize2, true );

    SizeType indexi = 0;
    for(SizeType i=GlobalSize1; i<MatSize1; i++)
    {
      SizeType indexj = 0;
      for(SizeType j=GlobalSize2; j<MatSize2; j++)
      {
        rLeftHandSideMatrix(i,j) = LocalLeftHandSideMatrix(indexi,indexj);
        indexj++;
      }
      indexi++;
    }

  }

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
								ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //set sizes to zero
  rRightHandSideVector.resize(0);

  MatrixType LocalLeftHandSideMatrix = Matrix();

  //Ask to the linked deformable element the LocalRightHandSide (only one element link)
  const SizeType inode = GetGeometry().PointsNumber()-1;
  WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyPointLinkSegregatedVCondition::COMPUTE_RHS_VECTOR);

    VectorType LocalRightHandSideVector;

    //get dofs size of the slave body element
    EquationIdVectorType EquationIds;
    ie->EquationIdVector(EquationIds,rCurrentProcessInfo);
    SizeType SlaveElementSize = EquationIds.size();

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
    SizeType GlobalSize1 = rRightHandSideVector.size();
    SizeType LocalSize1  = LocalRightHandSideVector.size();
    SizeType MatSize1    = GlobalSize1+LocalSize1;

    rRightHandSideVector.resize( MatSize1, true );

    SizeType indexi = 0;
    for(SizeType i=GlobalSize1; i<MatSize1; i++)
    {
      rRightHandSideVector[i]=LocalRightHandSideVector[indexi];
      indexi++;
    }

  }

  KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  rMassMatrix.resize(0, 0, false);

  KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  rDampingMatrix.resize(0, 0, false);

  KRATOS_CATCH("")
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::CalculateAndAddStiffness(MatrixType& rLeftHandSideMatrix,
							   MatrixType& rLinkedLeftHandSideMatrix,
							   GeneralVariables& rVariables)

{
  KRATOS_TRY

  //Set variables from the slave linked elements (deformable elements)
  //Get the stiffness tangent matrix from the linearization of the internal forces and the external forces
  //Predict the lagrange multiplier due to the link and add it to the stiffness components

  const SizeType dimension = GetGeometry().WorkingSpaceDimension();

  const SizeType dofs_size = this->GetDofsSize();

  // std::cout<<"  SLAVE NODE: "<< rVariables.SlaveNode <<std::endl;
  // std::cout<<"  Distance "<< rVariables.Distance <<" Skew "<< rVariables.SkewSymDistance<<std::endl;

  SizeType start_master = rLeftHandSideMatrix.size1() - dofs_size;
  SizeType start_slave  = rVariables.SlaveNode * dofs_size;
  SizeType end_slave    = start_slave + dofs_size;

  if(end_slave == rLinkedLeftHandSideMatrix.size1())
    end_slave = 0;

  // WriteMatrixInRows( "rLeftHandSideMatrix", rLeftHandSideMatrix );
  // WriteMatrixInRows( "rLinkedLeftHandSideMatrix", rLinkedLeftHandSideMatrix );

  for(SizeType i=0; i<dofs_size; i++)
  {
    for(SizeType j=0; j<dofs_size; j++)
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

  Matrix ForceMatrix(3,3);
  Matrix MomentMatrix(3,3);

  //WriteMatrixInRows( "rLeftHandSideMatrix S", rLeftHandSideMatrix );

  //******************************//

  //add column matrices due to the KINEMATIC LINK

  //column matrices related to the rotation constraint P (U1U2)
  noalias(ForceMatrix) = ZeroMatrix(3,3);

  for(SizeType i=0; i<dimension; i++)
  {
    for(SizeType j=0; j<dimension; j++)
    {
      ForceMatrix(i,j) = rLinkedLeftHandSideMatrix(end_slave+i,start_slave+j);
    }
  }

  noalias(MomentMatrix) = prod(ForceMatrix,rVariables.SkewSymDistance);

  for(SizeType i=0; i<dimension; i++)
  {
    for(SizeType j=0; j<dimension; j++)
    {
      rLeftHandSideMatrix(end_slave+i,start_master+dimension+j) -= MomentMatrix(i,j);
    }
  }


  //column matrices related to the rotation constraint H (U1O2)
  noalias(ForceMatrix) = ZeroMatrix(3,3);

  for(SizeType i=0; i<dimension; i++)
  {
    for(SizeType j=0; j<dimension; j++)
    {
      ForceMatrix(i,j) = rLinkedLeftHandSideMatrix(end_slave+dimension+i,start_slave+j);
    }
  }

  noalias(MomentMatrix) = prod(ForceMatrix,rVariables.SkewSymDistance);

  for(SizeType i=0; i<dimension; i++)
  {
    for(SizeType j=0; j<dimension; j++)
    {
      rLeftHandSideMatrix(end_slave+dimension+i,start_master+dimension+j) -= MomentMatrix(i,j);
    }
  }

  //******************************//

  //add matrices in the rigid body positions:

  //column matrices related to the rotation constraint P (U1U1)
  noalias(ForceMatrix) = ZeroMatrix(3,3);

  for(SizeType i=0; i<dimension; i++)
  {
    for(SizeType j=0; j<dimension; j++)
    {
      ForceMatrix(i,j) = rLinkedLeftHandSideMatrix(start_slave+i,start_slave+j);
    }
  }

  noalias(MomentMatrix) = prod(ForceMatrix,rVariables.SkewSymDistance);

  for(SizeType i=0; i<dimension; i++)
  {
    for(SizeType j=0; j<dimension; j++)
    {
      rLeftHandSideMatrix(start_master+i,start_master+dimension+j) -= MomentMatrix(i,j);
    }
  }

  //column matrices related to the rotation constraint H (U1O1)
  noalias(ForceMatrix) = ZeroMatrix(3,3);

  for(SizeType i=0; i<dimension; i++)
  {
    for(SizeType j=0; j<dimension; j++)
    {
      ForceMatrix(i,j) = rLinkedLeftHandSideMatrix(start_slave+dimension+i,start_slave+j);
    }
  }

  noalias(MomentMatrix) = prod(ForceMatrix,rVariables.SkewSymDistance);

  for(SizeType i=0; i<dimension; i++)
  {
    for(SizeType j=0; j<dimension; j++)
    {
      rLeftHandSideMatrix(start_master+dimension+i,start_master+dimension+j) -= MomentMatrix(i,j);
    }
  }

  //******************************//

  // //add row matrices due to the FORCE-MOMENTUM LINK

  // //row matrices related to the rotation constraint P (U1U2)
  // for(SizeType i=0; i<dimension; i++)
  //   {
  // 	for(SizeType j=0; j<dimension; j++)
  // 	  {
  // 	    ForceMatrix(i,j)= rLinkedLeftHandSideMatrix(start_slave+i,end_slave+j);
  // 	  }
  //   }

  // MomentMatrix = prod(rVariables.SkewSymDistance,ForceMatrix);

  // for(SizeType i=0; i<dimension; i++)
  //   {
  // 	for(SizeType j=0; j<dimension; j++)
  // 	  {
  // 	    rLeftHandSideMatrix(start_master+dimension+i,end_slave+j) += MomentMatrix(i,j);
  // 	  }
  //   }


  // //row matrices related to the rotation constraint H (U1O2)
  // ForceMatrix = ZeroMatrix(3,3);

  // for(SizeType i=0; i<dimension; i++)
  //   {
  // 	for(SizeType j=0; j<dimension; j++)
  // 	  {
  // 	    ForceMatrix(i,j)= rLinkedLeftHandSideMatrix(start_slave+i,end_slave+dimension+j);
  // 	  }
  //   }

  // MomentMatrix = prod(rVariables.SkewSymDistance,ForceMatrix);

  // for(SizeType i=0; i<dimension; i++)
  //   {
  // 	for(SizeType j=0; j<dimension; j++)
  // 	  {
  // 	    rLeftHandSideMatrix(start_master+dimension+i,end_slave+dimension+j) += MomentMatrix(i,j);
  // 	  }
  //   }

  // //******************************//

  // //add matrices in the rigid body positions:

  // //row matrices related to the rotation constraint P (U1U1)
  // ForceMatrix = ZeroMatrix(3,3);

  // for(SizeType i=0; i<dimension; i++)
  //   {
  // 	for(SizeType j=0; j<dimension; j++)
  // 	  {
  // 	    ForceMatrix(i,j)= rLinkedLeftHandSideMatrix(start_slave+i,start_slave+j);
  // 	  }
  //   }

  // MomentMatrix = prod(rVariables.SkewSymDistance,ForceMatrix);

  // for(SizeType i=0; i<dimension; i++)
  //   {
  // 	for(SizeType j=0; j<dimension; j++)
  // 	  {
  // 	    rLeftHandSideMatrix(start_master+dimension+i,start_master+j) += MomentMatrix(i,j);
  // 	  }
  //   }

  // //row matrices related to the rotation constraint H (U1O1)
  // ForceMatrix = ZeroMatrix(3,3);

  // for(SizeType i=0; i<dimension; i++)
  //   {
  // 	for(SizeType j=0; j<dimension; j++)
  // 	  {
  // 	    ForceMatrix(i,j)= rLinkedLeftHandSideMatrix(start_slave+i,start_slave+dimension+j);
  // 	  }
  //   }

  // MomentMatrix = prod(rVariables.SkewSymDistance,ForceMatrix);

  // for(SizeType i=0; i<dimension; i++)
  //   {
  // 	for(SizeType j=0; j<dimension; j++)
  // 	  {
  // 	    rLeftHandSideMatrix(start_master+dimension+i,start_master+dimension+j) += MomentMatrix(i,j);
  // 	  }
  //   }

  // //******************************//

  // //extra row matrices related to the rotation constraint P
  // ForceMatrix = ZeroMatrix(3,3);

  // for(SizeType i=0; i<dimension; i++)
  //   {
  // 	for(SizeType j=0; j<dimension; j++)
  // 	  {
  // 	    ForceMatrix(i,j)= rLinkedLeftHandSideMatrix(start_slave+i,start_slave+j);
  // 	  }
  //   }

  // MomentMatrix = prod(rVariables.SkewSymDistance,ForceMatrix);
  // MomentMatrix = prod(MomentMatrix, trans(rVariables.SkewSymDistance));

  // for(SizeType i=0; i<dimension; i++)
  //   {
  // 	for(SizeType j=0; j<dimension; j++)
  // 	  {
  // 	    rLeftHandSideMatrix(start_master+dimension+i,start_master+dimension+j) -= MomentMatrix(i,j);
  // 	  }
  //   }

  // WriteMatrixInRows( "rLeftHandSideMatrix S+R", rLeftHandSideMatrix );

  KRATOS_CATCH("")
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkSegregatedVCondition::CalculateAndAddForces(VectorType& rRightHandSideVector,
							VectorType& rLinkedRightHandSideVector,
							GeneralVariables& rVariables)

{
  KRATOS_TRY

  //Set variables from the slave linked elements (deformable elements)
  //Get the internal forces and the external forces, compute the residual
  //Predict the lagrange multiplier due to the link and add it to the forces

  const SizeType dofs_size = this->GetDofsSize();

  SizeType start_master = rRightHandSideVector.size() - dofs_size;
  SizeType start_slave  = rVariables.SlaveNode * dofs_size;


  //rLinkedRightHandSideVector.clear();

  for(SizeType i=0; i<dofs_size; i++)
  {
    rRightHandSideVector[start_master+i] += rLinkedRightHandSideVector[start_slave+i];
  }

  // Account for the moment of the shear force
  // Vector Force = ZeroVector(3);
  // for(SizeType i=0; i<dimension; i++)
  //   {
  //  	Force[i] = rLinkedRightHandSideVector[start_slave+i];
  //   }

  // Vector Moment = prod(rVariables.SkewSymDistance,Force);  //this is equal to (D x F)

  // for(SizeType i=0; i<dimension; i++)
  //   {
  // 	rRightHandSideVector[start_master+dimension+i] += Moment[i];
  //   }


  //std::cout<<" [LINK ID:"<<this->Id()<<"]  RHS Vector "<<rRightHandSideVector<<std::endl;

  //std::cout<<" [ID:"<<this->Id()<<"] [Link Force "<<Force<<"][Link Moment "<<Moment<<"][NodeId "<<GetGeometry()[0].Id()<<"]"<<std::endl;
  //std::cout<<" RHS Vector "<<rRightHandSideVector<<std::endl;


  //std::cout<<" LINK: "<<this->Id()<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

RigidBodyPointLinkSegregatedVCondition::SizeType RigidBodyPointLinkSegregatedVCondition::GetDofsSize()
{
  KRATOS_TRY

  const SizeType dimension       = GetGeometry().WorkingSpaceDimension();
  const SizeType number_of_nodes = GetGeometry().PointsNumber();

  SizeType size = number_of_nodes * (dimension * (dimension-1)); //usual size for displacement and rotation based elements

  return size;

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::VectorToSkewSymmetricTensor( const Vector& rVector,
							       Matrix& rSkewSymmetricTensor )
{
  KRATOS_TRY

  //Initialize Local Matrices
  if( rSkewSymmetricTensor.size1() != 3 )
    rSkewSymmetricTensor.resize(3, 3, false);

  noalias(rSkewSymmetricTensor) = ZeroMatrix(3,3);

  rSkewSymmetricTensor( 0, 1 ) = -rVector[2];
  rSkewSymmetricTensor( 0, 2 ) =  rVector[1];
  rSkewSymmetricTensor( 1, 2 ) = -rVector[0];

  rSkewSymmetricTensor( 1, 0 ) =  rVector[2];
  rSkewSymmetricTensor( 2, 0 ) = -rVector[1];
  rSkewSymmetricTensor( 2, 1 ) =  rVector[0];

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkSegregatedVCondition::WriteMatrixInRows( std::string MatrixName, const MatrixType& rMatrix )
{
  std::cout<<MatrixName<<" ["<<rMatrix.size1()<<","<<rMatrix.size2()<<"]:"<<std::endl;

  for(SizeType i=0; i<rMatrix.size1(); i++)
  {
    std::cout<<"(";
    for(SizeType j=0; j<rMatrix.size2()-1; j++)
    {
      std::cout<<rMatrix(i,j)<<", ";
    }
    std::cout<<rMatrix(i,rMatrix.size2()-1)<<")"<<std::endl;
  }
}

//***********************************************************************************
//***********************************************************************************

int RigidBodyPointLinkSegregatedVCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
  return 0;
}

} // Namespace Kratos.
