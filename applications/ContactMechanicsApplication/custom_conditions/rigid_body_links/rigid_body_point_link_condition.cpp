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
#include "custom_conditions/rigid_body_links/rigid_body_point_link_condition.hpp"

#include "contact_mechanics_application_variables.h"

namespace Kratos
{

const std::array<const VariableData,6> RigidBodyPointLinkCondition::mLinearDofs = {DISPLACEMENT_X,DISPLACEMENT_Y,DISPLACEMENT_Z,VELOCITY_X,VELOCITY_Y,VELOCITY_Z};
const std::array<const VariableData,3> RigidBodyPointLinkCondition::mAngularDofs = {ROTATION_X,ROTATION_Y,ROTATION_Z};


/**
 * Flags related to the condition computation
 */
KRATOS_CREATE_LOCAL_FLAG(RigidBodyPointLinkCondition, COMPUTE_RHS_VECTOR, 0);
KRATOS_CREATE_LOCAL_FLAG(RigidBodyPointLinkCondition, COMPUTE_LHS_MATRIX, 1);

//***********************************************************************************
//***********************************************************************************
RigidBodyPointLinkCondition::RigidBodyPointLinkCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
}

//***********************************************************************************
//***********************************************************************************
RigidBodyPointLinkCondition::RigidBodyPointLinkCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

//************************************************************************************
//************************************************************************************
RigidBodyPointLinkCondition::RigidBodyPointLinkCondition( RigidBodyPointLinkCondition const& rOther )
    : Condition(rOther)
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer RigidBodyPointLinkCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
  return Kratos::make_shared<RigidBodyPointLinkCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}


//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer RigidBodyPointLinkCondition::Clone(IndexType NewId, NodesArrayType const& rThisNodes) const
{
  return this->Create( NewId, rThisNodes, pGetProperties() );
}

//***********************************************************************************
//***********************************************************************************
RigidBodyPointLinkCondition::~RigidBodyPointLinkCondition()
{
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::GetDofList(DofsVectorType& rConditionDofList,
					     ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  rConditionDofList.resize(0);

  const SizeType inode = GetGeometry().PointsNumber()-1;

  WeakPointerVector<Element>& SlaveElements  = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

  for(WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {

    DofsVectorType SlaveDofList;
    ie->GetDofList(SlaveDofList, rCurrentProcessInfo);

    for(SizeType i=0; i<SlaveDofList.size(); i++)
      rConditionDofList.push_back(SlaveDofList[i]);
  }

  Element& MasterElement = (GetGeometry()[inode].GetValue(MASTER_ELEMENTS)).back();

  DofsVectorType MasterDofList;
  MasterElement.GetDofList(MasterDofList, rCurrentProcessInfo);

  for(SizeType i=0; i<MasterDofList.size(); i++)
    rConditionDofList.push_back(MasterDofList[i]);

  KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::EquationIdVector(EquationIdVectorType& rResult,
						   ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  rResult.resize( 0, false );

  const SizeType inode = GetGeometry().PointsNumber()-1;

  WeakPointerVector<Element>& SlaveElements  = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    EquationIdVectorType SlaveResult;
    ie->EquationIdVector(SlaveResult, rCurrentProcessInfo);

    for(SizeType i=0; i<SlaveResult.size(); i++)
      rResult.push_back(SlaveResult[i]);

  }

  Element& MasterElement = (GetGeometry()[inode].GetValue(MASTER_ELEMENTS)).back();

  EquationIdVectorType MasterResult;
  MasterElement.EquationIdVector(MasterResult, rCurrentProcessInfo);

  for(SizeType i=0; i<MasterResult.size(); i++)
    rResult.push_back(MasterResult[i]);

  KRATOS_CATCH("")
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::GetValuesVector(Vector& rValues, int Step)
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
  rValues.resize(sizei,true);

  for(SizeType i=0; i<MasterValues.size(); i++)
    rValues[indexi+i] = MasterValues[i];

  KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
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
  rValues.resize(sizei,true);

  for(SizeType i=0; i<MasterValues.size(); i++)
    rValues[indexi+i] = MasterValues[i];


  KRATOS_CATCH("")
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
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
  rValues.resize(sizei,true);

  for(SizeType i=0; i<MasterValues.size(); i++)
    rValues[indexi+i] = MasterValues[i];

  KRATOS_CATCH("")
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::Initialize()
{
  KRATOS_TRY

  //Fix linked deformable point dofs
  for ( SizeType i = 0; i < GetGeometry().size(); i++ )
  {
    DofsContainerType& rDofs = GetGeometry()[i].GetDofs();
    for(DofsContainerType::iterator it = rDofs.begin(); it != rDofs.end(); ++it)
    {
      if(it->GetVariable() != PRESSURE) // it must be some way to fix only kinematic dofs.
        it->FixDof();
    }
  }

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************
void RigidBodyPointLinkCondition::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
  KRATOS_TRY

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
{
  KRATOS_TRY

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
  KRATOS_TRY

  KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
                                                           VectorType& rRightHandSideVector,
                                                           const SizeType& rSlaveElementSize,
                                                           Flags& rCalculationFlags)

{
  //resizing as needed the LHS
  SizeType MatSize = rSlaveElementSize + this->GetDofsSize(); //slave elements + rigid-body element

  if ( rCalculationFlags.Is(RigidBodyPointLinkCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
  {
    if ( rLeftHandSideMatrix.size1() != MatSize )
      rLeftHandSideMatrix.resize( MatSize, MatSize, false );

    noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize); //resetting LHS
  }

  //resizing as needed the RHS
  if ( rCalculationFlags.Is(RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
  {
    if ( rRightHandSideVector.size() != MatSize )
      rRightHandSideVector.resize( MatSize, false );

    rRightHandSideVector = ZeroVector(MatSize); //resetting RHS
  }
}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::InitializeGeneralVariables(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  const SizeType inode = GetGeometry().PointsNumber()-1;

  // set the slave node
  GeometryType& SlaveGeometry = rVariables.pSlaveElement->GetGeometry();

  for(SizeType i=0; i<SlaveGeometry.size(); ++i)
  {
    if(SlaveGeometry[i].Id() == GetGeometry()[inode].Id()){
      rVariables.SlaveNode = i;
      //std::cout<<" Slave Node "<<i<<" Id "<<SlaveGeometry[i].Id()<<std::endl;
    }
    else{
      if(SlaveGeometry[i].Is(RIGID)){
        rVariables.RigidNodes.push_back(i);
        //std::cout<<" Rigid Node "<<i<<" Id "<<SlaveGeometry[i].Id()<<std::endl;
      }
      else{
        rVariables.DeformableNodes.push_back(i);
        //std::cout<<" Deformable Node "<<i<<" Id "<<SlaveGeometry[i].Id()<<std::endl;
      }
    }
  }

  DofsVectorType ElementalDofList;
  rVariables.pSlaveElement->GetDofList(ElementalDofList, rCurrentProcessInfo);
  rVariables.SlaveLinearBlockSize = 0;
  rVariables.SlaveAngularBlockSize = 0;
  for(const auto & elem_dof : ElementalDofList)
  {
    for(const auto& dof : mLinearDofs)
    {
      if( elem_dof->GetVariable() == dof )
        ++rVariables.SlaveLinearBlockSize;
    }
    for(const auto& dof : mAngularDofs)
    {
      if( elem_dof->GetVariable() == dof )
        ++rVariables.SlaveAngularBlockSize;
    }
  }

  rVariables.SlaveLinearBlockSize = SizeType(rVariables.SlaveLinearBlockSize/double(SlaveGeometry.size()));
  rVariables.SlaveAngularBlockSize = SizeType(rVariables.SlaveAngularBlockSize/double(SlaveGeometry.size()));

  //compute distance from the slave node to the master rigid body element node (center of gravity)
  Element& MasterElement = (GetGeometry()[inode].GetValue(MASTER_ELEMENTS)).back();
  MasterElement.GetDofList(ElementalDofList, rCurrentProcessInfo);
  rVariables.MasterLinearBlockSize = 0;
  rVariables.MasterAngularBlockSize = 0;
  for(const auto & elem_dof : ElementalDofList)
  {
    for(const auto& dof : mLinearDofs)
    {
      if( elem_dof->GetVariable() == dof )
        ++rVariables.MasterLinearBlockSize;
    }
    for(const auto& dof : mAngularDofs)
    {
      if( elem_dof->GetVariable() == dof )
        ++rVariables.MasterAngularBlockSize;
    }
  }

  if((rVariables.SlaveLinearBlockSize + rVariables.SlaveAngularBlockSize) != 0){
    if(rVariables.MasterLinearBlockSize != rVariables.SlaveLinearBlockSize)
      KRATOS_ERROR<<" Linear Deformable Dofs and Rigid Dofs not coincide "<<rVariables.SlaveLinearBlockSize<<" !+ "<<rVariables.MasterLinearBlockSize<<std::endl;

    if(rVariables.MasterAngularBlockSize < rVariables.SlaveAngularBlockSize)
      KRATOS_ERROR<<" Angular Deformable Dofs and Rigid Dofs not coincide "<<rVariables.SlaveAngularBlockSize<<" !+ "<<rVariables.MasterAngularBlockSize<<std::endl;
  }

  array_1d<double,3> Distance;
  for(SizeType i=0; i<rVariables.RigidNodes.size(); ++i)
  {
    Distance = SlaveGeometry[rVariables.RigidNodes[i]].Coordinates() - MasterElement.GetGeometry()[0].Coordinates();
    BeamMathUtilsType::VectorToSkewSymmetricTensor(Distance, rVariables.SlaveSkewSymDistance);
    rVariables.RigidSkewSymDistances.push_back(rVariables.SlaveSkewSymDistance);
    //std::cout<<" ["<<MasterElement.Id()<<"-"<<SlaveGeometry[rVariables.RigidNodes[i]].Id()<<"] "<<norm_2(Distance)<<std::endl;

  }

  Distance = SlaveGeometry[rVariables.SlaveNode].Coordinates() - MasterElement.GetGeometry()[0].Coordinates();
  BeamMathUtilsType::VectorToSkewSymmetricTensor(Distance, rVariables.SlaveSkewSymDistance);
  //std::cout<<" ["<<MasterElement.Id()<<"-"<<SlaveGeometry[rVariables.SlaveNode].Id()<<"] "<<norm_2(Distance)<<" SLV "<<std::endl;
  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateConditionSystem(LocalSystemComponents& rLocalSystem,
							   LocalSystemComponents& rLinkedSystem,
							   Element::Pointer& pSlaveElement,
							   ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

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

  KRATOS_CATCH("")
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, LocalSystemComponents& rLinkedSystem, GeneralVariables& rVariables)
{
  KRATOS_TRY

  MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();
  MatrixType& rLinkedLeftHandSideMatrix = rLinkedSystem.GetLeftHandSideMatrix();
  // MatrixType BeamLeftHandSideMatrix = rLeftHandSideMatrix;

  // operation performed: add Kg to the rLefsHandSideMatrix
  if(rLinkedLeftHandSideMatrix.size1()!=0)
    this->CalculateAndAddTangent(rLeftHandSideMatrix, rLinkedLeftHandSideMatrix, rVariables);

  if(rLinkedLeftHandSideMatrix.size1()!=0)
    this->CalculateAndAddTangentBeam(BeamLeftHandSideMatrix, rLinkedLeftHandSideMatrix, rVariables);

  // std::cout<<" Link "<<this->Id()<<std::endl;
  // WriteMatrixInRows( "rLinkHandSideMatrix", rLinkedLeftHandSideMatrix );
  // WriteMatrixInRows( "rLeftHandSideMatrix", rLeftHandSideMatrix );
  // WriteMatrixInRows( "BeamLHandSideMatrix", BeamLeftHandSideMatrix );
  // MatrixType ErrorMatrix = rLeftHandSideMatrix-BeamLeftHandSideMatrix;
  // WriteMatrixInRows( "ErrorMatrix", ErrorMatrix );

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, LocalSystemComponents& rLinkedSystem, GeneralVariables& rVariables)
{
  KRATOS_TRY

  //contribution of the internal and external forces
  VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
  VectorType& rLinkedRightHandSideVector = rLinkedSystem.GetRightHandSideVector();
  // VectorType  BeamRightHandSideVector = rRightHandSideVector;

  // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
  if(rLinkedRightHandSideVector.size()!=0)
    this->CalculateAndAddForces(rRightHandSideVector, rLinkedRightHandSideVector, rVariables);

  // if(rLinkedRightHandSideVector.size()!=0)
  //   this->CalculateAndAddForcesBeam(BeamRightHandSideVector, rLinkedRightHandSideVector, rVariables);

  // std::cout<<" Link "<<this->Id()<<std::endl;
  // KRATOS_WATCH( rLinkedRightHandSideVector )
  // KRATOS_WATCH( rRightHandSideVector )
  // KRATOS_WATCH( BeamRightHandSideVector )
  // KRATOS_WATCH( BeamRightHandSideVector-rRightHandSideVector )

  //rRightHandSideVector = BeamRightHandSideVector;

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
							VectorType& rRightHandSideVector,
							ProcessInfo& rCurrentProcessInfo )
{
  //Ask to the linked deformable element the LocalRightHandSide (only one element link)
  const SizeType inode = GetGeometry().PointsNumber()-1;
  WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

  // std::cout<<" LocalSystem [ID:"<<this->Id()<<"] [SlaveElements: "<<SlaveElements.size()<<"] [NodeId "<<GetGeometry()[0].Id()<<"]"<<std::endl;

  EquationIdVectorType SlaveResult;
  std::vector<SizeType> element_dofs;
  SizeType master_index = 0;
  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    ie->EquationIdVector(SlaveResult, rCurrentProcessInfo);
    element_dofs.push_back(SlaveResult.size());
    master_index += element_dofs.back();
  }

  SizeType total_size = master_index + this->GetDofsSize();

  //set sizes and initialize to zero
  rLeftHandSideMatrix.resize(total_size, total_size, false);
  noalias(rLeftHandSideMatrix) = ZeroMatrix(total_size, total_size);
  rRightHandSideVector.resize(total_size);
  noalias(rRightHandSideVector) = ZeroVector(total_size);

  SizeType counter = 0;
  SizeType local_index = 0;
  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    //std::cout<<" ["<<this->Id()<<"] SLAVE ELEMENT "<<ie->Id()<<" nodes "<<ie->GetGeometry().size()<<std::endl;

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_LHS_MATRIX);
    LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR);

    MatrixType LocalLeftHandSideMatrix;
    VectorType LocalRightHandSideVector;

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, element_dofs[counter], LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LocalLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(LocalRightHandSideVector);

    //create linked system components
    LocalSystemComponents LinkedSystem;

    MatrixType SlaveLeftHandSideMatrix;
    VectorType SlaveRightHandSideVector;
    ie->CalculateLocalSystem(SlaveLeftHandSideMatrix,SlaveRightHandSideVector,rCurrentProcessInfo);

    LinkedSystem.SetLeftHandSideMatrix(SlaveLeftHandSideMatrix);
    LinkedSystem.SetRightHandSideVector(SlaveRightHandSideVector);

    //std::cout<<"[LINK]: "<<ie->Id()<<std::endl;
    //Calculate condition system
    Element::Pointer pSlaveElement = Element::Pointer(*(ie.base()));
    this->CalculateConditionSystem( LocalSystem, LinkedSystem, pSlaveElement, rCurrentProcessInfo );

    //assemble the local system into the global system

    //std::cout<<" Assemble SLAVE start "<<ie->Id()<<std::endl;
    //std::cout<<" LHS "<<local_index<<" "<<element_dofs[counter]<<" "<<master_index<<" "<<total_size<<std::endl;

    this->AssembleLocalLHS(rLeftHandSideMatrix, LocalLeftHandSideMatrix, local_index, element_dofs[counter], master_index);

    //std::cout<<" RHS "<<local_index<<" "<<element_dofs[counter]<<" "<<master_index<<" "<<total_size<<std::endl;
    // std::cout<<" Local RighHandSide "<<LocalRightHandSideVector<<std::endl;

    this->AssembleLocalRHS(rRightHandSideVector, LocalRightHandSideVector, local_index, element_dofs[counter], master_index);

    local_index += element_dofs[counter];
    ++counter;

    //std::cout<<" Assemble SLAVE finish "<<ie->Id()<<std::endl;
  }

  //std::cout<<" LINK RighHandSide "<<rRightHandSideVector<<std::endl;
  //std::cout<<" LINK LeftHandSide "<<rLeftHandSideMatrix<<std::endl;
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::AssembleLocalLHS(MatrixType& rLeftHandSideMatrix,
                                                   const MatrixType& rLocalLeftHandSideMatrix,
                                                   const SizeType& local_index,
                                                   const SizeType& dofs_size,
                                                   const SizeType& master_index)
{
  SizeType total_size = rLeftHandSideMatrix.size1();
  SizeType local_size = local_index + dofs_size;
  SizeType indexi = 0;
  for(SizeType i=local_index; i<local_size; i++)
  {
    SizeType indexj = 0;
    for(SizeType j=local_index; j<local_size; j++)
    {
      rLeftHandSideMatrix(i,j) += rLocalLeftHandSideMatrix(indexi,indexj);
      indexj++;
    }

    indexj = dofs_size;
    for(SizeType j=master_index; j<total_size; j++)
    {
      rLeftHandSideMatrix(i,j) += rLocalLeftHandSideMatrix(indexi,indexj);
      rLeftHandSideMatrix(j,i) += rLocalLeftHandSideMatrix(indexj,indexi);
      indexj++;
    }

    indexi++;
  }

  indexi = dofs_size;
  for(SizeType i=master_index; i<total_size; i++)
  {
    SizeType indexj = dofs_size;
    for(SizeType j=master_index; j<total_size; j++)
    {
      rLeftHandSideMatrix(i,j) += rLocalLeftHandSideMatrix(indexi,indexj);
      indexj++;
    }
    indexi++;
  }

}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::AssembleLocalRHS(VectorType& rRightHandSideVector,
                                                   const VectorType& rLocalRightHandSideVector,
                                                   const SizeType& local_index,
                                                   const SizeType& dofs_size,
                                                   const SizeType& master_index)
{
  SizeType total_size = rRightHandSideVector.size();
  SizeType local_size = local_index + dofs_size;

  SizeType indexi = 0;
  for(SizeType i=local_index; i<local_size; i++)
  {
    rRightHandSideVector[i] += rLocalRightHandSideVector[indexi];
    indexi++;
  }
  indexi = dofs_size;
  for(SizeType i=master_index; i<total_size; i++)
  {
    rRightHandSideVector[i] += rLocalRightHandSideVector[indexi];
    indexi++;
  }

}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
									  VectorType& rRightHandSideVector,
									  ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //Ask to the linked deformable element the LocalRightHandSide (only one element link)
  const SizeType inode = GetGeometry().PointsNumber()-1;
  WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

  //std::cout<<" System2ndDerivatives [ID:"<<this->Id()<<"] [SlaveElements: "<<SlaveElements.size()<<"] [NodeId "<<GetGeometry()[0].Id()<<"]"<<std::endl;

  EquationIdVectorType SlaveResult;
  std::vector<SizeType> element_dofs;
  SizeType master_index = 0;
  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    ie->EquationIdVector(SlaveResult, rCurrentProcessInfo);
    element_dofs.push_back(SlaveResult.size());
    master_index += element_dofs.back();
  }

  SizeType total_size = master_index + this->GetDofsSize();

  //set sizes and initialize to zero
  rLeftHandSideMatrix.resize(total_size, total_size, false);
  noalias(rLeftHandSideMatrix) = ZeroMatrix(total_size, total_size);
  rRightHandSideVector.resize(total_size);
  noalias(rRightHandSideVector) = ZeroVector(total_size);

  SizeType counter = 0;
  SizeType local_index = 0;
  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    //std::cout<<" ["<<this->Id()<<"] 2nd SLAVE ELEMENT "<<ie->Id()<<" nodes "<<ie->GetGeometry().size()<<std::endl;

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_LHS_MATRIX);
    LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR);

    MatrixType LocalLeftHandSideMatrix;
    VectorType LocalRightHandSideVector;

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, element_dofs[counter], LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LocalLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(LocalRightHandSideVector);

    //create linked system components
    LocalSystemComponents LinkedSystem;

    MatrixType SlaveLeftHandSideMatrix;
    VectorType SlaveRightHandSideVector;
    //std::cout<<"[DERIVATIVE_LINK]: "<<ie->Id()<<std::endl;
    ie->CalculateSecondDerivativesContributions(SlaveLeftHandSideMatrix,SlaveRightHandSideVector,rCurrentProcessInfo);

    LinkedSystem.SetLeftHandSideMatrix(SlaveLeftHandSideMatrix);
    LinkedSystem.SetRightHandSideVector(SlaveRightHandSideVector);

    //Calculate condition system
    Element::Pointer pSlaveElement = Element::Pointer(*(ie.base()));
    this->CalculateConditionSystem( LocalSystem, LinkedSystem, pSlaveElement, rCurrentProcessInfo );

    //assemble the local system into the global system

    this->AssembleLocalLHS(rLeftHandSideMatrix, LocalLeftHandSideMatrix, local_index, element_dofs[counter], master_index);
    this->AssembleLocalRHS(rRightHandSideVector, LocalRightHandSideVector, local_index, element_dofs[counter], master_index);

    local_index += element_dofs[counter];
    ++counter;
  }

  //std::cout<<" LINK Dynamic RighHandSide "<<rRightHandSideVector<<std::endl;
  //std::cout<<" LINK Dynamic LeftHandSide "<<rLeftHandSideMatrix<<std::endl;

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                         ProcessInfo& rCurrentProcessInfo)
{
  //set sizes to zero
  MatrixType LocalLeftHandSideMatrix = Matrix();

  //Ask to the linked deformable element the LocalRightHandSide (only one element link)
  const SizeType inode = GetGeometry().PointsNumber()-1;
  WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

  EquationIdVectorType SlaveResult;
  std::vector<SizeType> element_dofs;
  SizeType master_index = 0;
  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    ie->EquationIdVector(SlaveResult, rCurrentProcessInfo);
    element_dofs.push_back(SlaveResult.size());
    master_index += element_dofs.back();
  }

  SizeType total_size = master_index + this->GetDofsSize();

  //set sizes and initialize to zero
  rRightHandSideVector.resize(total_size);
  noalias(rRightHandSideVector) = ZeroVector(total_size);

  SizeType counter = 0;
  SizeType local_index = 0;
  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR);

    VectorType LocalRightHandSideVector;

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, element_dofs[counter], LocalSystem.CalculationFlags );
    //Set Variables to Local system components
    LocalSystem.SetRightHandSideVector(LocalRightHandSideVector);

    //create linked system components
    LocalSystemComponents LinkedSystem;

    VectorType SlaveRightHandSideVector;
    ie->CalculateRightHandSide(SlaveRightHandSideVector,rCurrentProcessInfo);

    LinkedSystem.SetRightHandSideVector(SlaveRightHandSideVector);

    //Calculate condition system
    Element::Pointer pSlaveElement = ElementType::Pointer(*(ie.base()));
    this->CalculateConditionSystem( LocalSystem, LinkedSystem, pSlaveElement, rCurrentProcessInfo );

    //assemble the local system into the global system

    this->AssembleLocalRHS(rRightHandSideVector, LocalRightHandSideVector, local_index, element_dofs[counter], master_index);

    local_index += element_dofs[counter];
    ++counter;
  }
  //KRATOS_WATCH( rLeftHandSideMatrix )
  //KRATOS_WATCH( rRightHandSideVector )
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
								ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //set sizes to zero
  VectorType LocalRightHandSideVector = Vector(0);

  //Ask to the linked deformable element the LocalRightHandSide (only one element link)
  const SizeType inode = GetGeometry().PointsNumber()-1;
  WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

  EquationIdVectorType SlaveResult;
  std::vector<SizeType> element_dofs;
  SizeType master_index = 0;
  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    ie->EquationIdVector(SlaveResult, rCurrentProcessInfo);
    element_dofs.push_back(SlaveResult.size());
    master_index += element_dofs.back();
  }

  SizeType total_size = master_index + this->GetDofsSize();

  //set sizes and initialize to zero
  rLeftHandSideMatrix.resize(total_size, total_size, false);
  noalias(rLeftHandSideMatrix) = ZeroMatrix(total_size, total_size);

  SizeType counter = 0;
  SizeType local_index = 0;
  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_LHS_MATRIX);

    MatrixType LocalLeftHandSideMatrix;

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, element_dofs[counter], LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LocalLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(LocalRightHandSideVector);


    //create linked system components
    LocalSystemComponents LinkedSystem;

    MatrixType SlaveLeftHandSideMatrix;
    ie->CalculateSecondDerivativesLHS(SlaveLeftHandSideMatrix,rCurrentProcessInfo);

    LinkedSystem.SetLeftHandSideMatrix(SlaveLeftHandSideMatrix);

    //Calculate condition system
    Element::Pointer pSlaveElement = Element::Pointer(*(ie.base()));
    this->CalculateConditionSystem( LocalSystem, LinkedSystem, pSlaveElement, rCurrentProcessInfo );

    //assemble the local system into the global system

    this->AssembleLocalLHS(rLeftHandSideMatrix, LocalLeftHandSideMatrix, local_index, element_dofs[counter], master_index);

    local_index += element_dofs[counter];
    ++counter;
  }

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
								ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  //set sizes to zero
  MatrixType LocalLeftHandSideMatrix = Matrix();

  //Ask to the linked deformable element the LocalRightHandSide (only one element link)
  const SizeType inode = GetGeometry().PointsNumber()-1;
  WeakPointerVector< Element >& SlaveElements = (GetGeometry()[inode].GetValue(NEIGHBOUR_ELEMENTS));

  EquationIdVectorType SlaveResult;
  std::vector<SizeType> element_dofs;
  SizeType master_index = 0;
  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    ie->EquationIdVector(SlaveResult, rCurrentProcessInfo);
    element_dofs.push_back(SlaveResult.size());
    master_index += element_dofs.back();
  }

  SizeType total_size = master_index + this->GetDofsSize();

  //set sizes and initialize to zero
  rRightHandSideVector.resize(total_size);
  noalias(rRightHandSideVector) = ZeroVector(total_size);

  SizeType counter = 0;
  SizeType local_index = 0;
  for (WeakPointerVector<Element>::iterator ie= SlaveElements.begin(); ie!=SlaveElements.end(); ++ie)
  {
    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(RigidBodyPointLinkCondition::COMPUTE_RHS_VECTOR);

    VectorType LocalRightHandSideVector;

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LocalLeftHandSideMatrix, LocalRightHandSideVector, element_dofs[counter], LocalSystem.CalculationFlags );
    //Set Variables to Local system components
    LocalSystem.SetRightHandSideVector(LocalRightHandSideVector);

    //create linked system components
    LocalSystemComponents LinkedSystem;

    VectorType SlaveRightHandSideVector;
    ie->CalculateSecondDerivativesRHS(SlaveRightHandSideVector,rCurrentProcessInfo);

    LinkedSystem.SetRightHandSideVector(SlaveRightHandSideVector);

    //Calculate condition system
    Element::Pointer pSlaveElement = Element::Pointer(*(ie.base()));
    this->CalculateConditionSystem( LocalSystem, LinkedSystem, pSlaveElement, rCurrentProcessInfo );

    //assemble the local system into the global system

    this->AssembleLocalRHS(rRightHandSideVector, LocalRightHandSideVector, local_index, element_dofs[counter], master_index);

    local_index += element_dofs[counter];
    ++counter;
  }

  KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

      //std::cout<<" MassMatrix "<<this->Id()<<std::endl;

  rMassMatrix.resize(0, 0, false);

  KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

      //std::cout<<" DampingMatrix "<<this->Id()<<std::endl;

  rDampingMatrix.resize(0, 0, false);

  KRATOS_CATCH("")
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::CalculateAndAddTangent(MatrixType& rLeftHandSideMatrix,
                                                         MatrixType& rLinkedLeftHandSideMatrix,
                                                         GeneralVariables& rVariables)

{
  KRATOS_TRY

  //Set variables from the slave linked elements (deformable elements)
  //Get the stiffness tangent matrix from the linearization of the internal forces and the external forces
  //Predict the lagrange multiplier due to the link and add it to the stiffness components

  const SizeType dofs_size = this->GetDofsSize();

  SizeType SlaveBlockSize = rVariables.SlaveLinearBlockSize + rVariables.SlaveAngularBlockSize;

  //std::cout<<" Slave LHS matrix "<< rLinkedLeftHandSideMatrix <<std::endl;

  SizeType start_master = rLeftHandSideMatrix.size1() - dofs_size;
  SizeType start_slave = rVariables.SlaveNode * SlaveBlockSize;
  SizeType start_rotation = 3-rVariables.MasterAngularBlockSize;

  // Set slave counterpart
  Matrix ForceMatrix(3,3);
  Matrix MomentMatrix(3,3);
  Matrix MomentRowMatrix(3,3);
  Matrix MomentColumnMatrix(3,3);

  noalias(ForceMatrix) = ZeroMatrix(3,3);

  // slave rotation dofs
  bool add_row_moments = true;
  if( rVariables.SlaveAngularBlockSize != 0 )
    add_row_moments = false;

  for(SizeType i=0; i<rVariables.SlaveLinearBlockSize; i++)
  {
    for(SizeType j=0; j<rVariables.SlaveLinearBlockSize; j++)
    {
      //nodal block ii
      rLeftHandSideMatrix(start_master+i,start_master+j) += rLinkedLeftHandSideMatrix(start_slave+i,start_slave+j);

      ForceMatrix(i,j) = rLinkedLeftHandSideMatrix(start_slave+i,start_slave+j);
    }
  }

  noalias(MomentRowMatrix) = prod(rVariables.SlaveSkewSymDistance,ForceMatrix);
  noalias(MomentColumnMatrix) = prod(ForceMatrix,rVariables.SlaveSkewSymDistance);
  noalias(MomentMatrix) = prod(MomentRowMatrix,rVariables.SlaveSkewSymDistance);


  for(SizeType i=0; i<rVariables.SlaveAngularBlockSize; i++)
  {
    for(SizeType j=0; j<rVariables.SlaveLinearBlockSize; j++)
    {
      //row blocks
      rLeftHandSideMatrix(start_master+j,start_master+rVariables.MasterLinearBlockSize+i) += rLinkedLeftHandSideMatrix(start_slave+j,start_slave+rVariables.SlaveLinearBlockSize+i);
      rLeftHandSideMatrix(start_master+rVariables.MasterLinearBlockSize+i,start_master+j) += rLinkedLeftHandSideMatrix(start_slave+rVariables.SlaveLinearBlockSize+i,start_slave+j);
      ForceMatrix(i,j) = rLinkedLeftHandSideMatrix(start_slave+rVariables.SlaveLinearBlockSize+i,start_slave+j);
    }
  }
  // override Moment Matrix
  if(!add_row_moments)
    noalias(MomentRowMatrix) = prod(ForceMatrix,rVariables.SlaveSkewSymDistance);


  for(SizeType i=0; i<rVariables.MasterAngularBlockSize; i++)
  {
    for(SizeType j=0; j<rVariables.MasterLinearBlockSize; j++)
    {
      rLeftHandSideMatrix(start_master+j,start_master+rVariables.MasterLinearBlockSize+i) -= MomentColumnMatrix(j,start_rotation+i);
      if(add_row_moments)
        rLeftHandSideMatrix(start_master+rVariables.MasterLinearBlockSize+i,start_master+j) -= MomentRowMatrix(start_rotation+i,j);
    }
  }

  for(SizeType i=0; i<rVariables.MasterAngularBlockSize; i++)
  {
    for(SizeType j=0; j<rVariables.MasterAngularBlockSize; j++)
    {
      rLeftHandSideMatrix(start_master+rVariables.MasterLinearBlockSize+i,start_master+rVariables.MasterLinearBlockSize+j) += rLinkedLeftHandSideMatrix(start_slave+rVariables.SlaveLinearBlockSize+i,start_slave+rVariables.SlaveLinearBlockSize+j);
      if(add_row_moments)
        rLeftHandSideMatrix(start_master+rVariables.MasterLinearBlockSize+i,start_master+rVariables.MasterLinearBlockSize+j) += MomentMatrix(start_rotation+i,start_rotation+j);
      else
        rLeftHandSideMatrix(start_master+rVariables.MasterLinearBlockSize+i,start_master+rVariables.MasterLinearBlockSize+j) -= MomentRowMatrix(start_rotation+i,start_rotation+j);
    }
  }

  // Non rigid nodes of the linked element:
  for(SizeType d=0; d<rVariables.DeformableNodes.size(); ++d)
  {
    SizeType start_deformable = rVariables.DeformableNodes[d] * SlaveBlockSize;

    noalias(ForceMatrix) = ZeroMatrix(3,3);

    for(SizeType i=0; i<rVariables.SlaveLinearBlockSize; i++)
    {
      for(SizeType j=0; j<rVariables.SlaveLinearBlockSize; j++)
      {
        //nodal block ij
        rLeftHandSideMatrix(start_master+i,start_deformable+j) += rLinkedLeftHandSideMatrix(start_slave+i,start_deformable+j);
        rLeftHandSideMatrix(start_deformable+j,start_master+i) += rLinkedLeftHandSideMatrix(start_deformable+j,start_slave+i);
        ForceMatrix(i,j) = rLinkedLeftHandSideMatrix(start_slave+i,start_deformable+j);
      }
    }


    for(SizeType i=0; i<rVariables.SlaveAngularBlockSize; i++)
    {
      for(SizeType j=0; j<rVariables.SlaveAngularBlockSize; j++)
      {
        //nodal block ij
        rLeftHandSideMatrix(start_master+rVariables.SlaveLinearBlockSize+i,start_deformable+rVariables.SlaveLinearBlockSize+j) += rLinkedLeftHandSideMatrix(start_slave+rVariables.SlaveLinearBlockSize+i,start_deformable+rVariables.SlaveLinearBlockSize+j);
        //nodal block ji
        rLeftHandSideMatrix(start_deformable+rVariables.SlaveLinearBlockSize+j,start_master+rVariables.SlaveLinearBlockSize+i) += rLinkedLeftHandSideMatrix(start_deformable+rVariables.SlaveLinearBlockSize+j,start_slave+rVariables.SlaveLinearBlockSize+i);
      }

    }

    noalias(MomentRowMatrix) = prod(rVariables.SlaveSkewSymDistance,ForceMatrix);
    noalias(MomentColumnMatrix) = prod(ForceMatrix,rVariables.SlaveSkewSymDistance);


    //angular blocks
    for(SizeType i=0; i<rVariables.SlaveAngularBlockSize; i++)
    {
      for(SizeType j=0; j<rVariables.SlaveLinearBlockSize; j++)
      {
        //column blocks
        rLeftHandSideMatrix(start_deformable+i,start_master+rVariables.MasterLinearBlockSize+j) += rLinkedLeftHandSideMatrix(start_deformable+i,start_slave+rVariables.SlaveLinearBlockSize+j);
        rLeftHandSideMatrix(start_deformable+rVariables.MasterLinearBlockSize+i,start_master+j) += rLinkedLeftHandSideMatrix(start_deformable+rVariables.SlaveLinearBlockSize+i,start_slave+j);

        ForceMatrix(i,j) = rLinkedLeftHandSideMatrix(start_deformable+rVariables.SlaveLinearBlockSize+i,start_slave+j);

        //row blocks
        rLeftHandSideMatrix(start_master+i,start_deformable+rVariables.MasterLinearBlockSize+j) += rLinkedLeftHandSideMatrix(start_slave+i,start_deformable+rVariables.SlaveLinearBlockSize+j);
        rLeftHandSideMatrix(start_master+rVariables.MasterLinearBlockSize+i,start_deformable+j) += rLinkedLeftHandSideMatrix(start_slave+rVariables.SlaveLinearBlockSize+i,start_deformable+j);
      }
    }

    // override Moment Matrix
    if(!add_row_moments)
      noalias(MomentMatrix) = prod(ForceMatrix,rVariables.SlaveSkewSymDistance);

    for(SizeType i=0; i<rVariables.MasterAngularBlockSize; i++)
    {
      for(SizeType j=0; j<rVariables.MasterLinearBlockSize; j++)
      {
        rLeftHandSideMatrix(start_deformable+j,start_master+rVariables.MasterLinearBlockSize+i) -= MomentColumnMatrix(j,start_rotation+i);
        if(add_row_moments)
          rLeftHandSideMatrix(start_master+rVariables.MasterLinearBlockSize+i,start_deformable+j) -= MomentRowMatrix(start_rotation+i,j);
        else
          rLeftHandSideMatrix(start_deformable+rVariables.SlaveLinearBlockSize+j,start_master+rVariables.MasterLinearBlockSize+i) -= MomentMatrix(j,start_rotation+i);
      }
    }

  }

  // Other rigid nodes of the linked element:
  for(SizeType n=0; n<rVariables.RigidNodes.size(); ++n)
  {

    SizeType start_rigid = rVariables.RigidNodes[n] * SlaveBlockSize;

    noalias(ForceMatrix) = ZeroMatrix(3,3);

    for(SizeType i=0; i<rVariables.SlaveLinearBlockSize; i++)
    {
      for(SizeType j=0; j<rVariables.SlaveLinearBlockSize; j++)
      {
        //nodal block ij
        rLeftHandSideMatrix(start_master+i,start_master+j) += rLinkedLeftHandSideMatrix(start_slave+i,start_rigid+j);
        ForceMatrix(i,j) = rLinkedLeftHandSideMatrix(start_slave+i,start_rigid+j);
      }
    }

    for(SizeType i=0; i<rVariables.SlaveAngularBlockSize; i++)
    {
      for(SizeType j=0; j<rVariables.SlaveAngularBlockSize; j++)
      {
        //nodal block ij
        rLeftHandSideMatrix(start_master+rVariables.SlaveLinearBlockSize+i,start_rigid+rVariables.SlaveLinearBlockSize+j) += rLinkedLeftHandSideMatrix(start_slave+rVariables.SlaveLinearBlockSize+i,start_rigid+rVariables.SlaveLinearBlockSize+j);
      }
    }


    noalias(MomentRowMatrix) = prod(rVariables.SlaveSkewSymDistance,ForceMatrix);
    noalias(MomentColumnMatrix) = prod(ForceMatrix,rVariables.RigidSkewSymDistances[n]);
    noalias(MomentMatrix) = prod(MomentRowMatrix,rVariables.RigidSkewSymDistances[n]);

    for(SizeType i=0; i<rVariables.MasterAngularBlockSize; i++)
    {
      for(SizeType j=0; j<rVariables.MasterLinearBlockSize; j++)
      {
        rLeftHandSideMatrix(start_master+rVariables.MasterLinearBlockSize+i,start_master+j) -= MomentRowMatrix(start_rotation+i,j);
        rLeftHandSideMatrix(start_master+j,start_master+rVariables.MasterLinearBlockSize+i) -= MomentColumnMatrix(j,start_rotation+i);
      }
    }


    for(SizeType i=0; i<rVariables.MasterAngularBlockSize; i++)
    {
      for(SizeType j=0; j<rVariables.MasterAngularBlockSize; j++)
      {
        rLeftHandSideMatrix(start_master+rVariables.MasterLinearBlockSize+i,start_master+rVariables.MasterLinearBlockSize+j) += MomentMatrix(start_rotation+i,start_rotation+j);
      }
    }


    for(SizeType d=0; d<rVariables.DeformableNodes.size(); ++d)
    {
      SizeType start_deformable = rVariables.DeformableNodes[d] * SlaveBlockSize;

      noalias(ForceMatrix) = ZeroMatrix(3,3);

      for(SizeType i=0; i<rVariables.SlaveLinearBlockSize; i++)
      {
        for(SizeType j=0; j<rVariables.SlaveLinearBlockSize; j++)
        {
          //nodal block ij
          rLeftHandSideMatrix(start_deformable+j,start_master+i) += rLinkedLeftHandSideMatrix(start_deformable+j,start_slave+i);
          ForceMatrix(i,j) = rLinkedLeftHandSideMatrix(start_deformable+j,start_slave+i);
        }
      }

      noalias(MomentMatrix) = prod(ForceMatrix,rVariables.RigidSkewSymDistances[n]);

      for(SizeType i=0; i<rVariables.MasterAngularBlockSize; i++)
      {
        for(SizeType j=0; j<rVariables.MasterLinearBlockSize; j++)
        {
          rLeftHandSideMatrix(start_deformable+j,start_master+rVariables.MasterLinearBlockSize+i) -= MomentMatrix(j,start_rotation+i);
        }
      }
    }

  }

  // std::cout<<" Slave "<<rVariables.SlaveNode<<std::endl;
  // WriteMatrixInRows( "rLinkedLeftHandSideMatrix", rLinkedLeftHandSideMatrix );
  // WriteMatrixInRows( "rLeftHandSideMatrix", rLeftHandSideMatrix );

  KRATOS_CATCH("")
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

  const SizeType dofs_size = this->GetDofsSize();

  SizeType SlaveBlockSize = rVariables.SlaveLinearBlockSize + rVariables.SlaveAngularBlockSize;

  SizeType start_master = rRightHandSideVector.size() - dofs_size;
  SizeType start_slave  = rVariables.SlaveNode * SlaveBlockSize;

  //std::cout<<" RHS vector "<<rRightHandSideVector<<std::endl;

  if( rVariables.SlaveAngularBlockSize != 0 ){

    for(SizeType i=0; i<SlaveBlockSize; i++)
    {
      rRightHandSideVector[start_master+i] += rLinkedRightHandSideVector[start_slave+i];
    }

  }
  else{

    Vector ForceVector(3);
    noalias(ForceVector) = ZeroVector(3);

    //std::cout<<" RHS vector "<<rRightHandSideVector<<std::endl;

    for(SizeType i=0; i<rVariables.SlaveLinearBlockSize; i++)
    {
      rRightHandSideVector[start_master+i] += rLinkedRightHandSideVector[start_slave+i];
      ForceVector[i] = rLinkedRightHandSideVector[start_slave+i];
    }

    //std::cout<<" ForceVector "<<ForceVector<<" Moment "<<MomentVector<<" skew "<<rVariables.SlaveSkewSymDistance<<std::endl;

    SizeType start_rotation = 3-rVariables.MasterAngularBlockSize;

    Vector MomentVector(3);
    noalias(MomentVector) = prod(rVariables.SlaveSkewSymDistance,ForceVector);

    for(SizeType i=0; i<rVariables.MasterAngularBlockSize; i++)
    {
      rRightHandSideVector[start_master+rVariables.MasterLinearBlockSize+i] -= MomentVector[start_rotation+i];
    }

  }

  // std::cout<<" [LINK ID:"<<this->Id()<<"]  RHS Vector "<<rLinkedRightHandSideVector<<std::endl;
  // std::cout<<" [ID:"<<this->Id()<<"] [Link Force "<<ForceVector<<"][Link Moment "<<MomentVector<<"][NodeId "<<GetGeometry()[0].Id()<<"]"<<std::endl;
  //std::cout<<" RHS Vector "<<rRightHandSideVector<<std::endl;

  //std::cout<<" LINK: "<<this->Id()<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;

  KRATOS_CATCH("")
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointLinkCondition::CalculateAndAddTangentBeam(MatrixType& rLeftHandSideMatrix,
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

  if(end_slave == rLinkedLeftHandSideMatrix.size1())
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
  noalias(ForceMatrix) = ZeroMatrix(3,3);

  for(unsigned int i=0; i<dimension; i++)
  {
    for(unsigned int j=0; j<dimension; j++)
    {
      ForceMatrix(i,j) = rLinkedLeftHandSideMatrix(end_slave+i,start_slave+j);
    }
  }

  noalias(MomentMatrix) = prod(ForceMatrix,rVariables.SlaveSkewSymDistance);

  for(unsigned int i=0; i<dimension; i++)
  {
    for(unsigned int j=0; j<dimension; j++)
    {
      rLeftHandSideMatrix(end_slave+i,start_master+dimension+j) -= MomentMatrix(i,j);
    }
  }


  //column matrices related to the rotation constraint H (U1O2)
  noalias(ForceMatrix) = ZeroMatrix(3,3);

  for(unsigned int i=0; i<dimension; i++)
  {
    for(unsigned int j=0; j<dimension; j++)
    {
      ForceMatrix(i,j) = rLinkedLeftHandSideMatrix(end_slave+dimension+i,start_slave+j);
    }
  }

  noalias(MomentMatrix) = prod(ForceMatrix,rVariables.SlaveSkewSymDistance);

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
  noalias(ForceMatrix) = ZeroMatrix(3,3);

  for(unsigned int i=0; i<dimension; i++)
  {
    for(unsigned int j=0; j<dimension; j++)
    {
      ForceMatrix(i,j) = rLinkedLeftHandSideMatrix(start_slave+i,start_slave+j);
    }
  }

  noalias(MomentMatrix) = prod(ForceMatrix,rVariables.SlaveSkewSymDistance);

  for(unsigned int i=0; i<dimension; i++)
  {
    for(unsigned int j=0; j<dimension; j++)
    {
      rLeftHandSideMatrix(start_master+i,start_master+dimension+j) -= MomentMatrix(i,j);
    }
  }

  //column matrices related to the rotation constraint H (U1O1)
  noalias(ForceMatrix) = ZeroMatrix(3,3);

  for(unsigned int i=0; i<dimension; i++)
  {
    for(unsigned int j=0; j<dimension; j++)
    {
      ForceMatrix(i,j) = rLinkedLeftHandSideMatrix(start_slave+dimension+i,start_slave+j);
    }
  }

  noalias(MomentMatrix) = prod(ForceMatrix,rVariables.SlaveSkewSymDistance);

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

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::CalculateAndAddForcesBeam(VectorType& rRightHandSideVector,
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

  //std::cout<<" [ID:"<<this->Id()<<"] [Link Force "<<Force<<"][Link Moment "<<Moment<<"][NodeId "<<GetGeometry()[0].Id()<<"]"<<std::endl;
  //std::cout<<" RHS Vector "<<rRightHandSideVector<<std::endl;


  //std::cout<<" LINK: "<<this->Id()<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

RigidBodyPointLinkCondition::SizeType RigidBodyPointLinkCondition::GetDofsSize()
{
  KRATOS_TRY

  SizeType size = 0;
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();
  const SizeType number_of_nodes =  GetGeometry().size();
  bool rotation_dofs = GetGeometry()[0].HasDofFor(ROTATION_Z);

  if(rotation_dofs){
    if( dimension == 2 )
      size = number_of_nodes * 3;
    else
      size = number_of_nodes * 6;
  }
  else
    size = number_of_nodes * dimension;

  return size;

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void RigidBodyPointLinkCondition::WriteMatrixInRows( std::string MatrixName, const MatrixType& rMatrix )
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

int RigidBodyPointLinkCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
  return 0;
}

} // Namespace Kratos.
