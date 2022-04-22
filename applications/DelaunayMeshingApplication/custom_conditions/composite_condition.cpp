//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/composite_condition.hpp"

#include "delaunay_meshing_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

CompositeCondition::CompositeCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
  //DO NOT ADD DOFS HERE!!!
  this->Set(BOUNDARY);
  mInitializedChildren = false;
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

CompositeCondition::CompositeCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
  this->Set(BOUNDARY);
  mInitializedChildren = false;
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

CompositeCondition::CompositeCondition( CompositeCondition const& rOther)
    :Condition(rOther)
{

  //mChildConditions = rOther.mChildConditions;

  //deep copy
  ConditionsContainerType ChildConditionsTemporal;
  for (ConditionConstantIterator cn = rOther.mChildConditions.begin() ; cn != rOther.mChildConditions.end(); ++cn)
    {
      ChildConditionsTemporal.push_back(*(cn.base()));
    }

  mChildConditions.swap(ChildConditionsTemporal);

  this->Set(BOUNDARY);

  this->SetValue(CHILDREN_CONDITIONS, mChildConditions);

  mInitializedChildren = rOther.mInitializedChildren;

}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

CompositeCondition&  CompositeCondition::operator=(CompositeCondition const& rOther)
{
  Condition::operator=(rOther);

  //mChildConditions = rOther.mChildConditions;

  //deep copy
  ConditionsContainerType ChildConditionsTemporal;
  for (ConditionConstantIterator cn = rOther.mChildConditions.begin() ; cn != rOther.mChildConditions.end(); ++cn)
    {
      ChildConditionsTemporal.push_back(*(cn.base()));
    }

  mChildConditions.swap(ChildConditionsTemporal);

  this->Set(BOUNDARY);

  this->SetValue(CHILDREN_CONDITIONS, mChildConditions);

  mInitializedChildren = rOther.mInitializedChildren;

  return *this;
}


//*********************************CREATE*********************************************
//************************************************************************************

Condition::Pointer CompositeCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,
					       PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_intrusive<CompositeCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}



//*********************************CREATE*********************************************
//************************************************************************************

Condition::Pointer CompositeCondition::Create( IndexType NewId, GeometryType::Pointer pGeom,
					       PropertiesType::Pointer pProperties) const
{

  return Kratos::make_intrusive<CompositeCondition>( NewId, pGeom, pProperties );
}


//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer CompositeCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

  CompositeCondition NewCompositeCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

  if( this->mChildConditions.size() )
    NewCompositeCondition.mInitializedChildren = true;

  for (ConditionConstantIterator cn = this->mChildConditions.begin(); cn != this->mChildConditions.end(); ++cn)
    {
      Condition::Pointer pNewChildCondition;
      pNewChildCondition = cn->Clone(cn->Id(), rThisNodes);
      NewCompositeCondition.AddChild( pNewChildCondition );
      // std::cout<<" Add Child "<<std::endl;
      //NewCompositeCondition.AddChild(*(cn.base())); //problems in split, previous geometry is preserved, clone is needed.
    }

  NewCompositeCondition.SetData(this->GetData());
  NewCompositeCondition.SetFlags(this->GetFlags());

  return Kratos::make_intrusive< CompositeCondition > (NewCompositeCondition);
}



//*******************************DESTRUCTOR*******************************************
//************************************************************************************


CompositeCondition::~CompositeCondition()
{
}

//************************************************************************************
//************************************************************************************

void CompositeCondition::AddChild(ConditionType::Pointer pNewChildCondition)
{
  bool set = false;

  if( pNewChildCondition->Is(BOUNDARY) ){ //not add composite conditions as child conditions
    set = true;
  }
  else{
    for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
      {
	if(pNewChildCondition->Id() == cn->Id())
	  set=true;
      }
  }

  if(!set)
    mChildConditions.insert(mChildConditions.begin(), pNewChildCondition);
}

//************************************************************************************
//************************************************************************************

CompositeCondition::IntegrationMethod CompositeCondition::GetIntegrationMethod() const
{
  for (ConditionConstantIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      return cn->GetIntegrationMethod();
    }

  return Condition::GetIntegrationMethod();
}

//************************************************************************************
//************************************************************************************

void CompositeCondition::GetDofList( DofsVectorType& rConditionalDofList, const ProcessInfo& rCurrentProcessInfo ) const
{
  rConditionalDofList.resize(0);

  DofsVectorType LocalConditionalDofList;

  for (const auto& cn : mChildConditions)
    {
      if( IsActive(cn, rCurrentProcessInfo) ){

	cn.GetDofList(LocalConditionalDofList,rCurrentProcessInfo);

	for(unsigned int i=0; i<LocalConditionalDofList.size(); i++)
	  {
	    rConditionalDofList.push_back(LocalConditionalDofList[i]);
	  }
      }

    }

}

//************************************************************************************
//************************************************************************************

void CompositeCondition::EquationIdVector( EquationIdVectorType& rResult,
                                          const ProcessInfo& rCurrentProcessInfo ) const
{

  rResult.resize(0,false);

  EquationIdVectorType LocalResult;

  for (const auto& cn : mChildConditions)
    {
      if( IsActive(cn, rCurrentProcessInfo) ){

	cn.EquationIdVector(LocalResult,rCurrentProcessInfo);

	for(unsigned int i=0; i<LocalResult.size(); i++)
	  {
	    rResult.push_back(LocalResult[i]);
	  }
      }

    }

}


//*********************************SET VALUE TO CHILDREN******************************
//************************************************************************************

bool CompositeCondition::IsActive( const Condition& rChildCondition, const ProcessInfo& rCurrentProcessInfo ) const
{
  if(rCurrentProcessInfo.Is(THERMAL) && rChildCondition.Is(THERMAL))
    return true;

  if(rCurrentProcessInfo.IsNot(THERMAL) && rChildCondition.IsNot(THERMAL))
    return true;

  return false;
}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void CompositeCondition::GetValuesVector( Vector& rValues, int Step ) const
{
  Vector ChildValues;

  SizeType indexi = 0;
  SizeType sizei  = 0;
  for (const auto& cn : mChildConditions)
    {
      cn.GetValuesVector(ChildValues,Step);

      sizei += ChildValues.size();
      rValues.resize(sizei,true);

      for(SizeType i=0; i<ChildValues.size(); i++)
        rValues[indexi+i] = ChildValues[i];

      indexi += ChildValues.size();
  }
}


//************************************VELOCITY****************************************
//************************************************************************************

void CompositeCondition::GetFirstDerivativesVector( Vector& rValues, int Step ) const
{
  Vector ChildValues;

  SizeType indexi = 0;
  SizeType sizei  = 0;
  for (const auto& cn : mChildConditions)
    {
      cn.GetFirstDerivativesVector(ChildValues,Step);

      sizei += ChildValues.size();
      rValues.resize(sizei,true);

      for(SizeType i=0; i<ChildValues.size(); i++)
        rValues[indexi+i] = ChildValues[i];

      indexi += ChildValues.size();
  }
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void CompositeCondition::GetSecondDerivativesVector( Vector& rValues, int Step ) const
{
  Vector ChildValues;

  SizeType indexi = 0;
  SizeType sizei  = 0;
  for (const auto& cn : mChildConditions)
    {
      cn.GetSecondDerivativesVector(ChildValues,Step);

      sizei += ChildValues.size();
      rValues.resize(sizei,true);

      for(SizeType i=0; i<ChildValues.size(); i++)
        rValues[indexi+i] = ChildValues[i];

      indexi += ChildValues.size();
  }
}


//*********************************SET VECTOR VALUE***********************************
//************************************************************************************

void CompositeCondition::SetValuesOnIntegrationPoints( const Variable<Vector>& rVariable, const std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      cn->SetValuesOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
}


//*********************************SET MATRIX VALUE***********************************
//************************************************************************************

void CompositeCondition::SetValuesOnIntegrationPoints( const Variable<Matrix>& rVariable,
        const std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      cn->SetValuesOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
}

//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void CompositeCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
   KRATOS_TRY

   if(mInitializedChildren == false)
     this->InitializeChildren();

   for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      SetValueToChildren(MASTER_ELEMENTS);
      SetValueToChildren(MASTER_NODES);

      cn->Initialize(rCurrentProcessInfo);
    }


   KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void CompositeCondition::InitializeChildren()
{
  KRATOS_TRY

  if(!this->Has(CHILDREN_CONDITIONS)){

    this->SetValue(CHILDREN_CONDITIONS, mChildConditions);

  }
  else{

    if( mChildConditions.size() == 0 ){

      ConditionContainerType& ChildConditions = this->GetValue(CHILDREN_CONDITIONS);

      for (ConditionIterator cn = ChildConditions.begin() ; cn != ChildConditions.end(); ++cn)
	{
	  AddChild(*(cn.base()));
	}

      this->SetValue(CHILDREN_CONDITIONS, mChildConditions);
    }
    else{

      this->SetValue(CHILDREN_CONDITIONS, mChildConditions);
    }

  }

  mInitializedChildren = true;

  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void CompositeCondition::InitializeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{

  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      if( IsActive(*cn,rCurrentProcessInfo) ){

	SetValueToChildren(MASTER_ELEMENTS);
	SetValueToChildren(MASTER_NODES);

	bool set_geometry = false;  //check direct order
	for( unsigned int i=0; i<this->GetGeometry().size(); i++ )
	  {
 	    if( cn->GetGeometry()[i].Id() != this->GetGeometry()[i].Id() )
	      set_geometry = true;
	  }

	if( set_geometry ){ //check reverse order
	    set_geometry = false;
	    unsigned int size = this->GetGeometry().size()-1;
	    for( unsigned int i=0; i<this->GetGeometry().size(); i++ )
	    {
		if( cn->GetGeometry()[i].Id() != this->GetGeometry()[size].Id() )
		    set_geometry = true;
		--size;
	    }
	}

	if( set_geometry ){

	  std::cout<<" Set Geometry ( Something is wrong with children conditions ) "<<std::endl;

	  std::cout<<" Master "<<this->Id()<<" Geometry ["<<this->GetGeometry()[0].Id()<<", "<<this->GetGeometry()[1].Id()<<"] "<<std::endl;

	  std::cout<<" Pre Child "<<cn->Id()<<" Geometry ["<<cn->GetGeometry()[0].Id()<<", "<<cn->GetGeometry()[1].Id()<<"] "<<std::endl;

	  // GeometryType::PointsArrayType GeometryPoints = this->GetGeometry().Points();

	  // (cn->GetGeometry().Points()).swap(GeometryPoints);

	  // std::cout<<" Post Child "<<cn->Id()<<" Geometry ["<<cn->GetGeometry()[0].Id()<<", "<<cn->GetGeometry()[1].Id()<<"] "<<std::endl;

	}


	cn->InitializeSolutionStep(rCurrentProcessInfo);

      }
    }
}

//************************************************************************************
//************************************************************************************

void CompositeCondition::InitializeNonLinearIteration( const ProcessInfo& rCurrentProcessInfo )
{
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      if( IsActive(*cn,rCurrentProcessInfo) ){
	cn->InitializeNonLinearIteration(rCurrentProcessInfo);
    }
  }
}


//************************************************************************************
//************************************************************************************

void CompositeCondition::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      if( IsActive(*cn,rCurrentProcessInfo) ){
	cn->FinalizeSolutionStep(rCurrentProcessInfo);
      }
    }
}


//************************************************************************************
//************************************************************************************

void CompositeCondition::AddExplicitContribution(const VectorType& rRHS,
						 const Variable<VectorType>& rRHSVariable,
						 const Variable<array_1d<double,3> >& rDestinationVariable,
						 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      cn->AddExplicitContribution(rRHS, rRHSVariable, rDestinationVariable, rCurrentProcessInfo);
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

CompositeCondition::SizeType CompositeCondition::GetDofsSize(const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  SizeType size = 0;

  EquationIdVectorType  ChildResult;
  std::vector<SizeType> element_dofs;
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
  {
    cn->EquationIdVector(ChildResult, rCurrentProcessInfo);
    element_dofs.push_back(ChildResult.size());
    size += element_dofs.back();
  }

  return size;

  KRATOS_CATCH("")
}

//***********************COMPUTE LOCAL SYSTEM CONTRIBUTIONS***************************
//************************************************************************************


void CompositeCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY
  //std::cout<<" Calculate local system Skin "<<std::endl;
  SizeType size = this->GetDofsSize(rCurrentProcessInfo);

  if ( rLeftHandSideMatrix.size1() != size )
    rLeftHandSideMatrix.resize( size, size, false );

  noalias( rLeftHandSideMatrix ) = ZeroMatrix( size, size ); //resetting LHS

  if ( rRightHandSideVector.size() != size )
    rRightHandSideVector.resize( size, false );

  noalias(rRightHandSideVector) = ZeroVector( size ); //resetting RHS

  VectorType LocalRightHandSideVector;
  MatrixType LocalLeftHandSideMatrix;

  SizeType MatSize1 = 0;
  SizeType MatSize2 = 0;

  SizeType ChildSize1 = 0;
  SizeType ChildSize2 = 0;

  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
  {
    if( IsActive(*cn,rCurrentProcessInfo) ){

      cn->CalculateLocalSystem(LocalLeftHandSideMatrix,LocalRightHandSideVector,rCurrentProcessInfo);
      //std::cout<<" LocalRightHandSideVector "<<LocalRightHandSideVector<<std::endl;

      ChildSize1 = LocalLeftHandSideMatrix.size1();
      ChildSize2 = LocalLeftHandSideMatrix.size2();

      for(unsigned int i=0; i<ChildSize2; i++)
      {
        SizeType indexj = MatSize1;
        for(unsigned int j=0; j<ChildSize1; j++)
        {
          rLeftHandSideMatrix(indexj,MatSize2)=LocalLeftHandSideMatrix(j,i);
          ++indexj;
        }
        rRightHandSideVector[MatSize2]=LocalRightHandSideVector[i];
        ++MatSize2;
      }
      MatSize1+=ChildSize1;
    }
  }

  //std::cout<<" Skin "<<this->Id()<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;
  KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************


void CompositeCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
  //std::cout<<" Calculate local rhs system Skin "<<std::endl;

  SizeType size = this->GetDofsSize(rCurrentProcessInfo);

  if ( rRightHandSideVector.size() != size )
    rRightHandSideVector.resize( size, false );

  noalias(rRightHandSideVector) = ZeroVector( size ); //resetting RHS

  VectorType LocalRightHandSideVector;

  SizeType MatSize1 = 0;
  SizeType ChildSize1 = 0;

  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
  {
    if( IsActive(*cn,rCurrentProcessInfo) ){

      cn->CalculateRightHandSide(LocalRightHandSideVector,rCurrentProcessInfo);

      ChildSize1 = LocalRightHandSideVector.size();

      for(unsigned int i=0; i<ChildSize1; i++)
      {
        rRightHandSideVector[MatSize1]=LocalRightHandSideVector[i];
        ++MatSize1;
      }
    }
  }

  //std::cout<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;
}


//************************************************************************************
//************************************************************************************


void CompositeCondition::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
{

  SizeType size = this->GetDofsSize(rCurrentProcessInfo);

  if ( rLeftHandSideMatrix.size1() != size )
    rLeftHandSideMatrix.resize( size, size, false );

  noalias( rLeftHandSideMatrix ) = ZeroMatrix( size, size ); //resetting LHS

  MatrixType LocalLeftHandSideMatrix;

  SizeType MatSize1 = 0;
  SizeType MatSize2 = 0;

  SizeType ChildSize1 = 0;
  SizeType ChildSize2 = 0;

  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
  {
    if( IsActive(*cn,rCurrentProcessInfo) ){

      cn->CalculateLeftHandSide(LocalLeftHandSideMatrix,rCurrentProcessInfo);

      ChildSize1 = LocalLeftHandSideMatrix.size1();
      ChildSize2 = LocalLeftHandSideMatrix.size2();

      for(unsigned int i=0; i<ChildSize1; i++)
      {
        SizeType indexj = MatSize2;
        for(unsigned int j=0; j<ChildSize2; j++)
        {
          rLeftHandSideMatrix(MatSize1,indexj)=LocalLeftHandSideMatrix(i,j);
          ++indexj;
        }
        ++MatSize1;
      }
      MatSize2+=ChildSize2;
    }
  }
}


//************************************************************************************
//************************************************************************************

void CompositeCondition::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
  rMassMatrix.clear();
  rMassMatrix.resize(0,0,false);
}


//************************************************************************************
//************************************************************************************


void CompositeCondition::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
  rDampingMatrix.clear();
  rDampingMatrix.resize(0,0,false);
}

//************************************************************************************
//************************************************************************************

void CompositeCondition::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
  std::vector<double> LocalOutput;
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
  {
    cn->CalculateOnIntegrationPoints(rVariable,LocalOutput,rCurrentProcessInfo);

    if ( LocalOutput.size() != rOutput.size() )
      rOutput.resize(LocalOutput.size(),true);

    for(unsigned int i=0; i<LocalOutput.size(); i++)
    {
      rOutput[i]+=LocalOutput[i];
    }

  }
}


//************************************************************************************
//************************************************************************************

void CompositeCondition::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
  std::vector< Vector > LocalOutput;
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
  {
    cn->CalculateOnIntegrationPoints(rVariable,LocalOutput,rCurrentProcessInfo);

    if ( LocalOutput.size() != rOutput.size() )
      rOutput.resize(LocalOutput.size());

    for(unsigned int i=0; i<LocalOutput.size(); i++)
    {
      if ( LocalOutput[i].size() != rOutput[i].size() )
        rOutput[i].resize( LocalOutput[i].size(), false );
    }

    for(unsigned int i=0; i<LocalOutput.size(); i++)
    {
      rOutput[i]+=LocalOutput[i];
    }
  }
}

//************************************************************************************
//************************************************************************************

void CompositeCondition::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
  std::vector< Matrix > LocalOutput;
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
  {
    cn->CalculateOnIntegrationPoints(rVariable,LocalOutput,rCurrentProcessInfo);

    if ( LocalOutput.size() != rOutput.size() )
      rOutput.resize(LocalOutput.size());

    for(unsigned int i=0; i<LocalOutput.size(); i++)
    {
      if ( LocalOutput[i].size2() != rOutput[i].size2() )
        rOutput[i].resize( LocalOutput[i].size1(), LocalOutput[i].size2(), false );
    }

    for(unsigned int i=0; i<LocalOutput.size(); i++)
    {
      rOutput[i]+=LocalOutput[i];
    }
  }

}


//************************************************************************************
//************************************************************************************

void CompositeCondition::Calculate( const Variable<double>& rVariable, double& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
  double LocalOutput;
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
  {
    cn->Calculate(rVariable,LocalOutput,rCurrentProcessInfo);
    rOutput+=LocalOutput;
  }

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
int  CompositeCondition::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
  int check = 1;
  int child_check = 1;
  for (auto cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
  {
    child_check = cn->Check(rCurrentProcessInfo);
    if(child_check == 0)
      check = 0;
  }

  return check;
}


void CompositeCondition::save( Serializer& rSerializer ) const
{
  KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
  rSerializer.save("mChildConditions",mChildConditions);
  rSerializer.save("mInitializedChildren",mInitializedChildren);
}

void CompositeCondition::load( Serializer& rSerializer )
{
  KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
  rSerializer.load("mChildConditions",mChildConditions);
  rSerializer.load("mInitializedChildren",mInitializedChildren);
}



} // Namespace Kratos
