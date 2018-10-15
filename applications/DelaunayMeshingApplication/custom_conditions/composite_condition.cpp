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
  return Kratos::make_shared<CompositeCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}



//*********************************CREATE*********************************************
//************************************************************************************

Condition::Pointer CompositeCondition::Create( IndexType NewId, GeometryType::Pointer pGeom,
					       PropertiesType::Pointer pProperties) const
{

  return Kratos::make_shared<CompositeCondition>( NewId, pGeom, pProperties );
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

  return Kratos::make_shared< CompositeCondition > (NewCompositeCondition);
}



//*******************************DESTRUCTOR*******************************************
//************************************************************************************


CompositeCondition::~CompositeCondition()
{
}


//************* SETTING METHODS
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


//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

CompositeCondition::IntegrationMethod CompositeCondition::GetIntegrationMethod()
{
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      return cn->GetIntegrationMethod();
    }

  return Condition::GetIntegrationMethod();
}


//************************************************************************************
//************************************************************************************

void CompositeCondition::GetDofList( DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo )
{
  rConditionalDofList.clear();
  rConditionalDofList.resize(0);

  DofsVectorType LocalConditionalDofList;

  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      if( IsActive(cn, rCurrentProcessInfo) ){

	cn->GetDofList(LocalConditionalDofList,rCurrentProcessInfo);

	for(unsigned int i=0; i<LocalConditionalDofList.size(); i++)
	  {
	    rConditionalDofList.push_back(LocalConditionalDofList[i]);
	  }
      }

    }

}

//************************************************************************************
//************************************************************************************

void CompositeCondition::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{

  rResult.clear();
  rResult.resize(0);

  EquationIdVectorType LocalResult;

  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      if( IsActive(cn, rCurrentProcessInfo) ){

	cn->EquationIdVector(LocalResult,rCurrentProcessInfo);

	for(unsigned int i=0; i<LocalResult.size(); i++)
	  {
	    rResult.push_back(LocalResult[i]);
	  }
      }

    }

}


//*********************************SET VALUE TO CHILDREN******************************
//************************************************************************************

bool CompositeCondition::IsActive( ConditionIterator iChildCondition, const ProcessInfo& rCurrentProcessInfo )
{
  if(rCurrentProcessInfo.Is(THERMAL) && iChildCondition->Is(THERMAL))
    return true;

  if(rCurrentProcessInfo.IsNot(THERMAL) && iChildCondition->IsNot(THERMAL))
    return true;

  return false;
}



//*********************************DISPLACEMENT***************************************
//************************************************************************************

void CompositeCondition::GetValuesVector( Vector& rValues, int Step )
{
  Vector LocalValues;

  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      cn->GetValuesVector(LocalValues,Step);
      if(LocalValues.size() != rValues.size())
	rValues.resize(LocalValues.size(),false);

      rValues+=LocalValues;
  }
}


//************************************VELOCITY****************************************
//************************************************************************************

void CompositeCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
{
  Vector LocalValues;
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      cn->GetFirstDerivativesVector(LocalValues,Step);
      if(LocalValues.size() != rValues.size())
	rValues.resize(LocalValues.size(),false);

      rValues+=LocalValues;
    }
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void CompositeCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
{
  Vector LocalValues;
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      cn->GetSecondDerivativesVector(LocalValues,Step);
      if(LocalValues.size() != rValues.size())
	rValues.resize(LocalValues.size(),false);

      rValues+=LocalValues;
  }
}


//*********************************SET VECTOR VALUE***********************************
//************************************************************************************

void CompositeCondition::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      cn->SetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
}


//*********************************SET MATRIX VALUE***********************************
//************************************************************************************

void CompositeCondition::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      cn->SetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
}

//*********************************GET DOUBLE VALUE***********************************
//************************************************************************************

void CompositeCondition::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
							    std::vector<double>& rValues,
							    const ProcessInfo& rCurrentProcessInfo )
{

  std::vector< double > LocalValues;
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      cn->GetValueOnIntegrationPoints(rVariable,LocalValues,rCurrentProcessInfo);

    if ( LocalValues.size() != rValues.size() )
      rValues.resize( LocalValues.size(), false );

    for(unsigned int i=0; i<LocalValues.size(); i++)
      rValues[i] += LocalValues[i];
  }

}


//**********************************GET VECTOR VALUE**********************************
//************************************************************************************

void CompositeCondition::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
							    std::vector<Vector>& rValues,
							    const ProcessInfo& rCurrentProcessInfo )
{
  std::vector< Vector > LocalValues;
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      cn->GetValueOnIntegrationPoints(rVariable,LocalValues,rCurrentProcessInfo);

    if ( LocalValues.size() != rValues.size() )
      rValues.resize(LocalValues.size());

    for(unsigned int i=0; i<LocalValues.size(); i++)
      {
	if ( LocalValues[i].size() != rValues[i].size() )
	  rValues[i].resize( LocalValues[i].size(), false );
      }

    for(unsigned int i=0; i<LocalValues.size(); i++)
      rValues[i] += LocalValues[i];

  }
}

//***********************************GET MATRIX VALUE*********************************
//************************************************************************************

void CompositeCondition::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
							 std::vector<Matrix>& rValues,
							 const ProcessInfo& rCurrentProcessInfo )
{
  std::vector< Matrix > LocalValues;
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      cn->GetValueOnIntegrationPoints(rVariable,LocalValues,rCurrentProcessInfo);

    if ( LocalValues.size() != rValues.size() )
      rValues.resize(LocalValues.size());

    for(unsigned int i=0; i<LocalValues.size(); i++)
      {
	if ( LocalValues[i].size2() != rValues[i].size2() )
	  rValues[i].resize( LocalValues[i].size1(), LocalValues[i].size2(), false );
      }

    for(unsigned int i=0; i<LocalValues.size(); i++)
      rValues[i] += LocalValues[i];
  }
}



//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void CompositeCondition::Initialize()
{
   KRATOS_TRY

   if(mInitializedChildren == false)
     this->InitializeChildren();

   for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      SetValueToChildren(MASTER_ELEMENTS);
      SetValueToChildren(MASTER_NODES);

      cn->Initialize();
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

void CompositeCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{

  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      if( IsActive(cn,rCurrentProcessInfo) ){

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

void CompositeCondition::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      if( IsActive(cn,rCurrentProcessInfo) ){
	cn->InitializeNonLinearIteration(rCurrentProcessInfo);
    }
  }
}


//************************************************************************************
//************************************************************************************

void CompositeCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      if( IsActive(cn,rCurrentProcessInfo) ){
	cn->FinalizeSolutionStep(rCurrentProcessInfo);
      }
    }
}


//************************************************************************************
//************************************************************************************

void CompositeCondition::AddExplicitContribution(const VectorType& rRHS,
						 const Variable<VectorType>& rRHSVariable,
						 Variable<array_1d<double,3> >& rDestinationVariable,
						 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
      {
	cn->AddExplicitContribution(rRHS, rRHSVariable, rDestinationVariable, rCurrentProcessInfo);
      }

    KRATOS_CATCH( "" )
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//***********************COMPUTE LOCAL SYSTEM CONTRIBUTIONS***************************
//************************************************************************************


void CompositeCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
  rLeftHandSideMatrix.clear();
  rLeftHandSideMatrix.resize(0,0,false);

  rRightHandSideVector.clear();
  rRightHandSideVector.resize(0,false);

  //std::cout<<" Calculate local system Skin "<<std::endl;
  VectorType LocalRightHandSideVector;
  MatrixType LocalLeftHandSideMatrix;

  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      if( IsActive(cn,rCurrentProcessInfo) ){

	cn->CalculateLocalSystem(LocalLeftHandSideMatrix,LocalRightHandSideVector,rCurrentProcessInfo);
	//std::cout<<" LocalRightHandSideVector "<<LocalRightHandSideVector<<std::endl;

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
		rLeftHandSideMatrix(i,j)=LocalLeftHandSideMatrix(indexi,indexj);
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
    }

  //std::cout<<" Skin "<<this->Id()<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;

}

//************************************************************************************
//************************************************************************************


void CompositeCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
  //std::cout<<" Calculate local rhs system Skin "<<std::endl;
  rRightHandSideVector.clear();
  rRightHandSideVector.resize(0);


  VectorType LocalRightHandSideVector;
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      if( IsActive(cn,rCurrentProcessInfo) ){

	cn->CalculateRightHandSide(LocalRightHandSideVector,rCurrentProcessInfo);

	//resizing as needed the RHS
	unsigned int GlobalSize  = rRightHandSideVector.size();
	unsigned int LocalSize   = LocalRightHandSideVector.size();
	unsigned int MatSize     = GlobalSize+LocalSize;

	rRightHandSideVector.resize( MatSize, true );

	unsigned int indexi = 0;
	for(unsigned int i=GlobalSize; i<MatSize; i++)
	  {
	    rRightHandSideVector[i]=LocalRightHandSideVector[indexi];
	    indexi++;
	  }
      }
    }

  //std::cout<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;
}


//************************************************************************************
//************************************************************************************


void CompositeCondition::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
  rLeftHandSideMatrix.clear();
  rLeftHandSideMatrix.resize(0,0,false);

  MatrixType LocalLeftHandSideMatrix;
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
    {
      if( IsActive(cn,rCurrentProcessInfo) ){

	cn->CalculateLeftHandSide(LocalLeftHandSideMatrix,rCurrentProcessInfo);

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
		rLeftHandSideMatrix(i,j)=LocalLeftHandSideMatrix(indexi,indexj);
		indexj++;
	      }
	    indexi++;
	  }
      }
    }
}


//************************************************************************************
//************************************************************************************

void CompositeCondition::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
  rMassMatrix.clear();
  rMassMatrix.resize(0,0,false);

  // MatrixType LocalMassMatrix;
  // for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
  //    {
  //      cn->CalculateMassMatrix(LocalMassMatrix,rCurrentProcessInfo);
  //   rMassMatrix+=LocalMassMatrix;
  // }
}


//************************************************************************************
//************************************************************************************


void CompositeCondition::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
  rDampingMatrix.clear();
  rDampingMatrix.resize(0,0,false);

  // MatrixType LocalDampingMatrix;
  // for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
  //    {
  //      cn->DampingMatrix(LocalDampingMatrix,rCurrentProcessInfo);
  //   rDampingMatrix+=LocalDampingMatrix;
  // }
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
int  CompositeCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
  int check = 1;
  int child_check = 1;
  for (ConditionIterator cn = mChildConditions.begin() ; cn != mChildConditions.end(); ++cn)
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
