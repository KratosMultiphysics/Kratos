//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "custom_conditions/skin_multiple_condition.hpp"

#include "pfem_solid_mechanics_application.h"

namespace Kratos
{

  KRATOS_CREATE_LOCAL_FLAG( SkinMultipleCondition, SKIN, 5 );

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SkinMultipleCondition::SkinMultipleCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
  //DO NOT ADD DOFS HERE!!!
  this->Set(SkinMultipleCondition::SKIN);
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SkinMultipleCondition::SkinMultipleCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
  this->Set(SkinMultipleCondition::SKIN);
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

SkinMultipleCondition::SkinMultipleCondition( SkinMultipleCondition const& rOther)
    :Condition(rOther)
{
  ArrayPointerConditions.resize(rOther.ArrayPointerConditions.size());
  for(unsigned int cn=0; cn<rOther.ArrayPointerConditions.size(); cn++){
    ArrayPointerConditions[cn] = rOther.ArrayPointerConditions[cn];
  }

  this->Set(SkinMultipleCondition::SKIN);
}



//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

SkinMultipleCondition&  SkinMultipleCondition::operator=(SkinMultipleCondition const& rOther)
{
  Condition::operator=(rOther);

  for(unsigned int cn=0; cn<rOther.ArrayPointerConditions.size(); cn++){
    ArrayPointerConditions[cn] = rOther.ArrayPointerConditions[cn];
  }
     
  return *this;
}


//***************************SET MULTIPLE CONDITIONS**********************************
//************************************************************************************
void SkinMultipleCondition::SetCondition (Condition::Pointer pCondition)
{
  bool set = false;
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
    if(pCondition->Id() == ArrayPointerConditions[cn]->Id())
      set=true;
  }
  if(!set)
    ArrayPointerConditions.push_back(pCondition);
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Condition::Pointer SkinMultipleCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer(new SkinMultipleCondition( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************


SkinMultipleCondition::~SkinMultipleCondition()
{
}



//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

SkinMultipleCondition::IntegrationMethod SkinMultipleCondition::GetIntegrationMethod()
{
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
      return ArrayPointerConditions[cn]->GetIntegrationMethod();
  }

  return this->GetIntegrationMethod();
}


//************************************************************************************
//************************************************************************************

void SkinMultipleCondition::GetDofList( DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo )
{
  rConditionalDofList.clear();
  rConditionalDofList.resize(0);

  DofsVectorType LocalConditionalDofList;
   
  for (ConditionPointerVectorType::iterator cn = ArrayPointerConditions.begin() ; cn != ArrayPointerConditions.end(); ++cn)
    {
      if( IsActive((*cn), rCurrentProcessInfo) ){

	  (*cn)->GetDofList(LocalConditionalDofList,rCurrentProcessInfo);
	
	  for(unsigned int i=0; i<LocalConditionalDofList.size(); i++)
	    {
	      rConditionalDofList.push_back(LocalConditionalDofList[i]);
	    }
	}
      
    }
  
}

//************************************************************************************
//************************************************************************************

void SkinMultipleCondition::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
  
  rResult.clear();
  rResult.resize(0);

  EquationIdVectorType LocalResult;
   
  for (ConditionPointerVectorType::iterator cn = ArrayPointerConditions.begin() ; cn != ArrayPointerConditions.end(); ++cn)
    {
      if( IsActive((*cn),rCurrentProcessInfo) ){

	(*cn)->EquationIdVector(LocalResult,rCurrentProcessInfo);

	for(unsigned int i=0; i<LocalResult.size(); i++)
	  {
	    rResult.push_back(LocalResult[i]);
	  }
      }
      
    }
  
}


//*********************************SET VALUE TO CHILDREN******************************
//************************************************************************************

bool SkinMultipleCondition::IsActive(Condition::Pointer pCondition, const ProcessInfo& rCurrentProcessInfo )
{
  if(rCurrentProcessInfo.Is(THERMAL) && pCondition->Is(THERMAL))
    return true;

  if(rCurrentProcessInfo.IsNot(THERMAL) && pCondition->IsNot(THERMAL))
    return true;
  
  return false;
}



//*********************************DISPLACEMENT***************************************
//************************************************************************************

void SkinMultipleCondition::GetValuesVector( Vector& rValues, int Step )
{
  Vector LocalValues;
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
      ArrayPointerConditions[cn]->GetValuesVector(LocalValues,Step);
      if(LocalValues.size() != rValues.size())
	rValues.resize(LocalValues.size(),false);
      
      rValues+=LocalValues;
  }
}


//************************************VELOCITY****************************************
//************************************************************************************

void SkinMultipleCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
{
  Vector LocalValues;
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
      ArrayPointerConditions[cn]->GetFirstDerivativesVector(LocalValues,Step);
      if(LocalValues.size() != rValues.size())
	rValues.resize(LocalValues.size(),false);

      rValues+=LocalValues;
  }
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void SkinMultipleCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
{
  Vector LocalValues;
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
      ArrayPointerConditions[cn]->GetSecondDerivativesVector(LocalValues,Step);
      if(LocalValues.size() != rValues.size())
	rValues.resize(LocalValues.size(),false);

      rValues+=LocalValues;
  }
}


//*********************************SET VECTOR VALUE***********************************
//************************************************************************************

void SkinMultipleCondition::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
{
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
    ArrayPointerConditions[cn]->SetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
  }
}


//*********************************SET MATRIX VALUE***********************************
//************************************************************************************

void SkinMultipleCondition::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
    ArrayPointerConditions[cn]->SetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
  }
}

//*********************************GET DOUBLE VALUE***********************************
//************************************************************************************

void SkinMultipleCondition::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
							    std::vector<double>& rValues,
							    const ProcessInfo& rCurrentProcessInfo )
{

  std::vector< double > LocalValues;
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){

    ArrayPointerConditions[cn]->GetValueOnIntegrationPoints(rVariable,LocalValues,rCurrentProcessInfo);

    if ( LocalValues.size() != rValues.size() )
      rValues.resize( LocalValues.size(), false );

    for(unsigned int i=0; i<LocalValues.size(); i++)
      rValues[i] += LocalValues[i];
  }

}


//**********************************GET VECTOR VALUE**********************************
//************************************************************************************

void SkinMultipleCondition::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
							    std::vector<Vector>& rValues,
							    const ProcessInfo& rCurrentProcessInfo )
{
  std::vector< Vector > LocalValues;
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
    ArrayPointerConditions[cn]->GetValueOnIntegrationPoints(rVariable,LocalValues,rCurrentProcessInfo);

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

void SkinMultipleCondition::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
							 std::vector<Matrix>& rValues, 
							 const ProcessInfo& rCurrentProcessInfo )
{
  std::vector< Matrix > LocalValues;
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){

    ArrayPointerConditions[cn]->GetValueOnIntegrationPoints(rVariable,LocalValues,rCurrentProcessInfo);

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

void SkinMultipleCondition::Initialize()
{
    KRATOS_TRY
      
      for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){

	SetValueToChildren(MASTER_ELEMENTS);
	SetValueToChildren(MASTER_NODES);

	ArrayPointerConditions[cn]->Initialize();

      }
 

    KRATOS_CATCH( "" )
}



////************************************************************************************
////************************************************************************************

void SkinMultipleCondition::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){

    if( IsActive(ArrayPointerConditions[cn],rCurrentProcessInfo) ){

      ArrayPointerConditions[cn]->SetId(this->Id());

      SetValueToChildren(MASTER_ELEMENTS);
      SetValueToChildren(MASTER_NODES);

      ArrayPointerConditions[cn]->pGetGeometry() = this->pGetGeometry();

      ArrayPointerConditions[cn]->InitializeSolutionStep(rCurrentProcessInfo);

    }
  }
}

////************************************************************************************
////************************************************************************************

void SkinMultipleCondition::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
    if( IsActive(ArrayPointerConditions[cn],rCurrentProcessInfo) ){
      ArrayPointerConditions[cn]->InitializeNonLinearIteration(rCurrentProcessInfo);
    }
  }
}


//************************************************************************************
//************************************************************************************

void SkinMultipleCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
    if( IsActive(ArrayPointerConditions[cn],rCurrentProcessInfo) ){
      ArrayPointerConditions[cn]->FinalizeSolutionStep(rCurrentProcessInfo);
    }
  }
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//***********************COMPUTE LOCAL SYSTEM CONTRIBUTIONS***************************
//************************************************************************************


void SkinMultipleCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
  rLeftHandSideMatrix.clear();
  rLeftHandSideMatrix.resize(0,0);

  rRightHandSideVector.clear();
  rRightHandSideVector.resize(0);

  //std::cout<<" Calculate local system Skin "<<std::endl;
  VectorType LocalRightHandSideVector;
  MatrixType LocalLeftHandSideMatrix;

  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
    
    if( IsActive(ArrayPointerConditions[cn],rCurrentProcessInfo) ){

      ArrayPointerConditions[cn]->CalculateLocalSystem(LocalLeftHandSideMatrix,LocalRightHandSideVector,rCurrentProcessInfo);
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


void SkinMultipleCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
  //std::cout<<" Calculate local rhs system Skin "<<std::endl;
  rRightHandSideVector.clear();
  rRightHandSideVector.resize(0);


  VectorType LocalRightHandSideVector;
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){

    if( IsActive(ArrayPointerConditions[cn],rCurrentProcessInfo) ){
      ArrayPointerConditions[cn]->CalculateRightHandSide(LocalRightHandSideVector,rCurrentProcessInfo);

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


void SkinMultipleCondition::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
  rLeftHandSideMatrix.clear();
  rLeftHandSideMatrix.resize(0,0);

  MatrixType LocalLeftHandSideMatrix;
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){

    if( IsActive(ArrayPointerConditions[cn],rCurrentProcessInfo) ){

      ArrayPointerConditions[cn]->CalculateLeftHandSide(LocalLeftHandSideMatrix,rCurrentProcessInfo);

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


void SkinMultipleCondition:: MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
  // MatrixType LocalMassMatrix;
  // for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
  //   ArrayPointerConditions[cn]->MassMatrix(LocalMassMatrix,rCurrentProcessInfo);

  //   rMassMatrix+=LocalMassMatrix;
  // }
}


//************************************************************************************
//************************************************************************************


void SkinMultipleCondition:: DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
{
  // MatrixType LocalDampMatrix;
  // for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
  //   ArrayPointerConditions[cn]->DampMatrix(LocalDampMatrix,rCurrentProcessInfo);

  //   rDampMatrix+=LocalDampMatrix;
  // }
}

//************************************************************************************
//************************************************************************************

void SkinMultipleCondition::CalculateOnIntegrationPoints( const Variable<double>& rVariable, Vector& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
  Vector LocalOutput;
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
    ArrayPointerConditions[cn]->CalculateOnIntegrationPoints(rVariable,LocalOutput,rCurrentProcessInfo);
    rOutput+=LocalOutput;
  }
}


//************************************************************************************
//************************************************************************************

void SkinMultipleCondition::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
  std::vector< Vector > LocalOutput;
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
    ArrayPointerConditions[cn]->CalculateOnIntegrationPoints(rVariable,LocalOutput,rCurrentProcessInfo);

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

void SkinMultipleCondition::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
  std::vector< Matrix > LocalOutput;
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
    ArrayPointerConditions[cn]->CalculateOnIntegrationPoints(rVariable,LocalOutput,rCurrentProcessInfo);

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

void SkinMultipleCondition::Calculate( const Variable<double>& rVariable, double& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
  double LocalOutput;
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++){
    ArrayPointerConditions[cn]->Calculate(rVariable,LocalOutput,rCurrentProcessInfo);
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
int  SkinMultipleCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
  int check = 0;
  for(unsigned int cn=0; cn<ArrayPointerConditions.size(); cn++)
    ArrayPointerConditions[cn]->Check(rCurrentProcessInfo);
    
  return check;
}


void SkinMultipleCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    rSerializer.save("ArrayPointerConditions",ArrayPointerConditions);
}

void SkinMultipleCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    rSerializer.load("ArrayPointerConditions",ArrayPointerConditions);
}



} // Namespace Kratos


