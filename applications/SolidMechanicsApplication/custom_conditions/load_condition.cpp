//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/load_condition.hpp"


namespace Kratos
{

//***********************************************************************************
//***********************************************************************************
LoadCondition::LoadCondition()
    : BoundaryCondition()
{
  //DO NOT CALL IT: only needed for Register and Serialization!!!
}


//***********************************************************************************
//***********************************************************************************
LoadCondition::LoadCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : BoundaryCondition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//***********************************************************************************
//***********************************************************************************
LoadCondition::LoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : BoundaryCondition(NewId, pGeometry, pProperties)
{

    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
LoadCondition::LoadCondition( LoadCondition const& rOther )
    : BoundaryCondition(rOther)
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer LoadCondition::Create(IndexType NewId,
					 NodesArrayType const& ThisNodes,
					 PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new LoadCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}


//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer LoadCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  std::cout<<" Call base class LOAD CONDITION Clone "<<std::endl;
  
  LoadCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

  NewCondition.SetData(this->GetData());
  NewCondition.SetFlags(this->GetFlags());

  
  return Condition::Pointer( new LoadCondition(NewCondition) );
}


//***********************************************************************************
//***********************************************************************************
LoadCondition::~LoadCondition()
{
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//************************************************************************************
//************************************************************************************

void LoadCondition::InitializeConditionVariables(ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BoundaryCondition::InitializeConditionVariables(rVariables, rCurrentProcessInfo);
    
    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void LoadCondition::CalculateExternalLoad(ConditionVariables& rVariables)
{
    KRATOS_TRY
      
    KRATOS_ERROR << "calling the base class CalculateExternalLoad method for a load condition... " << std::endl;
    
    KRATOS_CATCH( "" )
}
  
//***********************************************************************************
//***********************************************************************************

void LoadCondition::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
						  ConditionVariables& rVariables,
						  double& rIntegrationWeight)

{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    int index = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        index = dimension * i;
	
        for ( unsigned int j = 0; j < dimension; j++ )
        {
	  rRightHandSideVector[index + j] += rVariables.N[i] * rVariables.ExternalVectorValue[j] * rIntegrationWeight;
        }

    }

    //std::cout<<" ExternalForces ["<<this->Id()<<"]"<<rRightHandSideVector<<std::endl;

    KRATOS_CATCH( "" )
}
  
//***********************************************************************************
//***********************************************************************************

double& LoadCondition::CalculateAndAddExternalEnergy(double& rEnergy,
						     ConditionVariables& rVariables,
						     double& rIntegrationWeight)

{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Energy Calculation:
    Vector CurrentValueVector(dimension);
    noalias(CurrentValueVector) = ZeroVector(dimension); 
    Vector Displacements(dimension);
    noalias(Displacements) = ZeroVector(dimension);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	//current displacements to compute energy
	CurrentValueVector = GetNodalCurrentValue( DISPLACEMENT, CurrentValueVector, i );

	Displacements += rVariables.N[i] * CurrentValueVector;
      }
    //------

    Vector ForceVector(dimension);
    noalias(ForceVector) = ZeroVector(dimension);
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
	ForceVector += rVariables.N[i] * rVariables.ExternalVectorValue * rIntegrationWeight;
    }

    rEnergy += inner_prod( ForceVector, Displacements );

    return rEnergy;

    KRATOS_CATCH( "" )
}



//***********************************************************************************
//***********************************************************************************

void LoadCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BoundaryCondition )
}

void LoadCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BoundaryCondition )
}


} // Namespace Kratos.
