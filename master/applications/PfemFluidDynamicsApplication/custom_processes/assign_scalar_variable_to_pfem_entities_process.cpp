//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                 August 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_processes/assign_scalar_variable_to_pfem_entities_process.hpp"

namespace Kratos
{

template< class TVarType, class TDataType >
void AssignScalarVariableToPfemEntitiesProcess::AssignValueToNodes(const TVarType& rVariable, const TDataType value)
{
  const int nnodes = mrModelPart.Nodes().size();

  typedef void (AssignScalarVariableToPfemEntitiesProcess::*AssignmentMethodPointer) (ModelPart::NodeType&, const TVarType&, const TDataType&);

  AssignmentMethodPointer AssignmentMethod = this->GetAssignmentMethod<AssignmentMethodPointer>();

  if(nnodes != 0)
  {
    ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh().NodesBegin();

#pragma omp parallel for
    for(int i = 0; i<nnodes; i++)
    {
      ModelPart::NodesContainerType::iterator it = it_begin + i;

      (this->*AssignmentMethod)(*it, rVariable, value);
    }
  }

}

template< class TVarType, class TDataType >
void AssignScalarVariableToPfemEntitiesProcess::AssignValueToConditions(const TVarType& rVariable, const TDataType value)
{
  const int nconditions = mrModelPart.GetMesh().Conditions().size();

  typedef void (AssignScalarVariableToPfemEntitiesProcess::*AssignmentMethodPointer) (ModelPart::ConditionType&, const TVarType&, const TDataType&);

  AssignmentMethodPointer AssignmentMethod = this->GetAssignmentMethod<AssignmentMethodPointer>();

  if(nconditions != 0)
  {
    ModelPart::ConditionsContainerType::iterator it_begin = mrModelPart.GetMesh().ConditionsBegin();

#pragma omp parallel for
    for(int i = 0; i<nconditions; i++)
    {
      ModelPart::ConditionsContainerType::iterator it = it_begin + i;

      (this->*AssignmentMethod)(*it, rVariable, value);
    }
  }
}

template< class TVarType, class TDataType >
void AssignScalarVariableToPfemEntitiesProcess::AssignValueToElements(const TVarType& rVariable, const TDataType value)
{
  const int nelements = mrModelPart.GetMesh().Elements().size();

  typedef void (AssignScalarVariableToPfemEntitiesProcess::*AssignmentMethodPointer) (ModelPart::ElementType&, const TVarType&, const TDataType&);

  AssignmentMethodPointer AssignmentMethod = this->GetAssignmentMethod<AssignmentMethodPointer>();

  if(nelements != 0)
  {
    ModelPart::ElementsContainerType::iterator it_begin = mrModelPart.GetMesh().ElementsBegin();

#pragma omp parallel for
    for(int i = 0; i<nelements; i++)
    {
      ModelPart::ElementsContainerType::iterator it = it_begin + i;

      (this->*AssignmentMethod)(*it, rVariable, value);
    }
  }
}


template<>
void AssignScalarVariableToPfemEntitiesProcess::AssignValueToNodes<Variable<bool>, bool>(const Variable<bool>& rVariable, const bool value)
{
  const int nnodes = mrModelPart.Nodes().size();

  if(nnodes != 0)
  {
    ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh().NodesBegin();

    #pragma omp parallel for
    for(int i = 0; i<nnodes; i++)
    {
      ModelPart::NodesContainerType::iterator it = it_begin + i;

      this->DirectAssignValue(*it, rVariable, value);
    }
  }
}

template<>
void AssignScalarVariableToPfemEntitiesProcess::AssignValueToConditions<Variable<bool>, bool>(const Variable<bool>& rVariable, const bool value)
{
  const int nconditions = mrModelPart.GetMesh().Conditions().size();

  if(nconditions != 0)
  {
    ModelPart::ConditionsContainerType::iterator it_begin = mrModelPart.GetMesh().ConditionsBegin();

    #pragma omp parallel for
    for(int i = 0; i<nconditions; i++)
    {
      ModelPart::ConditionsContainerType::iterator it = it_begin + i;

      this->DirectAssignValue(*it, rVariable, value);
    }
  }
}

template<>
void AssignScalarVariableToPfemEntitiesProcess::AssignValueToElements<Variable<bool>, bool>(const Variable<bool>& rVariable, const bool value)
{
  const int nelements = mrModelPart.GetMesh().Elements().size();

  if(nelements != 0)
  {
    ModelPart::ElementsContainerType::iterator it_begin = mrModelPart.GetMesh().ElementsBegin();

    #pragma omp parallel for
    for(int i = 0; i<nelements; i++)
    {
      ModelPart::ElementsContainerType::iterator it = it_begin + i;

      this->DirectAssignValue(*it, rVariable, value);
    }
  }
}

void AssignScalarVariableToPfemEntitiesProcess::Execute()
{

  KRATOS_TRY

  if( mEntity == EntityType::NODES ){

    if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
    {
      AssignValueToNodes<>(KratosComponents< Variable<double> >::Get(mvariable_name), mdouble_value);
    }
    else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
    {
      AssignValueToNodes<>(KratosComponents< Variable<int> >::Get(mvariable_name) , mint_value);
    }
    else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
    {
      AssignValueToNodes<>(KratosComponents< Variable<bool> >::Get(mvariable_name), mbool_value);
    }
    else
    {
      KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << mvariable_name << std::endl;
    }

  }
  else if( mEntity == EntityType::CONDITIONS ) {

    if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
    {
      AssignValueToConditions<>(KratosComponents< Variable<double> >::Get(mvariable_name), mdouble_value);
    }
    else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
    {
      AssignValueToConditions<>(KratosComponents< Variable<int> >::Get(mvariable_name) , mint_value);
    }
    else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
    {
      AssignValueToConditions<>(KratosComponents< Variable<bool> >::Get(mvariable_name), mbool_value);
    }
    else
    {
      KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << mvariable_name << std::endl;
    }

  }
  else if( mEntity == EntityType::ELEMENTS ) {

    if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
    {
      AssignValueToElements<>(KratosComponents< Variable<double> >::Get(mvariable_name), mdouble_value);
    }
    else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
    {
      AssignValueToElements<>(KratosComponents< Variable<int> >::Get(mvariable_name) , mint_value);
    }
    else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
    {
      AssignValueToElements<>(KratosComponents< Variable<bool> >::Get(mvariable_name), mbool_value);
    }
    else
    {
      KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << mvariable_name << std::endl;
    }

  }

  KRATOS_CATCH("");

}


void AssignScalarVariableToPfemEntitiesProcess::ExecuteFinalize()
{

  KRATOS_TRY

  if( mEntity == EntityType::CONDITIONS ){

    if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
    {
      double double_value = 0;
      AssignValueToConditions<>(KratosComponents< Variable<double> >::Get(mvariable_name), double_value);
    }
    else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
    {
      int int_value = 0;
      AssignValueToConditions<>(KratosComponents< Variable<int> >::Get(mvariable_name), int_value);
    }
    else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
    {
      bool bool_value = !mbool_value;
      AssignValueToConditions<>(KratosComponents< Variable<bool> >::Get(mvariable_name), bool_value);
    }
    else
    {
      KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << mvariable_name << std::endl;
    }
  }

  KRATOS_CATCH("")
}

}  // namespace Kratos.
