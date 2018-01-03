//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_strategies/time_integration_methods/bossak_step_rotation_method.hpp"
#include "utilities/beam_math_utilities.hpp"

namespace Kratos
{

  // specilization to array_1d
  
  template<>
  void BossakStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> >::UpdateStepRotationVariable(NodeType& rNode)
  {
      KRATOS_TRY

      // predict step variable from previous and current values
      array_1d<double,3>& CurrentStepVariable      = rNode.FastGetSolutionStepValue(*this->mpStepVariable,    0);
     
      array_1d<double,3>& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,        0);
      array_1d<double,3>& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,        1);

      // recover PreviousVariable:  ( previous iteration CurrentVariable = PreviousVariable + CurrentStepVariable ) **
      PreviousVariable += CurrentStepVariable;
      
      // update delta variable      
      array_1d<double,3> DeltaVariable;
      noalias(DeltaVariable) = CurrentVariable - PreviousVariable;
      
      Quaternion<double> DeltaVariableQuaternion = Quaternion<double>::FromRotationVector(DeltaVariable);

      // update step variable
      Quaternion<double> StepVariableQuaternion  = Quaternion<double>::FromRotationVector(CurrentStepVariable);

      StepVariableQuaternion = DeltaVariableQuaternion * StepVariableQuaternion;
            
      StepVariableQuaternion.ToRotationVector(CurrentStepVariable);
      
      // update variable:
      Quaternion<double> VariableQuaternion = Quaternion<double>::FromRotationVector(PreviousVariable);
      
      VariableQuaternion = DeltaVariableQuaternion * VariableQuaternion;

      VariableQuaternion.ToRotationVector( CurrentVariable );
      
      // update PreviousVariable: ( current iteration CurrentStepVariable = CurrentVariable - PreviousVariable ) **
      PreviousVariable = CurrentVariable - CurrentStepVariable;
       
      KRATOS_CATCH( "" )

   }

  
}  // namespace Kratos.


