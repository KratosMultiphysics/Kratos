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
#include "custom_solvers/time_integration_methods/static_step_rotation_method.hpp"
#include "utilities/beam_math_utilities.hpp"

namespace Kratos
{

  // other types
  template<class TVariableType, class TValueType>
  void StaticStepRotationMethod<TVariableType,TValueType>::Update(NodeType& rNode)
  {
      KRATOS_TRY

      KRATOS_ERROR << " Calling a non compatible type update for ROTATIONS in StaticStepRotationMethod " <<std::endl;

      KRATOS_CATCH( "" )
  }

  // specilization to array_1d
  template<>
  void StaticStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> >::Update(NodeType& rNode)
  {
      KRATOS_TRY

      // predict step variable from previous and current values
      array_1d<double,3>& CurrentStepVariable      = rNode.FastGetSolutionStepValue(*this->mpStepVariable,    0);

      array_1d<double,3>& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,        0);
      array_1d<double,3>& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,        1);

      // update delta variable
      array_1d<double,3> DeltaVariable;
      noalias(DeltaVariable) = CurrentVariable - PreviousVariable;

      Quaternion<double> DeltaVariableQuaternion = Quaternion<double>::FromRotationVector(DeltaVariable);

      // linear delta variable
      array_1d<double,3> LinearDeltaVariable;
      noalias(LinearDeltaVariable) = -CurrentStepVariable;

      // update step variable
      Quaternion<double> StepVariableQuaternion  = Quaternion<double>::FromRotationVector(CurrentStepVariable);

      StepVariableQuaternion = DeltaVariableQuaternion * StepVariableQuaternion;

      StepVariableQuaternion.ToRotationVector(CurrentStepVariable);

      LinearDeltaVariable += CurrentStepVariable;

      // update variable:
      Quaternion<double> VariableQuaternion = Quaternion<double>::FromRotationVector(PreviousVariable);

      VariableQuaternion = DeltaVariableQuaternion * VariableQuaternion;

      VariableQuaternion.ToRotationVector( CurrentVariable );

      // update variable previous iteration instead of previous step
      PreviousVariable     = CurrentVariable;

      // update linear delta variable:
      VariableQuaternion  = StepVariableQuaternion.conjugate() * VariableQuaternion;
      LinearDeltaVariable = BeamMathUtils<double>::MapToCurrentLocalFrame( VariableQuaternion, LinearDeltaVariable );


      KRATOS_CATCH( "" )
    }

    template class StaticStepRotationMethod< VariableComponent<VectorComponentAdaptor<array_1d<double,3>>>, double >;
    template class StaticStepRotationMethod< Variable<array_1d<double,3>>, array_1d<double,3> >;
    template class StaticStepRotationMethod< Variable<double>, double >;

}  // namespace Kratos.


