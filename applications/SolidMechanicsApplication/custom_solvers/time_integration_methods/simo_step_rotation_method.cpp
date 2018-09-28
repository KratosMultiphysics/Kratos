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
#include "custom_solvers/time_integration_methods/simo_step_rotation_method.hpp"
#include "utilities/beam_math_utilities.hpp"

namespace Kratos
{

  // specilization to array_1d
  template<class TVariableType, class TValueType>
  void SimoStepRotationMethod<TVariableType,TValueType>::Update(NodeType& rNode)
  {
      KRATOS_TRY

      KRATOS_ERROR << " Calling a non compatible type update for ROTATIONS in SimoStepRotationScheme " <<std::endl;

      KRATOS_CATCH( "" )
  }

  template<>
  void SimoStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> >::Update(NodeType& rNode)
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
      PreviousVariable   = CurrentVariable;

      // update linear delta variable:
      VariableQuaternion  = StepVariableQuaternion.conjugate() * VariableQuaternion;
      LinearDeltaVariable = BeamMathUtils<double>::MapToCurrentLocalFrame( VariableQuaternion, LinearDeltaVariable );

      // update first derivative
      array_1d<double,3>& CurrentFirstDerivative = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      noalias(CurrentFirstDerivative) += this->mNewmark.c1 * LinearDeltaVariable;

      // update second derivative
      array_1d<double,3>& CurrentSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative,  0);
      noalias(CurrentSecondDerivative) += this->mNewmark.c0 * LinearDeltaVariable;

      // std::cout<<*this->mpVariable<<" Update Node["<<rNode.Id()<<"]"<<CurrentVariable<<" "<<CurrentStepVariable<<" "<<CurrentFirstDerivative<<" "<<CurrentSecondDerivative<<std::endl;


      KRATOS_CATCH( "" )
    }

    template class KRATOS_API(SOLID_MECHANICS_APPLICATION) SimoStepRotationMethod< VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> ,  double >;
    template class KRATOS_API(SOLID_MECHANICS_APPLICATION) SimoStepRotationMethod< Variable<array_1d<double, 3> >,                                  array_1d<double,3>>;
    template class KRATOS_API(SOLID_MECHANICS_APPLICATION) SimoStepRotationMethod< Variable<double>, double >;

}  // namespace Kratos.


