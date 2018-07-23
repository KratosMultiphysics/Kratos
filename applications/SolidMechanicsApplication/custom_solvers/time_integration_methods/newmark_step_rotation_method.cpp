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
#include "custom_solvers/time_integration_methods/newmark_step_rotation_method.hpp"
#include "utilities/beam_math_utilities.hpp"

namespace Kratos
{

  // other types
  template<class TVariableType, class TValueType>
  void NewmarkStepRotationMethod<TVariableType,TValueType>::Update(NodeType& rNode)
  {
      KRATOS_TRY

      KRATOS_ERROR << " Calling a non compatible type update for ROTATIONS in NewmarkStepRotationMethod " <<std::endl;

      KRATOS_CATCH( "" )
  }

  // specilization to array_1d
  template<>
  void NewmarkStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> >::Update(NodeType& rNode)
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

      // update step variable
      Quaternion<double> StepVariableQuaternion  = Quaternion<double>::FromRotationVector(CurrentStepVariable);

      StepVariableQuaternion = DeltaVariableQuaternion * StepVariableQuaternion;

      StepVariableQuaternion.ToRotationVector(CurrentStepVariable);

      // update variable:
      Quaternion<double> VariableQuaternion = Quaternion<double>::FromRotationVector(PreviousVariable);

      VariableQuaternion = DeltaVariableQuaternion * VariableQuaternion;

      VariableQuaternion.ToRotationVector( CurrentVariable );

      // update variable previous iteration instead of previous step
      PreviousVariable = CurrentVariable;

      // update first derivative
      array_1d<double,3>& CurrentFirstDerivative = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const array_1d<double,3>& PreviousFirstDerivative = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
      const array_1d<double,3>& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative,  1);

      CurrentFirstDerivative = this->mNewmark.c1 * CurrentStepVariable - this->mNewmark.c4 * PreviousFirstDerivative - this->mNewmark.c5 * PreviousSecondDerivative;

      // update second derivative
      array_1d<double,3>& CurrentSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative,  0);
      CurrentSecondDerivative = mNewmark.c0 * CurrentStepVariable - mNewmark.c2 * PreviousFirstDerivative - mNewmark.c3 * PreviousSecondDerivative;

      //std::cout<<*this->mpVariable<<" Update Node["<<rNode.Id()<<"]"<<CurrentVariable<<" "<<CurrentStepVariable<<" "<<CurrentFirstDerivative<<" "<<CurrentSecondDerivative<<std::endl;


      KRATOS_CATCH( "" )

   }

    template class KRATOS_API(SOLID_MECHANICS_APPLICATION) NewmarkStepRotationMethod< VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>, double>;
    template class KRATOS_API(SOLID_MECHANICS_APPLICATION) NewmarkStepRotationMethod< Variable<array_1d<double, 3>>, array_1d<double,3>>;
    template class KRATOS_API(SOLID_MECHANICS_APPLICATION) NewmarkStepRotationMethod< Variable<double>, double >;

}  // namespace Kratos.


