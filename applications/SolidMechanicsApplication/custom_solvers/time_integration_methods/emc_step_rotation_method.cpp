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
#include "custom_solvers/time_integration_methods/emc_step_rotation_method.hpp"
#include "utilities/beam_math_utilities.hpp"

namespace Kratos
{

  // specilization to array_1d
  template<class TVariableType, class TValueType>
  void EmcStepRotationMethod<TVariableType,TValueType>::Update(NodeType& rNode)
  {
      KRATOS_TRY

      KRATOS_ERROR << " Calling a non compatible type update for ROTATIONS in EmcStepRotationMethod " <<std::endl;

      KRATOS_CATCH( "" )
  }

  template<>
  void EmcStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> >::Update(NodeType& rNode)
  {
    KRATOS_TRY

    // predict step variable from previous and current values
    array_1d<double,3>& CurrentStepVariable      = rNode.FastGetSolutionStepValue(*this->mpStepVariable,    0);
    // array_1d<double,3>& PreviousStepVariable     = rNode.FastGetSolutionStepValue(*this->mpStepVariable,    1);

    array_1d<double,3>& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,        0);
    array_1d<double,3>& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,        1);

    // update step variable previous iteration instead of previous step
    // PreviousStepVariable = CurrentStepVariable;

    // update delta variable
    array_1d<double,3> DeltaVariable;
    noalias(DeltaVariable) = CurrentVariable - PreviousVariable;

    Matrix CayleyDeltaVariable(3,3);
    noalias(CayleyDeltaVariable) = ZeroMatrix(3,3);
    BeamMathUtils<double>::CayleyTransform( DeltaVariable, CayleyDeltaVariable);


    Matrix CayleyStepVariable(3,3);
    noalias(CayleyStepVariable) = ZeroMatrix(3,3);
    BeamMathUtils<double>::CayleyTransform( CurrentStepVariable, CayleyStepVariable);

    // update step variable
    Matrix ComposedVariable(3,3);
    noalias(ComposedVariable) = prod(CayleyDeltaVariable,CayleyStepVariable);

    BeamMathUtils<double>::InverseCayleyTransform( ComposedVariable, CurrentStepVariable);

    // update variable:
    Matrix CayleyVariable(3,3);
    //noalias(CayleyVariable) = ZeroMatrix(3,3);
    //BeamMathUtils<double>::CayleyTransform( PreviousVariable, CayleyVariable);
    Quaternion<double> VariableQuaternion = Quaternion<double>::FromRotationVector(PreviousVariable);
    VariableQuaternion.ToRotationMatrix( CayleyVariable );

    noalias(ComposedVariable) = prod(CayleyDeltaVariable,CayleyVariable);

    //BeamMathUtils<double>::InverseCayleyTransform( ComposedVariable, CurrentVariable);
    VariableQuaternion = Quaternion<double>::FromRotationMatrix( ComposedVariable );
    VariableQuaternion.ToRotationVector( CurrentVariable );

    // update variable previous iteration instead of previous step
    PreviousVariable = CurrentVariable;

    // update first derivative
    array_1d<double,3>& CurrentFirstDerivative = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
    const array_1d<double,3>& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, 1);

    BeamMathUtils<double>::CayleyTransform( CurrentStepVariable, CayleyStepVariable);
    //here DeltaVariable is updated previous first derivative
    noalias(DeltaVariable) = prod(CayleyStepVariable,PreviousFirstDerivative);
    noalias(CurrentFirstDerivative) = this->mEmc.c0 * CurrentStepVariable - DeltaVariable;

    // update second derivative
    array_1d<double,3>& CurrentSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative,  0);
    noalias(CurrentSecondDerivative) = this->mEmc.c1 * (CurrentFirstDerivative - DeltaVariable);
    //noalias(CurrentSecondDerivative) = (this->mEmc.c1/this->mEmc.c0) * (CurrentFirstDerivative - DeltaVariable);

    //std::cout<<*this->mpVariable<<" Update Node["<<rNode.Id()<<"]"<<CurrentVariable<<" "<<CurrentStepVariable<<" "<<CurrentFirstDerivative<<" "<<CurrentSecondDerivative<<std::endl;

    KRATOS_CATCH( "" )
  }

  template class EmcStepRotationMethod< VariableComponent<VectorComponentAdaptor<array_1d<double,3>>>, double>;
  template class EmcStepRotationMethod< Variable<array_1d<double,3>>, array_1d<double,3>>;
  template class EmcStepRotationMethod< Variable<double>, double >;

}  // namespace Kratos.


