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
#include "custom_strategies/time_integration_methods/emc_step_rotation_method.hpp"
#include "utilities/beam_math_utilities.hpp"

namespace Kratos
{
 
  // specilization to array_1d
  
  template<>
  void EmcStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> >::Update(NodeType& rNode)
  {
    KRATOS_TRY

    // predict step variable from previous and current values
    array_1d<double,3>& CurrentStepVariable      = rNode.FastGetSolutionStepValue(*this->mpStepVariable,    0);
    array_1d<double,3>& PreviousStepVariable     = rNode.FastGetSolutionStepValue(*this->mpStepVariable,    1);
     
    array_1d<double,3>& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,        0);
    array_1d<double,3>& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,        1);

    // update step variable previous iteration instead of previous step
    PreviousStepVariable = CurrentStepVariable;
      
    // update delta variable      
    array_1d<double,3> DeltaVariable;
    noalias(DeltaVariable) = CurrentVariable - PreviousVariable;

    // update step variable
    Matrix CayleyDeltaVariable(3,3);
    noalias(CayleyDeltaVariable) = ZeroMatrix(3,3);
    BeamMathUtils<double>::CayleyTransform( DeltaVariable, CayleyDeltaVariable);

    Matrix CayleyStepVariable(3,3);
    noalias(CayleyStepVariable) = ZeroMatrix(3,3);
    BeamMathUtils<double>::CayleyTransform( CurrentStepVariable, CayleyStepVariable);

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
    noalias(CurrentFirstDerivative)  = -prod(CayleyStepVariable,PreviousFirstDerivative);
    noalias(CurrentFirstDerivative) += this->mEmc.c0 * CurrentStepVariable;
    
    // update second derivative
    array_1d<double,3>& CurrentSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative,  0);          

    noalias(CurrentSecondDerivative) = prod(CayleyStepVariable,PreviousFirstDerivative);
    CurrentSecondDerivative *= -this->mEmc.c1;

    noalias(CurrentSecondDerivative) += this->mEmc.c1 * CurrentFirstDerivative;
   
       
    KRATOS_CATCH( "" )
  }

  template<>
  void EmcStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> >::UpdateFirstDerivative(NodeType& rNode)
  {
    KRATOS_TRY

    // predict step variable from previous and current values
    const array_1d<double,3>& CurrentStepVariable      = rNode.FastGetSolutionStepValue(*this->mpStepVariable,    0);

    array_1d<double,3>& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,       0);
    
    const array_1d<double,3>& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, 1);
     
    Matrix CayleyStepVariable(3,3);
    noalias(CayleyStepVariable) = ZeroMatrix(3,3);
    BeamMathUtils<double>::CayleyTransform( CurrentStepVariable, CayleyStepVariable);

    noalias(CurrentFirstDerivative) = -prod(CayleyStepVariable,PreviousFirstDerivative);

    noalias(CurrentFirstDerivative) += this->mEmc.c0 * CurrentStepVariable;
    
    KRATOS_CATCH( "" )
  }

  template<>
  void EmcStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> >::UpdateSecondDerivative(NodeType& rNode)
  {
    KRATOS_TRY

    // predict step variable from previous and current values
    const array_1d<double,3>& CurrentStepVariable      = rNode.FastGetSolutionStepValue(*this->mpStepVariable,    0);
    const array_1d<double,3>& CurrentFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
    const array_1d<double,3>& PreviousFirstDerivative = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, 1);
    
    array_1d<double,3>& CurrentSecondDerivative  = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);
     
    Matrix CayleyStepVariable(3,3);
    noalias(CayleyStepVariable) = ZeroMatrix(3,3);
    BeamMathUtils<double>::CayleyTransform( CurrentStepVariable, CayleyStepVariable);

    noalias(CurrentSecondDerivative) = prod(CayleyStepVariable,PreviousFirstDerivative);
    CurrentSecondDerivative *= -this->mEmc.c1;

    noalias(CurrentSecondDerivative) += this->mEmc.c1 * CurrentFirstDerivative;
    
    KRATOS_CATCH( "" )
  }
  
  
}  // namespace Kratos.


