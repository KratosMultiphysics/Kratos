// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

#include "generalized_newmark_scheme.hpp"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class NewmarkQuasistaticUPwScheme : public GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(NewmarkQuasistaticUPwScheme);

    NewmarkQuasistaticUPwScheme(double beta, double gamma, double theta)
        : GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>(
              {FirstOrderScalarVariable(WATER_PRESSURE, DT_WATER_PRESSURE, DT_PRESSURE_COEFFICIENT)},
              {SecondOrderVectorVariable(DISPLACEMENT), SecondOrderVectorVariable(ROTATION)},
              beta,
              gamma,
              theta)
    {
    }


protected:
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            // static here means no velocities and accelerations for the displacement/rotation D.O.F.
            for (const auto& r_first_order_scalar_variable : this->GetFirstOrderScalarVariables()) {
                this->UpdateScalarTimeDerivative(rNode, r_first_order_scalar_variable.instance,
                                                 r_first_order_scalar_variable.first_time_derivative);
            }
        });

        KRATOS_CATCH("")
    }

}; // Class NewmarkQuasistaticUPwScheme

} // namespace Kratos