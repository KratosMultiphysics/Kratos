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

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_utilities/math_utilities.hpp"
#include "utilities/math_utils.h"
#include "utilities/variable_utils.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class PeriodicInterfaceProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(PeriodicInterfaceProcess);

    PeriodicInterfaceProcess( ModelPart& model_part,
                              Parameters rParameters ) : Process(Flags()) , mrModelPart(model_part)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "dimension": 2,
                "stress_limit": 100.0e6
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mDimension   = rParameters["dimension"].GetInt();
        mStressLimit = rParameters["stress_limit"].GetDouble();

        KRATOS_CATCH("");
    }

    PeriodicInterfaceProcess(const PeriodicInterfaceProcess&) = delete;
    PeriodicInterfaceProcess& operator=(const PeriodicInterfaceProcess&) = delete;
    ~PeriodicInterfaceProcess() override = default;

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY

        block_for_each(mrModelPart.Conditions(), [&](Condition& rCondition) {
            Condition::GeometryType& rGeom = rCondition.GetGeometry();

            rCondition.Set(PERIODIC,true);

            rGeom[0].FastGetSolutionStepValue(PERIODIC_PAIR_INDEX) = rGeom[1].Id();
            rGeom[1].FastGetSolutionStepValue(PERIODIC_PAIR_INDEX) = rGeom[0].Id();
        });
        VariableUtils{}.SetFlag(ACTIVE, false, mrModelPart.Elements());

        KRATOS_CATCH("")
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
        KRATOS_TRY

        if (mrModelPart.NumberOfConditions() > 0) {
            SizeType StressTensorSize = STRESS_TENSOR_SIZE_2D;
            if (mDimension == N_DIM_3D) StressTensorSize = STRESS_TENSOR_SIZE_3D;

            block_for_each(mrModelPart.Conditions(), [&StressTensorSize](Condition& rCondition) {
                Condition::GeometryType& rGeom = rCondition.GetGeometry();

                Matrix NodalStressMatrix(StressTensorSize, StressTensorSize);
                noalias(NodalStressMatrix) = 0.5 * ( rGeom[0].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR)
                                                   + rGeom[1].FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) );
                Vector PrincipalStresses(StressTensorSize);
                noalias(PrincipalStresses) = SolidMechanicsMathUtilities<double>::EigenValuesDirectMethod(NodalStressMatrix);

                // Check whether the principal stress S1 at the node is higher than the prescribed limit to activate the joints
                if (PrincipalStresses[0] >= mStressLimit) {
                    rCondition.Set(PERIODIC,false);
                    rGeom[0].FastGetSolutionStepValue(PERIODIC_PAIR_INDEX) = 0;
                    rGeom[1].FastGetSolutionStepValue(PERIODIC_PAIR_INDEX) = 0;

                    GlobalPointersVector<Element>& rE = rGeom[0].GetValue(NEIGHBOUR_ELEMENTS);
                    for (unsigned int ie = 0; ie < rE.size(); ie++) {
                        #pragma omp critical
                        {
                            rE[ie].Set(ACTIVE,true);
                        }
                    }
                }
            });
        }

        KRATOS_CATCH("")
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "PeriodicInterfaceProcess";
    }

private:
    ModelPart& mrModelPart;
    int mDimension;
    double mStressLimit;

}

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  PeriodicInterfaceProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PeriodicInterfaceProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}