// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#if !defined(KRATOS_GEO_APPLY_TIME_DEPENDENT_INTERPOLATE_LINE_PRESSURE_PROCESS)
#define  KRATOS_GEO_APPLY_TIME_DEPENDENT_INTERPOLATE_LINE_PRESSURE_PROCESS

#include <algorithm>
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"
#include "custom_processes/apply_constant_interpolate_line_pressure_process.hpp"

namespace Kratos
{

class ApplyTimeDependentInterpolateLinePressureProcess :
    public ApplyConstantInterpolateLinePressureProcess
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyTimeDependentInterpolateLinePressureProcess);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyTimeDependentInterpolateLinePressureProcess(ModelPart& model_part,
                                                    Parameters rParameters
                                                    ) : ApplyConstantInterpolateLinePressureProcess(model_part, rParameters)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyTimeDependentInterpolateLinePressureProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyTimeDependentInterpolateLinePressureProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY

        if (mrModelPart.NumberOfNodes() > 0) {
            const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);

            if (mIsSeepage) {
                block_for_each(mrModelPart.Nodes(), [&var, this](Node& rNode) {
                    const double pressure = CalculatePressure(rNode);

                    if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < 0.0) {
                        rNode.FastGetSolutionStepValue(var) = pressure;
                        if (mIsFixed) rNode.Fix(var);
                    } else {
                        if (mIsFixedProvided) rNode.Free(var);
                    }
                });
            } else {
                block_for_each(mrModelPart.Nodes(), [&var, this](Node& rNode) {
                    if (mIsFixed) rNode.Fix(var);
                    else if (mIsFixedProvided) rNode.Free(var);

                    const double pressure = CalculatePressure(rNode);

                    if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < mPressureTensionCutOff) {
                        rNode.FastGetSolutionStepValue(var) = pressure;
                    } else {
                        rNode.FastGetSolutionStepValue(var) = mPressureTensionCutOff;
                    }
                });
            }
        }

        KRATOS_CATCH("")
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyTimeDependentInterpolateLinePressureProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyTimeDependentInterpolateLinePressureProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ApplyTimeDependentInterpolateLinePressureProcess& operator=(ApplyTimeDependentInterpolateLinePressureProcess const& rOther);

    /// Copy constructor.
    //ApplyTimeDependentInterpolateLinePressureProcess(ApplyTimeDependentInterpolateLinePressureProcess const& rOther);

}; // Class ApplyTimeDependentInterpolateLinePressureProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyTimeDependentInterpolateLinePressureProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyTimeDependentInterpolateLinePressureProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_GEO_APPLY_TIME_DEPENDENT_INTERPOLATE_LINE_PRESSURE_PROCESS defined */
