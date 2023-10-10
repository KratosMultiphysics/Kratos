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
//                   Jonathan Nuttall

#include "apply_phreatic_multi_line_pressure_table_process.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

ApplyPhreaticMultiLinePressureTableProcess::ApplyPhreaticMultiLinePressureTableProcess(ModelPart& model_part,
Parameters rParameters
) : ApplyConstantPhreaticMultiLinePressureProcess(model_part, rParameters)
{
    KRATOS_TRY

    for (auto value : rParameters["table"].GetVector())
    {
        const auto TableId = static_cast<unsigned int>(value);
        if (TableId > 0)
        {
            auto pTable = model_part.pGetTable(TableId);
            mpTable.push_back(pTable);
        }
        else
        {
            mpTable.push_back(nullptr);
        }
    }

    mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

    KRATOS_CATCH("")
}

void ApplyPhreaticMultiLinePressureTableProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

        if (mrModelPart.NumberOfNodes() <= 0) {
            return;
        }

        const Variable<double> &var = KratosComponents< Variable<double> >::Get(VariableName());

        const double Time = mrModelPart.GetProcessInfo()[TIME]/mTimeUnitConverter;
        std::vector<double> deltaH;
        std::transform(mpTable.begin(), mpTable.end(), std::back_inserter(deltaH),
                       [Time](auto element){return element ? element->GetValue(Time) : 0.0;});

        block_for_each(mrModelPart.Nodes(), [&var, &deltaH, this](Node& rNode) {
            const double pressure = CalculatePressure(rNode, deltaH);
            if (IsSeepage()) {
                if (pressure < PORE_PRESSURE_SIGN_FACTOR * PressureTensionCutOff()) {
                    rNode.FastGetSolutionStepValue(var) = pressure;
                    if (IsFixed()) rNode.Fix(var);
                } else {
                    if (IsFixedProvided()) rNode.Free(var);
                }
            } else {
                if (IsFixed()) rNode.Fix(var);
                else if (IsFixedProvided()) rNode.Free(var);
                rNode.FastGetSolutionStepValue(var) = std::min(pressure, PORE_PRESSURE_SIGN_FACTOR * PressureTensionCutOff());
            }
        });

    KRATOS_CATCH("")
}

std::string ApplyPhreaticMultiLinePressureTableProcess::Info() const
{
    return "ApplyPhreaticMultiLinePressureTableProcess";
}

void ApplyPhreaticMultiLinePressureTableProcess::PrintInfo(std::ostream &rOStream) const
{
    rOStream << "ApplyPhreaticMultiLinePressureTableProcess";
}

}
