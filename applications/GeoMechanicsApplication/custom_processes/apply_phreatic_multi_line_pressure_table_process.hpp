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

#if !defined(KRATOS_GEO_APPLY_PHREATIC_MULTI_LINE_PRESSURE_TABLE_PROCESS )
#define  KRATOS_GEO_APPLY_PHREATIC_MULTI_LINE_PRESSURE_TABLE_PROCESS

#include "includes/table.h"

#include "custom_processes/apply_constant_phreatic_multi_line_pressure_process.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyPhreaticMultiLinePressureTableProcess : public ApplyConstantPhreaticMultiLinePressureProcess
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyPhreaticMultiLinePressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    using TableType = Table<double, double>;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyPhreaticMultiLinePressureTableProcess(ModelPart& model_part,
                                         Parameters rParameters
                                         ) : ApplyConstantPhreaticMultiLinePressureProcess(model_part, rParameters)
    {
        KRATOS_TRY

        for (auto value : rParameters["table"].GetVector())
        {
            const auto TableId = static_cast<unsigned int>(value);
            if (TableId > 0) {
                auto pTable = model_part.pGetTable(TableId);
                mpTable.push_back(pTable);
            } else {
                mpTable.push_back(nullptr);
            }
        }

        mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

        KRATOS_CATCH("")
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyPhreaticMultiLinePressureTableProcess() override = default;

    /// Assignment operator.
    ApplyPhreaticMultiLinePressureTableProcess& operator=(ApplyPhreaticMultiLinePressureTableProcess const&) = delete;

    /// Copy constructor.
    ApplyPhreaticMultiLinePressureTableProcess(ApplyPhreaticMultiLinePressureTableProcess const&) = delete;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
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

        if (IsSeepage()) {
            block_for_each(mrModelPart.Nodes(), [&var, &deltaH, this](Node& rNode) {
                const double pressure = CalculateTimeDependentPressure(rNode, deltaH);

                if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < 0) {
                    rNode.FastGetSolutionStepValue(var) = pressure;
                    if (IsFixed()) rNode.Fix(var);
                } else {
                    rNode.Free(var);
                }
            });
        } else {
            block_for_each(mrModelPart.Nodes(), [&var, &deltaH, this](Node& rNode) {
                const double pressure = CalculateTimeDependentPressure(rNode, deltaH);

                if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < PressureTensionCutOff()) {
                    rNode.FastGetSolutionStepValue(var) = pressure;
                } else {
                    rNode.FastGetSolutionStepValue(var) = PressureTensionCutOff();
                }
            });
        }

        KRATOS_CATCH("")
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyPhreaticMultiLinePressureTableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyPhreaticMultiLinePressureTableProcess";
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    double CalculateTimeDependentPressure(const Node &rNode, std::vector<double> &deltaH) const
    {

        double height = 0.0;
        int firstPointIndex;
        int secondPointIndex;
        double slope;
		array_1d<double, 2> y;

        // find nodes in horizontalDirectionCoordinates

        firstPointIndex = findIndex(rNode);
        secondPointIndex = firstPointIndex + 1;

        y[0] = deltaH[firstPointIndex] + GravityDirectionCoordinates()[firstPointIndex];
        y[1] = deltaH[secondPointIndex] + GravityDirectionCoordinates()[secondPointIndex];

        slope = (y[1] - y[0])
            / (HorizontalDirectionCoordinates()[secondPointIndex] - HorizontalDirectionCoordinates()[firstPointIndex]);


        height = slope * (rNode.Coordinates()[HorizontalDirection()] - HorizontalDirectionCoordinates()[firstPointIndex]) + y[0];

        const double distance = height - rNode.Coordinates()[GravityDirection()];
        const double pressure = - PORE_PRESSURE_SIGN_FACTOR * SpecificWeight() * distance;
        return pressure;
    }

private:
    std::vector<TableType::Pointer> mpTable;
    double mTimeUnitConverter;
}; // Class ApplyPhreaticMultiLinePressureTableProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyPhreaticMultiLinePressureTableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyPhreaticMultiLinePressureTableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


} // namespace Kratos.

#endif /* KRATOS_GEO_APPLY_PHREATIC_MULTI_LINE_PRESSURE_TABLE_PROCESS defined */