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

#if !defined(KRATOS_GEO_APPLY_HYDROSTATIC_PRESSURE_TABLE_PROCESS )
#define  KRATOS_GEO_APPLY_HYDROSTATIC_PRESSURE_TABLE_PROCESS

#include "includes/table.h"

#include "custom_processes/apply_constant_hydrostatic_pressure_process.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyHydrostaticPressureTableProcess : public ApplyConstantHydrostaticPressureProcess
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyHydrostaticPressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyHydrostaticPressureTableProcess(ModelPart& model_part,
                                         Parameters rParameters
                                         ) : ApplyConstantHydrostaticPressureProcess(model_part, rParameters)
    {
        KRATOS_TRY

        unsigned int TableId = rParameters["table"].GetInt();
        mpTable = model_part.pGetTable(TableId);
        mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

        KRATOS_CATCH("")
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyHydrostaticPressureTableProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyHydrostaticPressureTableProcess algorithms.
    void Execute() override
    {
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY

        if (mrModelPart.NumberOfNodes() > 0) {
            const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);
            const double Time = mrModelPart.GetProcessInfo()[TIME]/mTimeUnitConverter;
            const double deltaH = mpTable->GetValue(Time);

            if (mIsSeepage) {
                block_for_each(mrModelPart.Nodes(), [&var, &deltaH, this](Node& rNode) {
                    const double distance = mReferenceCoordinate - rNode.Coordinates()[mGravityDirection];
                    const double pressure = - PORE_PRESSURE_SIGN_FACTOR * mSpecificWeight * (distance + deltaH);

                    if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < 0.0) {
                        rNode.FastGetSolutionStepValue(var) = pressure;
                        if (mIsFixed) rNode.Fix(var);
                    } else {
                        if (mIsFixedProvided) rNode.Free(var);
                    }
                });
            } else {
                block_for_each(mrModelPart.Nodes(), [&var, &deltaH, this](Node& rNode) {
                    const double distance = mReferenceCoordinate - rNode.Coordinates()[mGravityDirection];
                    const double pressure = - PORE_PRESSURE_SIGN_FACTOR * mSpecificWeight * (distance + deltaH);

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
        return "ApplyHydrostaticPressureTableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyHydrostaticPressureTableProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    TableType::Pointer mpTable;
    double mTimeUnitConverter;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ApplyHydrostaticPressureTableProcess& operator=(ApplyHydrostaticPressureTableProcess const& rOther);

    /// Copy constructor.
    //ApplyHydrostaticPressureTableProcess(ApplyHydrostaticPressureTableProcess const& rOther);
}; // Class ApplyHydrostaticPressureTableProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyHydrostaticPressureTableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyHydrostaticPressureTableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


} // namespace Kratos.

#endif /* KRATOS_GEO_APPLY_HYDROSTATIC_PRESSURE_TABLE_PROCESS defined */