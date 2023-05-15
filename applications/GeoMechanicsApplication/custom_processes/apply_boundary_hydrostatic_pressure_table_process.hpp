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

#if !defined(KRATOS_GEO_APPLY_BOUNDARY_HYDROSTATIC_PRESSURE_TABLE_PROCESS )
#define  KRATOS_GEO_APPLY_BOUNDARY_HYDROSTATIC_PRESSURE_TABLE_PROCESS

#include "includes/table.h"

#include "custom_processes/apply_constant_boundary_hydrostatic_pressure_process.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyBoundaryHydrostaticPressureTableProcess : public ApplyConstantBoundaryHydrostaticPressureProcess
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyBoundaryHydrostaticPressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyBoundaryHydrostaticPressureTableProcess(ModelPart& model_part,
                                                 Parameters rParameters
                                                 ) : ApplyConstantBoundaryHydrostaticPressureProcess(model_part, rParameters)
    {
        KRATOS_TRY

        unsigned int TableId = rParameters["table"].GetInt();
        mpTable = model_part.pGetTable(TableId);
        mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];

        KRATOS_CATCH("")
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyBoundaryHydrostaticPressureTableProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyBoundaryHydrostaticPressureTableProcess algorithms.
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

            block_for_each(mrModelPart.Nodes(), [&deltaH, &var, this](Node& rNode){
                const double distance = mReferenceCoordinate - rNode.Coordinates()[mGravityDirection];
                const double pressure = mSpecificWeight * (distance + deltaH);

                if (pressure > 0.0) {
                    rNode.FastGetSolutionStepValue(var) = pressure;
                } else {
                    rNode.FastGetSolutionStepValue(var) = 0.0;
                }
            });
        }

        KRATOS_CATCH("")
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyBoundaryHydrostaticPressureTableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyBoundaryHydrostaticPressureTableProcess";
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
    ApplyBoundaryHydrostaticPressureTableProcess& operator=(ApplyBoundaryHydrostaticPressureTableProcess const& rOther);

    /// Copy constructor.
    //ApplyBoundaryHydrostaticPressureTableProcess(ApplyBoundaryHydrostaticPressureTableProcess const& rOther);
}; // Class ApplyBoundaryHydrostaticPressureTableProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyBoundaryHydrostaticPressureTableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyBoundaryHydrostaticPressureTableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


} // namespace Kratos.

#endif /* KRATOS_GEO_APPLY_BOUNDARY_BOUNDARY_HYDROSTATIC_PRESSURE_TABLE_PROCESS defined */