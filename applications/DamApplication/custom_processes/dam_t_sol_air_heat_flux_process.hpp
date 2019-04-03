//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//
//

#if !defined(KRATOS_DAM_T_SOL_AIR_HEAT_FLUX_PROCESS)
#define KRATOS_DAM_T_SOL_AIR_HEAT_FLUX_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"

namespace Kratos
{

class DamTSolAirHeatFluxProcess : public Process
{

  public:
    KRATOS_CLASS_POINTER_DEFINITION(DamTSolAirHeatFluxProcess);

    typedef Table<double, double> TableType;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamTSolAirHeatFluxProcess(ModelPart &rModelPart,
                              Parameters &rParameters) : Process(Flags()), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "h_0"                             : 0.0,
                "ambient_temperature"             : 0.0,
                "table_ambient_temperature"       : 0,
                "emisivity"                       : 0.0,
                "delta_R"                         : 0.0,
                "absorption_index"                : 0.0,
                "total_insolation"                : 0.0
            }  )");

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["h_0"];
        rParameters["delta_R"];
        rParameters["absorption_index"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mH0 = rParameters["h_0"].GetDouble();
        mAmbientTemperature = rParameters["ambient_temperature"].GetDouble();
        mEmisivity = rParameters["emisivity"].GetDouble();
        mDeltaR = rParameters["delta_R"].GetDouble();
        mAbsorption_index = rParameters["absorption_index"].GetDouble();
        mTotalInsolation = rParameters["total_insolation"].GetDouble();

        mTimeUnitConverter = mrModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];
        mTableId = rParameters["table_ambient_temperature"].GetInt();

        if (mTableId != 0)
            mpTable = mrModelPart.pGetTable(mTableId);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~DamTSolAirHeatFluxProcess() {}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute() override
    {
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


    void ExecuteInitialize() override
    {

        KRATOS_TRY;

        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();
        Variable<double> var = KratosComponents<Variable<double>>::Get(mVariableName);

        // Computing the t_soil_air according to t_sol_air criteria
        double t_sol_air = mAmbientTemperature + (mAbsorption_index * mTotalInsolation / mH0) - (mEmisivity * mDeltaR / mH0);

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

            #pragma omp parallel for
            for (int i = 0; i < nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                const double temp_current = it->FastGetSolutionStepValue(TEMPERATURE);
                const double heat_flux = mH0 * (t_sol_air - temp_current);
                it->FastGetSolutionStepValue(var) = heat_flux;
            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep() override
    {

        KRATOS_TRY;

        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();
        Variable<double> var = KratosComponents<Variable<double>>::Get(mVariableName);

        // Getting the values of table in case that it exist
        if (mTableId != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mAmbientTemperature = mpTable->GetValue(time);
        }

        // Computing the t_soil_air according to t_sol_air criteria
        double t_sol_air = mAmbientTemperature + (mAbsorption_index * mTotalInsolation / mH0) - (mEmisivity * mDeltaR / mH0);

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

            #pragma omp parallel for
            for (int i = 0; i < nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                const double temp_current = it->FastGetSolutionStepValue(TEMPERATURE);
                const double heat_flux = mH0 * (t_sol_air - temp_current);
                it->FastGetSolutionStepValue(var) = heat_flux;
            }
        }

        KRATOS_CATCH("");
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DamTSolAirHeatFluxProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "DamTSolAirHeatFluxProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  protected:
    /// Member Variables
    ModelPart &mrModelPart;
    std::string mVariableName;
    double mH0;
    double mAmbientTemperature;
    double mEmisivity;
    double mDeltaR;
    double mAbsorption_index;
    double mTotalInsolation;
    double mTimeUnitConverter;
    TableType::Pointer mpTable;
    int mTableId;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  private:
    /// Assignment operator.
    DamTSolAirHeatFluxProcess &operator=(DamTSolAirHeatFluxProcess const &rOther);

}; //Class

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DamTSolAirHeatFluxProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DamTSolAirHeatFluxProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_T_SOL_AIR_HEAT_FLUX_PROCESS defined */
