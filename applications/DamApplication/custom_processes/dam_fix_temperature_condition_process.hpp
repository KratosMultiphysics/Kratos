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

#if !defined(KRATOS_DAM_FIX_TEMPERATURE_CONDITION_PROCESS)
#define KRATOS_DAM_FIX_TEMPERATURE_CONDITION_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

class DamFixTemperatureConditionProcess : public Process
{

  public:
    KRATOS_CLASS_POINTER_DEFINITION(DamFixTemperatureConditionProcess);

    typedef Table<double, double> TableType;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamFixTemperatureConditionProcess(ModelPart &rModelPart,
                                      Parameters &rParameters) : Process(Flags()), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters(R"(
            {
                "model_part_name" : "PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name"   : "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed"        : false,
                "value"           : 0.0,
                "table"           : 0,
                "interval":[
                0.0,
                0.0
                ]
            }  )");

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mIsFixed = rParameters["is_fixed"].GetBool();
        mTemperature = rParameters["value"].GetDouble();

        mTimeUnitConverter = mrModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];
        mTableId = rParameters["table"].GetInt();

        if (mTableId != 0)
            mpTable = mrModelPart.pGetTable(mTableId);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~DamFixTemperatureConditionProcess() {}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    void Execute() override
    {
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitialize() override
    {

        KRATOS_TRY;

        const Variable<double>& var = KratosComponents<Variable<double>>::Get(mVariableName);
        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

            #pragma omp parallel for
            for (int i = 0; i < nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                if (mIsFixed)
                {
                    it->Fix(var);
                }

                it->FastGetSolutionStepValue(var) = mTemperature;
            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep() override
    {

        KRATOS_TRY;

        const Variable<double>& var = KratosComponents<Variable<double>>::Get(mVariableName);

        // Getting the values of table in case that it exist
        if (mTableId != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mTemperature = mpTable->GetValue(time);
        }

        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

            #pragma omp parallel for
            for (int i = 0; i < nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                if (mIsFixed)
                {
                    it->Fix(var);
                }

                it->FastGetSolutionStepValue(var) = mTemperature;
            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteFinalizeSolutionStep() override
    {

        KRATOS_TRY;

        const Variable<double>& var = KratosComponents<Variable<double>>::Get(mVariableName);

        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();

        if (nnodes != 0)
        {

            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

            #pragma omp parallel for
            for (int i = 0; i < nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->Free(var);
            }
        }

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "FixTemperatureConditionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "FixTemperatureConditionProcess";
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
    std::string mGravityDirection;
    bool mIsFixed;
    double mTemperature;
    double mTimeUnitConverter;
    TableType::Pointer mpTable;
    int mTableId;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  private:
    /// Assignment operator.
    DamFixTemperatureConditionProcess &operator=(DamFixTemperatureConditionProcess const &rOther);

}; //Class

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DamFixTemperatureConditionProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DamFixTemperatureConditionProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_FIX_TEMPERATURE_CONDITION_PROCESS defined */
