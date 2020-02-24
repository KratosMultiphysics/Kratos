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

#if !defined(KRATOS_DAM_RESERVOIR_CONSTANT_TEMPERATURE_PROCESS)
#define KRATOS_DAM_RESERVOIR_CONSTANT_TEMPERATURE_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"

namespace Kratos
{

class DamReservoirConstantTemperatureProcess : public Process
{

  public:
    KRATOS_CLASS_POINTER_DEFINITION(DamReservoirConstantTemperatureProcess);

    typedef Table<double, double> TableType;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamReservoirConstantTemperatureProcess(ModelPart &rModelPart,
                                           Parameters &rParameters) : Process(Flags()), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed"                                         : false,
                "Gravity_Direction"                                : "Y",
                "Reservoir_Bottom_Coordinate_in_Gravity_Direction" : 0.0,
                "Water_temp"                                       : 0.0,
                "Water_temp_Table"                                 : 0,
                "Water_level"                                      : 0.0,
                "Water_level_Table"                                : 0
            }  )");

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mIsFixed = rParameters["is_fixed"].GetBool();
        mGravityDirection = rParameters["Gravity_Direction"].GetString();
        mReferenceCoordinate = rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"].GetDouble();
        mWaterTemp = rParameters["Water_temp"].GetDouble();
        mWaterLevel = rParameters["Water_level"].GetDouble();

        mTimeUnitConverter = mrModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];
        mTableIdWaterTemp = rParameters["Water_temp_Table"].GetInt();
        mTableIdWater = rParameters["Water_level_Table"].GetInt();

        if (mTableIdWaterTemp != 0)
            mpTableWaterTemp = mrModelPart.pGetTable(mTableIdWaterTemp);

        if (mTableIdWater != 0)
            mpTableWater = mrModelPart.pGetTable(mTableIdWater);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~DamReservoirConstantTemperatureProcess() {}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute() override
    {
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitialize() override
    {

        KRATOS_TRY;

        Variable<double> var = KratosComponents<Variable<double>>::Get(mVariableName);
        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();
        int direction;

        if (mGravityDirection == "X")
            direction = 0;
        else if (mGravityDirection == "Y")
            direction = 1;
        else
            direction = 2;

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

            #pragma omp parallel for
            for (int i = 0; i < nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                double aux = (mReferenceCoordinate + mWaterLevel) - it->Coordinates()[direction];
                if (aux >= 0.0)
                {
                    if (mIsFixed)
                    {
                        it->Fix(var);
                    }
                    it->FastGetSolutionStepValue(var) = mWaterTemp;
                }
            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep() override
    {

        KRATOS_TRY;

        Variable<double> var = KratosComponents<Variable<double>>::Get(mVariableName);

        // Getting the values of table in case that it exist

        if (mTableIdWaterTemp != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mWaterTemp = mpTableWaterTemp->GetValue(time);
        }

        if (mTableIdWater != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mWaterLevel = mpTableWater->GetValue(time);
        }

        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();
        int direction;

        if (mGravityDirection == "X")
            direction = 0;
        else if (mGravityDirection == "Y")
            direction = 1;
        else
            direction = 2;

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

            #pragma omp parallel for
            for (int i = 0; i < nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                double aux = (mReferenceCoordinate + mWaterLevel) - it->Coordinates()[direction];
                if (aux >= 0.0)
                {
                    if (mIsFixed)
                    {
                        it->Fix(var);
                    }
                    it->FastGetSolutionStepValue(var) = mWaterTemp;
                }
            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteFinalizeSolutionStep() override
    {

        KRATOS_TRY;

        Variable<double> var = KratosComponents<Variable<double>>::Get(mVariableName);

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
        return "DamReservoirConstantTemperatureProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "DamReservoirConstantTemperatureProcess";
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
    double mReferenceCoordinate;
    double mWaterTemp;
    double mWaterLevel;
    double mTimeUnitConverter;
    TableType::Pointer mpTableWaterTemp;
    TableType::Pointer mpTableWater;
    int mTableIdWaterTemp;
    int mTableIdWater;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  private:
    /// Assignment operator.
    DamReservoirConstantTemperatureProcess &operator=(DamReservoirConstantTemperatureProcess const &rOther);

}; //Class

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DamReservoirConstantTemperatureProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DamReservoirConstantTemperatureProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_RESERVOIR_CONSTANT_TEMPERATURE_PROCESS defined */
