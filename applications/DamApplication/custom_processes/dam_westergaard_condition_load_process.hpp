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
#if !defined(KRATOS_DAM_WESTERGAARD_CONDITION_LOAD_PROCESS)
#define KRATOS_DAM_WESTERGAARD_CONDITION_LOAD_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"

namespace Kratos
{

class DamWestergaardConditionLoadProcess : public Process
{

  public:
    KRATOS_CLASS_POINTER_DEFINITION(DamWestergaardConditionLoadProcess);

    typedef Table<double, double> TableType;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamWestergaardConditionLoadProcess(ModelPart &rModelPart,
                                       Parameters &rParameters) : Process(Flags()), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "Modify"                                                : true,
                "Gravity_Direction"                                     : "Y",
                "Reservoir_Bottom_Coordinate_in_Gravity_Direction"      : 0.0,
                "Spe_weight"                                            : 0.0,
                "Water_level"                                           : 0.0,
                "Water_Table"                                           : 0,
                "Aceleration"                                           : 0.0,
                "Aceleration_Table"                                     : 0
            }  )");

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mGravityDirection = rParameters["Gravity_Direction"].GetString();
        mReferenceCoordinate = rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"].GetDouble();
        mSpecific = rParameters["Spe_weight"].GetDouble();
        mWaterLevel = rParameters["Water_level"].GetDouble();
        mAcceleration = rParameters["Aceleration"].GetDouble();

        mTimeUnitConverter = mrModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];
        mTableIdWater = rParameters["Water_Table"].GetInt();
        mTableIdAcceleration = rParameters["Aceleration_Table"].GetInt();

        if (mTableIdWater != 0)
            mpTableWater = mrModelPart.pGetTable(mTableIdWater);

        if (mTableIdAcceleration != 0)
            mpTableAcceleration = mrModelPart.pGetTable(mTableIdAcceleration);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~DamWestergaardConditionLoadProcess() {}

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
        double pressure;

        if (mGravityDirection == "X")
            direction = 0;
        else if (mGravityDirection == "Y")
            direction = 1;
        else
            direction = 2;

        double ref_coord = mReferenceCoordinate + mWaterLevel;
        double unit_acceleration = mAcceleration / 9.81;

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

#pragma omp parallel for
            for (int i = 0; i < nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                double y_water = ref_coord - (it->Coordinates()[direction]);

                if (y_water < 0.0)
                {
                    y_water = 0.0;
                }

                pressure = (mSpecific * (y_water)) + 0.875 * (unit_acceleration)*mSpecific * sqrt(y_water * mWaterLevel);
                it->FastGetSolutionStepValue(var) = pressure;
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
        if (mTableIdWater != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mWaterLevel = mpTableWater->GetValue(time);
        }

        if (mTableIdAcceleration != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mAcceleration = mpTableAcceleration->GetValue(time);
        }

        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();
        int direction;
        double pressure;

        if (mGravityDirection == "X")
            direction = 0;
        else if (mGravityDirection == "Y")
            direction = 1;
        else
            direction = 2;

        double ref_coord = mReferenceCoordinate + mWaterLevel;
        double unit_acceleration = mAcceleration / 9.81;

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

#pragma omp parallel for
            for (int i = 0; i < nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                double y_water = ref_coord - (it->Coordinates()[direction]);

                if (y_water < 0.0)
                {
                    y_water = 0.0;
                }

                // Hydrodynamics Westergaard effects just contribute when the acceleration goes in the upstream direction
                if (unit_acceleration < 0.0)
                {
                    pressure = (mSpecific * (y_water)) + 0.875 * (-1.0 * unit_acceleration) * mSpecific * sqrt(y_water * mWaterLevel);
                }
                else
                {
                    pressure = (mSpecific * (y_water));
                }

                if (pressure > 0.0)
                {
                    it->FastGetSolutionStepValue(var) = pressure;
                }
                else
                {
                    it->FastGetSolutionStepValue(var) = 0.0;
                }
            }
        }

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DamWestergaardConditionLoadProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "DamWestergaardConditionLoadProcess";
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
    double mReferenceCoordinate;
    double mSpecific;
    double mWaterLevel;
    double mAcceleration;
    double mTimeUnitConverter;
    TableType::Pointer mpTableWater;
    TableType::Pointer mpTableAcceleration;
    int mTableIdWater;
    int mTableIdAcceleration;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  private:
    /// Assignment operator.
    DamWestergaardConditionLoadProcess &operator=(DamWestergaardConditionLoadProcess const &rOther);

}; //Class

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DamWestergaardConditionLoadProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DamWestergaardConditionLoadProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_WESTERGAARD_CONDITION_LOAD_PROCESS defined */
