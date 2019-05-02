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

#if !defined(KRATOS_DAM_BOFANG_CONDITION_TEMPERATURE_PROCESS)
#define KRATOS_DAM_BOFANG_CONDITION_TEMPERATURE_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

class DamBofangConditionTemperatureProcess : public Process
{

  public:
    KRATOS_CLASS_POINTER_DEFINITION(DamBofangConditionTemperatureProcess);

    typedef Table<double, double> TableType;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamBofangConditionTemperatureProcess(ModelPart &rModelPart,
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
                "Surface_Temp"                                     : 0.0,
                "Bottom_Temp"                                      : 0.0,
                "Height_Dam"                                       : 0.0,
                "Temperature_Amplitude"                            : 0.0,
                "Day_Ambient_Temp"                                 : 1,
                "Water_level"                                      : 0.0,
                "Water_level_Table"                                : 0,
                "Month"                                            : 1.0,
                "Month_Table"                                      : 0
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
        mSurfaceTemp = rParameters["Surface_Temp"].GetDouble();
        mBottomTemp = rParameters["Bottom_Temp"].GetDouble();
        mHeight = rParameters["Height_Dam"].GetDouble();
        mAmplitude = rParameters["Temperature_Amplitude"].GetDouble();
        mDay = rParameters["Day_Ambient_Temp"].GetInt();
        mWaterLevel = rParameters["Water_level"].GetDouble();
        mMonth = rParameters["Month"].GetDouble();
        mFreq = 0.52323;

        mTimeUnitConverter = mrModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];
        mTableIdWater = rParameters["Water_level_Table"].GetInt();
        mTableIdMonth = rParameters["Month_Table"].GetInt();

        if (mTableIdWater != 0)
            mpTableWater = mrModelPart.pGetTable(mTableIdWater);

        if (mTableIdMonth != 0)
            mpTableMonth = mrModelPart.pGetTable(mTableIdMonth);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~DamBofangConditionTemperatureProcess() {}

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
                    double aux1 = ((mBottomTemp - (mSurfaceTemp * exp(-0.04 * mHeight))) / (1 - (exp(-0.04 * mHeight))));
                    double Temperature = (aux1 + ((mSurfaceTemp - aux1) * (exp(-0.04 * aux))) + (mAmplitude * (exp(-0.018 * aux)) * (cos(mFreq * (mMonth - (mDay / 30.0) - 2.15 + (1.30 * exp(-0.085 * aux)))))));

                    it->FastGetSolutionStepValue(var) = Temperature;
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
        if (mTableIdWater != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mWaterLevel = mpTableWater->GetValue(time);
        }

        if (mTableIdMonth != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mMonth = mpTableMonth->GetValue(time);
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
                    double aux1 = ((mBottomTemp - (mSurfaceTemp * exp(-0.04 * mHeight))) / (1 - (exp(-0.04 * mHeight))));
                    double Temperature = (aux1 + ((mSurfaceTemp - aux1) * (exp(-0.04 * aux))) + (mAmplitude * (exp(-0.018 * aux)) * (cos(mFreq * (mMonth - (mDay / 30.0) - 2.15 + (1.30 * exp(-0.085 * aux)))))));

                    it->FastGetSolutionStepValue(var) = Temperature;
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
        return "DamBofangConditionTemperatureProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "DamBofangConditionTemperatureProcess";
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
    double mSurfaceTemp;
    double mBottomTemp;
    double mHeight;
    double mAmplitude;
    int mDay;
    double mMonth;
    double mWaterLevel;
    double mFreq;
    double mTimeUnitConverter;
    TableType::Pointer mpTableWater;
    TableType::Pointer mpTableMonth;
    int mTableIdWater;
    int mTableIdMonth;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  private:
    /// Assignment operator.
    DamBofangConditionTemperatureProcess &operator=(DamBofangConditionTemperatureProcess const &rOther);

}; //Class

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DamBofangConditionTemperatureProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DamBofangConditionTemperatureProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_BOFANG_CONDITION_TEMPERATURE_PROCESS defined */
