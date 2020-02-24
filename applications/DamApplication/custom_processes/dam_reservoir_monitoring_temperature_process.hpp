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

#if !defined(KRATOS_DAM_RESERVOIR_MONITORING_TEMPERATURE_PROCESS)
#define KRATOS_DAM_RESERVOIR_MONITORING_TEMPERATURE_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

class DamReservoirMonitoringTemperatureProcess : public Process
{

  public:
    KRATOS_CLASS_POINTER_DEFINITION(DamReservoirMonitoringTemperatureProcess);

    typedef Table<double, double> TableType;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamReservoirMonitoringTemperatureProcess(ModelPart &rModelPart,
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
                "Height_Dam"                                       : 0.0,
                "Ambient_temp"                                     : 0.0,
                "Ambient_temp_Table"                               : 0,
                "Water_level"                                      : 0.0,
                "Water_level_Table"                                : 0,
                "Z_Coord_1"                                        : 0.0,
                "Water_temp_1"                                     : 0.0,
                "Water_temp_Table_1"                               : 0,
                "Z_Coord_2"                                        : 0.0,
                "Water_temp_2"                                     : 0.0,
                "Water_temp_Table_2"                               : 0,
                "Z_Coord_3"                                        : 0.0,
                "Water_temp_3"                                     : 0.0,
                "Water_temp_Table_3"                               : 0
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
        mHeight = rParameters["Height_Dam"].GetDouble();
        mAmbientTemp = rParameters["Ambient_temp"].GetDouble();
        mWaterLevel = rParameters["Water_level"].GetDouble();
        mZCoord1 = rParameters["Z_Coord_1"].GetDouble();
        mWaterTemp1 = rParameters["Water_temp_1"].GetDouble();
        mZCoord2 = rParameters["Z_Coord_2"].GetDouble();
        mWaterTemp2 = rParameters["Water_temp_2"].GetDouble();
        mZCoord3 = rParameters["Z_Coord_3"].GetDouble();
        mWaterTemp3 = rParameters["Water_temp_3"].GetDouble();

        mTimeUnitConverter = mrModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];
        mTableIdWater = rParameters["Water_level_Table"].GetInt();
        mTableIdAmbientTemp = rParameters["Ambient_temp_Table"].GetInt();
        mTableIdWaterTemp1 = rParameters["Water_temp_Table_1"].GetInt();
        mTableIdWaterTemp2 = rParameters["Water_temp_Table_2"].GetInt();
        mTableIdWaterTemp3 = rParameters["Water_temp_Table_3"].GetInt();

// ELIMINAR        mTableIdMonth = rParameters["Month_Table"].GetInt();

        if (mTableIdWater != 0)
            mpTableWater = mrModelPart.pGetTable(mTableIdWater);

        if (mTableIdAmbientTemp != 0)
            mpTableAmbientTemp = mrModelPart.pGetTable(mTableIdAmbientTemp);

        if (mTableIdWaterTemp1 != 0)
            mpTableWaterTemp1 = mrModelPart.pGetTable(mTableIdWaterTemp1);

        if (mTableIdWaterTemp2 != 0)
            mpTableWaterTemp2 = mrModelPart.pGetTable(mTableIdWaterTemp2);

        if (mTableIdWaterTemp3 != 0)
            mpTableWaterTemp3 = mrModelPart.pGetTable(mTableIdWaterTemp3);



        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~DamReservoirMonitoringTemperatureProcess() {}

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

                double aux = it->Coordinates()[direction];
                if (aux < (mReferenceCoordinate + mWaterLevel))
                {
                    if (mIsFixed)
                    {
                        it->Fix(var);
                    }
                    if (aux > mZCoord1)
                    {
                        double Temperature = ((mAmbientTemp - mWaterTemp1)/((mReferenceCoordinate + mWaterLevel) - mZCoord1)) * (aux - mZCoord1) + mWaterTemp1;
                        it->FastGetSolutionStepValue(var) = Temperature;
                    }
                    else if ((aux <= mZCoord1) && (aux > mZCoord2))
                    {
                        double Temperature = ((mWaterTemp1 - mWaterTemp2)/(mZCoord1 - mZCoord2)) * (aux - mZCoord2) + mWaterTemp2;
                        it->FastGetSolutionStepValue(var) = Temperature;
                    }
                    else if ((aux <= mZCoord2) && (aux > mZCoord3))
                    {
                        double Temperature = ((mWaterTemp2 - mWaterTemp3)/(mZCoord2 - mZCoord3)) * (aux - mZCoord3) + mWaterTemp3;
                        it->FastGetSolutionStepValue(var) = Temperature;
                    }
                    else if (aux <= mZCoord3)
                    {
                        double Temperature = ((mWaterTemp3 - mWaterTemp2)/(mZCoord3 - mZCoord2)) * (aux - mZCoord3) + mWaterTemp3;
                        it->FastGetSolutionStepValue(var) = Temperature;
                    }
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

        if (mTableIdAmbientTemp != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mAmbientTemp = mpTableAmbientTemp->GetValue(time);
        }

        if (mTableIdWaterTemp1 != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mWaterTemp1 = mpTableWaterTemp1->GetValue(time);
        }

        if (mTableIdWaterTemp2 != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mWaterTemp2 = mpTableWaterTemp2->GetValue(time);
        }

        if (mTableIdWaterTemp3 != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mWaterTemp3 = mpTableWaterTemp3->GetValue(time);
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

                double aux = it->Coordinates()[direction];
                if (aux < (mReferenceCoordinate + mWaterLevel))
                {
                    if (mIsFixed)
                    {
                        it->Fix(var);
                    }
                    if (aux > mZCoord1)
                    {
                        double Temperature = ((mAmbientTemp - mWaterTemp1)/((mReferenceCoordinate + mWaterLevel) - mZCoord1)) * (aux - mZCoord1) + mWaterTemp1;
                        it->FastGetSolutionStepValue(var) = Temperature;
                    }
                    else if ((aux <= mZCoord1) && (aux > mZCoord2))
                    {
                        double Temperature = ((mWaterTemp1 - mWaterTemp2)/(mZCoord1 - mZCoord2)) * (aux - mZCoord2) + mWaterTemp2;
                        it->FastGetSolutionStepValue(var) = Temperature;
                    }
                    else if ((aux <= mZCoord2) && (aux > mZCoord3))
                    {
                        double Temperature = ((mWaterTemp2 - mWaterTemp3)/(mZCoord2 - mZCoord3)) * (aux - mZCoord3) + mWaterTemp3;
                        it->FastGetSolutionStepValue(var) = Temperature;
                    }
                    else if (aux <= mZCoord3)
                    {
                        double Temperature = ((mWaterTemp3 - mWaterTemp2)/(mZCoord3 - mZCoord2)) * (aux - mZCoord3) + mWaterTemp3;
                        it->FastGetSolutionStepValue(var) = Temperature;
                    }
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
        return "DamReservoirMonitoringTemperatureProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "DamReservoirMonitoringTemperatureProcess";
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
    double mHeight;
    double mAmbientTemp;
    double mWaterLevel;
    double mZCoord1;
    double mWaterTemp1;
    double mZCoord2;
    double mWaterTemp2;
    double mZCoord3;
    double mWaterTemp3;
    double mTimeUnitConverter;
    TableType::Pointer mpTableWater;
    TableType::Pointer mpTableAmbientTemp;
    TableType::Pointer mpTableWaterTemp1;
    TableType::Pointer mpTableWaterTemp2;
    TableType::Pointer mpTableWaterTemp3;
    int mTableIdWater;
    int mTableIdAmbientTemp;
    int mTableIdWaterTemp1;
    int mTableIdWaterTemp2;
    int mTableIdWaterTemp3;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  private:
    /// Assignment operator.
    DamReservoirMonitoringTemperatureProcess &operator=(DamReservoirMonitoringTemperatureProcess const &rOther);

}; //Class

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DamReservoirMonitoringTemperatureProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DamReservoirMonitoringTemperatureProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_RESERVOIR_MONITORING_TEMPERATURE_PROCESS defined */
