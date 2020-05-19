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

#if !defined(KRATOS_DAM_TEMPERATURE_BY_DEVICE_PROCESS)
#define KRATOS_DAM_TEMPERATURE_BY_DEVICE_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"

namespace Kratos
{

class DamTemperaturebyDeviceProcess : public Process
{

  public:
    KRATOS_CLASS_POINTER_DEFINITION(DamTemperaturebyDeviceProcess);

    typedef Table<double, double> TableType;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamTemperaturebyDeviceProcess(ModelPart &rModelPart,
                                  Parameters &rParameters) : Process(Flags()), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed"          : false,
                "value"             : 0.0,
                "table"             : 0,
                "position"          : [0.0,0.0,0.0],
                "interval":[
                0.0,
                0.0
                ]
            }  )");

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["value"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mIsFixed = rParameters["is_fixed"].GetBool();
        mValue = rParameters["value"].GetDouble();

        unsigned int Dim = rModelPart.GetProcessInfo()[DOMAIN_SIZE];

        // Getting the values of the device coordinates
        mDeviceCoordinates.resize(Dim, false);
        mDeviceCoordinates[0] = rParameters["position"][0].GetDouble();
        mDeviceCoordinates[1] = rParameters["position"][1].GetDouble();
        if (Dim ==3) mDeviceCoordinates[2] = rParameters["position"][2].GetDouble();

        mTimeUnitConverter = mrModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];
        mTableId = rParameters["table"].GetInt();

        if (mTableId != 0)
            mpTable = mrModelPart.pGetTable(mTableId);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~DamTemperaturebyDeviceProcess() {}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep() override
    {

        KRATOS_TRY;

        const int nelements = mrModelPart.GetMesh(0).Elements().size();
        const Variable<double>& var = KratosComponents<Variable<double>>::Get(mVariableName);
        bool IsInside = false;
        array_1d<double, 3> LocalCoordinates;
        Element::Pointer pSelectedElement;

        // Getting the values of table in case that it exist
        if (mTableId != 0)
        {
            double time = mrModelPart.GetProcessInfo()[TIME];
            time = time / mTimeUnitConverter;
            mValue = mpTable->GetValue(time);
        }

        if (nelements != 0)
        {
            ModelPart::ElementsContainerType::iterator el_begin = mrModelPart.ElementsBegin();
            int PointsNumber = 0;

            for (int k = 0; k < nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;
                pSelectedElement = (*(it.base()));
                IsInside = pSelectedElement->GetGeometry().IsInside(mDeviceCoordinates, LocalCoordinates);
                if (IsInside)
                    break;
            }
            if (IsInside == false)
            {
                for (int k = 0; k < nelements; k++)
                {
                    ModelPart::ElementsContainerType::iterator it = el_begin + k;
                    pSelectedElement = (*(it.base()));
                    IsInside = pSelectedElement->GetGeometry().IsInside(mDeviceCoordinates, LocalCoordinates, 1.0e-5);
                    if (IsInside)
                        break;
                }
            }
            if (IsInside == false)
            {
                KRATOS_ERROR << "ERROR!!, PLEASE REPEAT THE SEARCH " << std::endl;
            }

            PointsNumber = pSelectedElement->GetGeometry().PointsNumber();

            for (int j = 0; j < PointsNumber; j++)
            {
                pSelectedElement->GetGeometry().GetPoint(j).FastGetSolutionStepValue(var) = mValue;
                pSelectedElement->GetGeometry().GetPoint(j).Fix(var);
            }
        }

        KRATOS_CATCH("");
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DamTemperaturebyDeviceProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "DamTemperaturebyDeviceProcess";
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
    bool mIsFixed;
    double mValue;
    array_1d<double, 3> mDeviceCoordinates;
    double mTimeUnitConverter;
    TableType::Pointer mpTable;
    int mTableId;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  private:
    /// Assignment operator.
    DamTemperaturebyDeviceProcess &operator=(DamTemperaturebyDeviceProcess const &rOther);

}; //Class

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DamTemperaturebyDeviceProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DamTemperaturebyDeviceProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_TEMPERATURE_BY_DEVICE_PROCESS defined */
