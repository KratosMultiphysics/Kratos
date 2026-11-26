//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    D.J. Vicente
//
//

#if !defined(KRATOS_DAM_NODAL_REFERENCE_TEMPERATURE_PROCESS )
#define  KRATOS_DAM_NODAL_REFERENCE_TEMPERATURE_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"
#include "custom_utilities/dam_nodal_field_utilities.h"

namespace Kratos
{

class DamNodalReferenceTemperatureProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(DamNodalReferenceTemperatureProcess);

    typedef Table<double,double> TableType;
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamNodalReferenceTemperatureProcess(ModelPart& rModelPart, TableType& Table,
                                Parameters& rParameters
                                ) : Process(Flags()) , mrModelPart(rModelPart) , mrTable(Table)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name"      : "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "initial_value"      : 0.0,
                "input_file_name"    : "",
                "interval":[
                0.0,
                0.0
                ]
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate against defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mInitialValue = rParameters["initial_value"].GetDouble();
        mInputFile = rParameters["input_file_name"].GetString();

        // The GiD interface writes "- No file" / "- Add new file" when no
        // per-node input file is selected. The substring test is the historical
        // behavior of this process and is kept as is.
        mbUseTable = !(mInputFile.empty() ||
                       mInputFile.find("- No file") != std::string::npos ||
                       mInputFile.find("- Add new file") != std::string::npos);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~DamNodalReferenceTemperatureProcess() {}


    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute() override
    {
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteBeforeSolutionLoop() override
    {

        KRATOS_TRY;

        AssignReferenceTemperature();

        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep() override
    {

        KRATOS_TRY;

        AssignReferenceTemperature();

        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DamNodalReferenceTemperatureProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DamNodalReferenceTemperatureProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mrModelPart;
    TableType& mrTable;
    std::string mVariableName;
    double mInitialValue;
    std::string mInputFile;
    bool mbUseTable = false;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    void AssignReferenceTemperature()
    {
        const Variable<double>& var = KratosComponents< Variable<double> >::Get(mVariableName);
        DamNodalFieldUtilities::SetHistoricalValuesFromTableOrConstant(
            mrModelPart, var, mInitialValue, mrTable, mbUseTable);
    }

    /// Assignment operator.
    DamNodalReferenceTemperatureProcess& operator=(DamNodalReferenceTemperatureProcess const& rOther);

};//Class


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                    DamNodalReferenceTemperatureProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DamNodalReferenceTemperatureProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_NODAL_REFERENCE_TEMPERATURE_PROCESS defined */

