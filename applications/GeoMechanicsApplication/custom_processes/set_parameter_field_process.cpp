//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aron Noordam
//

// System includes

// External includes

// Project includes
#include "set_parameter_field_process.hpp"

#include <utilities/function_parser_utility.h>
#include <utilities/mortar_utilities.h>

#include "utilities/interval_utility.h"
#include "geo_mechanics_application_variables.h"


namespace Kratos
{

    SetParameterFieldProcess::SetParameterFieldProcess(ModelPart& rModelPart,
                                                       const Parameters& rSettings)
                                                            : Process(), mrModelPart(rModelPart),
                                                            mParameters(rSettings)
{
    // function type: python, cpp, input
    const Parameters default_parameters(R"(
        {
            "help"              : "This process sets a parameter field on a model part, where each element can have different material properties.",
            "model_part_name"   : "please_specify_model_part_name",
            "variable_name"     : "CUSTOM",
            "func_type"         : "input",               
            "function"          : "0",
            "dataset"           : "dummy",
            "dataset_file_name" : "dummy"
        }  )"
    );

    mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

}


void SetParameterFieldProcess::SetValueAtElement(Element& rElement, const Variable<double>& rVar, const double Value)
{

    Properties& r_prop = rElement.GetProperties();

    // Copies properties
    Properties::Pointer p_new_prop = Kratos::make_shared<Properties>(r_prop);

    // Adds new properties to the element
    p_new_prop->SetValue(rVar, Value);
    rElement.SetProperties(p_new_prop);

}


void SetParameterFieldProcess::SetParameterFieldUsingInputFunction(const Variable<double>& rVar)
{
    auto parameter_function = BasicGenericFunctionUtility(mParameters["function"].GetString());
    const double current_time = this->mrModelPart.GetProcessInfo().GetValue(TIME);

    for (Element& r_element : mrModelPart.Elements()) {

        const auto& r_geom = r_element.GetGeometry();

        // calculate parameter value at current element
        const double val = parameter_function.CallFunction(r_geom.Center().X(), r_geom.Center().Y(), r_geom.Center().Z(), current_time, 0, 0, 0);

        SetValueAtElement(r_element, rVar, val);

    }
}


void SetParameterFieldProcess::SetParameterFieldUsingParametersClass(const Variable<double>& rVar, Parameters& rParameters)
{

    const Vector& r_data_vector = rParameters["values"].GetVector();

    KRATOS_ERROR_IF_NOT(r_data_vector.size() == mrModelPart.Elements().size()) << "The parameter field "
    	"does not have the same size as the amount of elements within the model part!" << std::endl;

    // set new data on the elements
    IndexType i = 0;
    for (Element& r_element : mrModelPart.Elements())
    {
        SetValueAtElement(r_element, rVar, r_data_vector[i]);
        ++i;
    }
}

void SetParameterFieldProcess::SetParameterFieldUsingJsonString(const  Variable<double>& rVar)
{
    // get new data from the data set
    const std::string& r_dataset = mParameters["dataset"].GetString();


    Parameters new_data{ r_dataset };
    this->SetParameterFieldUsingParametersClass(rVar, new_data);

}
void SetParameterFieldProcess::SetParameterFieldUsingJsonFile(const Variable<double>& rVar)
{
    // Read json string in field parameters file, create Parameters
    const std::string& field_file_name = mParameters["dataset_file_name"].GetString();
    KRATOS_ERROR_IF_NOT(std::filesystem::exists(field_file_name)) << "The parameter field file specified with name \"" << field_file_name << "\" does not exist!" << std::endl;

    std::ifstream ifs(field_file_name);
    Parameters new_data{ ifs };

    this->SetParameterFieldUsingParametersClass(rVar, new_data);

}




void SetParameterFieldProcess::ExecuteInitialize()
{
    if (mrModelPart.GetProcessInfo()[IS_RESTARTED]) {
        return;
    }

    KRATOS_TRY

    const auto& r_var = KratosComponents< Variable<double> >::Get(mParameters["variable_name"].GetString());

    // set parameter field from input function
    if (mParameters["func_type"].GetString() == "input")
    {
        this->SetParameterFieldUsingInputFunction(r_var);
    }
    // set parameter field with a json string
    else if (mParameters["func_type"].GetString() == "json_string")
    {
        this->SetParameterFieldUsingJsonString(r_var);
    }
    // set parameter field from a json input file
    else if (mParameters["func_type"].GetString() == "json_file")
    {
        this->SetParameterFieldUsingJsonFile(r_var);
    }

    KRATOS_CATCH("")
}


}  // namespace Kratos.
