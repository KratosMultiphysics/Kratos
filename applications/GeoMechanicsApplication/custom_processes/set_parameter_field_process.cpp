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
                                                            Parameters Settings)
                                                            : mrModelPart(rModelPart),
                                                            mParameters(Settings)
{
    // function type: python, cpp, input
    Parameters default_parameters(R"(
        {
            "help"            : "This process applies a moving load condition belonging to a modelpart. The load moves over line elements.",
            "model_part_name" : "please_specify_model_part_name",
            "variable_name"   : "CUSTOM",
            "func_type"       : "input",               
            "function"        : "0",
            "dataset"         : "dummy"
        }  )"
    );
    Parameters mParameters;

    mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

}


void SetParameterFieldProcess::SetValueAtElement(Element& rElement, const Variable<double>& rVar, double Value)
{


    Properties& r_prop = rElement.GetProperties();

    // Copies properties
    Properties::Pointer p_new_prop = Kratos::make_shared<Properties>(r_prop);

    // Adds new properties to the element
    p_new_prop->SetValue(rVar, Value);
    rElement.SetProperties(p_new_prop);

}


void SetParameterFieldProcess::ExecuteInitialize()
{
    KRATOS_TRY

    if (!this->mrModelPart.GetProcessInfo()[IS_RESTARTED]){

        const Variable<double>& r_var = KratosComponents< Variable<double> >::Get(mParameters["variable_name"].GetString());

        // set parameter field from input function
        if (mParameters["func_type"].GetString().compare("input") == 0)
        {

            BasicGenericFunctionUtility parameter_function = BasicGenericFunctionUtility(mParameters["function"].GetString());

            const double current_time = this->mrModelPart.GetProcessInfo().GetValue(TIME);

            for (Element& r_element : mrModelPart.Elements()) {

                auto& r_geom = r_element.GetGeometry();

                // calculate parameter value at current element
                const double val = parameter_function.CallFunction(r_geom.Center().X(), r_geom.Center().Y(), r_geom.Center().Z(), current_time, 0, 0, 0);
                

                this->SetValueAtElement(r_element, r_var, val);

            }
        }
        // set parameter field with a python function
        else if (mParameters["func_type"].GetString().compare("python") == 0)
        {

            // get new data from the data set
            const std::string dataset = mParameters["dataset"].GetString();
            Parameters new_data = Parameters(dataset);
            Vector data_vector = new_data["values"].GetVector();

            // set new data on the elements
            IndexType i = 0;
            for (Element& r_element : mrModelPart.Elements())
            {
                this->SetValueAtElement(r_element, r_var, data_vector[i]);
                i++;
            }
        }
        // set parameter field from an input json
        else if (mParameters["func_type"].GetString().compare("json") == 0)
        {


            // Read json string in field parameters file, create Parameters
            const std::string& field_file_name = mParameters["dataset"].GetString();

            KRATOS_ERROR_IF_NOT(std::filesystem::exists(field_file_name)) << "The parameter field file specified with name \"" << field_file_name << "\" does not exist!" << std::endl;

            std::ifstream ifs(field_file_name);
            Parameters new_data(ifs);
            Vector data_vector = new_data["values"].GetVector();

            KRATOS_ERROR_IF_NOT(data_vector.size() == mrModelPart.Elements().size()) << "The parameter field: \""
        	<< field_file_name << "\" does not have the same size as the amount of elements within the model part!" << std::endl;

            IndexType i = 0;
            for (Element& r_element : mrModelPart.Elements())
            {
                this->SetValueAtElement(r_element, r_var, data_vector[i]);
                i++;
            }
        }
    }
    KRATOS_CATCH("")
}


}  // namespace Kratos.
