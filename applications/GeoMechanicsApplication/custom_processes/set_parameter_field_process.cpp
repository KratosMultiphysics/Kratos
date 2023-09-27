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
        : Process(),
        mrModelPart(rModelPart),
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
            "dataset_file_name" : "dummy",
            "vector_variable_indices"      : []
        }  )"
        );

        mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
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

    void SetParameterFieldProcess::SetParameterFieldUsingInputFunction(const Variable<Vector>& rVar)
    {
        auto parameter_function = BasicGenericFunctionUtility(mParameters["function"].GetString());
        const double current_time = this->mrModelPart.GetProcessInfo().GetValue(TIME);
        const auto indices = GetVectorIndices();

        for (Element& r_element : mrModelPart.Elements()) {
            const auto& r_geom = r_element.GetGeometry();

            // calculate parameter value at current element
            const double val = parameter_function.CallFunction(r_geom.Center().X(), r_geom.Center().Y(), r_geom.Center().Z(), current_time, 0, 0, 0);
            // get properties per element
            Properties& r_prop = r_element.GetProperties();
            Vector& vector = r_prop.GetValue(rVar);

            // loop through the indexes that need to be set
            for (auto index : indices) {
                vector[index] = val;
            }

            SetValueAtElement(r_element, rVar, vector);
        }
    }

    void SetParameterFieldProcess::SetParameterFieldUsingParametersClass(const Variable<double>& rVar,
        const Parameters& rParameters)
    {
        const Vector& r_data_vector = rParameters["values"].GetVector();

        KRATOS_ERROR_IF_NOT(r_data_vector.size() == mrModelPart.Elements().size()) << "The parameter field "
            "does not have the same size as the amount of elements within the model part!" << std::endl;

        // set new data on the elements
        IndexType i = 0;
        for (Element& r_element : mrModelPart.Elements()) {
            SetValueAtElement(r_element, rVar, r_data_vector[i]);
            ++i;
        }
    }

    void SetParameterFieldProcess::SetParameterFieldUsingParametersClass(const Variable<Vector>& rVar,
        const Parameters& rParameters)
    {
        const Matrix& r_data_matrix = rParameters["values"].GetMatrix();

        KRATOS_ERROR_IF_NOT(r_data_matrix.size1() == mrModelPart.Elements().size()) << "The parameter field "
            "does not have the same size as the amount of elements within the model part!" << std::endl;

        // set new data on the elements
        const IndexType vector_size = r_data_matrix.size2();
        IndexType i = 0;
        for (Element& r_element : mrModelPart.Elements()) {
            Vector sub_vector;
            sub_vector.resize(vector_size);
            for (IndexType j = 0; j < vector_size; j++) {
                sub_vector[j] = r_data_matrix(i, j);
            }

            SetValueAtElement(r_element, rVar, sub_vector);
            ++i;
        }
    }

    void SetParameterFieldProcess::ExecuteInitialize()
    {
        if (mrModelPart.GetProcessInfo()[IS_RESTARTED]) {
            return;
        }

        KRATOS_TRY

            const auto variable_name = mParameters["variable_name"].GetString();
        if (KratosComponents<Variable<double>>::Has(variable_name)) {
            SetParameterFieldForVariableType(KratosComponents<Variable<double>>::Get(variable_name));
        }
        else if (KratosComponents<Variable<Vector>>::Has(variable_name)) {
            SetParameterFieldForVariableType(KratosComponents<Variable<Vector>>::Get(variable_name));
        }

        KRATOS_CATCH("")
    }

    std::vector<IndexType> SetParameterFieldProcess::GetVectorIndices() const
    {
        std::vector<IndexType> result;
        for (auto index : mParameters["vector_variable_indices"].GetVector()) {
            result.push_back(static_cast<IndexType>(index));
        }
        return result;
    }

} // namespace Kratos.