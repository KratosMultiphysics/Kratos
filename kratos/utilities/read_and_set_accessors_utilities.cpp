//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes
#include "utilities/read_and_set_accessors_utilities.h"

namespace Kratos {

/***********************************************************************************/
/***********************************************************************************/

void ReadAndSetAccessorsUtilities::ReadAndSetAccessors(
    const Parameters MaterialData,
    Properties& rProperty
    )
{
    if (MaterialData.Has("accessors")) {
        Parameters accessors = MaterialData["accessors"];

        // Loop over the accessors list
        for (auto iter = accessors.begin(); iter != accessors.end(); ++iter) {
            auto accessor_param = accessors.GetValue(iter.name());
            
            // Table Accessor
            if (accessor_param["accessor_type"].GetString() == "table_accessor") {
                // Independent Variable
                std::string input_var_name = accessor_param["properties"]["table_input_variable"].GetString();
                Variable<double> *p_input_var = static_cast<Variable<double> *>(KratosComponents<VariableData>::pGet(input_var_name));

                // Dependent Variable
                std::string output_var_name = accessor_param["properties"]["table_output_variable"].GetString();
                const auto& r_output_var  = KratosComponents<Variable<double>>().Get(output_var_name);

                // We set the variable type of the input variable (node_historical, node_non_historical and element)
                std::string input_var_type = accessor_param["properties"].Has("table_input_variable_type") ? accessor_param["properties"]["table_input_variable_type"].GetString() : "node_historical";

                KRATOS_ERROR_IF(rProperty.HasAccessor(r_output_var)) << "You are trying to add an TableAccessor between " << input_var_name << " and " << output_var_name << " which already exists..." << std::endl;

                rProperty.SetAccessor(r_output_var, (TableAccessor(*p_input_var, input_var_type)).Clone());
            } else {
                // No more accessors implemented currently
                KRATOS_ERROR << "This Accessor type is not available, only TableAccessor is ready for now" << std::endl;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

}  // namespace Kratos.
