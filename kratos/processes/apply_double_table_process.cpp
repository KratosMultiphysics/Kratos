//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "processes/apply_double_table_process.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class ApplyComponentTableProcess
 * @ingroup KratosCore
 * @brief This class assings a value to the BC or loads according to a table defined in the .mdpa
 * @author Ignasi de Pouplana
 * @author Alejandro Cornejo
*/
/// Constructor
ApplyDoubleTableProcess::ApplyDoubleTableProcess(
    ModelPart& model_part, 
    Parameters rParameters) 
    : ApplyComponentTableProcess(model_part, rParameters) 
{

}

/***********************************************************************************/
/***********************************************************************************/

/// this function is designed for being called at the beginning of the computations
/// right after reading the model and the groups
void ApplyDoubleTableProcess::ExecuteInitialize()
{
    KRATOS_TRY;
    
    Variable<double> var = KratosComponents<Variable<double>>::Get(mVariableName);
    const int nnodes = static_cast<int>(mrModelPart.Nodes().size());

    if (nnodes != 0) {
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        #pragma omp parallel for
        for (int i = 0; i < nnodes; i++) {
            ModelPart::NodesContainerType::iterator it = it_begin + i;
            if (mIsFixed) {
                it->Fix(var);
            }
            it->FastGetSolutionStepValue(var) = mInitialValue;
        }
    }
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

/// this function will be executed at every time step BEFORE performing the solve phase
void ApplyDoubleTableProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY;

    Variable<double> var = KratosComponents<Variable<double>>::Get(mVariableName);
	typedef Variable<double> double_var;
	double_var time_var_component = KratosComponents<double_var>::Get(mTimeVariableName);
    
    const double time = mrModelPart.GetProcessInfo()[time_var_component];
    const double value = mpTable->GetValue(time);
    
    const int nnodes = static_cast<int>(mrModelPart.Nodes().size());

    if (nnodes != 0) {
        auto& r_nodes_array = mrModelPart.Nodes();
        VariableUtils().SetScalarVar<double_var>(var, value, r_nodes_array);
    }
    KRATOS_CATCH("");
}

} // namespace Kratos.

