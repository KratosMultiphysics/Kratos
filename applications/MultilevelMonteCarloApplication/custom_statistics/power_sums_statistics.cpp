//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//
//

// System includes

// External includes

// Project includes

// Application includes
#include "custom_statistics/power_sums_statistics.h"


namespace Kratos
{
    /* Public functions *******************************************************/

    // Constructor
    PowerSumsStatistics::PowerSumsStatistics(
        ModelPart& rModelPart,
        Parameters ThisParameters):
        Process(),
        mrModelPart(rModelPart)
    {
        /**
         * Configure using the following parameters:
         * reference_variable_name: variable for which computing the power sums
         */
        Parameters default_parameters = Parameters(R"(
        {
            "reference_variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME"
        })"
        );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // Set parameters
        mReferenceVariable = ThisParameters["reference_variable_name"].GetString();
    }

    std::string PowerSumsStatistics::Info() const
    {
        return "PowerSumsStatistics";
    }

    void PowerSumsStatistics::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "PowerSumsStatistics";
    }

    void PowerSumsStatistics::PrintData(std::ostream& rOStream) const
    {
        this->PrintInfo(rOStream);
    }

    void PowerSumsStatistics::ExecuteInitialize()
    {
        KRATOS_TRY;

        auto& r_nodes_array = mrModelPart.Nodes();
        // Initialize power sums variables
        VariableUtils().SetNonHistoricalVariableToZero(POWER_SUM_1, r_nodes_array);
        VariableUtils().SetNonHistoricalVariableToZero(POWER_SUM_2, r_nodes_array);
        VariableUtils().SetNonHistoricalVariableToZero(POWER_SUM_3, r_nodes_array);
        VariableUtils().SetNonHistoricalVariableToZero(POWER_SUM_4, r_nodes_array);

        KRATOS_CATCH("");
    }

    void PowerSumsStatistics::ExecuteFinalizeSolutionStep()
    {
        KRATOS_TRY;

        // Check and set number of elements and check dimension
        KRATOS_ERROR_IF(mrModelPart.NumberOfNodes() == 0) << "The number of nodes in the domain is zero. The power sums statistic cannot be applied."<< std::endl;
        const unsigned int number_nodes = mrModelPart.NumberOfNodes();

        // Extract informations
        if (KratosComponents<Variable<double>>::Has(mReferenceVariable)) {
            const auto& r_reference_var = KratosComponents<Variable<double>>::Get(mReferenceVariable);
            // Iterate over the nodes
            #pragma omp parallel for
            for(int i_node = 0; i_node < static_cast<int>(number_nodes); ++i_node) {
                auto it_node = mrModelPart.NodesBegin() + i_node;
                // Retrieve current time step local variable
                double variable_current = it_node->GetSolutionStepValue(r_reference_var);
                // Retrieve power sums
                const double S1_old = it_node->GetValue(POWER_SUM_1);
                const double S2_old = it_node->GetValue(POWER_SUM_2);
                const double S3_old = it_node->GetValue(POWER_SUM_3);
                const double S4_old = it_node->GetValue(POWER_SUM_4);
                // Update power sums
                const double S1 = S1_old + std::pow(variable_current,1);
                const double S2 = S2_old + std::pow(variable_current,2);
                const double S3 = S3_old + std::pow(variable_current,3);
                const double S4 = S4_old + std::pow(variable_current,4);
                // Set power sums
                it_node->SetValue(POWER_SUM_1,S1);
                it_node->SetValue(POWER_SUM_2,S2);
                it_node->SetValue(POWER_SUM_3,S3);
                it_node->SetValue(POWER_SUM_4,S4);
            }
        }
        else {
            KRATOS_ERROR << "Variable " << mReferenceVariable << " not found in the scalar variables list of Kratos.\n";
        }

        KRATOS_CATCH("");
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const PowerSumsStatistics& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

}; // namespace Kratos
