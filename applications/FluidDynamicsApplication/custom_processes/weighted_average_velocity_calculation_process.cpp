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
#include "includes/checks.h"
#include "weighted_average_velocity_calculation_process.h"


namespace Kratos
{
    /* Public functions *******************************************************/

    // Constructor
    WeightedAverageVelocityCalculationProcess::WeightedAverageVelocityCalculationProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters):
        Process(),
        mrModelPart(rModelPart)
    {
        /**
         * We configure using the following parameters:
         * time_coefficient: Coefficient determining initial time for computing the average, i.e. TIME_START = time_coefficient * TIME_END
         */
        Parameters default_parameters = Parameters(R"(
        {
            "time_coefficient"  : 0.2
        })"
        );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // Set time coefficient: computations will be performed ONLY AFTER time_coefficient% of total time of the simulation
        mTimeCoefficient = ThisParameters["time_coefficient"].GetDouble();
    }

    std::string WeightedAverageVelocityCalculationProcess::Info() const
    {
        return "WeightedAverageVelocityCalculationProcess";
    }

    void WeightedAverageVelocityCalculationProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "WeightedAverageVelocityCalculationProcess";
    }

    void WeightedAverageVelocityCalculationProcess::PrintData(std::ostream& rOStream) const
    {
        this->PrintInfo(rOStream);
    }

    // Execution
    void WeightedAverageVelocityCalculationProcess::Execute()
    {
        KRATOS_TRY;

        // Extract time information
        const auto& r_current_process_info = mrModelPart.GetProcessInfo();
        const auto& r_previous_process_info = r_current_process_info.GetPreviousTimeStepInfo();
        const double& time_step_previous = r_previous_process_info[TIME];
        const double& final_time = r_current_process_info[END_TIME];
        const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];


        if (time_step_previous >= mTimeCoefficient * final_time) {
            // Check and set number of elements and check dimension
            KRATOS_ERROR_IF(dimension < 2 || dimension > 3) << "Inconsinstent dimension to execute the process. Dimension =" << dimension << " but needed dimension = 2 or dimension = 3." << std::endl;
            KRATOS_ERROR_IF(mrModelPart.NumberOfNodes() == 0) << "the number of nodes in the domain is zero. weighted divergence calculation cannot be applied"<< std::endl;
            const unsigned int number_nodes = mrModelPart.NumberOfNodes();

            const auto& r_node = *mrModelPart.NodesBegin();
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(AVERAGE_VELOCITY, r_node);

            // Iterate over the elements
            #pragma omp parallel for
            for(int i_node = 0; i_node < static_cast<int>(number_nodes); ++i_node) {
                auto it_node = mrModelPart.NodesBegin() + i_node;

                // Retrieve divergence from previous time step
                const auto variable_current = it_node->FastGetSolutionStepValue(VELOCITY);

                const auto average_old = it_node->FastGetSolutionStepValue(AVERAGE_VELOCITY);

                // Compute divergence weighted time average
                auto average_new = ComputeWeightedTimeAverage(average_old, variable_current);
                it_node->SetValue(AVERAGE_VELOCITY,average_new);

            }
        }

        KRATOS_CATCH("");
    }

    // Compute time average
    array_1d<double, 3> WeightedAverageVelocityCalculationProcess::ComputeWeightedTimeAverage(const array_1d<double, 3>& old_average, const array_1d<double, 3>& current_value)
    {
        // Extract time information
        const auto& r_current_process_info = mrModelPart.GetProcessInfo();
        const double& time_step_current  = r_current_process_info[TIME];
        const auto& r_previous_process_info = r_current_process_info.GetPreviousTimeStepInfo();
        const double& time_step_previous = r_previous_process_info[TIME];
        const double& final_time = r_current_process_info[END_TIME];

        array_1d<double, 3> new_average;
        for (std::size_t i_dim = 0; i_dim<3;i_dim++){
            new_average[i_dim] = std::sqrt(((time_step_previous-mTimeCoefficient*final_time) * std::pow(old_average[i_dim],2) + (time_step_current - time_step_previous) * current_value[i_dim]) / (time_step_current - mTimeCoefficient*final_time));
        }
        return new_average;
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const WeightedAverageVelocityCalculationProcess& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

}; // namespace Kratos
