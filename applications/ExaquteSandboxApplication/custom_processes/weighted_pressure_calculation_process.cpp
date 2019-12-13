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
#include "custom_processes/weighted_pressure_calculation_process.h"


namespace Kratos
{
    /* Public functions *******************************************************/

    // Constructor
    WeightedPressureCalculationProcess::WeightedPressureCalculationProcess(
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

    std::string WeightedPressureCalculationProcess::Info() const
    {
        return "WeightedPressureCalculationProcess";
    }

    void WeightedPressureCalculationProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "WeightedPressureCalculationProcess";
    }

    void WeightedPressureCalculationProcess::PrintData(std::ostream& rOStream) const
    {
        this->PrintInfo(rOStream);
    }

    void WeightedPressureCalculationProcess::ExecuteInitialize()
    {
        KRATOS_TRY;

        auto& r_nodes_array = mrModelPart.Nodes();
        // Initialize variable
        VariableUtils().SetNonHistoricalVariableToZero(PRESSURE_WEIGHTED, r_nodes_array);

        mIsStartingAverageTimeSet = false;

        KRATOS_CATCH("");
    }

    void WeightedPressureCalculationProcess::ExecuteFinalizeSolutionStep()
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
            KRATOS_ERROR_IF(mrModelPart.NumberOfNodes() == 0) << "the number of nodes in the domain is zero. weighted pressure calculation cannot be applied"<< std::endl;
            const unsigned int number_nodes = mrModelPart.NumberOfNodes();

            // Set InitialTimeAverage if not done before
            if (!mIsStartingAverageTimeSet) {
                mStartingAverageTime = time_step_previous;
                mIsStartingAverageTimeSet = true;
            }

            // Iterate over the nodes
            #pragma omp parallel for
            for(int i_node = 0; i_node < static_cast<int>(number_nodes); ++i_node) {
                auto it_node = mrModelPart.NodesBegin() + i_node;

                // Retrieve current time step local variable
                double pressure_current = it_node->GetSolutionStepValue(PRESSURE);

                // Retrieve variable from previous time step
                const double pressure_old = it_node->GetValue(PRESSURE_WEIGHTED);

                // Compute weighted time average
                auto pressure_current_avg = ComputeWeightedTimeAverage(pressure_old, pressure_current);
                it_node->SetValue(PRESSURE_WEIGHTED,pressure_current_avg);
            }
        }

        KRATOS_CATCH("");
    }

    // Compute time average
    double WeightedPressureCalculationProcess::ComputeWeightedTimeAverage(const double& old_average, const double& current_value)
    {
        // Extract time information
        const auto& r_current_process_info = mrModelPart.GetProcessInfo();
        const double& time_step_current  = r_current_process_info[TIME];
        const auto& r_previous_process_info = r_current_process_info.GetPreviousTimeStepInfo();
        const double& time_step_previous = r_previous_process_info[TIME];
        const double& final_time = r_current_process_info[END_TIME];

        // const double new_average = std::sqrt(((time_step_previous-mTimeCoefficient*final_time) * std::pow(old_average,2) + (time_step_current - time_step_previous) * current_value) / (time_step_current - mTimeCoefficient*final_time));
        const double new_average = ((time_step_previous-mStartingAverageTime) * old_average + (time_step_current - time_step_previous) * current_value) / (time_step_current - mStartingAverageTime);
        return new_average;
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const WeightedPressureCalculationProcess& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

}; // namespace Kratos
