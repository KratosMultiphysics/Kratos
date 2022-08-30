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
#include "custom_processes/weighted_divergence_calculation_process.h"


namespace Kratos
{
    /* Public functions *******************************************************/

    // Constructor
    WeightedDivergenceCalculationProcess::WeightedDivergenceCalculationProcess(
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

    std::string WeightedDivergenceCalculationProcess::Info() const
    {
        return "WeightedDivergenceCalculationProcess";
    }

    void WeightedDivergenceCalculationProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "WeightedDivergenceCalculationProcess";
    }

    void WeightedDivergenceCalculationProcess::PrintData(std::ostream& rOStream) const
    {
        this->PrintInfo(rOStream);
    }

    void WeightedDivergenceCalculationProcess::ExecuteInitialize()
    {
        KRATOS_TRY;

        auto& r_nodes_array = mrModelPart.Nodes();
        // Initialize variable
        VariableUtils().SetNonHistoricalVariableToZero(AVERAGED_DIVERGENCE, r_nodes_array);

        KRATOS_CATCH("");
    }

    void WeightedDivergenceCalculationProcess::ExecuteFinalizeSolutionStep()
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
            KRATOS_ERROR_IF(mrModelPart.NumberOfElements() == 0) << "the number of elements in the domain is zero. weighted divergence calculation cannot be applied"<< std::endl;
            const unsigned int number_elements = mrModelPart.NumberOfElements();

            // Auxiliar containers
            GeometryData::ShapeFunctionsGradientsType DN_DX;

            // Iterate over the elements
            #pragma omp parallel for firstprivate(DN_DX)
            for(int i_elem = 0; i_elem < static_cast<int>(number_elements); ++i_elem) {
                auto it_elem = mrModelPart.ElementsBegin() + i_elem;
                auto& r_geometry = it_elem->GetGeometry();

                // Current geometry information
                const std::size_t number_nodes_element = r_geometry.PointsNumber();

                // Build values vectors of the velocity
                Vector values_x(number_nodes_element);
                Vector values_y(number_nodes_element);
                Vector values_z(number_nodes_element);
                for(int i_node=0; i_node < static_cast<int>(number_nodes_element); ++i_node){
                    const auto &r_velocity = r_geometry[i_node].FastGetSolutionStepValue(VELOCITY);
                    values_x[i_node] = r_velocity[0];
                    values_y[i_node] = r_velocity[1];
                    if (dimension == 3) {
                        values_z[i_node] = r_velocity[2];
                    }
                }

                // Set integration points
                const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod(); // Default is 0
                const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
                const std::size_t number_of_integration_points = r_integration_points.size(); // Default is 1

                // Initialize auxiliary local variables
                double divergence_current = 0;
                double velocity_seminorm_current = 0;

                // Get gradient shape functions and detJ0 on integration points
                Vector detJ0;
                r_geometry.ShapeFunctionsIntegrationPointsGradients(DN_DX, detJ0, r_integration_method);

                // Loop over integration points
                for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ){

                    // Compute local gradient
                    Vector grad_x;
                    Vector grad_y;
                    Vector grad_z;
                    grad_x = prod(trans(DN_DX[point_number]), values_x);
                    grad_y = prod(trans(DN_DX[point_number]), values_y);
                    if (dimension == 3) {
                        grad_z = prod(trans(DN_DX[point_number]), values_z);
                    }

                    // Compute divergence and velocity seminorm
                    const double aux_current_divergence = ComputeAuxiliaryElementDivergence(grad_x, grad_y, grad_z);
                    const double aux_current_velocity_seminorm = ComputeAuxiliaryElementVelocitySeminorm(grad_x, grad_y, grad_z);
                    const double gauss_point_volume = r_integration_points[point_number].Weight() * detJ0[point_number];
                    divergence_current += std::pow(aux_current_divergence,2) * gauss_point_volume;
                    velocity_seminorm_current += aux_current_velocity_seminorm * gauss_point_volume;
                }

                // Retrieve divergence from previous time step
                const double divergence_old = it_elem->GetValue(AVERAGED_DIVERGENCE);
                const double velocity_seminorm_old = it_elem->GetValue(VELOCITY_H1_SEMINORM);

                // Compute divergence weighted time average
                auto divergence_current_avg = ComputeWeightedTimeAverage(divergence_old, divergence_current);
                it_elem->SetValue(AVERAGED_DIVERGENCE,divergence_current_avg);

                // Compute divergence_norm weighted time average
                auto velocity_seminorm_current_avg = ComputeWeightedTimeAverage(velocity_seminorm_old, velocity_seminorm_current);
                it_elem->SetValue(VELOCITY_H1_SEMINORM,velocity_seminorm_current_avg);
            }
        }

        KRATOS_CATCH("");
    }

    // Compute time average
    double WeightedDivergenceCalculationProcess::ComputeWeightedTimeAverage(const double& old_average, const double& current_value)
    {
        // Extract time information
        const auto& r_current_process_info = mrModelPart.GetProcessInfo();
        const double& time_step_current  = r_current_process_info[TIME];
        const auto& r_previous_process_info = r_current_process_info.GetPreviousTimeStepInfo();
        const double& time_step_previous = r_previous_process_info[TIME];
        const double& final_time = r_current_process_info[END_TIME];

        const double new_average = std::sqrt(((time_step_previous-mTimeCoefficient*final_time) * std::pow(old_average,2) + (time_step_current - time_step_previous) * current_value) / (time_step_current - mTimeCoefficient*final_time));
        return new_average;
    }

    // Compute local auxiliar divergence
    double WeightedDivergenceCalculationProcess::ComputeAuxiliaryElementDivergence(Vector& grad_x, Vector& grad_y, Vector& grad_z)
    {
        const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
        double aux_current_divergence;
        if (dimension == 2) {
            aux_current_divergence = grad_x[0] + grad_y[1];
        }
        else if (dimension == 3) {
            aux_current_divergence = grad_x[0] + grad_y[1] + grad_z[2];
        }
        return aux_current_divergence;
    }

    // Compute local auxiliar velocity seminorm
    double WeightedDivergenceCalculationProcess::ComputeAuxiliaryElementVelocitySeminorm(Vector& grad_x, Vector& grad_y, Vector& grad_z)
    {
        const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
        double aux_current_velocity_seminorm;
        if (dimension == 2) {
            aux_current_velocity_seminorm = inner_prod(grad_x, grad_x) + inner_prod(grad_y, grad_y);
        }
        else if (dimension == 3) {
            aux_current_velocity_seminorm = inner_prod(grad_x, grad_x) + inner_prod(grad_y, grad_y) + inner_prod(grad_z,grad_z);
        }
        return aux_current_velocity_seminorm;
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const WeightedDivergenceCalculationProcess& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

}; // namespace Kratos
