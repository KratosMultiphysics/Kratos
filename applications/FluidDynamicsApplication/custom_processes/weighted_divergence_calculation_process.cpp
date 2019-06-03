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
    WeightedDivergenceCalculationProcess::WeightedDivergenceCalculationProcess(ModelPart& rModelPart):
        Process(),
        mrModelPart(rModelPart)
    {

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

    // Execution
    void WeightedDivergenceCalculationProcess::Execute()
    {
        KRATOS_TRY;

        // Set time coefficient: computations will be performed ONLY AFTER (time_coefficient * END_TIME)
        const double time_coefficient = 0.2;

        // Extract time information
        const auto& r_current_process_info = mrModelPart.GetProcessInfo();
        const double& time_step_current  = r_current_process_info[TIME];
        const auto& r_previous_process_info = r_current_process_info.GetPreviousTimeStepInfo();
        const double& time_step_previous = r_previous_process_info[TIME];
        const double& final_time = r_current_process_info[END_TIME];

        if (time_step_current >= time_coefficient * final_time) {

            // Check and set number of elements
            KRATOS_ERROR_IF(mrModelPart.NumberOfElements() == 0) << "the number of elements in the domain is zero. weighted divergence calculation cannot be applied"<< std::endl;
            const unsigned int number_elements = mrModelPart.NumberOfElements();

            // Current domain size
            const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

            // Auxiliar containers
            Matrix DN_DX, J0;
            Vector N;

            // Iterate over the elements
            #pragma omp parallel for firstprivate(DN_DX,  N, J0)
            for(int i_elem = 0; i_elem < static_cast<int>(number_elements); ++i_elem) {
                auto it_elem = mrModelPart.ElementsBegin() + i_elem;
                auto& r_geometry = it_elem->GetGeometry();

                // Current geometry information
                const std::size_t local_space_dimension = r_geometry.LocalSpaceDimension();
                const std::size_t number_nodes_element = r_geometry.PointsNumber();

                // Resize if needed
                if (DN_DX.size1() != number_nodes_element || DN_DX.size2() != dimension) {
                    DN_DX.resize(number_nodes_element, dimension);
                }
                if (N.size() != number_nodes_element) {
                    N.resize(number_nodes_element);
                }
                if (J0.size1() != dimension || J0.size2() != local_space_dimension) {
                    J0.resize(dimension, local_space_dimension);
                }

                // Build values vectors of the velocity
                Vector values_x(number_nodes_element);
                Vector values_y(number_nodes_element);
                Vector values_z(number_nodes_element);
                for(int i_node=0; i_node < static_cast<int>(number_nodes_element); ++i_node){
                    const auto &r_velocity = r_geometry[i_node].FastGetSolutionStepValue(VELOCITY);
                    values_x[i_node] = r_velocity[0];
                    values_y[i_node] = r_velocity[1];
                    values_z[i_node] = r_velocity[2];
                }

                // Set integration points
                const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod(); // Default is 0
                const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
                const std::size_t number_of_integration_points = r_integration_points.size(); // Default is 1

                // Set containers of the shape functions and the local gradients
                const Matrix& rNcontainer = r_geometry.ShapeFunctionsValues(r_integration_method);
                const auto& rDN_DeContainer = r_geometry.ShapeFunctionsLocalGradients(r_integration_method);

                // Initialize auxiliary local variables
                double divergence_current = 0;
                double velocity_seminorm_current = 0;

                // Loop over integration points
                for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ){
                    // Getting the shape functions
                    noalias(N) = row(rNcontainer, point_number);

                    // Get the jacobians
                    GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_integration_points[point_number], J0);
                    double detJ0;
                    Matrix InvJ0;
                    MathUtils<double>::GeneralizedInvertMatrix(J0, InvJ0, detJ0);
                    const Matrix& rDN_De = rDN_DeContainer[point_number];
                    GeometryUtils::ShapeFunctionsGradients(rDN_De, InvJ0, DN_DX);

                    // Compute local gradient
                    const Vector grad_x = prod(trans(DN_DX), values_x);
                    const Vector grad_y = prod(trans(DN_DX), values_y);
                    const Vector grad_z = prod(trans(DN_DX), values_z);

                    // Compute divergence and velocity seminorm
                    const double aux_current_divergence = grad_x[0] + grad_y[1] + grad_z[2];
                    const double aux_current_velocity_seminorm = grad_x[0]*grad_x[0] + grad_x[1]*grad_x[1] + grad_x[2]*grad_x[2] + grad_y[0]*grad_y[0] + grad_y[1]*grad_y[1] + grad_y[2]*grad_y[2] + grad_z[0]*grad_z[0] + grad_z[1]*grad_z[1] + grad_z[2]*grad_z[2];
                    const double gauss_point_volume = r_integration_points[point_number].Weight() * detJ0;
                    divergence_current += std::pow(aux_current_divergence,2) * gauss_point_volume;
                    velocity_seminorm_current += aux_current_velocity_seminorm * gauss_point_volume;
                }
                // Retrieve divergence from previous time step
                auto divergence_old = it_elem->GetValue(DIVERGENCE);
                const double velocity_seminorm_old = it_elem->GetValue(VELOCITY_H1_SEMINORM);

                // Compute divergence weighted time average
                auto divergence_current_avg = std::sqrt(((time_step_previous-time_coefficient*final_time) * std::pow(divergence_old,2) + (time_step_current - time_step_previous) * divergence_current) /  (time_step_current-time_coefficient*final_time));
                it_elem->SetValue(DIVERGENCE,divergence_current_avg);

                // Compute divergence_norm weighted time average
                auto velocity_seminorm_current_avg = std::sqrt(((time_step_previous-time_coefficient*final_time) * std::pow(velocity_seminorm_old,2) + (time_step_current - time_step_previous) * velocity_seminorm_current) /  (time_step_current-time_coefficient*final_time));
                it_elem->SetValue(VELOCITY_H1_SEMINORM,velocity_seminorm_current_avg);

            }
        }

        KRATOS_CATCH("");
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
