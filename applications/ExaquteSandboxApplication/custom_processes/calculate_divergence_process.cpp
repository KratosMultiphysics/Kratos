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
#include "utilities/geometry_utilities.h"
#include "utilities/variable_utils.h"
#include "exaqute_sandbox_application_variables.h"

// Application includes
#include "custom_processes/calculate_divergence_process.h"


namespace Kratos
{
    /* Public functions *******************************************************/

    // Constructor
    CalculateDivergenceProcess::CalculateDivergenceProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters):
        Process(),
        mrModelPart(rModelPart)
    {}

    std::string CalculateDivergenceProcess::Info() const
    {
        return "CalculateDivergenceProcess";
    }

    void CalculateDivergenceProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "CalculateDivergenceProcess";
    }

    void CalculateDivergenceProcess::PrintData(std::ostream& rOStream) const
    {
        this->PrintInfo(rOStream);
    }

    void CalculateDivergenceProcess::ExecuteInitialize()
    {
        KRATOS_TRY;

        auto& r_elements_array = mrModelPart.Elements();
        // Initialize variables
        VariableUtils().SetNonHistoricalVariableToZero(DIVERGENCE, r_elements_array);
        VariableUtils().SetNonHistoricalVariableToZero(VELOCITY_H1_SEMINORM, r_elements_array);

        KRATOS_CATCH("");
    }

    void CalculateDivergenceProcess::ExecuteBeforeOutputStep()
    {
        KRATOS_TRY;

        // Extract information
        const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

        // Check and set number of elements and check dimension
        KRATOS_ERROR_IF(dimension < 2 || dimension > 3) << "Inconsinstent dimension to execute the process. Dimension =" << dimension << " but needed dimension = 2 or dimension = 3." << std::endl;
        KRATOS_ERROR_IF(mrModelPart.NumberOfElements() == 0) << "The number of elements in the domain is zero. The process can not be applied."<< std::endl;
        const unsigned int number_elements = mrModelPart.NumberOfElements();

        // Auxiliar containers
        GeometryData::ShapeFunctionsGradientsType DN_DX;
        Vector grad_x;
        Vector grad_y;
        Vector grad_z;

        // Iterate over the elements
        #pragma omp parallel for firstprivate(DN_DX,grad_x,grad_y,grad_z)
        for(int i_elem = 0; i_elem < static_cast<int>(number_elements); ++i_elem) {
            auto it_elem = mrModelPart.ElementsBegin() + i_elem;
            const auto& r_geometry = it_elem->GetGeometry();

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
                values_z[i_node] = r_velocity[2];
            }

            // Set integration points
            const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod(); // Default is 0
            const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
            const std::size_t number_of_integration_points = r_integration_points.size(); // Default is 1

            // Initialize auxiliary local variables
            double aux_divergence = 0;
            double aux_velocity_seminorm = 0;

            // Get gradient shape functions and detJ0 on integration points
            Vector detJ0;
            r_geometry.ShapeFunctionsIntegrationPointsGradients(DN_DX, detJ0, r_integration_method);

            // Loop over integration points
            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ){

                // Compute local gradient
                grad_x = prod(trans(DN_DX[point_number]), values_x);
                grad_y = prod(trans(DN_DX[point_number]), values_y);
                grad_z = prod(trans(DN_DX[point_number]), values_z);

                // Compute divergence and velocity seminorm
                const double aux_current_divergence = ComputeAuxiliaryElementDivergence(grad_x, grad_y, grad_z);
                const double aux_current_velocity_seminorm = ComputeAuxiliaryElementVelocitySeminorm(grad_x, grad_y, grad_z);
                const double gauss_point_volume = r_integration_points[point_number].Weight() * detJ0[point_number];
                aux_divergence += std::pow(aux_current_divergence,2) * gauss_point_volume;
                aux_velocity_seminorm += aux_current_velocity_seminorm * gauss_point_volume;
            }

            // Set values to corresponding variables
            it_elem->SetValue(DIVERGENCE,aux_divergence);
            it_elem->SetValue(VELOCITY_H1_SEMINORM,aux_velocity_seminorm);
        }

        KRATOS_CATCH("");
    }

    // Compute local auxiliar divergence
    double CalculateDivergenceProcess::ComputeAuxiliaryElementDivergence(Vector& grad_x, Vector& grad_y, Vector& grad_z)
    {
        double aux_current_divergence;
        aux_current_divergence = grad_x[0] + grad_y[1] + grad_z[2];
        return aux_current_divergence;
    }

    // Compute local auxiliar velocity seminorm
    double CalculateDivergenceProcess::ComputeAuxiliaryElementVelocitySeminorm(Vector& grad_x, Vector& grad_y, Vector& grad_z)
    {
        double aux_current_velocity_seminorm;
        aux_current_velocity_seminorm = inner_prod(grad_x, grad_x) + inner_prod(grad_y, grad_y) + inner_prod(grad_z,grad_z);
        return aux_current_velocity_seminorm;
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const CalculateDivergenceProcess& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

}; // namespace Kratos
