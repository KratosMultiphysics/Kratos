//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//                   Riccardo Tosi

#ifndef KRATOS_WEIGHTED_DIVERGENCE_CALCULATION_PROCESS_H
#define KRATOS_WEIGHTED_DIVERGENCE_CALCULATION_PROCESS_H


#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "processes/process.h"
#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "processes/compute_nodal_gradient_process.h"
#include "utilities/variable_utils.h"

#include <string>
#include <iostream>
#include <sstream>

#include <boost/functional/hash.hpp> //TODO: remove this dependence when Kratos has en internal one
#include <unordered_map> //TODO: remove this dependence when Kratos has en internal one
#include <utility>

namespace Kratos
{

class WeightedDivergenceCalculationProcess: public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(WeightedDivergenceCalculationProcess);

    typedef ModelPart::ElementType ElementType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for WeightedDivergenceCalculationProcess Process
//     WeightedDivergenceCalculationProcess(ModelPart& rModelPart,
//                      KratosParameters& parameters
//                     ):
//         Process(),
//         mrModelPart(rModelPart),
//         mrOptions(Flags()),
//         mrParameters(parameters)
//     {
//     }
    /// Constructor for WeightedDivergenceCalculationProcess Process
    WeightedDivergenceCalculationProcess(ModelPart& rModelPart
                    ):
        Process(),
        mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    ~WeightedDivergenceCalculationProcess() override {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    void Execute() override
    {
        KRATOS_TRY;

        // Set time coefficient: computations will be performed ONLY AFTER time_coefficient * END_TIME
        const double time_coefficient = 0.2;
        // Extract time information
        const ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
    	const double& time_step_current  = rCurrentProcessInfo[TIME];
	    const ProcessInfo& rPreviousProcessInfo = rCurrentProcessInfo.GetPreviousTimeStepInfo();
        const double& time_step_previous = rPreviousProcessInfo[TIME];
        const double& final_time = rCurrentProcessInfo[END_TIME];

        if (time_step_current >= time_coefficient * final_time) {

            // Check and set number of elements
            int number_elements;
            KRATOS_ERROR_IF(mrModelPart.Elements().size() == 0) << "the number of elements in the domain is zero. weighted divergence calculation cannot be applied"<< std::endl;
            number_elements = mrModelPart.Elements().size();

            // Check and set number of nodes
            int number_nodes;
            KRATOS_ERROR_IF(mrModelPart.Nodes().size() == 0) << "the number of nodes in the domain is zero. weighted divergence calculation cannot be applied" << std::endl;
            number_nodes = mrModelPart.Nodes().size();

            // Current domain size
            const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

            // Auxiliar containers
            Matrix DN_DX, J0;
            Vector N;

            // Iterate over the elements
            #pragma omp parallel for firstprivate(DN_DX,  N, J0)
            for(int i_elem=0; i_elem<number_elements; ++i_elem) {
                auto it_elem = mrModelPart.ElementsBegin() + i_elem;
                auto& r_geometry = it_elem->GetGeometry();

                // Current geometry information
                const std::size_t local_space_dimension = r_geometry.LocalSpaceDimension();
                // KRATOS_WATCH(local_space_dimension);
                const std::size_t number_of_nodes = r_geometry.PointsNumber();
                // KRATOS_WATCH(number_of_nodes);

                // Resize if needed
                if (DN_DX.size1() != number_of_nodes || DN_DX.size2() != dimension)
                    DN_DX.resize(number_of_nodes, dimension);
                if (N.size() != number_of_nodes)
                    N.resize(number_of_nodes);
                if (J0.size1() != dimension || J0.size2() != local_space_dimension)
                    J0.resize(dimension, local_space_dimension);

                // Build values vectors of the velocity
                Vector values_x(number_of_nodes);
                Vector values_y(number_of_nodes);
                Vector values_z(number_of_nodes);
                for(int i_node=0; i_node<number_of_nodes; ++i_node){
                    values_x[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_X);
                    values_y[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_Y);
                    values_z[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_Z);
                }
                // KRATOS_WATCH(values_x);
                // KRATOS_WATCH(values_y);
                // KRATOS_WATCH(values_z);

                // Set integration points
                const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod(); // Default is 0
                const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method); // Default is [3 dimensional integration point(0.333333 , 0.333333 , 0), weight = 0.5]
                const std::size_t number_of_integration_points = r_integration_points.size(); // Default is 1

                // Set containers of the shape functions and the local gradients
                const Matrix& rNcontainer = r_geometry.ShapeFunctionsValues(r_integration_method); // [1,3]((0.333333,0.333333,0.333333))
                const auto& rDN_DeContainer = r_geometry.ShapeFunctionsLocalGradients(r_integration_method); // [1]([3,2]((-1,-1),(1,0),(0,1)))
                // KRATOS_WATCH(rNcontainer);
                // KRATOS_WATCH(rDN_DeContainer);

                // Loop over integration points
                double divergence_current = 0;
                for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ){
                    // Getting the shape functions
                    noalias(N) = row(rNcontainer, point_number); // [3](0.333333,0.333333,0.333333)

                    // Get the jacobians
                    GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_integration_points[point_number], J0);
                    double detJ0;
                    Matrix InvJ0;
                    MathUtils<double>::GeneralizedInvertMatrix(J0, InvJ0, detJ0);
                    const Matrix& rDN_De = rDN_DeContainer[point_number];
                    GeometryUtils::ShapeFunctionsGradients(rDN_De, InvJ0, DN_DX);

                    // Compute local gradient
                    const Vector grad_x = prod(trans(DN_DX), values_x);
                    // KRATOS_WATCH(grad_x);
                    const Vector grad_y = prod(trans(DN_DX), values_y);
                    // KRATOS_WATCH(grad_y);
                    const Vector grad_z = prod(trans(DN_DX), values_z);
                    // KRATOS_WATCH(grad_z);
                    const double aux_current_divergence = grad_x[0] + grad_y[1] + grad_z[2];
                    const double gauss_point_volume = r_integration_points[point_number].Weight() * detJ0;
                    divergence_current += aux_current_divergence * gauss_point_volume;
                    // KRATOS_WATCH(divergence_current);
                }

                // Retrieve divergence from previous time step
                auto divergence_old_avg = it_elem->GetValue(DIVERGENCE);

                // Compute weighetd in time divergence average
                // TODO: coefficient indipendent!
                auto divergence_current_avg = std::sqrt(((time_step_previous-time_coefficient*final_time) * std::pow(divergence_old_avg,2) + (time_step_current - time_step_previous) * std::pow(divergence_current,2)) /  (time_step_current-time_coefficient*final_time));
                it_elem->SetValue(DIVERGENCE,divergence_current_avg);
                // if (i_elem == 3) {
                //     KRATOS_WATCH(it_elem->GetValue(DIVERGENCE));
                // }
            }
        }


        // // Initialize auxiliar variables
        // array_1d<double, 3> aux_zero_vector = ZeroVector(3);
        // for(int i_node=0; i_node < number_nodes; ++i_node) {
        //     auto it_node = mrModelPart.NodesBegin() + i_node;
        //     it_node->SetValue(NODAL_AREA, 0.0);
        //     it_node->SetValue(AUXILIAR_GRADIENT_X, aux_zero_vector);
        //     it_node->SetValue(AUXILIAR_GRADIENT_Y, aux_zero_vector);
        //     it_node->SetValue(AUXILIAR_GRADIENT_Z, aux_zero_vector);
        // }

        // // Compute gradient VELOCITY_X
        // auto gradient_process_x = ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>(mrModelPart, VELOCITY_X, AUXILIAR_GRADIENT_X, NODAL_AREA);
        // gradient_process_x.Execute();

        // // Compute gradient VELOCITY_Y
        // for(int i_node = 0; i_node < number_nodes; ++i_node) {
        //     auto it_node = mrModelPart.NodesBegin() + i_node;
        //     it_node->SetValue(NODAL_AREA, 0.0);
        // }
        // auto gradient_process_y = ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>(mrModelPart, VELOCITY_Y, AUXILIAR_GRADIENT_Y, NODAL_AREA);
        // gradient_process_y.Execute();

        // // Compute gradient VELOCITY_Z
        // for(int i_node = 0; i_node < number_nodes; ++i_node) {
        //     auto it_node = mrModelPart.NodesBegin() + i_node;
        //     it_node->SetValue(NODAL_AREA, 0.0);
        // }
        // auto gradient_process_z = ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>(mrModelPart, VELOCITY_Z, AUXILIAR_GRADIENT_Z, NODAL_AREA);
        // gradient_process_z.Execute();

        // for(int i_node = 0; i_node < number_nodes; ++i_node){
        //     auto it_node = mrModelPart.NodesBegin() + i_node;
        //     auto aux_gradient_x = it_node->GetValue(AUXILIAR_GRADIENT_X);
        //     auto aux_gradient_y = it_node->GetValue(AUXILIAR_GRADIENT_Y);
        //     auto aux_gradient_z = it_node->GetValue(AUXILIAR_GRADIENT_Z);
        //     // Retrieve divergence from previous time step
        //     auto divergence_old_avg = it_node->GetValue(DIVERGENCE);
        //     if (i_node == 1){
        //         // KRATOS_WATCH(time_step_current);
        //         // KRATOS_WATCH(time_step_previous);
        //         // KRATOS_WATCH(divergence_old_avg);
        //         // KRATOS_WATCH(it_node->GetValue(DIVERGENCE));
        //         // KRATOS_WATCH(aux_gradient_x);
        //         // KRATOS_WATCH(aux_gradient_y);
        //         // KRATOS_WATCH(aux_gradient_z);
        //     }
        //     // Compute divergence for current time step
        //     auto divergence_current = aux_gradient_x[0] + aux_gradient_y[1] + aux_gradient_z[2];
        //     // Compute weighetd in time divergence average
        //     auto divergence_current_avg = (time_step_previous * divergence_old_avg + (time_step_current - time_step_previous) * divergence_current) /  (time_step_current);
        //     it_node->SetValue(DIVERGENCE,divergence_current_avg);
        //     if (i_node == 1){
        //         // KRATOS_WATCH(divergence_current);
        //         // KRATOS_WATCH(it_node->GetValue(DIVERGENCE));
        //     }
        // }

        KRATOS_CATCH("");
    }



    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "WeightedDivergenceCalculationProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "WeightedDivergenceCalculationProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}


private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    Flags mrOptions;


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    WeightedDivergenceCalculationProcess& operator=(WeightedDivergenceCalculationProcess const& rOther);

    /// Copy constructor.
    WeightedDivergenceCalculationProcess(WeightedDivergenceCalculationProcess const& rOther);


    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  WeightedDivergenceCalculationProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const WeightedDivergenceCalculationProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



} // namespace Kratos


#endif // KRATOS_WEIGHTED_DIVERGENCE_CALCULATION_PROCESS_H
