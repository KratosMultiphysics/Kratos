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

    /**
     * This process computes the element average in time of the divergence and of the seminorm of the divergence of the velocity field.
     * We define the seminorm of the divergence as: \left \| \nabla \cdot u_{h} \right \|_{L^2(K)}^2 ,
     * where u is the velocity field and K an element of the domain \Omega.
     * The time average does not consider the transient first 20% part of the simulation.
     * The process requires a model part as input.
     */

    void Execute() override
    {
        KRATOS_TRY;

        // Set time coefficient: computations will be performed ONLY AFTER (time_coefficient * END_TIME)
        const double time_coefficient = 0.2;
        // Set copute maximums boolean
        const bool compute_maximums = false;
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
            for(int i_elem=0; i_elem < number_elements; ++i_elem) {
                auto it_elem = mrModelPart.ElementsBegin() + i_elem;
                auto& r_geometry = it_elem->GetGeometry();

                // Current geometry information
                const std::size_t local_space_dimension = r_geometry.LocalSpaceDimension();
                // KRATOS_WATCH(local_space_dimension);
                const std::size_t number_nodes_element = r_geometry.PointsNumber();
                // KRATOS_WATCH(number_nodes_element);

                // Resize if needed
                if (DN_DX.size1() != number_nodes_element || DN_DX.size2() != dimension)
                    DN_DX.resize(number_nodes_element, dimension);
                if (N.size() != number_nodes_element)
                    N.resize(number_nodes_element);
                if (J0.size1() != dimension || J0.size2() != local_space_dimension)
                    J0.resize(dimension, local_space_dimension);

                // Build values vectors of the velocity
                Vector values_x(number_nodes_element);
                Vector values_y(number_nodes_element);
                Vector values_z(number_nodes_element);
                for(int i_node=0; i_node < number_nodes_element; ++i_node){
                    values_x[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_X);
                    values_y[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_Y);
                    values_z[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VELOCITY_Z);
                }

                // Set integration points
                const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod(); // Default is 0
                const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method); // Default is [3 dimensional integration point(0.333333 , 0.333333 , 0), weight = 0.5]
                const std::size_t number_of_integration_points = r_integration_points.size(); // Default is 1

                // Set containers of the shape functions and the local gradients
                const Matrix& rNcontainer = r_geometry.ShapeFunctionsValues(r_integration_method); // [1,3]((0.333333,0.333333,0.333333))
                const auto& rDN_DeContainer = r_geometry.ShapeFunctionsLocalGradients(r_integration_method); // [1]([3,2]((-1,-1),(1,0),(0,1)))

                // Initialize auxiliary local variables
                double divergence_current = 0;
                double norm_divergence_current = 0;

                // Loop over integration points
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
                    const Vector grad_y = prod(trans(DN_DX), values_y);
                    const Vector grad_z = prod(trans(DN_DX), values_z);

                    // Compute divergence and
                    const double aux_current_divergence = grad_x[0] + grad_y[1] + grad_z[2];
                    const double aux_current_divergence_norm = grad_x[0]*grad_x[0] + grad_x[1]*grad_x[1] + grad_x[2]*grad_x[2] + grad_y[0]*grad_y[0] + grad_y[1]*grad_y[1] + grad_y[2]*grad_y[2] + grad_z[0]*grad_z[0] + grad_z[1]*grad_z[1] + grad_z[2]*grad_z[2];
                    const double gauss_point_volume = r_integration_points[point_number].Weight() * detJ0;
                    divergence_current += std::pow(aux_current_divergence,2) * gauss_point_volume;
                    norm_divergence_current += aux_current_divergence_norm * gauss_point_volume;
                }

                // Retrieve divergence from previous time step
                auto divergence_old = it_elem->GetValue(DIVERGENCE);
                auto norm_divergence_old = it_elem->GetValue(DIVERGENCE_H1SEMINORM);

                // Save element volume
                it_elem->SetValue(AUX_VOLUME,it_elem->GetGeometry().Area());

                // Compute divergence weighted time average
                auto divergence_current_avg = std::sqrt(((time_step_previous-time_coefficient*final_time) * std::pow(divergence_old,2) + (time_step_current - time_step_previous) * divergence_current) /  (time_step_current-time_coefficient*final_time));
                it_elem->SetValue(DIVERGENCE,divergence_current_avg);

                // Compute divergence_norm weighted time average
                auto norm_divergence_current_avg = std::sqrt(((time_step_previous-time_coefficient*final_time) * std::pow(norm_divergence_old,2) + (time_step_current - time_step_previous) * norm_divergence_current) /  (time_step_current-time_coefficient*final_time));
                it_elem->SetValue(DIVERGENCE_H1SEMINORM,norm_divergence_current_avg);

                // Compute maximums
                if (compute_maximums) {
                    auto divergence_current_max = std::sqrt(divergence_current);
                    auto norm_divergence_current_max = std::sqrt(norm_divergence_current);
                    if (divergence_old > divergence_current_max) divergence_current_max = divergence_old;
                    if (norm_divergence_old > norm_divergence_current_max) norm_divergence_current_max = norm_divergence_old;
                    it_elem->SetValue(DIVERGENCE,divergence_current_max);
                    it_elem->SetValue(DIVERGENCE_H1SEMINORM,norm_divergence_current_max);
                }

            }
        }

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
