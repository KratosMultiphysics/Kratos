//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//

#ifndef KRATOS_DIVERGENCE_CALCULATION_PROCESS_H
#define KRATOS_DIVERGENCE_CALCULATION_PROCESS_H


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

class DivergenceCalculationProcess: public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(DivergenceCalculationProcess);

    typedef ModelPart::ElementType ElementType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for DivergenceCalculationProcess Process
//     DivergenceCalculationProcess(ModelPart& rModelPart,
//                      KratosParameters& parameters
//                     ):
//         Process(),
//         mrModelPart(rModelPart),
//         mrOptions(Flags()),
//         mrParameters(parameters)
//     {
//     }
    /// Constructor for DivergenceCalculationProcess Process
    DivergenceCalculationProcess(ModelPart& rModelPart
                    ):
        Process(),
        mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    ~DivergenceCalculationProcess() override {}


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

        // Extract time information
        const ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
    	const double& time_step_current  = rCurrentProcessInfo[TIME];
	    const ProcessInfo& rPreviousProcessInfo = rCurrentProcessInfo.GetPreviousTimeStepInfo();
        const double& time_step_previous = rPreviousProcessInfo[TIME];

        // Check and set number of elements
        int number_elements;
        KRATOS_ERROR_IF(mrModelPart.Elements().size() == 0) << "the number of elements in the domain is zero. divergence calculation cannot be applied"<< std::endl;
        number_elements = mrModelPart.Elements().size();

        // Check and set number of nodes
        int number_nodes;
        KRATOS_ERROR_IF(mrModelPart.Nodes().size() == 0) << "the number of nodes in the domain is zero. divergence calculation cannot be applied" << std::endl;
        number_nodes = mrModelPart.Nodes().size();

        // Geometry information
        const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

        // Initialize auxiliar variables
        array_1d<double, 3> aux_zero_vector = ZeroVector(3);
        for(int i_node=0; i_node < number_nodes; ++i_node) {
            auto it_node = mrModelPart.NodesBegin() + i_node;
            it_node->SetValue(NODAL_AREA, 0.0);
            it_node->SetValue(AUXILIAR_GRADIENT_X, aux_zero_vector);
            it_node->SetValue(AUXILIAR_GRADIENT_Y, aux_zero_vector);
            it_node->SetValue(AUXILIAR_GRADIENT_Z, aux_zero_vector);
        }

        // Compute gradient VELOCITY_X
        auto gradient_process_x = ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>(mrModelPart, VELOCITY_X, AUXILIAR_GRADIENT_X, NODAL_AREA);
        gradient_process_x.Execute();

        // Compute gradient VELOCITY_Y
        for(int i_node = 0; i_node < number_nodes; ++i_node) {
            auto it_node = mrModelPart.NodesBegin() + i_node;
            it_node->SetValue(NODAL_AREA, 0.0);
        }
        auto gradient_process_y = ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>(mrModelPart, VELOCITY_Y, AUXILIAR_GRADIENT_Y, NODAL_AREA);
        gradient_process_y.Execute();

        // Compute gradient VELOCITY_Z
        for(int i_node = 0; i_node < number_nodes; ++i_node) {
            auto it_node = mrModelPart.NodesBegin() + i_node;
            it_node->SetValue(NODAL_AREA, 0.0);
        }
        auto gradient_process_z = ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>(mrModelPart, VELOCITY_Z, AUXILIAR_GRADIENT_Z, NODAL_AREA);
        gradient_process_z.Execute();

        for(int i_node = 0; i_node < number_nodes; ++i_node){
            auto it_node = mrModelPart.NodesBegin() + i_node;
            auto aux_gradient_x = it_node->GetValue(AUXILIAR_GRADIENT_X);
            auto aux_gradient_y = it_node->GetValue(AUXILIAR_GRADIENT_Y);
            auto aux_gradient_z = it_node->GetValue(AUXILIAR_GRADIENT_Z);
            auto divergence_old_avg = it_node->GetValue(DIVERGENCE);
            // KRATOS_WATCH(time_step_current);
            // KRATOS_WATCH(time_step_previous);
            // KRATOS_WATCH(it_node->GetValue(DIVERGENCE));
            auto divergence_current = aux_gradient_x[0] + aux_gradient_y[1] + aux_gradient_z[2];
            // Compute divergence weighetd in time average
            auto divergence_current_avg = (time_step_previous*divergence_old_avg +(time_step_current-time_step_previous)*divergence_current)/(time_step_current);
            it_node->SetValue(DIVERGENCE,divergence_current_avg);
            // KRATOS_WATCH(divergence_current);
            // KRATOS_WATCH(it_node->GetValue(DIVERGENCE));
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
        return "DivergenceCalculationProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DivergenceCalculationProcess";
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
    DivergenceCalculationProcess& operator=(DivergenceCalculationProcess const& rOther);

    /// Copy constructor.
    DivergenceCalculationProcess(DivergenceCalculationProcess const& rOther);


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
                                  DivergenceCalculationProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DivergenceCalculationProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



} // namespace Kratos


#endif // KRATOS_DIVERGENCE_CALCULATION_PROCESS_H
