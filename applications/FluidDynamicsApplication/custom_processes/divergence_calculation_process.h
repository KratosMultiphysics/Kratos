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

        // Declare auxiliar doubles
        double divergence_avg;
        int number_elements;
        int number_nodes;
        double time_step_previous;
        double time_step_current;

        // Check and set number of elements
        KRATOS_ERROR_IF(mrModelPart.Elements().size() == 0) << "the number of elements in the domain is zero. divergence calculation cannot be applied"<< std::endl;
        number_elements = mrModelPart.Elements().size();
        KRATOS_WATCH(number_elements);

        // Check and set number of nodes
        KRATOS_ERROR_IF(mrModelPart.Nodes().size() == 0) << "the number of nodes in the domain is zero. divergence calculation cannot be applied" << std::endl;
        number_nodes = mrModelPart.Nodes().size();
        KRATOS_WATCH(number_nodes);

        // Geometry information
        const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
        KRATOS_WATCH(dimension)

        // Initialize auxiliar variables
        array_1d<double, 3> aux_zero_vector = ZeroVector(3);
        for(int i_node = 0; i_node < number_nodes; ++i_node) {
            auto it_node = mrModelPart.NodesBegin() + i_node;
            it_node->SetValue(NODAL_AREA, 0.0);
            it_node->SetValue(AUXILIAR_GRADIENT, aux_zero_vector);
        }

        // Compute gradient VELOCITY_X
        auto gradient_process = ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>(mrModelPart, VELOCITY_X, AUXILIAR_GRADIENT, NODAL_AREA);
        gradient_process.Execute();
        // Store gradient into auxiliary array

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
