//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#ifndef KRATOS_EMBEDDED_NODES_INITIALIZATION_PROCESS_H
#define KRATOS_EMBEDDED_NODES_INITIALIZATION_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

// Application includes


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Utility to initialize velocity and pressure in those nodes whose
/// domain has changed because of the distance function movement
class EmbeddedNodesInitializationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EmbeddedNodesInitializationProcess
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedNodesInitializationProcess);

    typedef Node                     NodeType;
    typedef NodeType::Pointer    NodePointerType;
    typedef Geometry<NodeType>      GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EmbeddedNodesInitializationProcess(ModelPart& rModelPart, unsigned int MaxIterations = 10) : mrModelPart(rModelPart)
    {
        mMaxIterations = MaxIterations;
    }

    /// Destructor.
    ~EmbeddedNodesInitializationProcess() override{}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    void ExecuteInitializeSolutionStep() override
    {
        const unsigned int BufferSize = mrModelPart.GetBufferSize();
        ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();
        ModelPart::ElementsContainerType& rElements = mrModelPart.Elements();

        // Simple check
        if( mrModelPart.NodesBegin()->SolutionStepsDataHas( DISTANCE ) == false )
            KRATOS_ERROR << "Nodes do not have DISTANCE variable!";
        if( mrModelPart.NodesBegin()->SolutionStepsDataHas( PRESSURE ) == false )
            KRATOS_ERROR << "Nodes do not have PRESSURE variable!";
        if( mrModelPart.NodesBegin()->SolutionStepsDataHas( VELOCITY ) == false )
            KRATOS_ERROR << "Nodes do not have VELOCITY variable!";

        // Mark the nodes that have switched from structure to fluid during the last step
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(rNodes.size()); ++k)
        {
            ModelPart::NodesContainerType::iterator itNode = rNodes.begin() + k;
            const double& d  = itNode->FastGetSolutionStepValue(DISTANCE);      // Current step distance value
            const double& dn = itNode->FastGetSolutionStepValue(DISTANCE,1);    // Previous step distance value

            if ((d>0) && (dn<=0))
                itNode->Set(SELECTED, true);
        }

        for (unsigned int it=0; it<mMaxIterations; ++it)
        {
            // Loop along the elements to find which ones have a unique selected node
            for (int k = 0; k < static_cast<int>(rElements.size()); ++k)
            {
                unsigned int NewNodes = 0;
                ModelPart::ElementsContainerType::iterator itElement = rElements.begin() + k;
                const GeometryType& rGeometry = itElement->GetGeometry();
                const unsigned int ElemNumNodes = rGeometry.PointsNumber();

                // Get the number of nodes that have switched from structure to fluid in the current element
                for (unsigned int j=0; j<ElemNumNodes; ++j)
                {
                    if (rGeometry[j].Is(SELECTED))
                        NewNodes++;
                }

                // If there is only one unique "new" node in the element, initialize it.
                // Otherwise it remains to be initialized in a further iteration.
                if (NewNodes == 1)
                {
                    NodeType::Pointer pNode;
                    double p_avg = 0.0;
                    array_1d<double, 3> v_avg = ZeroVector(3);

                    for (unsigned int j=0; j<ElemNumNodes; ++j)
                    {
                        if (rGeometry[j].IsNot(SELECTED))
                        {
                            // Compute the average velocity in the non selected nodes
                            p_avg += rGeometry[j].FastGetSolutionStepValue(PRESSURE);
                            v_avg += rGeometry[j].FastGetSolutionStepValue(VELOCITY);
                        }
                        else
                        {
                            // Get a pointer to the unique SELECTED node
                            pNode = rGeometry(j);
                        }
                    }

                    // Compute the non-SELECTED nodes average values
                    p_avg /= (ElemNumNodes-1);
                    v_avg /= (ElemNumNodes-1);

                    // Historical values initialization
                    pNode->Set(SELECTED, false);                                    // Once a node has been initialized it is marked as non SELECTED
                    for (unsigned int step=0; step<BufferSize; ++step)              // Fill the velocity and pressure buffer
                    {
                        pNode->FastGetSolutionStepValue(PRESSURE, step) = p_avg;
                        pNode->FastGetSolutionStepValue(VELOCITY, step) = v_avg;
                    }
                }
            }
        }

        // If there still exist some remaining nodes, initialize their values to zero
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(rNodes.size()); ++k)
        {
            ModelPart::NodesContainerType::iterator itNode = rNodes.begin() + k;

            if (itNode->Is(SELECTED))
            {
                itNode->Set(SELECTED, false);                                   // Once a node has been initialized it is marked as non SELECTED
                for (unsigned int step=0; step<BufferSize; ++step)              // Fill the velocity and pressure buffer
                {
                    itNode->FastGetSolutionStepValue(PRESSURE, step) = 0.0;
                    itNode->FastGetSolutionStepValue(VELOCITY, step) = ZeroVector(3);
                }
            }
        }
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "EmbeddedNodesInitializationProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "EmbeddedNodesInitializationProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ModelPart&                                 mrModelPart;
    unsigned int                            mMaxIterations;

    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Default constructor.
    EmbeddedNodesInitializationProcess() = delete;

    /// Assignment operator.
    EmbeddedNodesInitializationProcess& operator=(EmbeddedNodesInitializationProcess const& rOther) = delete;

    /// Copy constructor.
    EmbeddedNodesInitializationProcess(EmbeddedNodesInitializationProcess const& rOther) = delete;


    ///@}

}; // Class EmbeddedNodesInitializationProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_EMBEDDED_NODES_INITIALIZATION_PROCESS_H
