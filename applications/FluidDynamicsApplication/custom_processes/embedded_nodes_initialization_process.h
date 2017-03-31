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
// #include "processes/find_nodal_h_process.h"

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

    typedef Node<3>                     NodeType;
    typedef NodeType::Pointer    NodePointerType;
    typedef Geometry<NodeType>      GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EmbeddedNodesInitializationProcess(ModelPart& rModelPart, unsigned int MaxIterations = 10)
    {
        mrModelPart = rModelPart;
        mMaxIterations = MaxIterations;
    }

    /// Destructor.
    virtual ~EmbeddedNodesInitializationProcess(){}

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
        // Mark the nodes that have switched from structure to fluid during the last step
        for (ModelPart::NodeIterator itNode=mrModelPart.NodesBegin(); itNode!=mrModelPart.NodesEnd(); ++itNode)
        {
            const double& d  = itNode->FastGetSolutionStepValue(DISTANCE);      // Current step distance value
            const double& dn = itNode->FastGetSolutionStepValue(DISTANCE,1);    // Previous step distance value

            if ((d>0) && (dn<=0))
            {
                itNode->Set(SELECTED, true);
            }
        }

        for (unsigned int it=0; it<mMaxIterations; ++it)
        {
            // Loop along the elements to find which ones have a unique selected node
            for (ModelPart::ElementIterator itElement=mrModelPart.ElementsBegin(); itElement<mrModelPart.ElementsEnd(); ++itElement)
            {

                unsigned int NewNodes = 0;
                GeometryType& rGeometry = itElement->GetGeometry();
                const unsigned int ElemNumNodes = rGeometry.PointsNumber();

                // Get the number of nodes that have switch from structure to fluid in the current element
                for (unsigned int j=0; j<ElemNumNodes; ++j)
                {
                    if (rGeometry[j].Is(SELECTED))
                    {
                        NewNodes++;
                    }
                }

                // std::cout << "Element number " << itElement->Id() << " has " << NewNodes << " nodes switching from structure to fluid." << std::endl;

                // If there is only one unique "new" node in the element, initialize it
                if (NewNodes == 1)
                {
                    NodeType::Pointer pNode;
                    double p_avg = 0.0;
                    array_1d<double, 3> v_avg;
                    v_avg[0] = 0.0;
                    v_avg[1] = 0.0;
                    v_avg[2] = 0.0;

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
                            NodeType::Pointer pAux = rGeometry(j);
                            std::swap(pAux, pNode);
                        }
                    }

                    p_avg /= (ElemNumNodes-1);
                    v_avg /= (ElemNumNodes-1);

                    // std::cout << "Structure to fluid initialization of node " << NewId << " in iteration number " << it << "." << std::endl;
                    // Once a node is initialized it is marked as non SELECTED,
                    // disregarding the fact that it can be shared by other elements
                    pNode->Set(SELECTED, false);
                    // Historical values initialization
                    pNode->FastGetSolutionStepValue(PRESSURE, 0) = p_avg;
                    pNode->FastGetSolutionStepValue(PRESSURE, 1) = p_avg;
                    pNode->FastGetSolutionStepValue(PRESSURE, 2) = p_avg;
                    pNode->FastGetSolutionStepValue(VELOCITY, 0) = v_avg;
                    pNode->FastGetSolutionStepValue(VELOCITY, 1) = v_avg;
                    pNode->FastGetSolutionStepValue(VELOCITY, 2) = v_avg;

                }
            }
        }

        array_1d<double, 3> v_null;
        v_null[0] = 0.0;
        v_null[1] = 0.0;
        v_null[2] = 0.0;

        // If there still exist some remaining nodes, initialize their values to zero
        for (ModelPart::NodeIterator itNode=mrModelPart.NodesBegin(); itNode!=mrModelPart.NodesEnd(); ++itNode)
        {
            if (itNode->Is(SELECTED))
            {
                itNode->Set(SELECTED, false);

                // std::cout << "Node " << itNode->Id() << " remained without structure to fluid initialization. Initializing it to 0." << std::endl;

                itNode->FastGetSolutionStepValue(PRESSURE, 0) = 0.0;
                itNode->FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
                itNode->FastGetSolutionStepValue(PRESSURE, 2) = 0.0;
                itNode->FastGetSolutionStepValue(VELOCITY, 0) = v_null;
                itNode->FastGetSolutionStepValue(VELOCITY, 1) = v_null;
                itNode->FastGetSolutionStepValue(VELOCITY, 2) = v_null;
            }
        }
    }

    // void ExecuteFinalize() override
    void ExecuteFinalizeSolutionStep() override
    {
        const array_1d<double, 3> aux_zero = ZeroVector(3);

        // Store the positive distance nodes VELOCITY in EMBEDDED_WET_VELOCITY variable. The negative distance
        // nodes EMBEDDED_WET_VELOCITY is set to zero for visualization purposes
        // The same is done for the PRESSURE using the EMBEDDED_WET_PRESSURE variable
        for (ModelPart::NodeIterator itNode=mrModelPart.NodesBegin(); itNode!=mrModelPart.NodesEnd(); ++itNode)
        {
            const double dist = itNode->FastGetSolutionStepValue(DISTANCE);
            double& emb_wet_pres = itNode->FastGetSolutionStepValue(EMBEDDED_WET_PRESSURE);
            array_1d<double, 3>& emb_wet_vel = itNode->FastGetSolutionStepValue(EMBEDDED_WET_VELOCITY);

            if (dist >= 0.0)
            {
                emb_wet_pres = itNode->FastGetSolutionStepValue(PRESSURE);
                emb_wet_vel = itNode->FastGetSolutionStepValue(VELOCITY);
            }
            else
            {
                emb_wet_pres = 0.0;
                emb_wet_vel = aux_zero;
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
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "EmbeddedNodesInitializationProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "EmbeddedNodesInitializationProcess";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}


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

    ModelPart                                  mrModelPart;
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
    EmbeddedNodesInitializationProcess(){}

    /// Assignment operator.
    EmbeddedNodesInitializationProcess& operator=(EmbeddedNodesInitializationProcess const& rOther){return *this;}

    /// Copy constructor.
    EmbeddedNodesInitializationProcess(EmbeddedNodesInitializationProcess const& rOther){}


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
