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

#ifndef KRATOS_DISTANCE_MODIFICATION_PROCESS_H
#define KRATOS_DISTANCE_MODIFICATION_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/cfd_variables.h"
#include "processes/find_nodal_h_process.h"

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

/// Utility to modify the distances of an embedded object in order to avoid bad intersections
/// Besides, it also deactivate the full negative distance elements
class DistanceModificationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DistanceModificationProcess
    KRATOS_CLASS_POINTER_DEFINITION(DistanceModificationProcess);

    typedef Node<3>                     NodeType;
    typedef Geometry<NodeType>      GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DistanceModificationProcess(ModelPart& rModelPart, const bool& rCheckAtEachStep)
    {
        mFactorCoeff = 2.0;
        mrModelPart = rModelPart;
        mrCheckAtEachStep = rCheckAtEachStep;
    }

    /// Destructor.
    virtual ~DistanceModificationProcess(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    void ExecuteInitialize() override
    {
        KRATOS_TRY;

        // Obtain NODAL_H values
        FindNodalHProcess NodalHCalculator(mrModelPart);
        NodalHCalculator.Execute();

        KRATOS_CATCH("");
    }


    void ExecuteBeforeSolutionLoop() override
    {
        KRATOS_TRY;

        unsigned int counter = 1;
        unsigned int bad_cuts = 1;
        double factor = 0.01;

        // Modify the nodal distance values until there is no bad intersections
        while (bad_cuts > 0)
        {
            this->ModifyDistance(factor, bad_cuts);
            std::cout << "Distance modification iteration: " << counter << " Total bad cuts: " << bad_cuts << " Factor: " << factor << std::endl;
            factor /= mFactorCoeff;
            counter++;
        }

        this->DeactivateFullNegativeElements();

        KRATOS_CATCH("");
    }


    void ExecuteInitializeSolutionStep() override
    {
        if(mrCheckAtEachStep == true)
        {
            ExecuteBeforeSolutionLoop();
        }
    }


    void ExecuteFinalize() override
    {
        if(mrCheckAtEachStep == true)
        {
            RecoverOriginalDistance();
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
        buffer << "DistanceModificationProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "DistanceModificationProcess";}

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
    double                                    mFactorCoeff;
    bool                                 mrCheckAtEachStep;
    std::vector<unsigned int>        mModifiedDistancesIDs;
    std::vector<double>           mModifiedDistancesValues;

    ///@}
    ///@name Protected Operators
    ///@{

    void ModifyDistance(const double& factor,
                        unsigned int& bad_cuts)
    {
        double tol_d;
        bad_cuts = 0;

        // Simple check
        if( mrModelPart.NodesBegin()->SolutionStepsDataHas( DISTANCE ) == false )
            KRATOS_ERROR << "Nodes do not have DISTANCE variable!";
        if( mrModelPart.NodesBegin()->SolutionStepsDataHas( NODAL_H ) == false )
            KRATOS_ERROR << "Nodes do not have NODAL_H variable!";

        // Distance modification
        if (mrCheckAtEachStep == false) // Case in where the original distance does not need to be recomputed (e.g. CFD)
        {
            for (auto itNode=mrModelPart.NodesBegin(); itNode!=mrModelPart.NodesEnd(); itNode++)
            {
                const double h = itNode->FastGetSolutionStepValue(NODAL_H);
                double& d = itNode->FastGetSolutionStepValue(DISTANCE);
                tol_d = factor*h;

                if((d >= 0.0) && (d < tol_d))
                {
                    // Modify the distance to avoid almost empty fluid elements
                    std::cout << "Node: " << itNode->Id() << " distance " << d;
                    d = -0.001*tol_d;
                    std::cout << " modified to " << d << std::endl;
                }
            }

            // Syncronize data between partitions (the modified distance has always a lower value)
            mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(DISTANCE);
        }
        else // Case in where the original distance needs to be kept to track the interface (e.g. FSI)
        {
            for (auto itNode=mrModelPart.NodesBegin(); itNode!=mrModelPart.NodesEnd(); itNode++)
            {
                const double h = itNode->FastGetSolutionStepValue(NODAL_H);
                double& d = itNode->FastGetSolutionStepValue(DISTANCE);
                tol_d = factor*h;

                if((d >= 0.0) && (d < tol_d))
                {
                    // Store the original distance to be recovered at the end of the step
                    mModifiedDistancesIDs.push_back(itNode->Id());
                    mModifiedDistancesValues.push_back(d);

                    // Modify the distance to avoid almost empty fluid elements
                    std::cout << "Node: " << itNode->Id() << " distance " << d;
                    d = -0.001*tol_d;
                    std::cout << " modified to " << d << std::endl;
                }
            }

            // Syncronize data between partitions (the modified distance has always a lower value)
            mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(DISTANCE);
        }

        // Check if there still exist bad cuts
        for (auto itElement=mrModelPart.ElementsBegin(); itElement!=mrModelPart.ElementsEnd(); itElement++)
        {
            unsigned int npos = 0;
            unsigned int nneg = 0;

            GeometryType& rGeometry = itElement->GetGeometry();

            for (unsigned int itNode=0; itNode<rGeometry.size(); itNode++)
            {
                double d = rGeometry[itNode].FastGetSolutionStepValue(DISTANCE);

                if(d > 0.0)
                {
                    npos++;
                }
                else
                {
                    nneg++;
                }
            }

            if((nneg > 0) && (npos > 0)) // The element is cut
            {
                for(unsigned int itNode=0; itNode<rGeometry.size(); itNode++)
                {
                    const double h = rGeometry[itNode].GetValue(NODAL_H);
                    double d = rGeometry[itNode].FastGetSolutionStepValue(DISTANCE);
                    tol_d = (factor*mFactorCoeff)*h;

                    if((d >= 0.0) && (d < tol_d))
                    {
                        bad_cuts++;
                        break;
                    }
                }
            }
        }
    }


    void RecoverOriginalDistance()
    {
        for(unsigned int i=0; i<mModifiedDistancesIDs.size(); i++)
        {
            unsigned int nodeId = mModifiedDistancesIDs[i];
            mrModelPart.GetNode(nodeId).FastGetSolutionStepValue(DISTANCE) = mModifiedDistancesValues[i];
        }

        // Empty the modified distance vectors
        mModifiedDistancesIDs.resize(0);
        mModifiedDistancesValues.resize(0);
        mModifiedDistancesIDs.shrink_to_fit();
        mModifiedDistancesValues.shrink_to_fit();

    }


    void DeactivateFullNegativeElements()
    {
        // Deactivate the full negative distance elements
        for (auto itElement=mrModelPart.ElementsBegin(); itElement!=mrModelPart.ElementsEnd(); itElement++)
        {
            unsigned int fixed = 0;
            unsigned int inside = 0;

            GeometryType& rGeometry = itElement->GetGeometry();
            const unsigned int NumNodes = rGeometry.size();

            // Check the distance function sign and fixity at the element nodes
            if (rGeometry.Dimension() == 2)
            {
                for (unsigned int itNode=0; itNode<NumNodes; itNode++)
                {
                    if (rGeometry[itNode].GetSolutionStepValue(DISTANCE)<0.0)
                        inside++;
                    if (rGeometry[itNode].IsFixed(VELOCITY_X) && rGeometry[itNode].IsFixed(VELOCITY_Y))
                        fixed++;
                }
            } else {
                for (unsigned int itNode=0; itNode<NumNodes; itNode++)
                {
                    if (rGeometry[itNode].GetSolutionStepValue(DISTANCE)<0.0)
                        inside++;
                    if (rGeometry[itNode].IsFixed(VELOCITY_X) && rGeometry[itNode].IsFixed(VELOCITY_Y) && rGeometry[itNode].IsFixed(VELOCITY_Z))
                        fixed++;
                }
            }

            // If all the nodes have negative distance value (non-fluid domain) deactivate the element
            // If the sum of inside nodes (negative distance) and fixed nodes equals the number of nodes,
            // deactivate the element as well. In this way non-well defined elements are avoided.
            if ((inside == NumNodes) || (inside + fixed == NumNodes))
            {
                itElement->Set(ACTIVE, false);

                // If deactivated, set element DOFs to zero
                const array_1d<double,3> auxVec = ZeroVector(3);
                for (unsigned int i=0; i<NumNodes; ++i)
                {
                    double& pres = rGeometry[i].GetSolutionStepValue(PRESSURE);
                    array_1d<double,3>& vel = rGeometry[i].GetSolutionStepValue(VELOCITY);

                    pres = 0.0;
                    vel = auxVec;
                }
            }
            // Otherwise, activate the element (it might have been deactivated in a previous time step)
            else
            {
                itElement->Set(ACTIVE, true);
            }
        }
    }

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
    DistanceModificationProcess(){}

    /// Assignment operator.
    DistanceModificationProcess& operator=(DistanceModificationProcess const& rOther){return *this;}

    /// Copy constructor.
    DistanceModificationProcess(DistanceModificationProcess const& rOther){}


    ///@}

}; // Class DistanceModificationProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_DISTANCE_MODIFICATION_PROCESS_H
