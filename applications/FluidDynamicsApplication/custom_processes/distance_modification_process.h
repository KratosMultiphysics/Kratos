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
// #include "includes/define.h"
#include "processes/process.h"
#include "includes/cfd_variables.h"
#include "processes/find_nodal_h_process.h"
#include "utilities/openmp_utils.h"
#include "includes/kratos_parameters.h"

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
    DistanceModificationProcess(
        ModelPart& rModelPart, 
        const bool CheckAtEachStep, 
        const bool NegElemDeactivation,
        const bool RecoverOriginalDistance)
        : Process(), mrModelPart(rModelPart) {

        mFactorCoeff = 2.0;
        mCheckAtEachStep = CheckAtEachStep;
        mNegElemDeactivation = NegElemDeactivation;
        mRecoverOriginalDistance = RecoverOriginalDistance;
    }

    DistanceModificationProcess(
        ModelPart& rModelPart,
        Parameters& rParameters)
        : Process(), mrModelPart(rModelPart) {

        Parameters default_parameters( R"(
        {
            "mesh_id"                                : 0,
            "model_part_name"                        : "default_model_part_name",
            "distance_factor"                        : 2.0,
            "check_at_each_time_step"                : false,
            "deactivate_full_negative_elements"      : true,
            "recover_original_distance_at_each_step" : false
        }  )" );

        rParameters.ValidateAndAssignDefaults(default_parameters);

        mFactorCoeff = rParameters["distance_factor"].GetDouble();
        mCheckAtEachStep = rParameters["check_at_each_time_step"].GetBool();
        mNegElemDeactivation = rParameters["deactivate_full_negative_elements"].GetBool();
        mRecoverOriginalDistance = rParameters["recover_original_distance_at_each_step"].GetBool();
    }

    /// Destructor.
    ~DistanceModificationProcess() override{}

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
            factor /= mFactorCoeff;
            counter++;
        }

        if (mNegElemDeactivation) {
            this->DeactivateFullNegativeElements();
        }

        KRATOS_CATCH("");
    }


    void ExecuteInitializeSolutionStep() override
    {
        if(mCheckAtEachStep == true)
        {
            ExecuteBeforeSolutionLoop();
        }
    }


    void ExecuteFinalizeSolutionStep() override
    {
        if(mRecoverOriginalDistance == true)
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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "DistanceModificationProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "DistanceModificationProcess";}

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

    ModelPart&                                              mrModelPart;
    double                                                 mFactorCoeff;
    bool                                               mCheckAtEachStep;
    bool                                           mNegElemDeactivation;
    bool                                       mRecoverOriginalDistance;
    std::vector<std::vector<unsigned int>>        mModifiedDistancesIDs;
    std::vector<std::vector<double>>           mModifiedDistancesValues;

    ///@}
    ///@name Protected Operators
    ///@{

    void ModifyDistance(const double& factor,
                        unsigned int& bad_cuts)
    {
        ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();
        ModelPart::ElementsContainerType& rElements = mrModelPart.Elements();

        // Simple check
        if( mrModelPart.NodesBegin()->SolutionStepsDataHas( DISTANCE ) == false )
            KRATOS_ERROR << "Nodes do not have DISTANCE variable!";
        if( mrModelPart.NodesBegin()->SolutionStepsDataHas( NODAL_H ) == false )
            KRATOS_ERROR << "Nodes do not have NODAL_H variable!";

        // Distance modification
        if (mRecoverOriginalDistance == false) // Case in where the original distance does not need to be preserved (e.g. CFD)
        {
            #pragma omp parallel for
            for (int k = 0; k < static_cast<int>(rNodes.size()); ++k)
            {
                ModelPart::NodesContainerType::iterator itNode = rNodes.begin() + k;
                const double h = itNode->FastGetSolutionStepValue(NODAL_H);
                const double tol_d = factor*h;
                double& d = itNode->FastGetSolutionStepValue(DISTANCE);

                // Modify the distance to avoid almost empty fluid elements
                if((d >= 0.0) && (d < tol_d))
                    d = -0.001*tol_d;
            }
        }
        else // Case in where the original distance needs to be kept to track the interface (e.g. FSI)
        {
            const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
            std::vector<std::vector<unsigned int>> AuxModifiedDistancesIDs(NumThreads);
            std::vector<std::vector<double>> AuxModifiedDistancesValues(NumThreads);

            #pragma omp parallel shared(AuxModifiedDistancesIDs, AuxModifiedDistancesValues)
            {
                const int ThreadId = OpenMPUtils::ThisThread();             // Get the thread id
                std::vector<unsigned int>   LocalModifiedDistancesIDs;      // Local modified distances nodes id vector
                std::vector<double>      LocalModifiedDistancesValues;      // Local modified distances original values vector

                #pragma omp for
                for (int k = 0; k < static_cast<int>(rNodes.size()); ++k)
                {
                    ModelPart::NodesContainerType::iterator itNode = rNodes.begin() + k;
                    const double h = itNode->FastGetSolutionStepValue(NODAL_H);
                    const double tol_d = factor*h;
                    double& d = itNode->FastGetSolutionStepValue(DISTANCE);

                    if((d >= 0.0) && (d < tol_d))
                    {
                        // Store the original distance to be recovered at the end of the step
                        LocalModifiedDistancesIDs.push_back(d);
                        LocalModifiedDistancesValues.push_back(itNode->Id());

                        // Modify the distance to avoid almost empty fluid elements
                        d = -0.001*tol_d;
                    }
                }

                AuxModifiedDistancesIDs[ThreadId] = LocalModifiedDistancesIDs;
                AuxModifiedDistancesValues[ThreadId] = LocalModifiedDistancesValues;
            }

            mModifiedDistancesIDs = AuxModifiedDistancesIDs;
            mModifiedDistancesValues = AuxModifiedDistancesValues;
        }

        // Syncronize data between partitions (the modified distance has always a lower value)
        mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(DISTANCE);

        // Check if there still exist bad cuts
        unsigned int num_bad_cuts = 0;
        /* Note: I'm defining a temporary variable because 'num_bad_cuts'
        *  instead of writing directly into input argument 'bad_cuts'
        *  because using a reference variable in a reduction pragma does
        *  not compile in MSVC 2015 nor in clang-3.8 (Is it even allowed by omp?)
        */
        #pragma omp parallel for reduction(+ : num_bad_cuts)
        for (int k = 0; k < static_cast<int>(rElements.size()); ++k)
        {
            unsigned int npos = 0;
            unsigned int nneg = 0;

            ModelPart::ElementsContainerType::iterator itElement = rElements.begin() + k;
            GeometryType& rGeometry = itElement->GetGeometry();

            for (unsigned int iNode=0; iNode<rGeometry.size(); iNode++)
            {
                const double d = rGeometry[iNode].FastGetSolutionStepValue(DISTANCE);
                (d > 0.0) ? npos++ : nneg++;
            }

            if((nneg > 0) && (npos > 0)) // The element is cut
            {
                for(unsigned int iNode=0; iNode<rGeometry.size(); iNode++)
                {
                    const Node<3> &rConstNode = rGeometry[iNode];
                    const double h = rConstNode.GetValue(NODAL_H);
                    const double tol_d = (factor*mFactorCoeff)*h;
                    const double d = rConstNode.FastGetSolutionStepValue(DISTANCE);

                    if((d >= 0.0) && (d < tol_d))
                    {
                        num_bad_cuts++;
                        break;
                    }
                }
            }
        }

        bad_cuts = num_bad_cuts;
    }


    void RecoverOriginalDistance()
    {
        #pragma omp parallel
        {
            const int ThreadId = OpenMPUtils::ThisThread();
            const std::vector<unsigned int> LocalModifiedDistancesIDs = mModifiedDistancesIDs[ThreadId];
            const std::vector<double> LocalModifiedDistancesValues = mModifiedDistancesValues[ThreadId];

            for(unsigned int i=0; i<LocalModifiedDistancesIDs.size(); ++i)
            {
                const unsigned int nodeId = LocalModifiedDistancesIDs[i];
                mrModelPart.GetNode(nodeId).FastGetSolutionStepValue(DISTANCE) = LocalModifiedDistancesValues[i];
            }
        }

        // Syncronize data between partitions (the modified distance has always a lower value)
        mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(DISTANCE);

        // Empty the modified distance vectors
        mModifiedDistancesIDs.resize(0);
        mModifiedDistancesValues.resize(0);
        mModifiedDistancesIDs.shrink_to_fit();
        mModifiedDistancesValues.shrink_to_fit();

    }


    void DeactivateFullNegativeElements()
    {
        ModelPart::ElementsContainerType& rElements = mrModelPart.Elements();

        // Deactivate those elements whose fixed nodes and negative distance nodes summation is equal (or larger) to their number of nodes
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(rElements.size()); ++k)
        {
            unsigned int fixed = 0;
            unsigned int inside = 0;
            ModelPart::ElementsContainerType::iterator itElement = rElements.begin() + k;
            GeometryType& rGeometry = itElement->GetGeometry();

            // Check the distance function sign at the element nodes
            for (unsigned int itNode=0; itNode<rGeometry.size(); itNode++)
            {
                if (rGeometry[itNode].GetSolutionStepValue(DISTANCE)<0.0)
                    inside++;
                if (rGeometry[itNode].IsFixed(VELOCITY_X) && rGeometry[itNode].IsFixed(VELOCITY_Y) && rGeometry[itNode].IsFixed(VELOCITY_Z))
                    fixed++;
            }

            (inside+fixed >= rGeometry.size()) ? itElement->Set(ACTIVE, false) : itElement->Set(ACTIVE, true);
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
    DistanceModificationProcess() = delete;

    /// Assignment operator.
    DistanceModificationProcess& operator=(DistanceModificationProcess const& rOther) = delete;

    /// Copy constructor.
    DistanceModificationProcess(DistanceModificationProcess const& rOther) = delete;


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
