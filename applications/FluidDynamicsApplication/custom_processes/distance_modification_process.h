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
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) DistanceModificationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DistanceModificationProcess
    KRATOS_CLASS_POINTER_DEFINITION(DistanceModificationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DistanceModificationProcess(
        ModelPart& rModelPart, 
        const double FactorCoeff, //TODO: Remove it (here for legacy reasons)
        const double DistanceThreshold,
        const bool CheckAtEachStep, 
        const bool NegElemDeactivation,
        const bool RecoverOriginalDistance);

    /// Constructor with Kratos parameters.
    DistanceModificationProcess(
        ModelPart& rModelPart,
        Parameters& rParameters);

    /// Destructor.
    ~DistanceModificationProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    void ExecuteInitialize() override;

    void ExecuteBeforeSolutionLoop() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteFinalizeSolutionStep() override;

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

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart&                                       mrModelPart;
    double                                    mDistanceThreshold;
    bool                                             mIsModified;
    bool                                     mContinuousDistance;
    bool                                        mCheckAtEachStep;
    bool                                    mNegElemDeactivation;
    bool                               mAvoidAlmostEmptyElements;
    bool                                mRecoverOriginalDistance;
    std::vector<unsigned int>              mModifiedDistancesIDs;
    std::vector<double>                 mModifiedDistancesValues;
    std::vector<Vector>        mModifiedElementalDistancesValues;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void ModifyDistance();

    void ModifyDiscontinuousDistance();

    void RecoverDeactivationPreviousState();

    void RecoverOriginalDistance();
    
    void RecoverOriginalDiscontinuousDistance();

    void DeactivateFullNegativeElements();

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
