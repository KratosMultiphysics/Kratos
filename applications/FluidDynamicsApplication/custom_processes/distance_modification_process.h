//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
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
    typedef Variable<double> ComponentType;

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

    /// Constructor with Kratos model
    DistanceModificationProcess(
        Model& rModel,
        Parameters& rParameters);

    /// Destructor.
    ~DistanceModificationProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    void Execute() override;

    void ExecuteInitialize() override;

    void ExecuteBeforeSolutionLoop() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteFinalizeSolutionStep() override;

    ///@}
    ///@name Operations
    ///@{

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
    std::vector<const Variable<double>*>    mDoubleVariablesList;
    std::vector<const ComponentType*>    mComponentVariablesList;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CheckDefaultsAndProcessSettings(Parameters &rParameters);

    /**
     * @brief Initialize the EMBEDDED_IS_ACTIVE variable
     * This method initializes the non historical variable EMBEDDED_IS_ACTIVE.
     * It needs to be called in the constructor to do a threadsafe initialization
     * of such nodal variable before any other operation is done.
     */
    void InitializeEmbeddedIsActive();

    void ModifyDistance();

    void ModifyDiscontinuousDistance();

    void RecoverDeactivationPreviousState();

    void RecoverOriginalDistance();

    void RecoverOriginalDiscontinuousDistance();

    void DeactivateFullNegativeElements();

    template<class TDistancesVectorType>
    void SetElementToSplitFlag(
        Element &rElem,
        const TDistancesVectorType& rDistancesVector)
    {
        unsigned int n_pos = 0;
        unsigned int n_neg = 0;
        for (double i_dist : rDistancesVector) {
            if (i_dist < 0.0) {
                n_neg++;
            } else {
                n_pos++;
            }
        }
        if (n_neg != 0 && n_pos != 0) {
            rElem.Set(TO_SPLIT, true);
        } else {
            rElem.Set(TO_SPLIT, false);
        }
    }

    void SetContinuousDistanceToSplitFlag();

    void SetDiscontinuousDistanceToSplitFlag();

    /**
     * @brief Reads the variables list specified in the Parameters to be fixed in the elements
     * that are fully negative, storing them in mDoubleVariablesList and mComponentVariablesList.
     * It also checks that the variables and the DOFs are defined in the rmModelPart.
     * @param rVariableStringArray Array containing the variables to be fixed in the full negative elements
    */
    void CheckAndStoreVariablesList(const std::vector<std::string>& rVariableStringArray);

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
