//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aron Noordam
//

#if !defined(KRATOS_SET_MOVING_LOAD_PROCESS_H_INCLUDED )
#define  KRATOS_SET_MOVING_LOAD_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@name Kratos Classes
///@{

/// Process to create the animated Eigenvectors
/** This process distributes a load on surface load conditions belonging to a modelpart.
 *  The load is distributed according to the surface area.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SetMovingLoadProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of SetMovingLoadProcess
    KRATOS_CLASS_POINTER_DEFINITION(SetMovingLoadProcess);

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    SetMovingLoadProcess(ModelPart& rModelPart,
                                  Parameters Parameters);

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteFinalizeSolutionStep() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override {
        return "SetMovingLoadProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "SetMovingLoadProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    Parameters mParameters;

    std::vector<Condition> mSortedConditions;
    array_1d<double,3> mLoad;
    double mLoadVelocity;
    double mCurrentDistance;

    ///@}
    ///@name Private Operations
    ///@{

    ///@}

}; // Class SetMovingLoadProcess

///@}

}  // namespace Kratos.

#endif // KRATOS_SET_MOVING_LOAD_PROCESS_H_INCLUDED  defined
