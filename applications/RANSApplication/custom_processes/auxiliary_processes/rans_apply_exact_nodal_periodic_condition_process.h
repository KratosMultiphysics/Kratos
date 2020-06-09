//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_APPLY_EXACT_NODAL_PERIODIC_CONDITION_PROCESS_H_INCLUDED)
#define KRATOS_RANS_APPLY_EXACT_NODAL_PERIODIC_CONDITION_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"

namespace Kratos
{
///@addtogroup RANSApplication
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

/**
 * @brief A process to create periodic conditions between two boundaries
 *
 * This process does node to node matching to create periodic boundary
 * conditions between the master and slave model parts, to be used in ResidualBasedBlockBuilderAndSolverPeriodic.
 *
 * This process requires two boundaries to have maching meshes.
 * Make sure to avoid calling "ReplaceElementsAndConditions" method call, and directly impose element
 * and condition name in the "*.mdpa" file in order to correctly use this process alongside with
 * ResidualBasedBlockBuilderAndSolverPeriodic
 *
 * @see class ResidualBasedBlockBuilderAndSolverPeriodic
 */
class KRATOS_API(RANS_APPLICATION) RansApplyExactNodalPeriodicConditionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansApplyExactNodalPeriodicConditionProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansApplyExactNodalPeriodicConditionProcess);

    using NodeType = ModelPart::NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansApplyExactNodalPeriodicConditionProcess(Model& rModel, Parameters rParameters);

    /// Destructor.
    ~RansApplyExactNodalPeriodicConditionProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitialize() override;

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
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

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

    ///@}
    ///@name Protected Operators
    ///@{

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

    Model& mrModel;
    Parameters mrParameters;
    int mEchoLevel;

    std::string mBaseModelPartName;
    std::string mMasterModelPartName;
    std::string mSlaveModelPartName;

    std::vector<std::string> mVariablesList;
    double mTolerance;

    // translation settings
    array_1d<double, 3> mTranslationDirection;
    double mTranslationMagnitude;

    // rotation settings
    array_1d<double, 3> mRotationAxis;
    array_1d<double, 3> mRotationCenter;
    double mRotationAngle;

    bool mReorder;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CreatePeriodicConditions();

    array_1d<double, 3> CalculateRotatedPosition(const array_1d<double, 3>& rInitialPosition) const;

    void CalculateRotationMatrix(BoundedMatrix<double, 3, 3>& rOutput) const;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    RansApplyExactNodalPeriodicConditionProcess& operator=(
        RansApplyExactNodalPeriodicConditionProcess const& rOther);

    /// Copy constructor.
    RansApplyExactNodalPeriodicConditionProcess(RansApplyExactNodalPeriodicConditionProcess const& rOther);

    ///@}

}; // Class RansApplyExactNodalPeriodicConditionProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansApplyExactNodalPeriodicConditionProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_APPLY_EXACT_NODAL_PERIODIC_CONDITION_PROCESS_H_INCLUDED defined
