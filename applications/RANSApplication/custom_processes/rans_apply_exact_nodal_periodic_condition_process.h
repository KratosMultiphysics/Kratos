//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
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

///@name Kratos Classes
///@{

/**
 * @brief A process to create periodic conditions between two boundaries
 *
 * This process does node to node matching to create periodic boundary
 * conditions between the master and slave model parts, to be used in ResidualBasedBlockBuilderAndSolverPeriodic.
 *
 * This process requires two boundaries to have maching meshes.
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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansApplyExactNodalPeriodicConditionProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~RansApplyExactNodalPeriodicConditionProcess() override = default;

    /// Assignment operator.
    RansApplyExactNodalPeriodicConditionProcess& operator=(RansApplyExactNodalPeriodicConditionProcess const& rOther) = delete;

    /// Copy constructor.
    RansApplyExactNodalPeriodicConditionProcess(RansApplyExactNodalPeriodicConditionProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override;

    const Parameters GetDefaultParameters() const override;

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

private:
    ///@name Member Variables
    ///@{

    Model& mrModel;

    int mEchoLevel;

    std::string mMasterModelPartName;
    std::string mSlaveModelPartName;

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
    ///@name Private Operations
    ///@{

    /**
     * @brief Create a Periodic Conditions
     *
     * This method iterates through master and slave model parts nodes, and find matching
     * slave node in slave model part for its master node in master model part. If translation
     * is given with magnitude greater than zero, then first master node is translated with direction and magnitude
     * then if rotation is given this translated coordinates are rotated to locate slave node in slave model part.
     *
     * Afterwards, it creates two node Periodic conditions (in 2D and 3D) in the base model part.
     *
     */
    void CreatePeriodicConditions();

    /**
     * @brief The method calculates rotated position
     *
     * This method returns rotated coordinates of rInitialPosition. It rotates around
     * mRotationCenter in the mRotationAxis for given mRotationAngle
     *
     * @param rInitialPosition       Initial position to be rotated
     * @return array_1d<double, 3>   Rotated initial position
     */
    array_1d<double, 3> CalculateRotatedPosition(
        const array_1d<double, 3>& rInitialPosition) const;

    /**
     * @brief Calculates rotation matrix
     *
     * Calculates rotation matrix based on mRotationAxis and mRotationAngle
     *
     * @param rOutput               Bounded rotation matrix
     */
    void CalculateRotationMatrix(
        BoundedMatrix<double, 3, 3>& rOutput) const;

    ///@}

}; // Class RansApplyExactNodalPeriodicConditionProcess

///@}

///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansApplyExactNodalPeriodicConditionProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_APPLY_EXACT_NODAL_PERIODIC_CONDITION_PROCESS_H_INCLUDED defined
