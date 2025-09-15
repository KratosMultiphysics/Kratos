//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
//
#if !defined(KRATOS_APPLY_CHIMERA_FRACTIONALSTEP_H_INCLUDED)
#define KRATOS_APPLY_CHIMERA_FRACTIONALSTEP_H_INCLUDED

// System includes

// External includes

// Project includes

// Application includes
#include "apply_chimera_process.h"

namespace Kratos {

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
 * @class ApplyChimeraProcessFractionalStep
 *
 * @ingroup ChimeraApplication
 *
 * @brief This class extends ApplyChimera base class and overwrites the function ApplyContinuityWithMpcs to use different containers for storing pressure and velocity constraints.
 *
*/

template <int TDim>
class KRATOS_API(CHIMERA_APPLICATION) ApplyChimeraProcessFractionalStep
    : public ApplyChimera<TDim> {
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ApplyChimeraProcessFractionalStep
    KRATOS_CLASS_POINTER_DEFINITION(ApplyChimeraProcessFractionalStep);
    typedef ApplyChimera<TDim> BaseType;
    typedef typename BaseType::MasterSlaveConstraintContainerType MasterSlaveConstraintContainerType;
    typedef typename BaseType::PointLocatorType PointLocatorType;
    typedef typename BaseType::PointLocatorPointerType PointLocatorPointerType;
    typedef typename BaseType::MasterSlaveContainerVectorType MasterSlaveContainerVectorType;
    typedef typename BaseType::ConstraintIdsVectorType ConstraintIdsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{
    ApplyChimeraProcessFractionalStep(ModelPart& rMainModelPart, Parameters iParameters);

    /// Destructor.
    ~ApplyChimeraProcessFractionalStep() = default;

    virtual void ExecuteFinalizeSolutionStep() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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

    /**
     * @brief Applies the continuity between the boundary modelpart and the
     * background.
     * @param rBoundaryModelPart The boundary modelpart for which the continuity
     * is to be enforced.
     * @param rBinLocator The bin based locator formulated on the background.
     * This is used to locate nodes on rBoundaryModelPart.
     */
    void ApplyContinuityWithMpcs(ModelPart& rBoundaryModelPart,
                                 PointLocatorType& rBinLocator) override;

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

    /// Assignment operator.
    ApplyChimeraProcessFractionalStep& operator=(ApplyChimeraProcessFractionalStep const& rOther);

    ///@}

}; // Class ApplyChimeraProcessFractionalStep

} // namespace Kratos.

#endif // KRATOS_APPLY_CHIMERA_FRACTIONALSTEP_H_INCLUDED
