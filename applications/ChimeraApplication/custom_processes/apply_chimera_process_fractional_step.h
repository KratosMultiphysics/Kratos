// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
// 					 Navaneeth K Narayanan
//					 Rishith Ellath Meethal
// ==============================================================================
//
#if !defined(KRATOS_APPLY_CHIMERA_FRACTIONALSTEP_H_INCLUDED)
#define KRATOS_APPLY_CHIMERA_FRACTIONALSTEP_H_INCLUDED

// System includes

// External includes

// Project includes

// Application includes
#include "apply_chimera_process_monolithic.h"

namespace Kratos
{

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

template <int TDim, class TDistanceCalculatorType>
class ApplyChimeraProcessFractionalStep : public ApplyChimeraProcessMonolithic<TDim, TDistanceCalculatorType>
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ApplyChimeraProcessFractionalStep
    KRATOS_CLASS_POINTER_DEFINITION(ApplyChimeraProcessFractionalStep);
    typedef ApplyChimeraProcessMonolithic<TDim, TDistanceCalculatorType>            BaseType;
    typedef typename BaseType::MasterSlaveConstraintContainerType                   MasterSlaveConstraintContainerType;
    typedef typename BaseType::PointLocatorType                                     PointLocatorType;
    typedef typename BaseType::PointLocatorPointerType                              PointLocatorPointerType;
    typedef typename BaseType::MasterSlaveContainerVectorType                       MasterSlaveContainerVectorType;
    typedef typename BaseType::ConstraintIdsVectorType                              ConstraintIdsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{
    ApplyChimeraProcessFractionalStep(ModelPart &rMainModelPart, Parameters iParameters)
        : BaseType(rMainModelPart, iParameters)
    {
    }


    /// Destructor.
    virtual ~ApplyChimeraProcessFractionalStep()
    {
    }
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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
     * @brief Applies the continuity between the boundary modelpart and the background.
     * @param rBoundaryModelPart The boundary modelpart for which the continuity is to be enforced.
     * @param pBinLocator The bin based locator formulated on the background. This is used to locate nodes on rBoundaryModelPart.
     */
    void ApplyContinuityWithMpcs(ModelPart &rBoundaryModelPart, PointLocatorPointerType &pBinLocator) override;

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
    ApplyChimeraProcessFractionalStep &operator=(ApplyChimeraProcessFractionalStep const &rOther);

    ///@}




}; // Class ApplyChimeraProcessFractionalStep

} // namespace Kratos.

#endif // KRATOS_APPLY_CHIMERA_FRACTIONALSTEP_H_INCLUDED
