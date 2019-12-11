// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//  license:     structural_mechanics_application/license.txt
//
//  Main authors: Mahmoud Zidan
//

#if !defined(KRATOS_FIBER_BEAM_COLUMN_SECTION_H_INCLUDED )
#define  KRATOS_FIBER_BEAM_COLUMN_SECTION_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/element.h"
#include "custom_conditions/fiber_beam_column_uniaxial_fiber.hpp"

namespace Kratos
{

///@}
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
///@name  Kratos Classes
///@{

/**
 * @class FiberBeamColumnSection
 *
 * @brief A 3D section for fiber beam-column element
 * @details
 *
 * @author Mahmoud Zidan
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) FiberBeamColumnSection
{

public:

    ///@name Type Definitions
    ///@{

    typedef Element                                     BaseType;
    typedef BaseType::PropertiesType              PropertiesType;
    typedef BaseType::IndexType                        IndexType;
    typedef BaseType::SizeType                          SizeType;
    typedef BaseType::MatrixType                      MatrixType;
    typedef BaseType::VectorType                      VectorType;

    // by default: integration point is 3-dimensional
    typedef IntegrationPoint<3>             IntegrationPointType;

    ///@}
    ///@name Pointer Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor
    FiberBeamColumnSection(IndexType NewId = 0);

    /// Constructor using an array of nodes
    FiberBeamColumnSection(IndexType NewId, IntegrationPointType integrationPoint, PropertiesType::Pointer pProperties);

    // /// Copy Constructor
    // FiberBeamColumnSection(FiberBeamColumnSection const& rOther);

    // /// Destructor
    // ~FiberBeamColumnSection(){}

    void CalculateLocalFlexibilityMatrix();
    void CalculateBMatrix();
    Matrix GetGlobalFlexibilityMatrix();
    Vector GetGlobalDeformationResiduals();

    bool StateDetermination(const Vector& rElementForceIncrements);
    void ResetResidual();

    void FinalizeSolutionStep();

    // ///@}
    // ///@name Operators
    // ///@{

    // /// Assignment Operator
    // FiberBeamColumnSection & operator=(FiberBeamColumnSection const& rOther);

    ///@}
    ///@name Operations
    ///@{

    void Initialize();

    ///@}
    ///@name Access
    ///@{

    /// number of fibers
    SizeType Size() { return mFibers.size(); }
    double GetWeight() { return mWeight; }
    void SetFibers (std::vector<FiberBeamColumnUniaxialFiber> Fibers) { mFibers = std::move(Fibers); }
    std::vector<FiberBeamColumnUniaxialFiber> const& GetFibers() const { return mFibers; }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

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

    IndexType mId;
    double mPosition;
    double mWeight;
    std::vector<FiberBeamColumnUniaxialFiber> mFibers;
    double mTolerance;

    Matrix mBMatrix = ZeroMatrix(3, 5);
    Matrix mLocalFlexibilityMatrix = ZeroMatrix(3, 3);
    Vector mForces = ZeroVector(3);
    Vector mUnbalanceForces = ZeroVector(3);
    Vector mDeformationResiduals = ZeroVector(3);

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    IndexType Id() const { return mId; }

    ///@}
    ///@name Serialization
    ///@{

    friend Serializer;
    void save(Serializer& rSerializer) const ;//override;
    void load(Serializer& rSerializer) ;//override;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

};  // class FiberBeamColumnSection

/// output stream
inline std::ostream & operator <<(std::ostream& rOStream, const FiberBeamColumnSection& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}  // namespace Kratos


#endif
