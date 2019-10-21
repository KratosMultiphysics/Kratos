// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
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
 * @class FiberBeamColumnElement3D2N
 *
 * @brief A 3D-2node fiber beam-column element for reinforced concrete modeling
 *
 * @author Mahmoud Zidan
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) FiberBeamColumnSection
{

public:

    ///@name Type Definitions
    ///@{

    typedef Element                                     BaseType;
    typedef BaseType::GeometryType                  GeometryType;
    typedef BaseType::NodesArrayType              NodesArrayType;
    typedef BaseType::PropertiesType              PropertiesType;
    typedef BaseType::IndexType                        IndexType;
    typedef BaseType::SizeType                          SizeType;
    typedef BaseType::MatrixType                      MatrixType;
    typedef BaseType::VectorType                      VectorType;
    typedef BaseType::EquationIdVectorType  EquationIdVectorType;
    typedef BaseType::DofsVectorType              DofsVectorType;

    ///@}
    ///@name Pointer Definitions
    ///@brief Pointer definition of FiberBeamColumnElement3D2N
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FiberBeamColumnSection);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    FiberBeamColumnSection(IndexType NewId = 0);

    /// Constructor using an array of nodes
    // FiberBeamColumnSection(IndexType NewId, IntegrationPoint<1>& rIntegrationPoint);
    // FiberBeamColumnSection(IndexType NewId, double position, double weight, PropertiesType::Pointer pProperties);
    FiberBeamColumnSection(IndexType, double, double, int, int, double, double);

    // /// Copy Constructor
    // FiberBeamColumnSection(FiberBeamColumnSection const& rOther);

    // /// Destructor
    // ~FiberBeamColumnSection() ;//override;

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

    std::vector<FiberBeamColumnUniaxialFiber> mFibers;

    ///@}
    ///@name Protected member Variables
    ///@{

    /// FIXME: we need the constitutive law on the fiber level ...
    // ConstitutiveLaw::Pointer mpConstitutiveLaw = nullptr;

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
