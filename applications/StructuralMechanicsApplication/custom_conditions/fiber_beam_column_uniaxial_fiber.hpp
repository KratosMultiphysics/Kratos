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

#if !defined(KRATOS_FIBER_BEAM_COLUMN_UNIAXIAL_FIBER_H_INCLUDED )
#define  KRATOS_FIBER_BEAM_COLUMN_UNIAXIAL_FIBER_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/element.h"

#include "custom_constitutive/uniaxial_fiber_beam_column_material_law.hpp"

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
 * @class FiberBeamColumnUniaxialFiber
 *
 * @brief A 3D unaxial fiber for the beam-column element for reinforced concrete modeling
 *
 * @author Mahmoud Zidan
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) FiberBeamColumnUniaxialFiber
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
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FiberBeamColumnUniaxialFiber);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor
    FiberBeamColumnUniaxialFiber(IndexType NewId = 0);

    /// Constructor using an array of nodes
    FiberBeamColumnUniaxialFiber(
        IndexType NewId, double Y, double Z, double Area, UniaxialFiberBeamColumnMaterialLaw::Pointer pMaterial);

    Matrix CreateGlobalFiberStiffnessMatrix();

    void StateDetermination(const Vector& rSectionDeformationIncrements);
    Vector CreateGlobalFiberInternalForces();
    void FinalizeSolutionStep();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize();

    ///@}
    ///@name Access
    ///@{

    double const& GetStress() const {return mStress;}
    double const& GetStrain() const {return mStrain;}

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

    IndexType mId;
    Vector mTransformationVector = ZeroVector(3);
    double mArea = 0;
    UniaxialFiberBeamColumnMaterialLaw::Pointer mpMaterial = nullptr;
    PropertiesType::Pointer mpProperties = nullptr;
    double mStrain = 0.0;
    double mStress = 0.0;

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
    void save(Serializer& rSerializer) const;
    void load(Serializer& rSerializer);

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

};  // class FiberBeamColumnUniaxialFiber

/// output stream
inline std::ostream & operator <<(std::ostream& rOStream, const FiberBeamColumnUniaxialFiber& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos

#endif