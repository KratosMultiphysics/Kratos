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

#if !defined(KRATOS_FIBER_BEAM_COLUMN_Uniaxial_Fiber_H_INCLUDED )
#define  KRATOS_FIBER_BEAM_COLUMN_Uniaxial_Fiber_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/element.h"

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
    ///@brief Pointer definition of FiberBeamColumnElement3D2N
    ///@{

    // KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(FiberBeamColumnSection); //FIXME: This gives an error

    ///@}
    ///@name Life Cycle
    ///@{

    // /// Constructor
    // FiberBeamColumnSection(IndexType NewId = 0);

    /// Constructor using an array of nodes
    // FiberBeamColumnSection(IndexType NewId, IntegrationPoint<1>& rIntegrationPoint);
    FiberBeamColumnUniaxialFiber(IndexType NewId, double Y, double Z, double Area, ConstitutiveLaw Material);

}; // class FiberBeamColumnUniaxialFiber

} // namespace Kratos

#endif