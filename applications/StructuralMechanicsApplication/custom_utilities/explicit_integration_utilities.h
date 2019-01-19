// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B Sautter
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_EXPLICIT_LAW_UTILITIES)
#define KRATOS_EXPLICIT_LAW_UTILITIES

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

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
/**
 * @class ExplicitIntegrationUtilities
 * @ingroup StructuralMechanicsApplication
 * @brief This class includes several utilities necessaries for the computation of the explicit integration
 * @details The methods are static, so it can be called without constructing the class
 * @author Klaus B Sautter
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ExplicitIntegrationUtilities
{
  public:
    ///@name Type definitions
    ///@{

    /// The size type definition
    typedef std::size_t SizeType;

    /// The index type definition
    typedef std::size_t IndexType;

    /// Node type definition
    typedef Node<3> NodeType;

    /// Geometry definitions
    typedef Geometry<NodeType> GeometryType;

    /// The arrays of elements and nodes
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::NodesContainerType NodesArrayType;

    /// The zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method computes the necessry delta time to avoid numerical instabilities
     * @param rModelPart The model of the problem to solve
     * @param PredictionLevel The prediction level
     * @param Maximum The maximum delta time to be considered
     * @param SafetyFactor The factor to not consider exactly the theoretical value
     */
    static void CalculateDeltaTime(
        ModelPart& rModelPart,
        const double PredictionLevel = 2.0,
        const double Maximum = 1.0e-3,
        const double SafetyFactor = 0.5
        );

}; // class ExplicitIntegrationUtilities
} // namespace Kratos
#endif /* KRATOS_EXPLICIT_LAW_UTILITIES defined */
