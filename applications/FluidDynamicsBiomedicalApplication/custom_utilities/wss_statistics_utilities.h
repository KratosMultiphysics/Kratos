//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla,
//                   Eduardo Soudah
//

#if !defined(KRATOS_WSS_STATISTICS_UTILITIES_H )
#define  KRATOS_WSS_STATISTICS_UTILITIES_H

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// A set of functions to compute the Wall Shear Stress (WSS)
class KRATOS_API(FLUID_DYNAMICS_BIOMEDICAL_APPLICATION) WssStatisticsUtilities {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of WssStatisticsUtilities
    KRATOS_CLASS_POINTER_DEFINITION(WssStatisticsUtilities);

    /// Node definition
    using NodeType = ModelPart::NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Deleted default constructor
    WssStatisticsUtilities() = delete;

    /// Deleted copy constructor.
    WssStatisticsUtilities(WssStatisticsUtilities const& rOther) = delete;

    /// Destructor.
    ~WssStatisticsUtilities() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    WssStatisticsUtilities& operator=(WssStatisticsUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initialize the WSS variables
     * This method initializes all the Wall Shear Stress (WSS) related variables.
     * Note that it is required to be called once before all the other WSS methods.
     * @param rModelPart Model part in which the WSS magnitudes are to be computed
     */
    static void InitializeWSSVariables(ModelPart &rModelPart);

    /**
     * @brief Calculate the Wall Shear Stress (WSS)
     * This method computes the Wall Shear Stress (WSS).
     * @param rModelPart Model part in where the WSS is computed
     * @param rNormalVariable Variable storing the wall nodal normal
     * @param IsNormalHistorical Bool variable indicating wether the historical or non-historical nodal database is used
     */
    static void CalculateWSS(
        ModelPart &rModelPart,
        const Variable<array_1d<double,3>>& rNormalVariable,
        const bool IsNormalHistorical);

    /**
     * @brief Calculate the Oscillatory Shear Index (OSI)
     * This method computes the Oscillatory Shear Index (OSI)
     * @param rModelPart Model part in where the OSI is computed
     */
    static void CalculateOSI(ModelPart &rModelPart);

    ///@}

};  // Class WssStatisticsUtilities

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_WSS_STATISTICS_UTILITIES_H  defined
