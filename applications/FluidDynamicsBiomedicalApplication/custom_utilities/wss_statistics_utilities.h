//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduardo Soudah,
//                   Ruben Zorrilla,
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
class WssStatisticsUtilities {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of WssStatisticsUtilities
    KRATOS_CLASS_POINTER_DEFINITION(WssStatisticsUtilities);

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
     * @brief Calculate the Wall Shear Stress (WSS)
     * This method computes the Wall Shear Stress (WSS).
     * @param rModelPart Model part in where the WSS is computed
     */
    static void CalculateWSS(ModelPart &rModelPart);

    /**
     * @brief Calculate the Oscillatory Shear Index (OSI)
     * This method computes the Oscillatory Shear Index (OSI)
     * @param rModelPart Model part in where the OSI is computed
     */
    static void CalculateOSI(ModelPart &rModelPart);

    /**
     * @brief Calulate the Temporal Wall Shear Stress (TWSS)
     * This method computes the Temporal Wall Shear Stress (TWSS)
     * @param rModelPart Model part in where the TWSS is computed
     */
    static void CalculateTWSS(ModelPart &rModelPart);

    ///@}

};  // Class WssStatisticsUtilities

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_WSS_STATISTICS_UTILITIES_H  defined
