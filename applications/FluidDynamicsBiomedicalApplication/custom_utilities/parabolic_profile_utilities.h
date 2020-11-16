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

#if !defined(KRATOS_PARABOLIC_PROFILE_UTILITIES_H )
#define  KRATOS_PARABOLIC_PROFILE_UTILITIES_H

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
class KRATOS_API(FLUID_DYNAMICS_BIOMEDICAL_APPLICATION) ParabolicProfileUtilities {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ParabolicProfileUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ParabolicProfileUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    // Deleted default constructor
    ParabolicProfileUtilities() = delete;

    /// Deleted copy constructor.
    ParabolicProfileUtilities(ParabolicProfileUtilities const& rOther) = delete;

    /// Destructor.
    ~ParabolicProfileUtilities() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    ParabolicProfileUtilities& operator=(ParabolicProfileUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Impose parabolic profile based on the distance map
     * This method compute the parabolic profile.
     * @param rModelPart Model part in where the parabolic profile is computed
     * @param rSkinModelPart Model part in where the parabolic profile is computed
     */
    static void ComputeMaxDist(ModelPart &rModelPart);
    //, ModelPart &rSkinModelPart);

    /**
     * @brief Impose parabolic profile based on the distance map
     * This method compute the parabolic profile.
     * @param rModelPart Model part in where the parabolic profile is computed
     * @param rSkinModelPart Model part in where the parabolic profile is computed
     */
    static void ImposeParabolic(ModelPart &rModelPart);
    //, ModelPart &rSkinModelPart);


    // /**
    //  * @brief Calculate the Wall Shear Stress (WSS) using Gauss point
    //  * This method computes the Wall Shear Stress (WSS).
    //  * @param rModelPart Model part in where the WSS is computed
    //  */
    // static void CalculateWSSGauss(ModelPart &rModelPart);


    ///@}

};  // Class ParabolicProfileUtilities

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PARABOLIC_PROFILE_UTILITIES_H  defined
