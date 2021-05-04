//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Suneth Warnakulasuriya
//

#if !defined(KRATOS_FLUID_AUXILIARY_UTILITIES_H)
#define KRATOS_FLUID_AUXILIARY_UTILITIES_H

// System includes
#include <tuple>
#include <type_traits>

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"

// Application includes


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidAuxiliaryUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    ///@}
    ///@name Static Operations
    ///@{

    /**
     * @brief Calculate the fluid volume
     * For the given model part, this function iterates its elements to calculate the total volume (or area in 2D)
     * @param rModelPart The model part to calculate the volume of
     * @return double Fluid volume
     */
    static double CalculateFluidVolume(ModelPart& rModelPart);


    static double CalculateFluidPositiveVolume(ModelPart& rModelPart);

    static double CalculateFluidNegativeVolume(ModelPart& rModelPart);

    ///@}
private:
    ///@name Private Operations
    ///@{

    // template<bool IsPositiveSide>
    // static double

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_AUXILIARY_UTILITIES_H
