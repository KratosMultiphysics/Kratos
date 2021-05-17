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

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "modified_shape_functions/modified_shape_functions.h"

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

    using ModifiedShapeFunctionsFactoryType = std::function<ModifiedShapeFunctions::UniquePointer(const GeometryType::Pointer, const Vector&)>;

    ///@}
    ///@name Static Operations
    ///@{

    /**
     * @brief Checks if an element is split
     * From the given vector containing the distance in each node, this method checks if the element is split
     * @param rElementDistancesVector Vector containing the distance values at each node
     * @return true The element is split
     * @return false The element is not split
     */
    static inline bool IsSplit(const Vector& rElementDistancesVector);

    /**
     * @brief Checks if an element is positive
     * From the given vector containing the distance in each node, this method checks if the element is positive
     * @param rElementDistancesVector Vector containing the distance values at each node
     * @return true The element is positive
     * @return false The element is not positive
     */
    static inline bool IsPositive(const Vector &rElementDistancesVector);

    /**
     * @brief Checks if an element is negative
     * From the given vector containing the distance in each node, this method checks if the element is negative
     * @param rElementDistancesVector Vector containing the distance values at each node
     * @return true The element is negative
     * @return false The element is not negative
     */
    static inline bool IsNegative(const Vector &rElementDistancesVector);

    /**
     * @brief Calculate the fluid volume
     * For the given model part, this function iterates its elements to calculate the total volume (or area in 2D)
     * @param rModelPart The model part to calculate the volume of
     * @return double Fluid volume
     */
    static double CalculateFluidVolume(const ModelPart& rModelPart);

    /**
     * @brief Calculate the fluid positive volume
     * For the given model part, this function calculates the total positive fluid volume fraction.
     * It is assumed that a unique element geometry type it is present in the mesh.
     * Only simplex geometries (linear triangle and tetrahedron) are supported.
     * The positive fraction must be described in terms of a continuous level set function stored
     * in the variable DISTANCE of the historical database
     * @param rModelPart The model part to calculate the positive volume
     * @return double Fluid positive volume
     */
    static double CalculateFluidPositiveVolume(const ModelPart& rModelPart);

    /**
     * @brief Calculate the fluid negative volume
     * For the given model part, this function calculates the total negative fluid volume fraction.
     * It is assumed that a unique element geometry type it is present in the mesh.
     * Only simplex geometries (linear triangle and tetrahedron) are supported.
     * The negative fraction must be described in terms of a continuous level set function stored
     * in the variable DISTANCE of the historical database
     * @param rModelPart The model part to calculate the negative volume
     * @return double Fluid negative volume
     */
    static double CalculateFluidNegativeVolume(const ModelPart& rModelPart);

    /**
     * @brief Calculate the flow rate through the given model part conditions
     * This method calculates the flow rate throught the given model part conditions
     * It is assumed that only linear elements are employed for the discretization of velocity field
     * @param rModelPart The model part to calculate the flow rate
     * @return double Flow rate
     */
    static double CalculateFlowRate(const ModelPart& rModelPart);

    /**
     * @brief Get the Standard Modified Shape Functions Factory object
     * For the given geometry, this methods returns a factory function to create the corresponding
     * standard modified shape functions. The factory returns a unique pointer and expects a
     * pointer to the geometry as well as the nodal distances vector.
     * @param rGeometry Geometry reference to check the geometry type
     * @return ModifiedShapeFunctionsFactoryType Factory function to create the standard modified shape functions
     */
    static ModifiedShapeFunctionsFactoryType GetStandardModifiedShapeFunctionsFactory(const GeometryType& rGeometry);

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_AUXILIARY_UTILITIES_H
