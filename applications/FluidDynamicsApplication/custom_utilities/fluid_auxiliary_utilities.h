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

#pragma once

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

    using NodeType = Node;

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
    static inline bool IsSplit(const Vector& rElementDistancesVector)
    {
        std::size_t n_pos(0);
        std::size_t n_neg(0);
        const std::size_t pts_number = rElementDistancesVector.size();
        for (std::size_t i_node = 0; i_node < pts_number; ++i_node){
            if (rElementDistancesVector[i_node] > 0.0)
                n_pos++;
            else
                n_neg++;
        }
        return (n_pos > 0 && n_neg > 0) ? true : false;
    }

    /**
     * @brief Checks if an element is positive
     * From the given vector containing the distance in each node, this method checks if the element is positive
     * @param rElementDistancesVector Vector containing the distance values at each node
     * @return true The element is positive
     * @return false The element is not positive
     */
    static inline bool IsPositive(const Vector &rElementDistancesVector)
    {
        std::size_t n_pos (0);
        const std::size_t pts_number = rElementDistancesVector.size();
        for (std::size_t i_node = 0; i_node < pts_number; ++i_node){
            if (rElementDistancesVector[i_node] > 0.0)
                n_pos++;
        }
        return (n_pos == pts_number) ? true : false;
    }

    /**
     * @brief Checks if an element is negative
     * From the given vector containing the distance in each node, this method checks if the element is negative
     * @param rElementDistancesVector Vector containing the distance values at each node
     * @return true The element is negative
     * @return false The element is not negative
     */
    static inline bool IsNegative(const Vector &rElementDistancesVector)
    {
        std::size_t n_neg (0);
        const std::size_t pts_number = rElementDistancesVector.size();
        for (std::size_t i_node = 0; i_node < pts_number; ++i_node){
            if (rElementDistancesVector[i_node] < 0.0)
                n_neg++;
        }
        return n_neg == pts_number;
    }

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
     * @brief Calculate the fluid negative volume in cut elements
     * For the given model part, this function calculates the total negative fluid volume fraction in the cut elements.
     * It is assumed that a unique element geometry type it is present in the mesh.
     * Only simplex geometries (linear triangle and tetrahedron) are supported.
     * The negative fraction must be described in terms of a continuous level set function stored
     * in the variable DISTANCE of the historical database
     * @param rModelPart The model part to calculate the negative volume
     * @return double Fluid negative volume over the cut elements
     */
    static double CalculateFluidCutElementNegativeVolume(const ModelPart &rModelPart);

    /**
     * @brief Calculate the fluid positive volume in cut elements
     * For the given model part, this function calculates the total positive fluid volume fraction in the cut elements.
     * It is assumed that a unique element geometry type it is present in the mesh.
     * Only simplex geometries (linear triangle and tetrahedron) are supported.
     * The positive fraction must be described in terms of a continuous level set function stored
     * in the variable DISTANCE of the historical database
     * @param rModelPart The model part to calculate the positive volume
     * @return double Fluid positive volume over the cut elements
     */
    static double CalculateFluidCutElementsPositiveVolume(const ModelPart &rModelPart);

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
     * @brief Calculate the flow rate through the given model part conditions (positive subdomain)
     * This method calculates the flow rate throught the positive part of given model part conditions
     * It is assumed that only linear elements are employed for the discretization of velocity field
     * @param rModelPart The model part to calculate the flow rate
     * @return double Flow rate
     */
    static double CalculateFlowRatePositiveSkin(const ModelPart& rModelPart);

    /**
     * @brief Calculate the flow rate through the given model part conditions (negative subdomain)
     * This method calculates the flow rate throught the negative part of given model part conditions
     * It is assumed that only linear elements are employed for the discretization of velocity field
     * @param rModelPart The model part to calculate the flow rate
     * @return double Flow rate
     */
    static double CalculateFlowRateNegativeSkin(const ModelPart& rModelPart);

    /**
     * @brief Calculate the flow rate through the given model part conditions (positive subdomain)
     * This method calculates the flow rate throught the positive part of given model part conditions
     * It is assumed that only linear elements are employed for the discretization of velocity field
     * @param rModelPart The model part to calculate the flow rate
     * @param rSkinFlag Flag that marks the conditions to be included in the calculation
     * @return double Flow rate
     */
    static double CalculateFlowRatePositiveSkin(
        const ModelPart& rModelPart,
        const Flags& rSkinFlag);

    /**
     * @brief Calculate the flow rate through the given model part conditions (negative subdomain)
     * This method calculates the flow rate throught the negative part of given model part conditions
     * It is assumed that only linear elements are employed for the discretization of velocity field
     * @param rModelPart The model part to calculate the flow rate
     * @param rSkinFlag Flag that marks the conditions to be included in the calculation
     * @return double Flow rate
     */
    static double CalculateFlowRateNegativeSkin(
        const ModelPart& rModelPart,
        const Flags& rSkinFlag);

    /**
     * @brief Get the Standard Modified Shape Functions Factory object
     * For the given geometry, this methods returns a factory function to create the corresponding
     * standard modified shape functions. The factory returns a unique pointer and expects a
     * pointer to the geometry as well as the nodal distances vector.
     * @param rGeometry Geometry reference to check the geometry type
     * @return ModifiedShapeFunctionsFactoryType Factory function to create the standard modified shape functions
     */
    static ModifiedShapeFunctionsFactoryType GetStandardModifiedShapeFunctionsFactory(const GeometryType& rGeometry);

    /**
     * @brief This function maps the velocity of the nodes "on the skin" by the use of RBF interpolation
     * to the variable EMBEDDED_VELOCITY of the nodes which correspond to the cut cells in the volume.
     * It is designed for use together with the "embedded solver"
     * @param rVolumeModelPart is the destination domain on which EMBEDDED_VELOCITY will be calculated
     * @param rSkinModelPart is the skin of the object from which the velocity will be taken
     * @param SearchRadius is the radius which will be used in searching the neighbours. It needs to be sufficiently large otherwise the method will fail in the calculation of the RBF basis
     * NOTE: historical variable VELOCITY is assumed to be present on the nodes of the rSkinModelPArt
     * NOTE: non historical variable EMBEDDED_VELOCITY is assumed to be present in rVolumeModelPart prior to calling the function
     */
    static void MapVelocityFromSkinToVolumeRBF(
        ModelPart& rVolumeModelPart,
        ModelPart& rSkinModelPart,
        const double SearchRadius);

    /**
     * @brief Find the maximum edge length
     * This function finds and returns the maximum edge length in the given model part from the nodal neighbours
     * @param rModelPart The model part to find the maximum edge length
     * @param CalculateNodalNeighbours Indicates if the nodal neighbours calculation is required (true by default)
     * @return double The maximum edge length
     */
    static double FindMaximumEdgeLength(
        ModelPart& rModelPart,
        const bool CalculateNodalNeighbours = true);

    /**
     * @brief Postprocess the midpoint nodes pressure in P2P1 elements
     * This function takes the edges' midpoint nodes in P2P1 elements and postprocess the pressure, which
     * is assumed to be stored in PRESSURE historical variable, from the edges' endpoint values.
     * Note that the nodal flag VISITED is used to mark the nodes which pressure has been already set.
     * @param rModelPart The model part to which the pressure is to be postprocessed
     */
    static void PostprocessP2P1ContinuousPressure(ModelPart& rModelPart);

    ///@}
private:

    /**
     * @brief Calculate one condition flow rate
     * This method calculates the flow rate through the given condition geometry
     * @param rGeometry Condition geometry
     * @return double Flow rate
     */
    static double CalculateConditionFlowRate(const GeometryType& rGeometry);

    /**
     * @brief Auxiliary function to bypass the positive and negative standard methods
     * This auxiliary function shouldn't be called from outside and serves to
     * avoid the reimplementation for positive and negative subdomains
     * @tparam IsPositiveSubdomain Positive for positive subdomain and viceversa
     * @param rModelPart The model part to calculate the flow rate
     * @param rSkinFlag Flag that marks the conditions to be included in the calculation
     * @return double Flow rate
     */
    template<bool IsPositiveSubdomain, bool CheckConditionFlag>
    static double CalculateFlowRateAuxiliary(
        const ModelPart& rModelPart,
        const Flags& rSkinFlag = Flags());

    /**
     * @brief
     *
     * @tparam CheckConditionFlag
     * @param rCondition
     * @param rSkinFlag
     * @return true
     * @return false
     */
    template<bool CheckConditionFlag>
    static bool CheckConditionFlagAuxiliary(
        const Condition& rCondition,
        const Flags& rSkinFlag);

    /**
     * @brief Auxiliary function to bypass the positive and negative subdomain check
     * This auxiliary function shouldn't be called from outside and serves to
     * standarize the implementation of the methods that check for an either positive or negative subdomain
     * @tparam IsPositiveSubdomain Positive for positive subdomain and viceversa
     * @param rElementDistancesVector Vector containing the distance values at each node
     * @return true If agrees the template argument
     * @return false If not agrees the template argument
     */
    template<bool IsPositiveSubdomain>
    static inline bool CheckNonSplitConditionSubdomain(const Vector &rElementDistancesVector);

    /**
     * @brief Auxiliary function to bypass the calculation of modified shape functions
     * This auxiliary function shouldn't be called from outside and serves to
     * avoid checking for positive and negative methods of the modified shape functions
     * @tparam IsPositiveSubdomain Positive for positive subdomain and viceversa
     * @param rpModShapeFunc Pointer to the modified shape functions utility of the condition parent
     * @param FaceId Parent geometry face id corresponding the condition of interest
     * @param rShapeFunctions Matrix container to store the computed shape functions
     * @param rNormals Vector containing the outwards normals to the face of interest
     * @param rWeights Vector containing the face of interest weights
     */
    template<bool IsPositiveSubdomain>
    static void CalculateSplitConditionGeometryData(
        const ModifiedShapeFunctions::UniquePointer& rpModShapeFunc,
        const std::size_t FaceId,
        Matrix& rShapeFunctions,
        std::vector<array_1d<double,3>>& rNormals,
        Vector& rWeights);

    /**
     * @brief Auxilary function to postprocess one P2P1 edge pressure
     * This function postprocesses the PRESSURE in a P2P1 element edge midpoint.
     * Once the pressure value is set, the node is marked as VISITED.
     * @param rGeometry Reference to current element geometry
     * @param PostNodeLocalId Local id of the node to which the pressure is to be set
     * @param EdgeNodeLocalIdI Local id of the i-node of the edge to which previous node belongs
     * @param EdgeNodeLocalIdJ Local id of the j-node of the edge to which previous node belongs
     */
    static void PostprocessP2P1NodePressure(
        GeometryType& rGeometry,
        const std::size_t PostNodeLocalId,
        const std::size_t EdgeNodeLocalIdI,
        const std::size_t EdgeNodeLocalIdJ)
    {
        if (rGeometry[PostNodeLocalId].IsNot(VISITED)) {
            rGeometry[PostNodeLocalId].SetLock();
            const double p_i = rGeometry[EdgeNodeLocalIdI].FastGetSolutionStepValue(PRESSURE);
            const double p_j = rGeometry[EdgeNodeLocalIdJ].FastGetSolutionStepValue(PRESSURE);
            rGeometry[PostNodeLocalId].FastGetSolutionStepValue(PRESSURE) = 0.5 * (p_i + p_j);
            rGeometry[PostNodeLocalId].Set(VISITED, true);
            rGeometry[PostNodeLocalId].UnSetLock();
        }
    }

};

///@}

} // namespace Kratos
