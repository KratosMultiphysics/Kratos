//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_SHALLOW_WATER_UTILITIES_H_INCLUDED
#define KRATOS_SHALLOW_WATER_UTILITIES_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

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
 * @ingroup ShallowWaterApplication
 * @class ShallowWaterUtilities
 * @brief This class is a wrapper of useful utilities for shallow water computations
 */
class KRATOS_API(SHALLOW_WATER_APPLICATION) ShallowWaterUtilities
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node NodeType;

    typedef Geometry<NodeType> GeometryType;

    typedef ModelPart::NodesContainerType NodesContainerType;

    ///@}
    ///@name Pointer definition
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ShallowWaterUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void ComputeFreeSurfaceElevation(ModelPart& rModelPart);

    void ComputeHeightFromFreeSurface(ModelPart& rModelPart);

    void ComputeVelocity(ModelPart& rModelPart, bool PerformProjection = false);

    void ComputeSmoothVelocity(ModelPart& rModelPart);

    void ComputeMomentum(ModelPart& rModelPart);

    void ComputeLinearizedMomentum(ModelPart& rModelPart);

    template<bool THistorical>
    void ComputeFroude(ModelPart& rModelPart, const double Epsilon);

    template<bool THistorical>
    void ComputeEnergy(ModelPart& rModelPart);

    void FlipScalarVariable(Variable<double>& rOriginVariable, Variable<double>& rDestinationVariable, ModelPart& rModelPart);

    void IdentifySolidBoundary(ModelPart& rModelPart, double SeaWaterLevel, Flags SolidBoundaryFlag);

    void FlagWetElements(ModelPart& rModelPart, Flags WetFlag, double RelativeDryHeight = -1.0);

    void ExtrapolateElementalFlagToNodes(ModelPart& rModelPart, Flags Flag);

    void NormalizeVector(ModelPart& rModelPart, const Variable<array_1d<double,3>>& rVariable);

    template<class TDataType, class TVarType = Variable<TDataType>>
    void SmoothHistoricalVariable(
        const TVarType& rVariable,
        NodesContainerType& rNodes,
        const double ElapsedTime,
        const double SemiPeriod)
    {
        const double smooth = -std::expm1(-ElapsedTime / SemiPeriod);
        block_for_each(rNodes, [&](NodeType& rNode){
            TDataType& initial = rNode.FastGetSolutionStepValue(rVariable, 1);
            TDataType& current = rNode.FastGetSolutionStepValue(rVariable);
            TDataType increment = current - initial;
            current = initial + smooth * increment;
        });
    }

    template<class TVarType>
    void CopyVariableToPreviousTimeStep(ModelPart& rModelPart, const TVarType& rVariable)
    {
        block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
            rNode.FastGetSolutionStepValue(rVariable, 1) = rNode.FastGetSolutionStepValue(rVariable);
        });
    }

    void SetMinimumValue(ModelPart& rModelPart, const Variable<double>& rVariable, double MinValue);

    /**
     * @brief Set the z-coordinate of the mesh to zero
     */
    void SetMeshZCoordinateToZero(ModelPart& rModelPart);

    /**
     * @brief Set the z0-coordinate of the mesh to zero
     */
    void SetMeshZ0CoordinateToZero(ModelPart& rModelPart);

    /**
     * @brief Move the z-coordinate of the mesh according to a variable
     */
    void SetMeshZCoordinate(ModelPart& rModelPart, const Variable<double>& rVariable);

    /**
     * @brief Move the z-coordinate of the mesh according to a variable
     */
    void OffsetMeshZCoordinate(ModelPart& rModelPart, const double Increment);

    /**
     * @brief Swap the Y and Z coordinates of the nodes
     */
    void SwapYZCoordinates(ModelPart& rModelPart);

    /**
     * @brief Swap the Y and Z coordinates of the nodes
     */
    void SwapY0Z0Coordinates(ModelPart& rModelPart);

    /**
     * @brief Store a double variable as NonHistorical and set the value to no-data if the node is dry
     */
    void StoreNonHistoricalGiDNoDataIfDry(ModelPart& rModelPart, const Variable<double>& rVariable);

    /**
     * @brief Swap the Y and Z components of a vector variable
     */
    void SwapYZComponents(const Variable<array_1d<double,3>>& rVariable, NodesContainerType& rNodes)
    {
        block_for_each(rNodes, [&](NodeType& rNode){
            array_1d<double,3>& r_value = rNode.FastGetSolutionStepValue(rVariable);
            std::swap(r_value[1], r_value[2]);
        });
    }

    /**
     * @brief Swap the Y and Z components of a vector variable
     */
    template<class TContainerType>
    void SwapYZComponentsNonHistorical(const Variable<array_1d<double,3>>& rVariable, TContainerType& rContainer)
    {
        block_for_each(rContainer, [&](typename TContainerType::value_type& rEntity){
            array_1d<double,3>& r_value = rEntity.GetValue(rVariable);
            std::swap(r_value[1], r_value[2]);
        });
    }

    /**
     * @brief Offset the ids of the given container for visualization purpose in GiD
     */
    template<class TContainerType>
    void OffsetIds(TContainerType& rContainer, const double Offset)
    {
        block_for_each(rContainer, [&](typename TContainerType::value_type& rEntity){
            rEntity.SetId(rEntity.Id() + Offset);
        });
    }

    /**
     * @brief Offset the ids of the given container for visualization purpose in GiD
     */
    template<class TContainerType>
    void OffsetIds(TContainerType& rContainer)
    {
        const std::size_t offset = rContainer.size();
        OffsetIds(rContainer, offset);
    }

    /**
     * @brief Compute the L-2 norm for the given double variable
     */
    template<bool THistorical>
    double ComputeL2Norm(ModelPart& rModelPart, const Variable<double>& rVariable);

    /**
     * @brief Compute the L-2 norm for the given double variable inside an axis-aligned bounding box
     */
    template<bool THistorical>
    double ComputeL2NormAABB(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        Point& rLow,
        Point& rHigh);

    /**
     * @brief Compute the horizontal hydrostatic pressures
     */
    template<class TContainerType>
    array_1d<double,3> ComputeHydrostaticForces(
        TContainerType& rContainer,
        const ProcessInfo& rProcessInfo,
        const double RelativeDryHeight = -1.0)
    {
        KRATOS_ERROR_IF_NOT(rProcessInfo.Has(GRAVITY)) << "ShallowWaterUtilities::ComputeHydrostaticForces : GRAVITY is not defined in the ProcessInfo" << std::endl;
        if (rContainer.size() > 0) {
            const auto& r_prop = rContainer.begin()->GetProperties();
            KRATOS_ERROR_IF_NOT(r_prop.Has(DENSITY)) << "ShallowWaterUtilities::ComputeHydrostaticForces : DENSITY is not defined in the Properties" << std::endl;
        }

        array_1d<double,3> forces = ZeroVector(3);
        forces = block_for_each<SumReduction<array_1d<double,3>>>(
            rContainer, [&](typename TContainerType::value_type& rEntity){
                array_1d<double,3> local_force = ZeroVector(3);
                if (RelativeDryHeight >= 0.0) {
                    if (IsWet(rEntity.GetGeometry(), RelativeDryHeight)) {
                        rEntity.Calculate(FORCE, local_force, rProcessInfo);
                    }
                } else {
                    rEntity.Calculate(FORCE, local_force, rProcessInfo);
                }
                return local_force;
            }
        );
        return forces;
    }

    ///@}

private:

    ///@name Operations
    ///@{

    void CalculateMassMatrix(Matrix& rMassMatrix, const GeometryType& rGeometry);

    template<bool THistorical>
    double& GetValue(NodeType& rNode, const Variable<double>& rVariable);

    bool IsWet(const GeometryType& rGeometry, const double RelativeDryHeight);

    bool IsWet(const GeometryType& rGeometry, const double Height, const double RelativeDryHeight);

    bool IsWet(const double Height, const double DryHeight);

    ///@}

}; // Class ShallowWaterUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SHALLOW_WATER_UTILITIES_H_INCLUDED  defined
