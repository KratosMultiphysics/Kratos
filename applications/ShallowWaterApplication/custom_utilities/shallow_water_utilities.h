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

    typedef Node<3> NodeType;

    typedef Geometry<NodeType> GeometryType;

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

    template<bool THistorical>
    void ComputeFroude(ModelPart& rModelPart, const double Epsilon);

    template<bool THistorical>
    void ComputeEnergy(ModelPart& rModelPart);

    double InverseHeight(const double Height, const double Epsilon);

    double WetFraction(double Height, double Epsilon);

    void FlipScalarVariable(Variable<double>& rOriginVariable, Variable<double>& rDestinationVariable, ModelPart& rModelPart);

    void IdentifySolidBoundary(ModelPart& rModelPart, double SeaWaterLevel, Flags SolidBoundaryFlag);

    void IdentifyWetDomain(ModelPart& rModelPart, Flags WetFlag, double RelativeDryHeight = 0.1);

    template<class TContainerType>
    void CopyFlag(Flags OriginFlag, Flags DestinationFlag, TContainerType& rContainer)
    {
        block_for_each(rContainer, [&](typename TContainerType::value_type& rEntity){
            rEntity.Set(DestinationFlag, rEntity.Is(OriginFlag));
        });
    }

    void NormalizeVector(ModelPart& rModelPart, Variable<array_1d<double,3>>& rVariable);

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
     * @brief Store a double variable as NonHistorical and set the value to no-data if the node is dry
     */
    void StoreNonHistoricalGiDNoDataIfDry(ModelPart& rModelPart, const Variable<double>& rVariable);

    /**
     * @brief Offset the ids of the given container for visualization purpose in GiD
     */
    template<class TContainerType>
    void OffsetIds(TContainerType& rContainer)
    {
        const std::size_t offset = rContainer.size();
        block_for_each(rContainer, [&](typename TContainerType::value_type& rEntity){
            rEntity.SetId(rEntity.Id() + offset);
        });
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
        KRATOS_ERROR_IF_NOT(rProcessInfo.Has(DENSITY)) << "ShallowWaterUtilities::ComputeHydrostaticForces : DENSITY is not defined in the ProcessInfo" << std::endl;
        const double gravity = rProcessInfo.GetValue(GRAVITY_Z);
        const double density = rProcessInfo.GetValue(DENSITY);

        array_1d<double,3> forces = ZeroVector(3);
        forces = block_for_each<SumReduction<array_1d<double,3>>>(
            rContainer, [&](typename TContainerType::value_type& rEntity){
                const auto& r_geom = rEntity.GetGeometry();
                const double area = r_geom.Area();
                const array_1d<double,3> normal = r_geom.UnitNormal(r_geom[0]); // At the first Point
                double height = 0.0;
                for (auto& r_node : r_geom) {
                    height += r_node.FastGetSolutionStepValue(HEIGHT);
                }
                height /= r_geom.size();
                array_1d<double,3> local_force;
                if (RelativeDryHeight >= 0.0) {
                    if (IsWet(r_geom, height, RelativeDryHeight)) {
                        local_force = EvaluateHydrostaticForce<TContainerType>(density, gravity, height, area, normal);
                    } else {
                        local_force = ZeroVector(3);
                    }
                } else {
                    local_force = EvaluateHydrostaticForce<TContainerType>(density, gravity, height, area, normal);
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

    template<class TContainerType>
    void IdentifyWetEntities(TContainerType& rContainer, Flags WetFlag, double RelativeDryHeight)
    {
        block_for_each(rContainer, [&](typename TContainerType::value_type& rEntity){
            const auto& r_geom = rEntity.GetGeometry();
            const bool is_wet = IsWet(r_geom, RelativeDryHeight);
            rEntity.Set(WetFlag, is_wet);
            for (auto& r_node : r_geom)
            {
                if (is_wet)
                {
                    if (r_node.IsNot(WetFlag))
                    {
                        r_node.SetLock();
                        r_node.Set(WetFlag);
                        r_node.UnSetLock();
                    }
                }
            }
        });
    }

    bool IsWet(const GeometryType& rGeometry, const double RelativeDryHeight);

    bool IsWet(const GeometryType& rGeometry, const double Height, const double RelativeDryHeight);

    bool IsWet(const double Height, const double DryHeight);

    template<class TContainerType>
    array_1d<double,3> EvaluateHydrostaticForce(
        const double Density,
        const double Gravity,
        const double Height,
        const double Area,
        const array_1d<double,3>& rNormal);

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
