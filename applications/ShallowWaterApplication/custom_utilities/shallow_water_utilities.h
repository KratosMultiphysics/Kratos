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

    /// Pointer definition of ShallowWaterUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ShallowWaterUtilities);

    typedef Node<3> NodeType;

    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    /// Destructor.

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void ComputeFreeSurfaceElevation(ModelPart& rModelPart);

    void ComputeHeightFromFreeSurface(ModelPart& rModelPart);

    void ComputeVelocity(ModelPart& rModelPart, bool PerformProjection = false);

    void ComputeSmoothVelocity(ModelPart& rModelPart);

    void ComputeMomentum(ModelPart& rModelPart);

    void ComputeEnergy(ModelPart& rModelPart);

    void ComputeAccelerations(ModelPart& rModelPart);

    double InverseHeight(const double Height, const double Epsilon);

    double WetFraction(double Height, double Epsilon);

    void FlipScalarVariable(Variable<double>& rOriginVariable, Variable<double>& rDestinationVariable, ModelPart& rModelPart);

    void IdentifySolidBoundary(ModelPart& rModelPart, double SeaWaterLevel, Flags SolidBoundaryFlag);

    void IdentifyWetDomain(ModelPart& rModelPart, Flags WetFlag, double Thickness = 0.0);

    void ResetDryDomain(ModelPart& rModelPart, double Thickness = 0.0);

    template<class TContainerType>
    void CopyFlag(Flags OriginFlag, Flags DestinationFlag, TContainerType& rContainer)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rContainer.size()); ++i)
        {
            auto it = rContainer.begin() + i;
            it->Set(DestinationFlag, it->Is(OriginFlag));
        }
    }

    void NormalizeVector(ModelPart& rModelPart, Variable<array_1d<double,3>>& rVariable);

    template<class TVarType>
    void CopyVariableToPreviousTimeStep(ModelPart& rModelPart, const TVarType& rVariable)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rModelPart.NumberOfNodes()); ++i)
        {
            auto const it_node = rModelPart.NodesBegin() + i;
            it_node->FastGetSolutionStepValue(rVariable,1) = it_node->FastGetSolutionStepValue(rVariable);
        }
    }

    void SetMinimumValue(ModelPart& rModelPart, const Variable<double>& rVariable, double MinValue);

    /*
     * @brief Set the z-coordinate of the mesh to zero
     */
    void SetMeshZCoordinateToZero(ModelPart& rModelPart);

    /*
     * @brief Set the z0-coordinate of the mesh to zero
     */
    void SetMeshZ0CoordinateToZero(ModelPart& rModelPart);

    /*
     * @brief Move the z-coordinate of the mesh according to a variable
     */
    void SetMeshZCoordinate(ModelPart& rModelPart, const Variable<double>& rVariable);

    /*
     *@brief Compute the L-2 norm for the given double variable
     */
    template<bool THistorical>
    double ComputeL2Norm(ModelPart& rModelPart, const Variable<double>& rVariable)
    {
        double l2_norm = block_for_each<SumReduction<double>>(rModelPart.Elements(), [&](Element& rElem){
            double partial_l2_norm = 0.0;
            for (auto& r_node : rElem.GetGeometry()) {
                partial_l2_norm += std::pow(GetValue<THistorical>(r_node, rVariable), 2);
            }
            partial_l2_norm *= rElem.GetGeometry().Area();
            partial_l2_norm /= rElem.GetGeometry().size();
            return partial_l2_norm;
        });
        return std::sqrt(l2_norm);
    }

    /*
     *@brief Compute the L-2 norm for the given double variable inside an axis-aligned bounding box
     */
    template<bool THistorical>
    double ComputeL2NormAABB(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        Point& rLow,
        Point& rHigh)
    {
        double l2_norm = block_for_each<SumReduction<double>>(rModelPart.Elements(), [&](Element& rElem){
            double partial_l2_norm = 0.0;
            if (rElem.GetGeometry().HasIntersection(rLow, rHigh)) {
                for (auto& r_node : rElem.GetGeometry()) {
                    partial_l2_norm += std::pow(GetValue<THistorical>(r_node, rVariable), 2);
                }
                partial_l2_norm *= rElem.GetGeometry().Area();
                partial_l2_norm /= rElem.GetGeometry().size();
            }
            return partial_l2_norm;
        });
        return std::sqrt(l2_norm);
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}
private:
    ///@name Operations
    ///@{

    void CalculateMassMatrix(Matrix& rMassMatrix, const GeometryType& rGeometry);

    template<bool THistorical>
    double GetValue(NodeType& rNode, const Variable<double>& rVariable);

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
