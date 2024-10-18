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

#ifndef KRATOS_DEPTH_INTEGRATION_PROCESS_H_INCLUDED
#define KRATOS_DEPTH_INTEGRATION_PROCESS_H_INCLUDED


// System includes


// External includes


// Project includes
#include "containers/model.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"
#include "utilities/binbased_fast_point_locator.h"

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
 * @brief Forward declaration of ModelPart
 */
class ModelPart;

/**
 * @class DepthIntegrationProcess
 * @ingroup ShallowWaterApplication
 * @brief Calculate the minimum distance from all the nodes to a boundary condition in 2D
 * @details The boundary conditions are assumed to be contained in a line
 * @author Miguel Maso Sotomayor
 */
template<std::size_t TDim>
class KRATOS_API(SHALLOW_WATER_APPLICATION) DepthIntegrationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DepthIntegrationProcess
    KRATOS_CLASS_POINTER_DEFINITION(DepthIntegrationProcess);

    /// Definition of the node type
    using NodeType = Node;

    /// Definition of the geometry type
    using GeometryType = Geometry<NodeType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @details Removed
     */
    DepthIntegrationProcess() = delete;

    /**
     * @brief Constructor with Model and Parameters
     */
    DepthIntegrationProcess(Model& rModel, Parameters ThisParameters = Parameters());

    /**
     * @brief Destructor
     */
    ~DepthIntegrationProcess() override = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    int Check() override;

    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override {
        std::stringstream buffer;
        buffer << "DepthIntegrationProcess";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrVolumeModelPart;
    ModelPart& mrInterfaceModelPart;
    array_1d<double,3> mDirection;
    bool mStoreHistorical;

    // Option to print a debug file
    bool mPrintVelocityProfile;

    // Option to substitute the boundaries by the neighbor
    bool mExtrapolateBoundaries;
    NodeType::Pointer mpFirstBoundaryNode;
    NodeType::Pointer mpSecondBoundaryNode;
    NodeType::Pointer mpFirstBoundaryNeighbor;
    NodeType::Pointer mpSecondBoundaryNeighbor;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void GetBoundingVolumeLimits(double& rMin, double& rMax);

    void Integrate(
        NodeType& rNode,
        const double Bottom,
        const double Top,
        BinBasedFastPointLocator<TDim>& rLocator,
        typename BinBasedFastPointLocator<TDim>::ResultContainerType& rResults,
        Vector& rShapeFunctionsValues);

    array_1d<double,3> InterpolateVelocity(
        const Element::Pointer ElementId,
        const Vector& rShapeFunctionValues) const;

    template<class TDataType, class TVarType = Variable<TDataType>>
    void SetValue(NodeType& rNode, const TVarType& rVariable, TDataType rValue)
    {
        if (mStoreHistorical)
            rNode.FastGetSolutionStepValue(rVariable) = rValue;
        else
            rNode.GetValue(rVariable) = rValue;
    }

    template<class TDataType, class TVarType = Variable<TDataType>>
    TDataType GetValue(const NodeType& rNode, const TVarType& rVariable)
    {
        if (mStoreHistorical)
            return rNode.FastGetSolutionStepValue(rVariable);
        else
            return rNode.GetValue(rVariable);
    }

    void FindBoundaryNeighbors();

    void CopyValues(const NodeType& rOriginNode, NodeType& rDestinationNode);

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    DepthIntegrationProcess& operator=(DepthIntegrationProcess const& rOther) = delete;

    /// Copy constructor.
    DepthIntegrationProcess(DepthIntegrationProcess const& rOther) = delete;

    ///@}

}; // Class DepthIntegrationProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_DEPTH_INTEGRATION_PROCESS_H_INCLUDED defined
