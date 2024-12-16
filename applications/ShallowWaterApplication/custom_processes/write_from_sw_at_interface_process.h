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

#pragma once

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
 * @class WriteFromSwAtInterfaceProcess
 * @ingroup ShallowWaterApplication
 * @brief Calculate the minimum distance from all the nodes to a boundary condition in 2D
 * @details The boundary conditions are assumed to be contained in a line
 * @author Miguel Maso Sotomayor/Andrea Montanino
 */
template<std::size_t TDim>
class KRATOS_API(SHALLOW_WATER_APPLICATION) WriteFromSwAtInterfaceProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of WriteFromSwAtInterfaceProcess
    KRATOS_CLASS_POINTER_DEFINITION(WriteFromSwAtInterfaceProcess);

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
    WriteFromSwAtInterfaceProcess() = delete;

    /**
     * @brief Constructor with Model and Parameters
     */
    WriteFromSwAtInterfaceProcess(Model& rModel, Parameters ThisParameters = Parameters());

    /**
     * @brief Destructor
     */
    ~WriteFromSwAtInterfaceProcess() override = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    struct locator_tls {
        Vector N;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results;
        locator_tls(const int max_results = 10000) {
            N.resize(TDim+1);
            results.resize(max_results);
        }
    };

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
        buffer << "WriteFromSwAtInterfaceProcess";
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

    // void GetBoundingVolumeLimits(double& rMin, double& rMax);

    void ReadAndSetValues(
        NodeType& rNode,
        BinBasedFastPointLocator<TDim>& rLocator,
        typename BinBasedFastPointLocator<TDim>::ResultContainerType& rResults);

    // array_1d<double,3> InterpolateVelocity(
    //     const Element::Pointer ElementId,
    //     const Vector& rShapeFunctionValues) const;

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

    // void FindBoundaryNeighbors();

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
    WriteFromSwAtInterfaceProcess& operator=(WriteFromSwAtInterfaceProcess const& rOther) = delete;

    /// Copy constructor.
    WriteFromSwAtInterfaceProcess(WriteFromSwAtInterfaceProcess const& rOther) = delete;

    ///@}

}; // Class WriteFromSwAtInterfaceProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

