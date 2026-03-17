// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#pragma once

#include "geo_mechanics_application_variables.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include <string>

namespace Kratos
{
class ModelPart;
class Node;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyConstantInterpolateLinePressureProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantInterpolateLinePressureProcess);

    ApplyConstantInterpolateLinePressureProcess(ModelPart& model_part, Parameters rParameters);
    ~ApplyConstantInterpolateLinePressureProcess() override = default;
    ApplyConstantInterpolateLinePressureProcess(const ApplyConstantInterpolateLinePressureProcess&) = delete;
    ApplyConstantInterpolateLinePressureProcess& operator=(const ApplyConstantInterpolateLinePressureProcess&) = delete;
    ApplyConstantInterpolateLinePressureProcess(ApplyConstantInterpolateLinePressureProcess&&) = delete;
    ApplyConstantInterpolateLinePressureProcess& operator=(ApplyConstantInterpolateLinePressureProcess&&) = delete;

    /// this function is designed for being called at the first solution step
    /// Note that it depends on the values left by other processes at the nodes that it is interpolating
    void ExecuteInitializeSolutionStep() override;

    /// Turn back information as a string.
    std::string Info() const override;

private:
    /// Member Variables
    ModelPart&         mrModelPart;
    std::string        mVariableName;
    bool               mIsFixed;
    bool               mIsFixedProvided;
    bool               mIsSeepage;
    bool               mIsInitialized = false;
    unsigned int       mGravityDirection;
    unsigned int       mOutOfPlaneDirection;
    unsigned int       mHorizontalDirection;
    std::vector<Node*> mBoundaryNodes;
    double             mPressureTensionCutOff;

    double CalculatePressure(const Node& rNode);

    void CalculateBoundaryPressure(const Node&               rNode,
                                   const std::vector<Node*>& BoundaryNodes,
                                   double&                   pressure,
                                   double&                   coordinate,
                                   bool                      isBottom = false);
    void InterpolateBoundaryPressureWithOneContainer(const Node&               rNode,
                                                     const std::vector<Node*>& rBoundaryNodes,
                                                     double&                   rPressure,
                                                     double&                   rCoordinate) const;

    void InterpolateBoundaryPressure(const Node& rNode,
                                     const Node* pLeftNode,
                                     const Node* pRightNode,
                                     double&     rPressure,
                                     double&     rCoordinate) const;

    void FindTwoClosestNodeOnBoundaryNodes(const Node&               rNode,
                                           const std::vector<Node*>& rBoundaryNodes,
                                           std::vector<Node*>&       rFoundNodes) const;

    Node* FindClosestNodeOnBoundaryNodes(const Node& rNode, const std::vector<Node*>& BoundaryNodes, const bool isBottom);

    void FindTopBoundaryNodes(const Node& rNode, std::vector<Node*>& TopBoundaryNodes) const;

    void FindBottomBoundaryNodes(const Node& rNode, std::vector<Node*>& BottomBoundaryNodes) const;

    void FindLeftBoundaryNodes(const Node&               rNode,
                               const std::vector<Node*>& rBoundaryNodes,
                               std::vector<Node*>&       rLeftBoundaryNodes) const;

    void FindRightBoundaryNodes(const Node&               rNode,
                                const std::vector<Node*>& rBoundaryNodes,
                                std::vector<Node*>&       rRightBoundaryNodes) const;

    int GetMaxNodeID();

    void FindBoundaryNodes();

    void FillListOfBoundaryNodesFast(std::vector<int>& BoundaryNodes);

    bool IsMoreThanOneElementWithThisEdgeFast(const std::vector<int>&              rFaceIDs,
                                              const std::vector<std::vector<int>>& rELementsOfNodes,
                                              const std::vector<int>& rELementsOfNodesSize) const;

}; // Class ApplyConstantInterpolateLinePressureProcess

} // namespace Kratos