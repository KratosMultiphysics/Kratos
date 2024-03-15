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

#include <algorithm>
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyConstantInterpolateLinePressureProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantInterpolateLinePressureProcess);

    ApplyConstantInterpolateLinePressureProcess(ModelPart& model_part,
                                                Parameters rParameters
                                                ) : Process(Flags()) , mrModelPart(model_part)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed": false,
                "is_seepage": false,
                "gravity_direction": 1,
                "out_of_plane_direction": 2,
                "pressure_tension_cut_off" : 0.0,
                "table" : 1
            }  )" );

        // Some values need to be mandatory prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["variable_name"];
        rParameters["model_part_name"];

        mIsFixedProvided = rParameters.Has("is_fixed");

        // Now validate against defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();

        FindBoundaryNodes();

        mIsFixed             = rParameters["is_fixed"].GetBool();
        mIsSeepage           = rParameters["is_seepage"].GetBool();
        mGravityDirection    = rParameters["gravity_direction"].GetInt();
        mOutOfPlaneDirection = rParameters["out_of_plane_direction"].GetInt();
        if (mGravityDirection == mOutOfPlaneDirection)
            KRATOS_ERROR << "Gravity direction cannot be the same as Out-of-Plane directions"
                         << rParameters
                         << std::endl;

        mHorizontalDirection = 0;
        for (unsigned int i=0; i<N_DIM_3D; ++i)
           if (i!=mGravityDirection && i!=mOutOfPlaneDirection) mHorizontalDirection = i;

        mPressureTensionCutOff = rParameters["pressure_tension_cut_off"].GetDouble();

        KRATOS_CATCH("")
    }

    ~ApplyConstantInterpolateLinePressureProcess() override = default;
    ApplyConstantInterpolateLinePressureProcess(const ApplyConstantInterpolateLinePressureProcess&) = delete;
    ApplyConstantInterpolateLinePressureProcess& operator=(const ApplyConstantInterpolateLinePressureProcess&) = delete;
    ApplyConstantInterpolateLinePressureProcess(ApplyConstantInterpolateLinePressureProcess&&) = delete;
    ApplyConstantInterpolateLinePressureProcess& operator=(ApplyConstantInterpolateLinePressureProcess&&) = delete;

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY

        const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);

        block_for_each(mrModelPart.Nodes(), [&var, this](Node& rNode) {
            const double pressure = CalculatePressure(rNode);
            if (mIsSeepage) {
                if (pressure < PORE_PRESSURE_SIGN_FACTOR * mPressureTensionCutOff) {
                    rNode.FastGetSolutionStepValue(var) = pressure;
                    if (mIsFixed) rNode.Fix(var);
                } else {
                    rNode.Free(var);
                }
            } else {
                rNode.FastGetSolutionStepValue(var) = std::min(pressure,PORE_PRESSURE_SIGN_FACTOR * mPressureTensionCutOff);
                if (mIsFixed) rNode.Fix(var);
                else if (mIsFixedProvided) rNode.Free(var);
            }
        });

        KRATOS_CATCH("")
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyConstantInterpolateLinePressureProcess";
    }

private:
    /// Member Variables
    ModelPart& mrModelPart;
    std::string mVariableName;
    bool mIsFixed;
    bool mIsFixedProvided;
    bool mIsSeepage;
    unsigned int mGravityDirection;
    unsigned int mOutOfPlaneDirection;
    unsigned int mHorizontalDirection;
    std::vector< Node * > mBoundaryNodes;
    double mPressureTensionCutOff;

    double CalculatePressure(const Node &rNode)
    {
        // find top boundary
        std::vector< Node* > TopBoundaryNodes;
        FindTopBoundaryNodes(rNode, TopBoundaryNodes);
        double PressureTop;
        double CoordinateTop;
        CalculateBoundaryPressure(rNode, TopBoundaryNodes, PressureTop, CoordinateTop);

        // find bottom boundary
        std::vector< Node* > BottomBoundaryNodes;
        FindBottomBoundaryNodes(rNode, BottomBoundaryNodes);
        double PressureBottom;
        double CoordinateBottom;
        CalculateBoundaryPressure(rNode, BottomBoundaryNodes, PressureBottom, CoordinateBottom, true);

        // calculate pressure
        if (std::abs(CoordinateTop - CoordinateBottom) > TINY) {
            const double slopeP   = (PressureTop - PressureBottom) / (CoordinateTop - CoordinateBottom);
            return slopeP * (rNode.Coordinates()[mGravityDirection] - CoordinateBottom ) + PressureBottom;
        } else {
            return PressureBottom;
        }
    }

    void CalculateBoundaryPressure( const Node &rNode,
                                    const std::vector< Node*> &BoundaryNodes,
                                    double &pressure,
                                    double &coordinate,
                                    bool isBottom=false )
    {
        // find top boundary
        std::vector< Node*> LeftBoundaryNodes;
        FindLeftBoundaryNodes(rNode, BoundaryNodes, LeftBoundaryNodes);

        std::vector< Node*> RightBoundaryNodes;
        FindRightBoundaryNodes(rNode, BoundaryNodes, RightBoundaryNodes);

        if (!LeftBoundaryNodes.empty() && !RightBoundaryNodes.empty()) {
            const Node *LeftNode  = FindClosestNodeOnBoundaryNodes(rNode, LeftBoundaryNodes,  isBottom);
            const Node *RightNode = FindClosestNodeOnBoundaryNodes(rNode, RightBoundaryNodes, isBottom);

            InterpolateBoundaryPressure(rNode, LeftNode, RightNode, pressure, coordinate);
            return;

        } else if (!LeftBoundaryNodes.empty()) {
            InterpolateBoundaryPressureWithOneContainer(rNode, LeftBoundaryNodes, pressure, coordinate);
            return;

        } else if (!RightBoundaryNodes.empty()) {
            InterpolateBoundaryPressureWithOneContainer(rNode, RightBoundaryNodes, pressure, coordinate);
            return;

        } else {
            KRATOS_ERROR << "There is not enough points around interpolation, node Id" << rNode.Id() << std::endl;
        }

    }

    void InterpolateBoundaryPressureWithOneContainer(const Node&               rNode,
                                                     const std::vector<Node*>& rBoundaryNodes,
                                                     double&                   rPressure,
                                                     double&                   rCoordinate) const
    {
        std::vector< Node*> FoundNodes;
        FindTwoClosestNodeOnBoundaryNodes(rNode, rBoundaryNodes, FoundNodes);

        const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);

        const double &pressureLeft = FoundNodes[0]->FastGetSolutionStepValue(var);
        Vector3 CoordinatesLeft;
        noalias(CoordinatesLeft) = FoundNodes[0]->Coordinates();

        const double &pressureRight = FoundNodes[1]->FastGetSolutionStepValue(var);
        Vector3 CoordinatesRight;
        noalias(CoordinatesRight) = FoundNodes[1]->Coordinates();

        // calculate pressure
        if (std::abs(CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) > TINY) {
            const auto slope_p =
                (pressureRight - pressureLeft) /
                (CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]);
            rPressure = slope_p * (rNode.Coordinates()[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) +
                        pressureLeft;

            const auto slope_y =
                (CoordinatesRight[mGravityDirection] - CoordinatesLeft[mGravityDirection]) /
                (CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]);
            rCoordinate = slope_y * (rNode.Coordinates()[mHorizontalDirection] -
                                     CoordinatesLeft[mHorizontalDirection]) +
                          CoordinatesLeft[mGravityDirection];
        } else {
            rPressure   = pressureLeft;
            rCoordinate = CoordinatesLeft[mGravityDirection];
        }
    }

    void InterpolateBoundaryPressure(const Node& rNode,
                                     const Node* pLeftNode,
                                     const Node* pRightNode,
                                     double&     rPressure,
                                     double&     rCoordinate) const
    {
        const auto& r_variable = KratosComponents<Variable<double>>::Get(mVariableName);

        const auto pressure_left = pLeftNode->FastGetSolutionStepValue(r_variable);
        Vector3    coordinates_left;
        noalias(coordinates_left) = pLeftNode->Coordinates();

        const auto pressure_right = pRightNode->FastGetSolutionStepValue(r_variable);
        Vector3    coordinates_right;
        noalias(coordinates_right) = pRightNode->Coordinates();

        // calculate pressure
        if (std::abs(coordinates_right[mHorizontalDirection] - coordinates_left[mHorizontalDirection]) > TINY) {
            const auto slope_p =
                (pressure_right - pressure_left) /
                (coordinates_right[mHorizontalDirection] - coordinates_left[mHorizontalDirection]);
            rPressure = slope_p * (rNode.Coordinates()[mHorizontalDirection] -
                                   coordinates_left[mHorizontalDirection]) +
                        pressure_left;

            const auto slope_y =
                (coordinates_right[mGravityDirection] - coordinates_left[mGravityDirection]) /
                (coordinates_right[mHorizontalDirection] - coordinates_left[mHorizontalDirection]);
            rCoordinate = slope_y * (rNode.Coordinates()[mHorizontalDirection] -
                                     coordinates_left[mHorizontalDirection]) +
                          coordinates_left[mGravityDirection];
        } else {
            rPressure   = pressure_left;
            rCoordinate = coordinates_left[mGravityDirection];
        }
    }

    void FindTwoClosestNodeOnBoundaryNodes(const Node&               rNode,
                                           const std::vector<Node*>& rBoundaryNodes,
                                           std::vector<Node*>&       rFoundNodes) const
    {
        const double HorizontalCoordinate = rNode.Coordinates()[mHorizontalDirection];
        rFoundNodes.resize(2);

        unsigned int nFound = 0;
        double horizontalDistanceClosest_1 = LARGE;
        for (auto boundary_node : rBoundaryNodes) {
            Vector3 CoordinatesBoundary;
            noalias(CoordinatesBoundary) = boundary_node->Coordinates();

            if (std::abs(CoordinatesBoundary[mHorizontalDirection] - HorizontalCoordinate) <= horizontalDistanceClosest_1) {
                horizontalDistanceClosest_1 = std::abs(CoordinatesBoundary[mHorizontalDirection] - HorizontalCoordinate);
                rFoundNodes[0] = boundary_node;
                nFound++;
            }
        }

        double horizontalDistanceClosest_2 = LARGE;
        for (auto boundary_node : rBoundaryNodes) {
            Vector3 CoordinatesBoundary;
            noalias(CoordinatesBoundary) = boundary_node->Coordinates();

            if (std::abs(CoordinatesBoundary[mHorizontalDirection] - HorizontalCoordinate) <= horizontalDistanceClosest_2 &&
                std::abs(CoordinatesBoundary[mHorizontalDirection] - HorizontalCoordinate) > horizontalDistanceClosest_1) {
                horizontalDistanceClosest_2 = std::abs(CoordinatesBoundary[mHorizontalDirection] - HorizontalCoordinate);
                rFoundNodes[1] = boundary_node;
                nFound++;
            }
        }

        KRATOS_ERROR_IF(nFound < 2) << "Not enough points for interpolation: Coordinates"<< rNode.Coordinates() << std::endl;
    }


    Node* FindClosestNodeOnBoundaryNodes(const Node &rNode,
                                         const std::vector< Node* > &BoundaryNodes,
                                         const bool isBottom)
    {
        const double HorizontalCoordinate = rNode.Coordinates()[mHorizontalDirection];
        std::vector< Node*> FoundNodes;

        double horizontalDistance = LARGE;
        for (auto boundary_node : BoundaryNodes) {
            if (std::abs(boundary_node->Coordinates()[mHorizontalDirection] - HorizontalCoordinate) <= horizontalDistance) {
                horizontalDistance = std::abs(boundary_node->Coordinates()[mHorizontalDirection] - HorizontalCoordinate);
                FoundNodes.push_back(boundary_node);
            }
        }

        Node* p_result = nullptr;
        if (isBottom) {
            double height = LARGE;
            for (auto found_node : FoundNodes) {
                if (found_node->Coordinates()[mGravityDirection] < height) {
                    p_result = found_node;
                    height   = found_node->Coordinates()[mGravityDirection];
                }
            }
        } else {
            double height = -LARGE;
            for (auto found_node : FoundNodes) {
                if (found_node->Coordinates()[mGravityDirection] > height) {
                    p_result = found_node;
                    height   = found_node->Coordinates()[mGravityDirection];
                }
            }
        }

        return p_result;
    }

    void FindTopBoundaryNodes(const Node &rNode,
                              std::vector< Node* > &TopBoundaryNodes)
    {
        for (unsigned int i = 0; i < mBoundaryNodes.size(); ++i) {
            if (mBoundaryNodes[i]->Coordinates()[mGravityDirection] >= rNode.Coordinates()[mGravityDirection]) {
                // node is on top boundary
                TopBoundaryNodes.push_back(mBoundaryNodes[i]);
            }
        }
    }

    void FindBottomBoundaryNodes(const Node &rNode,
                                 std::vector< Node*> &BottomBoundaryNodes)
    {
        for (unsigned int i = 0; i < mBoundaryNodes.size(); ++i) {
            if (mBoundaryNodes[i]->Coordinates()[mGravityDirection] <= rNode.Coordinates()[mGravityDirection]) {
                // node is on top boundary
                BottomBoundaryNodes.push_back(mBoundaryNodes[i]);
            }
        }
    }

    void FindLeftBoundaryNodes(const Node&               rNode,
                               const std::vector<Node*>& rBoundaryNodes,
                               std::vector<Node*>&       rLeftBoundaryNodes) const
    {
        for (unsigned int i = 0; i < rBoundaryNodes.size(); ++i) {
            if (rBoundaryNodes[i]->Coordinates()[mHorizontalDirection] <= rNode.Coordinates()[mHorizontalDirection]) {
                // node is on top boundary
                rLeftBoundaryNodes.push_back(rBoundaryNodes[i]);
            }
        }
    }

    void FindRightBoundaryNodes(const Node&               rNode,
                                const std::vector<Node*>& rBoundaryNodes,
                                std::vector<Node*>&       rRightBoundaryNodes) const
    {
        for (unsigned int i = 0; i < rBoundaryNodes.size(); ++i) {
            if (rBoundaryNodes[i]->Coordinates()[mHorizontalDirection] >= rNode.Coordinates()[mHorizontalDirection]) {
                // node is on top boundary
                rRightBoundaryNodes.push_back(rBoundaryNodes[i]);
            }
        }
    }

    int GetMaxNodeID()
    {
        KRATOS_TRY

        int MaxNodeID = -1;
        block_for_each(mrModelPart.Nodes(), [&MaxNodeID](const Node& rNode) {
            #pragma omp critical
            MaxNodeID = std::max<int>(MaxNodeID, rNode.Id());
        });

        return MaxNodeID;

        KRATOS_CATCH("")
    }

    void FindBoundaryNodes()
    {
        KRATOS_TRY

        std::vector<int> BoundaryNodes;

        FillListOfBoundaryNodesFast(BoundaryNodes);
        mBoundaryNodes.resize(BoundaryNodes.size());

        unsigned int iPosition = 0;
        block_for_each(mrModelPart.Nodes(), [&iPosition, &BoundaryNodes, this](Node& rNode) {
            const int Id = rNode.Id();
            for (unsigned int j = 0; j < BoundaryNodes.size(); ++j) {
                if (Id == BoundaryNodes[j]) {
                    mBoundaryNodes[iPosition++] = &rNode;
                }
            }
        });

        KRATOS_CATCH("")
    }

    void FillListOfBoundaryNodesFast(std::vector<int> &BoundaryNodes)
    {
        const int ID_UNDEFINED = -1;
        const int N_ELEMENT = 10;

        std::vector<std::vector<int>> ELementsOfNodes;
        std::vector<int> ELementsOfNodesSize;

        int MaxNodeID = GetMaxNodeID();

        ELementsOfNodes.resize(MaxNodeID);
        ELementsOfNodesSize.resize(MaxNodeID);

        for (unsigned int i=0; i < ELementsOfNodes.size(); ++i) {
            ELementsOfNodes[i].resize(N_ELEMENT);
            ELementsOfNodesSize[i] = 0;
            std::fill(ELementsOfNodes[i].begin(), ELementsOfNodes[i].end(), ID_UNDEFINED);
        }

        const unsigned int nElements = mrModelPart.NumberOfElements();
        ModelPart::ElementsContainerType::iterator it_begin_elements = mrModelPart.ElementsBegin();

        for (unsigned int i=0; i < nElements; ++i) {
            ModelPart::ElementsContainerType::iterator pElemIt = it_begin_elements + i;
            for (unsigned int iPoint=0; iPoint < pElemIt->GetGeometry().PointsNumber(); ++iPoint) {
               int NodeID = pElemIt->GetGeometry()[iPoint].Id();
               int ElementId = pElemIt->Id();

               int index = NodeID-1;
               ELementsOfNodesSize[index]++;
               if (ELementsOfNodesSize[index] > N_ELEMENT-1) {
                   ELementsOfNodes[index].push_back(ElementId);
               } else {
                   ELementsOfNodes[index][ELementsOfNodesSize[index]-1] = ElementId;
               }
            }
        }

        for (unsigned int i=0; i < nElements; ++i) {
            ModelPart::ElementsContainerType::iterator pElemIt = it_begin_elements + i;

            int nEdges = pElemIt->GetGeometry().EdgesNumber();
            for (int iEdge = 0; iEdge < nEdges; ++iEdge) {
                const unsigned int nPoints = pElemIt->GetGeometry().GenerateEdges()[iEdge].PointsNumber();
                std::vector<int> FaceID(nPoints);
                for (unsigned int iPoint = 0; iPoint < nPoints; ++iPoint) {
                    FaceID[iPoint] = pElemIt->GetGeometry().GenerateEdges()[iEdge].GetPoint(iPoint).Id();
                }

                if (!IsMoreThanOneElementWithThisEdgeFast(FaceID, ELementsOfNodes, ELementsOfNodesSize)) {
                    // boundary nodes:
                    for (unsigned int iPoint = 0; iPoint < nPoints; ++iPoint) {
                        auto it = std::find(BoundaryNodes.begin(), BoundaryNodes.end(), FaceID[iPoint]);
                        if (it == BoundaryNodes.end()) {
                            BoundaryNodes.push_back(FaceID[iPoint]);
                        }
                    }
                }
            }
        }

        KRATOS_ERROR_IF(BoundaryNodes.empty())
            << "No boundary node is found for interpolate line pressure process" << std::endl;

    }

    bool IsMoreThanOneElementWithThisEdgeFast(const std::vector<int>&              rFaceIDs,
                                              const std::vector<std::vector<int>>& rELementsOfNodes,
                                              const std::vector<int>& rELementsOfNodesSize) const

    {
        const int ID_UNDEFINED = -1;
        int nMaxElements = 0;
        for (unsigned int iPoint = 0; iPoint < rFaceIDs.size(); ++iPoint) {
            int NodeID = rFaceIDs[iPoint];
            int index = NodeID-1;
            nMaxElements += rELementsOfNodesSize[index];
        }

        if (nMaxElements > 0) {
            std::vector<vector<int>> ElementIDs;
            ElementIDs.resize(rFaceIDs.size());
            for (unsigned int i=0; i<ElementIDs.size(); ++i) {
                ElementIDs[i].resize(nMaxElements);
                std::fill(ElementIDs[i].begin(), ElementIDs[i].end(), ID_UNDEFINED);
            }

            for (unsigned int iPoint = 0; iPoint < rFaceIDs.size(); ++iPoint) {
                int NodeID = rFaceIDs[iPoint];
                int index = NodeID-1;
                for (int i=0; i < rELementsOfNodesSize[index]; ++i) {
                    int iElementID = rELementsOfNodes[index][i];
                    ElementIDs[iPoint][i] = iElementID;
                }
            }

            std::vector<int> SharedElementIDs;
            for (unsigned int iPoint = 0; iPoint < rFaceIDs.size(); ++iPoint) {
                for (unsigned int i=0; i < ElementIDs[iPoint].size(); ++i) {
                    int iElementID = ElementIDs[iPoint][i];
                    bool found = false;
                    if (iElementID !=ID_UNDEFINED) {
                        for (unsigned int iPointInner = 0; iPointInner < rFaceIDs.size(); ++iPointInner) {
                            if (iPointInner != iPoint) {
                                //std::any_of followed by breaking out of 2 for loops
                                for (unsigned int j = 0; j < ElementIDs[iPointInner].size(); ++j) {
                                    if (ElementIDs[iPointInner][j]==iElementID) found = true;
                                }
                            }
                        }
                    }

                    if (found) {
                        auto it = std::find(SharedElementIDs.begin(), SharedElementIDs.end(), iElementID);
                        if (it == SharedElementIDs.end()) {
                            SharedElementIDs.push_back(iElementID);
                        }
                    }
                }
            }

            return SharedElementIDs.size() > 1;
        }

        return false;
    }

}; // Class ApplyConstantInterpolateLinePressureProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyConstantInterpolateLinePressureProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyConstantInterpolateLinePressureProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}