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

#include "apply_constant_interpolate_line_pressure_process.h"
#include "geo_mechanics_application_variables.h"
#include "includes/kratos_flags.h"
#include "includes/model_part.h"

#include <algorithm>

namespace Kratos
{
using namespace std::string_literals;

ApplyConstantInterpolateLinePressureProcess ::ApplyConstantInterpolateLinePressureProcess(ModelPart& model_part,
                                                                                          Parameters rParameters)
    : Process(Flags()), mrModelPart(model_part)
{
    KRATOS_TRY

    // only include validation with c++11 since raw_literals do not exist in c++03
    Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed": false,
                "is_seepage": false,
                "gravity_direction": 1,
                "out_of_plane_direction": 2,
                "pressure_tension_cut_off" : 0.0,
                "table" : 1
            }  )");

    // Some values need to be mandatory prescribed since no meaningful default value exist. For
    // this reason try accessing to them So that an error is thrown if they don't exist
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
                     << rParameters << std::endl;

    mHorizontalDirection = 0;
    for (unsigned int i = 0; i < N_DIM_3D; ++i)
        if (i != mGravityDirection && i != mOutOfPlaneDirection) mHorizontalDirection = i;

    mPressureTensionCutOff = rParameters["pressure_tension_cut_off"].GetDouble();

    KRATOS_CATCH("")
}

void ApplyConstantInterpolateLinePressureProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY
    if (mIsInitialized) return;
    const Variable<double>& var = KratosComponents<Variable<double>>::Get(mVariableName);

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
            rNode.FastGetSolutionStepValue(var) =
                std::min(pressure, PORE_PRESSURE_SIGN_FACTOR * mPressureTensionCutOff);
            if (mIsFixed) rNode.Fix(var);
            else if (mIsFixedProvided) rNode.Free(var);
        }
    });
    mIsInitialized = true;
    KRATOS_CATCH("")
}

std::string ApplyConstantInterpolateLinePressureProcess::Info() const
{
    return "ApplyConstantInterpolateLinePressureProcess"s;
}

double ApplyConstantInterpolateLinePressureProcess::CalculatePressure(const Node& rNode)
{
    // find top boundary
    std::vector<Node*> TopBoundaryNodes;
    FindTopBoundaryNodes(rNode, TopBoundaryNodes);
    double PressureTop;
    double CoordinateTop;
    CalculateBoundaryPressure(rNode, TopBoundaryNodes, PressureTop, CoordinateTop);

    // find bottom boundary
    std::vector<Node*> BottomBoundaryNodes;
    FindBottomBoundaryNodes(rNode, BottomBoundaryNodes);
    double PressureBottom;
    double CoordinateBottom;
    CalculateBoundaryPressure(rNode, BottomBoundaryNodes, PressureBottom, CoordinateBottom, true);

    // calculate pressure
    if (std::abs(CoordinateTop - CoordinateBottom) > TINY) {
        const double slopeP = (PressureTop - PressureBottom) / (CoordinateTop - CoordinateBottom);
        return slopeP * (rNode.Coordinates()[mGravityDirection] - CoordinateBottom) + PressureBottom;
    } else {
        return PressureBottom;
    }
}

void ApplyConstantInterpolateLinePressureProcess::CalculateBoundaryPressure(
    const Node& rNode, const std::vector<Node*>& BoundaryNodes, double& pressure, double& coordinate, bool isBottom)
{
    // find top boundary
    std::vector<Node*> LeftBoundaryNodes;
    FindLeftBoundaryNodes(rNode, BoundaryNodes, LeftBoundaryNodes);

    std::vector<Node*> RightBoundaryNodes;
    FindRightBoundaryNodes(rNode, BoundaryNodes, RightBoundaryNodes);

    if (!LeftBoundaryNodes.empty() && !RightBoundaryNodes.empty()) {
        const Node* LeftNode  = FindClosestNodeOnBoundaryNodes(rNode, LeftBoundaryNodes, isBottom);
        const Node* RightNode = FindClosestNodeOnBoundaryNodes(rNode, RightBoundaryNodes, isBottom);

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

void ApplyConstantInterpolateLinePressureProcess::InterpolateBoundaryPressureWithOneContainer(
    const Node& rNode, const std::vector<Node*>& rBoundaryNodes, double& rPressure, double& rCoordinate) const
{
    std::vector<Node*> FoundNodes;
    FindTwoClosestNodeOnBoundaryNodes(rNode, rBoundaryNodes, FoundNodes);

    const Variable<double>& var = KratosComponents<Variable<double>>::Get(mVariableName);

    const double& pressureLeft = FoundNodes[0]->FastGetSolutionStepValue(var);
    Vector3       CoordinatesLeft;
    noalias(CoordinatesLeft) = FoundNodes[0]->Coordinates();

    const double& pressureRight = FoundNodes[1]->FastGetSolutionStepValue(var);
    Vector3       CoordinatesRight;
    noalias(CoordinatesRight) = FoundNodes[1]->Coordinates();

    // calculate pressure
    if (std::abs(CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) > TINY) {
        const auto slope_p = (pressureRight - pressureLeft) / (CoordinatesRight[mHorizontalDirection] -
                                                               CoordinatesLeft[mHorizontalDirection]);
        rPressure =
            slope_p * (rNode.Coordinates()[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) + pressureLeft;

        const auto slope_y =
            (CoordinatesRight[mGravityDirection] - CoordinatesLeft[mGravityDirection]) /
            (CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]);
        rCoordinate =
            slope_y * (rNode.Coordinates()[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) +
            CoordinatesLeft[mGravityDirection];
    } else {
        rPressure   = pressureLeft;
        rCoordinate = CoordinatesLeft[mGravityDirection];
    }
}

void ApplyConstantInterpolateLinePressureProcess::InterpolateBoundaryPressure(const Node& rNode,
                                                                              const Node* pLeftNode,
                                                                              const Node* pRightNode,
                                                                              double& rPressure,
                                                                              double& rCoordinate) const
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
        rPressure = slope_p * (rNode.Coordinates()[mHorizontalDirection] - coordinates_left[mHorizontalDirection]) +
                    pressure_left;

        const auto slope_y =
            (coordinates_right[mGravityDirection] - coordinates_left[mGravityDirection]) /
            (coordinates_right[mHorizontalDirection] - coordinates_left[mHorizontalDirection]);
        rCoordinate =
            slope_y * (rNode.Coordinates()[mHorizontalDirection] - coordinates_left[mHorizontalDirection]) +
            coordinates_left[mGravityDirection];
    } else {
        rPressure   = pressure_left;
        rCoordinate = coordinates_left[mGravityDirection];
    }
}

void ApplyConstantInterpolateLinePressureProcess::FindTwoClosestNodeOnBoundaryNodes(
    const Node& rNode, const std::vector<Node*>& rBoundaryNodes, std::vector<Node*>& rFoundNodes) const
{
    const double HorizontalCoordinate = rNode.Coordinates()[mHorizontalDirection];
    rFoundNodes.resize(2);

    unsigned int nFound                      = 0;
    double       horizontalDistanceClosest_1 = LARGE;
    for (auto boundary_node : rBoundaryNodes) {
        Vector3 CoordinatesBoundary;
        noalias(CoordinatesBoundary) = boundary_node->Coordinates();

        if (std::abs(CoordinatesBoundary[mHorizontalDirection] - HorizontalCoordinate) <= horizontalDistanceClosest_1) {
            horizontalDistanceClosest_1 =
                std::abs(CoordinatesBoundary[mHorizontalDirection] - HorizontalCoordinate);
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
            horizontalDistanceClosest_2 =
                std::abs(CoordinatesBoundary[mHorizontalDirection] - HorizontalCoordinate);
            rFoundNodes[1] = boundary_node;
            nFound++;
        }
    }

    KRATOS_ERROR_IF(nFound < 2) << "Not enough points for interpolation: Coordinates"
                                << rNode.Coordinates() << std::endl;
}

Node* ApplyConstantInterpolateLinePressureProcess::FindClosestNodeOnBoundaryNodes(
    const Node& rNode, const std::vector<Node*>& BoundaryNodes, const bool isBottom)
{
    const double       HorizontalCoordinate = rNode.Coordinates()[mHorizontalDirection];
    std::vector<Node*> FoundNodes;

    double horizontalDistance = LARGE;
    for (auto boundary_node : BoundaryNodes) {
        if (std::abs(boundary_node->Coordinates()[mHorizontalDirection] - HorizontalCoordinate) <= horizontalDistance) {
            horizontalDistance =
                std::abs(boundary_node->Coordinates()[mHorizontalDirection] - HorizontalCoordinate);
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

void ApplyConstantInterpolateLinePressureProcess::FindTopBoundaryNodes(const Node& rNode,
                                                                       std::vector<Node*>& TopBoundaryNodes) const
{
    for (const auto& p_boundary_node : mBoundaryNodes) {
        if (p_boundary_node->Coordinates()[mGravityDirection] >= rNode.Coordinates()[mGravityDirection]) {
            TopBoundaryNodes.push_back(p_boundary_node);
        }
    }
}

void ApplyConstantInterpolateLinePressureProcess::FindBottomBoundaryNodes(const Node& rNode,
                                                                          std::vector<Node*>& BottomBoundaryNodes) const
{
    for (const auto& p_boundary_node : mBoundaryNodes) {
        if (p_boundary_node->Coordinates()[mGravityDirection] <= rNode.Coordinates()[mGravityDirection]) {
            BottomBoundaryNodes.push_back(p_boundary_node);
        }
    }
}

void ApplyConstantInterpolateLinePressureProcess::FindLeftBoundaryNodes(const Node& rNode,
                                                                        const std::vector<Node*>& rBoundaryNodes,
                                                                        std::vector<Node*>& rLeftBoundaryNodes) const
{
    for (const auto& p_boundary_node : rBoundaryNodes) {
        if (p_boundary_node->Coordinates()[mHorizontalDirection] <= rNode.Coordinates()[mHorizontalDirection]) {
            rLeftBoundaryNodes.push_back(p_boundary_node);
        }
    }
}

void ApplyConstantInterpolateLinePressureProcess::FindRightBoundaryNodes(const Node& rNode,
                                                                         const std::vector<Node*>& rBoundaryNodes,
                                                                         std::vector<Node*>& rRightBoundaryNodes) const
{
    for (const auto& p_boundary_node : rBoundaryNodes) {
        if (p_boundary_node->Coordinates()[mHorizontalDirection] >= rNode.Coordinates()[mHorizontalDirection]) {
            rRightBoundaryNodes.push_back(p_boundary_node);
        }
    }
}

int ApplyConstantInterpolateLinePressureProcess::GetMaxNodeID()
{
    KRATOS_TRY

    int MaxNodeID = -1;
    block_for_each(mrModelPart.Nodes(), [&MaxNodeID](const Node& rNode) {
#pragma omp critical
        MaxNodeID = std::max(MaxNodeID, static_cast<int>(rNode.Id()));
    });

    return MaxNodeID;

    KRATOS_CATCH("")
}

void ApplyConstantInterpolateLinePressureProcess::FindBoundaryNodes()
{
    KRATOS_TRY

    std::vector<int> BoundaryNodes;

    FillListOfBoundaryNodesFast(BoundaryNodes);
    mBoundaryNodes.resize(BoundaryNodes.size());

    unsigned int iPosition = 0;
    block_for_each(mrModelPart.Nodes(), [&iPosition, &BoundaryNodes, this](Node& rNode) {
        const auto Id = static_cast<int>(rNode.Id());
        for (const auto& r_boundary_node : BoundaryNodes) {
            if (Id == r_boundary_node) {
                mBoundaryNodes[iPosition] = &rNode;
                iPosition++;
            }
        }
    });

    KRATOS_CATCH("")
}

void ApplyConstantInterpolateLinePressureProcess::FillListOfBoundaryNodesFast(std::vector<int>& BoundaryNodes)
{
    const int ID_UNDEFINED = -1;
    const int N_ELEMENT    = 10;

    std::vector<std::vector<int>> ELementsOfNodes;
    std::vector<int>              ELementsOfNodesSize;

    int MaxNodeID = GetMaxNodeID();

    ELementsOfNodes.resize(MaxNodeID);
    ELementsOfNodesSize.resize(MaxNodeID);

    for (unsigned int i = 0; i < ELementsOfNodes.size(); ++i) {
        ELementsOfNodes[i].resize(N_ELEMENT);
        ELementsOfNodesSize[i] = 0;
        std::fill(ELementsOfNodes[i].begin(), ELementsOfNodes[i].end(), ID_UNDEFINED);
    }

    const auto nElements = static_cast<unsigned int>(mrModelPart.NumberOfElements());
    ModelPart::ElementsContainerType::iterator it_begin_elements = mrModelPart.ElementsBegin();

    for (unsigned int i = 0; i < nElements; ++i) {
        ModelPart::ElementsContainerType::iterator pElemIt = it_begin_elements + i;
        for (unsigned int iPoint = 0; iPoint < pElemIt->GetGeometry().PointsNumber(); ++iPoint) {
            auto NodeID    = static_cast<int>(pElemIt->GetGeometry()[iPoint].Id());
            auto ElementId = static_cast<int>(pElemIt->Id());

            int index = NodeID - 1;
            ELementsOfNodesSize[index]++;
            if (ELementsOfNodesSize[index] > N_ELEMENT - 1) {
                ELementsOfNodes[index].push_back(ElementId);
            } else {
                ELementsOfNodes[index][ELementsOfNodesSize[index] - 1] = ElementId;
            }
        }
    }

    for (unsigned int i = 0; i < nElements; ++i) {
        ModelPart::ElementsContainerType::iterator pElemIt = it_begin_elements + i;

        auto nEdges = static_cast<int>(pElemIt->GetGeometry().EdgesNumber());
        for (int iEdge = 0; iEdge < nEdges; ++iEdge) {
            const auto nPoints =
                static_cast<unsigned int>(pElemIt->GetGeometry().GenerateEdges()[iEdge].PointsNumber());
            std::vector<int> FaceID(nPoints);
            for (unsigned int iPoint = 0; iPoint < nPoints; ++iPoint) {
                FaceID[iPoint] =
                    static_cast<int>(pElemIt->GetGeometry().GenerateEdges()[iEdge].GetPoint(iPoint).Id());
            }

            if (IsMoreThanOneElementWithThisEdgeFast(FaceID, ELementsOfNodes, ELementsOfNodesSize))
                continue;
            auto add_boundary_node_if_missing = [&BoundaryNodes](int node_id) {
                if (std::ranges::find(BoundaryNodes, node_id) == BoundaryNodes.end()) {
                    BoundaryNodes.push_back(node_id);
                }
            };
            for (unsigned int iPoint = 0; iPoint < nPoints; ++iPoint) {
                add_boundary_node_if_missing(FaceID[iPoint]);
            }
        }
    }

    KRATOS_ERROR_IF(BoundaryNodes.empty())
        << "No boundary node is found for interpolate line pressure process" << std::endl;
}

bool ApplyConstantInterpolateLinePressureProcess::IsMoreThanOneElementWithThisEdgeFast(
    const std::vector<int>&              rFaceIDs,
    const std::vector<std::vector<int>>& rELementsOfNodes,
    const std::vector<int>&              rELementsOfNodesSize) const

{
    int nMaxElements = 0;
    for (auto node_id : rFaceIDs) {
        nMaxElements += rELementsOfNodesSize[node_id - 1];
    }

    if (nMaxElements == 0) return false;

    constexpr int ID_UNDEFINED = -1;

    std::vector<vector<int>> ElementIDs;
    ElementIDs.resize(rFaceIDs.size());
    for (auto& r_element_id : ElementIDs) {
        r_element_id.resize(nMaxElements);
        std::ranges::fill(r_element_id, ID_UNDEFINED);
    }

    for (unsigned int iPoint = 0; iPoint < rFaceIDs.size(); ++iPoint) {
        int NodeID = rFaceIDs[iPoint];
        int index  = NodeID - 1;
        for (int i = 0; i < rELementsOfNodesSize[index]; ++i) {
            int iElementID        = rELementsOfNodes[index][i];
            ElementIDs[iPoint][i] = iElementID;
        }
    }

    std::vector<int> SharedElementIDs;
    for (unsigned int iPoint = 0; iPoint < rFaceIDs.size(); ++iPoint) {
        for (const auto element_id : ElementIDs[iPoint]) {
            bool found = false;
            if (element_id == ID_UNDEFINED) continue;
            for (unsigned int iPointInner = 0; iPointInner < rFaceIDs.size(); ++iPointInner) {
                if (iPointInner == iPoint) continue;
                // std::any_of followed by breaking out of 2 for loops
                for (const auto element_id_loop : ElementIDs[iPointInner]) {
                    if (element_id_loop == element_id) found = true;
                }
            }

            if (!found) continue;
            auto it = std::find(SharedElementIDs.begin(), SharedElementIDs.end(), element_id);
            if (it == SharedElementIDs.end()) {
                SharedElementIDs.push_back(element_id);
            }
        }
    }

    return SharedElementIDs.size() > 1;
}

} // namespace Kratos