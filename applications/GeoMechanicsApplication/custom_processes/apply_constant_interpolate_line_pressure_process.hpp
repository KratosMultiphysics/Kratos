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

#if !defined(KRATOS_GEO_APPLY_CONSTANT_INTERPOLATE_LINE_PRESSURE_PROCESS )
#define  KRATOS_GEO_APPLY_CONSTANT_INTERPOLATE_LINE_PRESSURE_PROCESS

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

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
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
                "gravity_direction": 1,
                "out_of_plane_direction": 2,
                "pressure_tension_cut_off" : 0.0,
                "table" : 1
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();

        FindBoundaryNodes();

        mIsFixed = rParameters["is_fixed"].GetBool();
        mGravityDirection = rParameters["gravity_direction"].GetInt();
        mOutOfPlaneDirection = rParameters["out_of_plane_direction"].GetInt();
        if (mGravityDirection == mOutOfPlaneDirection)
            KRATOS_ERROR << "Gravity direction cannot be the same as Out-of-Plane directions"
                         << rParameters
                         << std::endl;

        mHorizontalDirection = 0;
        for (unsigned int i=0; i<N_DIM_3D; ++i)
           if (i!=mGravityDirection && i!=mOutOfPlaneDirection) mHorizontalDirection = i;

        if (rParameters.Has("pressure_tension_cut_off"))
          mPressureTensionCutOff = rParameters["pressure_tension_cut_off"].GetDouble();
        else
          mPressureTensionCutOff = 0.0;

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyConstantInterpolateLinePressureProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyConstantInterpolateLinePressureProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;

        const int nNodes = static_cast<int>(mrModelPart.Nodes().size());

        if (nNodes > 0) {
            const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);

            block_for_each(mrModelPart.Nodes(), [&var, this](Node<3>& rNode) {
                if (mIsFixed) rNode.Fix(var);
                else          rNode.Free(var);

                double pressure = CalculatePressure(rNode);

                if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < mPressureTensionCutOff) {
                    rNode.FastGetSolutionStepValue(var) = pressure;
                } else {
                    rNode.FastGetSolutionStepValue(var) = mPressureTensionCutOff;
                }
            });
        }

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyConstantInterpolateLinePressureProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyConstantInterpolateLinePressureProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mrModelPart;
    std::string mVariableName;
    bool mIsFixed;
    unsigned int mGravityDirection;
    unsigned int mOutOfPlaneDirection;
    unsigned int mHorizontalDirection;
    std::vector< Node<3> * > mBoundaryNodes;
    double mPressureTensionCutOff;

    double CalculatePressure(const Node<3> &rNode)
    {
        // find top boundary
        std::vector< Node<3>* > TopBoundaryNodes;
        FindTopBoundaryNodes(rNode, TopBoundaryNodes);
        double PressureTop;
        double CoordinateTop;
        CalculateBoundaryPressure(rNode, TopBoundaryNodes, PressureTop, CoordinateTop);

        // find bottom boundary
        std::vector< Node<3>* > BottomBoundaryNodes;
        FindBottomBoundaryNodes(rNode, BottomBoundaryNodes);
        double PressureBottom;
        double CoordinateBottom;
        CalculateBoundaryPressure(rNode, BottomBoundaryNodes, PressureBottom, CoordinateBottom, true);

        // calculate pressure
        if (std::abs(CoordinateTop - CoordinateBottom) > TINY) {
            const double slopeP = (PressureTop - PressureBottom) / (CoordinateTop - CoordinateBottom);
            const double pressure = slopeP * (rNode.Coordinates()[mGravityDirection] - CoordinateBottom ) + PressureBottom;
            return pressure;
        } else {
            return PressureBottom;
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ApplyConstantInterpolateLinePressureProcess& operator=(ApplyConstantInterpolateLinePressureProcess const& rOther);

    /// Copy constructor.
    //ApplyConstantInterpolateLinePressureProcess(ApplyConstantInterpolateLinePressureProcess const& rOther);

    void CalculateBoundaryPressure( const Node<3> &rNode,
                                    const std::vector< Node<3>*> &BoundaryNodes,
                                    double &pressure,
                                    double &coordinate,
                                    bool isBottom=false )
    {
        // find top boundary
        std::vector< Node<3>*> LeftBoundaryNodes;
        FindLeftBoundaryNodes(rNode, BoundaryNodes, LeftBoundaryNodes);

        std::vector< Node<3>*> RightBoundaryNodes;
        FindRightBoundaryNodes(rNode, BoundaryNodes, RightBoundaryNodes);

        if (LeftBoundaryNodes.size() > 0 && RightBoundaryNodes.size() > 0) {
            Node<3> *LeftNode;
            LeftNode = FindClosestNodeOnBoundaryNodes(rNode, LeftBoundaryNodes, isBottom);
            Node<3> *RightNode;
            RightNode = FindClosestNodeOnBoundaryNodes(rNode, RightBoundaryNodes, isBottom);

            InterpolateBoundaryPressure(rNode, LeftNode, RightNode, pressure, coordinate);
            return;

        } else if (LeftBoundaryNodes.size() > 0) {
            InterpolateBoundaryPressureWithOneContainer(rNode, LeftBoundaryNodes, pressure, coordinate);
            return;

        } else if (RightBoundaryNodes.size() > 0) {
            InterpolateBoundaryPressureWithOneContainer(rNode, RightBoundaryNodes, pressure, coordinate);
            return;

        } else {
            KRATOS_ERROR << "There is not enough points around interpolation, node Id" << rNode.Id() << std::endl;
        }

    }

    void InterpolateBoundaryPressureWithOneContainer(const Node<3> &rNode,
                                                     const std::vector< Node<3>*> &BoundaryNodes,
                                                     double &pressure,
                                                     double &coordinate )
    {
        std::vector< Node<3>*> FoundNodes;
        FindTwoClosestNodeOnBoundaryNodes(rNode, BoundaryNodes, FoundNodes);

        const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);

        const double &pressureLeft = FoundNodes[0]->FastGetSolutionStepValue(var);
        Vector3 CoordinatesLeft;
        noalias(CoordinatesLeft) = FoundNodes[0]->Coordinates();

        const double &pressureRight = FoundNodes[1]->FastGetSolutionStepValue(var);
        Vector3 CoordinatesRight;
        noalias(CoordinatesRight) = FoundNodes[1]->Coordinates();

        // calculate pressure
        if (std::abs(CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) > TINY) {
            const double slopeP = (pressureRight - pressureLeft) / (CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]);
            pressure = slopeP * (rNode.Coordinates()[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) + pressureLeft;

            const double slopeY = (CoordinatesRight[mGravityDirection] - CoordinatesLeft[mGravityDirection]) / (CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]);
            coordinate = slopeY * (rNode.Coordinates()[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) + CoordinatesLeft[mGravityDirection];
        } else {
            pressure   = pressureLeft;
            coordinate = CoordinatesLeft[mGravityDirection];
        }
    }

    void InterpolateBoundaryPressure(const Node<3> &rNode,
                                     const Node<3> *LeftNode,
                                     const Node<3> *RightNode,
                                     double &pressure,
                                     double &coordinate )
    {
        const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);

        const double &pressureLeft = LeftNode->FastGetSolutionStepValue(var);
        Vector3 CoordinatesLeft;
        noalias(CoordinatesLeft) = LeftNode->Coordinates();

        const double &pressureRight = RightNode->FastGetSolutionStepValue(var);
        Vector3 CoordinatesRight;
        noalias(CoordinatesRight) = RightNode->Coordinates();

        // calculate pressure
        if (std::abs(CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) > TINY) {
            const double slopeP = (pressureRight - pressureLeft) / (CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]);
            pressure = slopeP * (rNode.Coordinates()[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) + pressureLeft;

            const double slopeY = (CoordinatesRight[mGravityDirection] - CoordinatesLeft[mGravityDirection]) / (CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]);
            coordinate = slopeY * (rNode.Coordinates()[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) + CoordinatesLeft[mGravityDirection];
        } else {
            pressure   = pressureLeft;
            coordinate = CoordinatesLeft[mGravityDirection];
        }
    }

    void FindTwoClosestNodeOnBoundaryNodes(const Node<3> &rNode,
                                           const std::vector< Node<3>*> &BoundaryNodes,
                                           std::vector< Node<3>*> &FoundNodes)
    {
        const double HorizontalCoordiante = rNode.Coordinates()[mHorizontalDirection];
        FoundNodes.resize(2);

        unsigned int nFound = 0;
        double horizontalDistanceClosest_1 = LARGE;
        for (unsigned int i = 0; i < BoundaryNodes.size(); ++i) {
            Vector3 CoordinatesBoundary;
            noalias(CoordinatesBoundary) = BoundaryNodes[i]->Coordinates();

            if (std::abs(CoordinatesBoundary[mHorizontalDirection] - HorizontalCoordiante) <= horizontalDistanceClosest_1) {
                horizontalDistanceClosest_1 = std::abs(CoordinatesBoundary[mHorizontalDirection] - HorizontalCoordiante);
                FoundNodes[0] = BoundaryNodes[i];
                nFound++;
            }
        }

        double horizontalDistanceClosest_2 = LARGE;
        for (unsigned int i = 0; i < BoundaryNodes.size(); ++i) {
            Vector3 CoordinatesBoundary;
            noalias(CoordinatesBoundary) = BoundaryNodes[i]->Coordinates();

            if (std::abs(CoordinatesBoundary[mHorizontalDirection] - HorizontalCoordiante) <= horizontalDistanceClosest_2 &&
                std::abs(CoordinatesBoundary[mHorizontalDirection] - HorizontalCoordiante) > horizontalDistanceClosest_1) {
                horizontalDistanceClosest_2 = std::abs(CoordinatesBoundary[mHorizontalDirection] - HorizontalCoordiante);
                FoundNodes[1] = BoundaryNodes[i];
                nFound++;
            }
        }

        KRATOS_ERROR_IF(nFound < 2) << "Not enough points for interpolation: Coordinates"<< rNode.Coordinates() << std::endl;
    }


    Node<3>* 
        FindClosestNodeOnBoundaryNodes(const Node<3> &rNode,
                                       const std::vector< Node<3>* > &BoundaryNodes,
                                       const bool isBottom)
    {
        const double HorizontalCoordiante = rNode.Coordinates()[mHorizontalDirection];
        Node<3> *pNode;
        std::vector< Node<3>*> FoundNodes;

        double horizontalDistance = LARGE;
        for (unsigned int i = 0; i < BoundaryNodes.size(); ++i) {
            if (std::abs(BoundaryNodes[i]->Coordinates()[mHorizontalDirection] - HorizontalCoordiante) <= horizontalDistance) {
                horizontalDistance = std::abs(BoundaryNodes[i]->Coordinates()[mHorizontalDirection] - HorizontalCoordiante);
                FoundNodes.push_back(BoundaryNodes[i]);
            }
        }

        if (isBottom) {
            double height = LARGE;
            for (unsigned int i = 0; i < FoundNodes.size(); ++i) {
                if (FoundNodes[i]->Coordinates()[mGravityDirection] < height) {
                    pNode = FoundNodes[i];
                    height = FoundNodes[i]->Coordinates()[mGravityDirection];
                }
            }
        } else {
            double height = -LARGE;
            for (unsigned int i = 0; i < FoundNodes.size(); ++i) {
                if (FoundNodes[i]->Coordinates()[mGravityDirection] > height) {
                    pNode = FoundNodes[i];
                    height = FoundNodes[i]->Coordinates()[mGravityDirection];
                }
            }
        }

        return pNode;
    }

    void FindTopBoundaryNodes(const Node<3> &rNode,
                              std::vector< Node<3>* > &TopBoundaryNodes)
    {
        for (unsigned int i = 0; i < mBoundaryNodes.size(); ++i) {
            if (mBoundaryNodes[i]->Coordinates()[mGravityDirection] >= rNode.Coordinates()[mGravityDirection]) {
                // node is on top boundary
                TopBoundaryNodes.push_back(mBoundaryNodes[i]);
            }
        }
    }

    void FindBottomBoundaryNodes(const Node<3> &rNode,
                                 std::vector< Node<3>*> &BottomBoundaryNodes)
    {
        for (unsigned int i = 0; i < mBoundaryNodes.size(); ++i) {
            if (mBoundaryNodes[i]->Coordinates()[mGravityDirection] <= rNode.Coordinates()[mGravityDirection]) {
                // node is on top boundary
                BottomBoundaryNodes.push_back(mBoundaryNodes[i]);
            }
        }
    }

    void FindLeftBoundaryNodes(const Node<3> &rNode,
                               const std::vector< Node<3>*> &BoundaryNodes,
                               std::vector< Node<3>*> &LeftBoundaryNodes)
    {
        for (unsigned int i = 0; i < BoundaryNodes.size(); ++i) {
            if (BoundaryNodes[i]->Coordinates()[mHorizontalDirection] <= rNode.Coordinates()[mHorizontalDirection]) {
                // node is on top boundary
                LeftBoundaryNodes.push_back(BoundaryNodes[i]);
            }
        }
    }

    void FindRightBoundaryNodes(const Node<3> &rNode,
                                const std::vector< Node<3>*> &BoundaryNodes,
                                std::vector< Node<3>*> &RightBoundaryNodes)
    {
        for (unsigned int i = 0; i < BoundaryNodes.size(); ++i) {
            if (BoundaryNodes[i]->Coordinates()[mHorizontalDirection] >= rNode.Coordinates()[mHorizontalDirection]) {
                // node is on top boundary
                RightBoundaryNodes.push_back(BoundaryNodes[i]);
            }
        }
    }

    int GetMaxNodeID()
    {
        KRATOS_TRY

        int MaxNodeID = -1;
        block_for_each(mrModelPart.Nodes(), [&MaxNodeID](Node<3>& rNode) {
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
        block_for_each(mrModelPart.Nodes(), [&iPosition, &BoundaryNodes, this](Node<3>& rNode) {
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
                        std::vector<int>::iterator it = std::find(BoundaryNodes.begin(), BoundaryNodes.end(), FaceID[iPoint]);
                        if (it == BoundaryNodes.end()) {
                            BoundaryNodes.push_back(FaceID[iPoint]);
                        }
                    }
                }
            }
        }
        // KRATOS_INFO("FillListOfBoundaryNodesFast:BoundaryNodes") << BoundaryNodes << std::endl;

        KRATOS_ERROR_IF(BoundaryNodes.size()==0)
            << "No boundary node is found for interpolate line pressure process" << std::endl;

    }

    bool IsMoreThanOneElementWithThisEdgeFast(const std::vector<int> &FaceID,
                                              const std::vector<std::vector<int>> &ELementsOfNodes,
                                              const std::vector<int> &ELementsOfNodesSize)

    {
        const int ID_UNDEFINED = -1;
        int nMaxElements = 0;
        for (unsigned int iPoint = 0; iPoint < FaceID.size(); ++iPoint) {
            int NodeID = FaceID[iPoint];
            int index = NodeID-1;
            nMaxElements += ELementsOfNodesSize[index];
        }

        if (nMaxElements > 0) {
            std::vector<vector<int>> ElementIDs;
            ElementIDs.resize(FaceID.size());
            for (unsigned int i=0; i<ElementIDs.size(); ++i) {
                ElementIDs[i].resize(nMaxElements);
                std::fill(ElementIDs[i].begin(), ElementIDs[i].end(), ID_UNDEFINED);
            }
            

            for (unsigned int iPoint = 0; iPoint < FaceID.size(); ++iPoint) {
                int NodeID = FaceID[iPoint];
                int index = NodeID-1;
                for (int i=0; i < ELementsOfNodesSize[index]; ++i) {
                    int iElementID = ELementsOfNodes[index][i];
                    ElementIDs[iPoint][i] = iElementID;
                }
            }


            std::vector<int> SharedElementIDs;
            for (unsigned int iPoint = 0; iPoint < FaceID.size(); ++iPoint) {
                for (unsigned int i=0; i < ElementIDs[iPoint].size(); ++i) {
                    int iElementID = ElementIDs[iPoint][i];
                    bool found = false;
                    if (iElementID !=ID_UNDEFINED) {
                        for (unsigned int iPointInner = 0; iPointInner < FaceID.size(); ++iPointInner) {
                            if (iPointInner != iPoint) {
                                for (unsigned int j = 0; j < ElementIDs[iPointInner].size(); ++j) {
                                    if (ElementIDs[iPointInner][j]==iElementID) found = true;
                                }
                            }
                        }
                    }

                    if (found) {
                        std::vector<int>::iterator it = std::find(SharedElementIDs.begin(), SharedElementIDs.end(), iElementID);
                        if (it == SharedElementIDs.end()) {
                            SharedElementIDs.push_back(iElementID);
                        }
                    }
                }
            }

            if (SharedElementIDs.size() > 1)
                return true;
            else
                return false;
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

} // namespace Kratos.

#endif /* KRATOS_GEO_APPLY_CONSTANT_INTERPOLATE_LINE_PRESSURE_PROCESS defined */
