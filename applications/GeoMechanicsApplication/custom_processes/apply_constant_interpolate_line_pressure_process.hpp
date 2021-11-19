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

        if (nNodes != 0)
        {
            const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();


            #pragma omp parallel for
            for (int i = 0; i<nNodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                if (mIsFixed) it->Fix(var);
                else          it->Free(var);

                double pressure= CalculatePressure(it);

                if ((PORE_PRESSURE_SIGN_FACTOR * pressure) < mPressureTensionCutOff)
                {
                    it->FastGetSolutionStepValue(var) = pressure;
                }
                else
                {
                    it->FastGetSolutionStepValue(var) = mPressureTensionCutOff;
                }
            }
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
    std::vector< ModelPart::NodesContainerType::iterator> iterBoundaryNodes;
    double mPressureTensionCutOff;

    double CalculatePressure(const ModelPart::NodesContainerType::iterator &it)
    {
        Vector3 Coordinates;
        noalias(Coordinates) = it->Coordinates();

        // find top boundary
        std::vector< ModelPart::NodesContainerType::iterator> iterTopBoundaryNodes;
        FindTopBoundaryNodes(it, iterTopBoundaryNodes);
        double PressureTop;
        double CoordinateTop;
        CalculateBoundaryPressure(it, iterTopBoundaryNodes, PressureTop, CoordinateTop);

        // find bottom boundary
        std::vector< ModelPart::NodesContainerType::iterator> iterBottomBoundaryNodes;
        FindBottomBoundaryNodes(it, iterBottomBoundaryNodes);
        double PressureBottom;
        double CoordinateBottom;
        CalculateBoundaryPressure(it, iterBottomBoundaryNodes, PressureBottom, CoordinateBottom, true);

        // calculate pressure
        if (std::abs(CoordinateTop - CoordinateBottom) > TINY)
        {
            double slopeP = (PressureTop - PressureBottom) / (CoordinateTop - CoordinateBottom);
            double pressure = slopeP * (Coordinates[mGravityDirection] - CoordinateBottom ) + PressureBottom;
            return pressure;
        }
        else
        {
            return PressureBottom;
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ApplyConstantInterpolateLinePressureProcess& operator=(ApplyConstantInterpolateLinePressureProcess const& rOther);

    /// Copy constructor.
    //ApplyConstantInterpolateLinePressureProcess(ApplyConstantInterpolateLinePressureProcess const& rOther);

    void CalculateBoundaryPressure( const ModelPart::NodesContainerType::iterator &it,
                                    const std::vector< ModelPart::NodesContainerType::iterator> &iterBoundaryNodes,
                                    double& pressure,
                                    double& coordinate,
                                    bool isBottom=false )
    {
        Vector3 Coordinates;
        noalias(Coordinates) = it->Coordinates();

        // find top boundary
        std::vector< ModelPart::NodesContainerType::iterator> iterLeftBoundaryNodes;
        FindLeftBoundaryNodes(Coordinates, iterBoundaryNodes, iterLeftBoundaryNodes);

        std::vector< ModelPart::NodesContainerType::iterator> iterRightBoundaryNodes;
        FindRightBoundaryNodes(Coordinates, iterBoundaryNodes, iterRightBoundaryNodes);

        if (iterLeftBoundaryNodes.size() > 0 && iterRightBoundaryNodes.size() > 0)
        {
            ModelPart::NodesContainerType::iterator itLeft;
            itLeft = FindClosestNodeOnBoundaryNodes(Coordinates, iterLeftBoundaryNodes, isBottom);

            ModelPart::NodesContainerType::iterator itRight;
            itRight = FindClosestNodeOnBoundaryNodes(Coordinates, iterRightBoundaryNodes, isBottom);
            InterpolateBoundaryPressure(Coordinates, itLeft, itRight, pressure, coordinate);
            return;
        }
        else if (iterLeftBoundaryNodes.size() > 0)
        {
            InterpolateBoundaryPressureWithOneContainer(Coordinates, iterLeftBoundaryNodes, pressure, coordinate);
            return;
        }
        else if (iterRightBoundaryNodes.size() > 0)
        {
            InterpolateBoundaryPressureWithOneContainer(Coordinates, iterRightBoundaryNodes, pressure, coordinate);
            return;
        }
        else
        {
            KRATOS_INFO("CalculateBoundaryPressure:There is not enough points around interpolation, node Id") << it->Id() << std::endl;
            KRATOS_ERROR << "There is not enough points around interpolation, node Id" << it->Id() << std::endl;
        }

    }

    void InterpolateBoundaryPressureWithOneContainer(const Vector3& Coordinates,
                                                     const std::vector< ModelPart::NodesContainerType::iterator> &iterBoundaryNodes,
                                                     double &pressure,
                                                     double &coordinate )
    {
        std::vector< ModelPart::NodesContainerType::iterator> iterFoundNodes;
        FindTwoClosestNodeOnBoundaryNodes(Coordinates, iterBoundaryNodes, iterFoundNodes);

        const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);

        double pressureLeft = iterFoundNodes[0]->FastGetSolutionStepValue(var);
        Vector3 CoordinatesLeft;
        noalias(CoordinatesLeft) = iterFoundNodes[0]->Coordinates();

        double pressureRight = iterFoundNodes[1]->FastGetSolutionStepValue(var);
        Vector3 CoordinatesRight;
        noalias(CoordinatesRight) = iterFoundNodes[1]->Coordinates();

        // calculate pressure
        if (std::abs(CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) > TINY)
        {
            double slopeP = (pressureRight - pressureLeft) / (CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]);
            pressure = slopeP * (Coordinates[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) + pressureLeft;

            double slopeY = (CoordinatesRight[mGravityDirection] - CoordinatesLeft[mGravityDirection]) / (CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]);
            coordinate = slopeY * (Coordinates[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) + CoordinatesLeft[mGravityDirection];
        }
        else
        {
            pressure   = pressureLeft;
            coordinate = CoordinatesLeft[mGravityDirection];
        }
    }

    void InterpolateBoundaryPressure(const Vector3& Coordinates,
                                     const ModelPart::NodesContainerType::iterator &itLeft,
                                     const ModelPart::NodesContainerType::iterator &itRight,
                                     double &pressure,
                                     double &coordinate )
    {
        const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);

        double pressureLeft = itLeft->FastGetSolutionStepValue(var);
        Vector3 CoordinatesLeft;
        noalias(CoordinatesLeft) = itLeft->Coordinates();

        double pressureRight = itRight->FastGetSolutionStepValue(var);
        Vector3 CoordinatesRight;
        noalias(CoordinatesRight) = itRight->Coordinates();

        // calculate pressure
        if (std::abs(CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) > TINY)
        {
            double slopeP = (pressureRight - pressureLeft) / (CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]);
            pressure = slopeP * (Coordinates[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) + pressureLeft;

            double slopeY = (CoordinatesRight[mGravityDirection] - CoordinatesLeft[mGravityDirection]) / (CoordinatesRight[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]);
            coordinate = slopeY * (Coordinates[mHorizontalDirection] - CoordinatesLeft[mHorizontalDirection]) + CoordinatesLeft[mGravityDirection];
        }
        else
        {
            pressure   = pressureLeft;
            coordinate = CoordinatesLeft[mGravityDirection];
        }
    }

    void FindTwoClosestNodeOnBoundaryNodes(const Vector3 &Coordinates,
                                           const std::vector< ModelPart::NodesContainerType::iterator> &iterBoundaryNodes,
                                           std::vector< ModelPart::NodesContainerType::iterator> &iterFoundNodes)
    {
        const double Coordinate = Coordinates[mHorizontalDirection];
        iterFoundNodes.resize(2);

        unsigned int nFound = 0;
        double horizontalDistanceClosest_1 = LARGE;
        for (unsigned int i = 0; i < iterBoundaryNodes.size(); ++i)
        {
            Vector3 CoordinatesBoundary;
            noalias(CoordinatesBoundary) = iterBoundaryNodes[i]->Coordinates();

            if (std::abs(CoordinatesBoundary[mHorizontalDirection] - Coordinate) <= horizontalDistanceClosest_1)
            {
                horizontalDistanceClosest_1 = std::abs(CoordinatesBoundary[mHorizontalDirection] - Coordinate);
                iterFoundNodes[0] = iterBoundaryNodes[i];
                nFound++;
            }
        }

        double horizontalDistanceClosest_2 = LARGE;
        for (unsigned int i = 0; i < iterBoundaryNodes.size(); ++i)
        {
            Vector3 CoordinatesBoundary;
            noalias(CoordinatesBoundary) = iterBoundaryNodes[i]->Coordinates();

            if (std::abs(CoordinatesBoundary[mHorizontalDirection] - Coordinate) <= horizontalDistanceClosest_2 &&
                std::abs(CoordinatesBoundary[mHorizontalDirection] - Coordinate) > horizontalDistanceClosest_1)
            {
                horizontalDistanceClosest_2 = std::abs(CoordinatesBoundary[mHorizontalDirection] - Coordinate);
                iterFoundNodes[1] = iterBoundaryNodes[i];
                nFound++;
            }
        }

        if (nFound < 2)
        {
            KRATOS_INFO("FindTwoClosestNodeOnBoundaryNodes:There is not enough points around interpolation, Coordinates") << Coordinates << std::endl;
            KRATOS_ERROR << "Not enough points for interpolation: Coordinates"<< Coordinates << std::endl;
        }
    }


    ModelPart::NodesContainerType::iterator 
        FindClosestNodeOnBoundaryNodes(const Vector3 &Coordinates,
                                       const std::vector< ModelPart::NodesContainerType::iterator> &iterBoundaryNodes,
                                       const bool isBottom)
    {
        const double Coordinate = Coordinates[mHorizontalDirection];
        ModelPart::NodesContainerType::iterator it;
        double horizontalDistance = LARGE;
        std::vector< ModelPart::NodesContainerType::iterator> iterFoundNodes;

        for (unsigned int i = 0; i < iterBoundaryNodes.size(); ++i)
        {
            Vector3 CoordinatesBoundary;
            noalias(CoordinatesBoundary) = iterBoundaryNodes[i]->Coordinates();

            if (std::abs(CoordinatesBoundary[mHorizontalDirection] - Coordinate) <= horizontalDistance)
            {
                horizontalDistance = std::abs(CoordinatesBoundary[mHorizontalDirection] - Coordinate);
                iterFoundNodes.push_back(iterBoundaryNodes[i]);
            }
        }

        if (isBottom)
        {
            double height = LARGE;
            for (unsigned int i = 0; i < iterFoundNodes.size(); ++i)
            {
                Vector3 CoordinatesBoundary;
                noalias(CoordinatesBoundary) = iterFoundNodes[i]->Coordinates();
                if (CoordinatesBoundary[mGravityDirection] < height)
                {
                    it = iterFoundNodes[i];
                    height = CoordinatesBoundary[mGravityDirection];
                }
            }

        }
        else
        {
            double height = -LARGE;
            for (unsigned int i = 0; i < iterFoundNodes.size(); ++i)
            {
                Vector3 CoordinatesBoundary;
                noalias(CoordinatesBoundary) = iterFoundNodes[i]->Coordinates();
                if (CoordinatesBoundary[mGravityDirection] > height)
                {
                    it = iterFoundNodes[i];
                    height = CoordinatesBoundary[mGravityDirection];
                }
            }
        }

        return it;
    }

    void FindTopBoundaryNodes(const ModelPart::NodesContainerType::iterator &it,
                              std::vector< ModelPart::NodesContainerType::iterator> &iterTopBoundaryNodes)
    {
        Vector3 Coordinates;
        noalias(Coordinates) = it->Coordinates();

        for (unsigned int i = 0; i < iterBoundaryNodes.size(); ++i)
        {
            Vector3 CoordinatesBoundary;
            noalias(CoordinatesBoundary) = iterBoundaryNodes[i]->Coordinates();
            if (CoordinatesBoundary[mGravityDirection] >= Coordinates[mGravityDirection])
            {
                // it is on top boundary
                iterTopBoundaryNodes.push_back(iterBoundaryNodes[i]);
            }
        }
    }

    void FindBottomBoundaryNodes(const ModelPart::NodesContainerType::iterator &it,
                                 std::vector< ModelPart::NodesContainerType::iterator> &iterBottomBoundaryNodes)
    {
        Vector3 Coordinates;
        noalias(Coordinates) = it->Coordinates();

        for (unsigned int i = 0; i < iterBoundaryNodes.size(); ++i)
        {
            Vector3 CoordinatesBoundary;
            noalias(CoordinatesBoundary) = iterBoundaryNodes[i]->Coordinates();
            if (CoordinatesBoundary[mGravityDirection] <= Coordinates[mGravityDirection])
            {
                // it is on top boundary
                iterBottomBoundaryNodes.push_back(iterBoundaryNodes[i]);
            }
        }
    }

    void FindLeftBoundaryNodes(const Vector3 &Coordinates,
                               const std::vector< ModelPart::NodesContainerType::iterator> &iterBoundaryNodes,
                               std::vector< ModelPart::NodesContainerType::iterator> &iterLeftBoundaryNodes)
    {
        for (unsigned int i = 0; i < iterBoundaryNodes.size(); ++i)
        {
            Vector3 CoordinatesBoundary;
            noalias(CoordinatesBoundary) = iterBoundaryNodes[i]->Coordinates();
            if (CoordinatesBoundary[mHorizontalDirection] <= Coordinates[mHorizontalDirection])
            {
                // it is on top boundary
                iterLeftBoundaryNodes.push_back(iterBoundaryNodes[i]);
            }
        }
    }

    void FindRightBoundaryNodes(const Vector3 &Coordinates,
                                const std::vector< ModelPart::NodesContainerType::iterator> &iterBoundaryNodes,
                                std::vector< ModelPart::NodesContainerType::iterator> &iterRightBoundaryNodes)
    {
        for (unsigned int i = 0; i < iterBoundaryNodes.size(); ++i)
        {
            Vector3 CoordinatesBoundary;
            noalias(CoordinatesBoundary) = iterBoundaryNodes[i]->Coordinates();
            if (CoordinatesBoundary[mHorizontalDirection] >= Coordinates[mHorizontalDirection])
            {
                // it is on top boundary
                iterRightBoundaryNodes.push_back(iterBoundaryNodes[i]);
            }
        }
    }

    int GetMaxNodeID()
    {
        KRATOS_TRY;
        int MaxNodeID = -1;
        const int nNodes = mrModelPart.NumberOfNodes();
        ModelPart::NodesContainerType::iterator it_begin_nodes = mrModelPart.NodesBegin();
        #pragma omp parallel for
        for (int i = 0; i < nNodes; i++)
        {
            ModelPart::NodesContainerType::iterator it = it_begin_nodes + i;
            int id = it->Id();
            MaxNodeID = std::max(MaxNodeID, id);
        }
        return MaxNodeID;
        KRATOS_CATCH("");
    }

    void FindBoundaryNodes()
    {
        KRATOS_TRY;
        std::vector<int> BoundaryNodes;

        //FillListOfBoundaryNodes(BoundaryNodes);
        FillListOfBoundaryNodesFast(BoundaryNodes);
        iterBoundaryNodes.resize(BoundaryNodes.size());

        const int nNodes = mrModelPart.NumberOfNodes();
        ModelPart::NodesContainerType::iterator it_begin_nodes = mrModelPart.NodesBegin();

        unsigned int iPosition = -1;
        #pragma omp parallel for
        for (int i = 0; i < nNodes; i++)
        {
            ModelPart::NodesContainerType::iterator it = it_begin_nodes + i;
            int id = it->Id();
            for (unsigned int j = 0; j < BoundaryNodes.size(); ++j)
            {
                if (id == BoundaryNodes[j])
                {
                    iPosition++;
                    iterBoundaryNodes[iPosition] = it;
                }
            }
        }
        KRATOS_CATCH("");
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

        for (unsigned int i=0; i < ELementsOfNodes.size(); ++i)
        {
            ELementsOfNodes[i].resize(N_ELEMENT);
            ELementsOfNodesSize[i] = 0;
            std::fill(ELementsOfNodes[i].begin(), ELementsOfNodes[i].end(), ID_UNDEFINED);
        }

        const unsigned int nElements = mrModelPart.NumberOfElements();
        ModelPart::ElementsContainerType::iterator it_begin_elements = mrModelPart.ElementsBegin();

        for (unsigned int i=0; i < nElements; ++i)
        {
            ModelPart::ElementsContainerType::iterator pElemIt = it_begin_elements + i;
            for (unsigned int iPoint=0; iPoint < pElemIt->GetGeometry().PointsNumber(); ++iPoint)
            {
               int NodeID = pElemIt->GetGeometry()[iPoint].Id();
               int ElementId = pElemIt->Id();

               int index = NodeID-1;
               ELementsOfNodesSize[index]++;
               if (ELementsOfNodesSize[index] > N_ELEMENT-1)
               {
                   ELementsOfNodes[index].push_back(ElementId);
               }
               else
               {
                   ELementsOfNodes[index][ELementsOfNodesSize[index]-1] = ElementId;
               }
            }
        }

        for (unsigned int i=0; i < nElements; ++i)
        {
            ModelPart::ElementsContainerType::iterator pElemIt = it_begin_elements + i;

            int nEdges = pElemIt->GetGeometry().EdgesNumber();
            for (int iEdge = 0; iEdge < nEdges; ++iEdge)
            {
                const unsigned int nPoints = pElemIt->GetGeometry().GenerateEdges()[iEdge].PointsNumber();
                std::vector<int> FaceID(nPoints);
                for (unsigned int iPoint = 0; iPoint < nPoints; ++iPoint)
                {
                    FaceID[iPoint] = pElemIt->GetGeometry().GenerateEdges()[iEdge].GetPoint(iPoint).Id();
                }
                
                if (!IsMoreThanOneElementWithThisEdgeFast(FaceID, ELementsOfNodes, ELementsOfNodesSize))
                {
                    // boundary nodes:
                    for (unsigned int iPoint = 0; iPoint < nPoints; ++iPoint)
                    {
                        std::vector<int>::iterator it = std::find(BoundaryNodes.begin(), BoundaryNodes.end(), FaceID[iPoint]);
                        if (it == BoundaryNodes.end())
                        {
                            BoundaryNodes.push_back(FaceID[iPoint]);
                        }
                    }
                }
            }
        }
        // KRATOS_INFO("FillListOfBoundaryNodesFast:BoundaryNodes") << BoundaryNodes << std::endl;

        if (BoundaryNodes.size()==0)
        {
            KRATOS_INFO("No boundary node is found for interpolate line pressure process") << std::endl;
            KRATOS_ERROR << "No boundary node is found for interpolate line pressure process" << std::endl;
        }

    }

    bool IsMoreThanOneElementWithThisEdgeFast(const std::vector<int> &FaceID,
                                              const std::vector<std::vector<int>> &ELementsOfNodes,
                                              const std::vector<int> &ELementsOfNodesSize)

    {
        const int ID_UNDEFINED = -1;
        int nMaxElements = 0;
        for (unsigned int iPoint = 0; iPoint < FaceID.size(); ++iPoint)
        {
            int NodeID = FaceID[iPoint];
            int index = NodeID-1;
            nMaxElements += ELementsOfNodesSize[index];
        }

        if (nMaxElements > 0)
        {
            std::vector<vector<int>> ElementIDs;
            ElementIDs.resize(FaceID.size());
            for (unsigned int i=0; i<ElementIDs.size(); ++i)
            {
                ElementIDs[i].resize(nMaxElements);
                std::fill(ElementIDs[i].begin(), ElementIDs[i].end(), ID_UNDEFINED);
            }
            

            for (unsigned int iPoint = 0; iPoint < FaceID.size(); ++iPoint)
            {
                int NodeID = FaceID[iPoint];
                int index = NodeID-1;
                for (int i=0; i < ELementsOfNodesSize[index]; ++i)
                {
                    int iElementID = ELementsOfNodes[index][i];
                    ElementIDs[iPoint][i] = iElementID;
                }
            }


            std::vector<int> SharedElementIDs;
            for (unsigned int iPoint = 0; iPoint < FaceID.size(); ++iPoint)
            {
                for (unsigned int i=0; i < ElementIDs[iPoint].size(); ++i)
                {
                    int iElementID = ElementIDs[iPoint][i];
                    bool found = false;
                    if (iElementID !=ID_UNDEFINED)
                    {
                        for (unsigned int iPointInner = 0; iPointInner < FaceID.size(); ++iPointInner)
                        {
                            if (iPointInner != iPoint)
                            {
                                for (unsigned int j = 0; j < ElementIDs[iPointInner].size(); ++j)
                                {
                                    if (ElementIDs[iPointInner][j]==iElementID) found = true;
                                }
                            }
                        }
                    }
                    if (found)
                    {
                        std::vector<int>::iterator it = std::find(SharedElementIDs.begin(), SharedElementIDs.end(), iElementID);
                        if (it == SharedElementIDs.end())
                        {
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
