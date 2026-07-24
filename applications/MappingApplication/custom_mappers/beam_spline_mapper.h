//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
// Main authors:
// Contributor:

// This mapper is created with RBF technique with Cubic Kernel
// Cubic Kernel: phi(r) = |r|^3 in which r is the distance between the point to be projected and the point of the base geometry(beam).

#pragma once

// System includes
#include <limits>
#include <array>
#include <string>
#include <unordered_map>

// External includes

// Project includes
#include "mappers/mapper.h"
#include "custom_searching/interface_communicator.h"
#include "custom_utilities/interface_vector_container.h"
#include "custom_utilities/mapper_local_system.h"

#include "custom_utilities/projection_utilities.h"
#include "custom_utilities/beam_mapper_utilities.h"
#include "utilities/geometrical_projection_utilities.h"

namespace Kratos
{
class KRATOS_API(MAPPING_APPLICATION) BeamSplineMapperInterfaceInfo : public MapperInterfaceInfo
{
public:
    using InterfaceObjectPointerType = Kratos::shared_ptr<InterfaceObject>;
    using MatrixType = Matrix;
    using VectorType = Vector;
    using MapperInterfaceInfo::GetValue;

    explicit BeamSplineMapperInterfaceInfo(const double LocalCoordTol = 0.0)
        : mLocalCoordTol(LocalCoordTol)
    {
    }

    explicit BeamSplineMapperInterfaceInfo(
        const CoordinatesArrayType& rCoordinates,
        const IndexType SourceLocalSystemIndex,
        const IndexType SourceRank,
        const double LocalCoordTol = 0.0)
        : MapperInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank),
          mLocalCoordTol(LocalCoordTol)
    {
    }

    MapperInterfaceInfo::Pointer Create() const override
    {
        return Kratos::make_shared<BeamSplineMapperInterfaceInfo>(mLocalCoordTol);
    }

    MapperInterfaceInfo::Pointer Create(
        const CoordinatesArrayType& rCoordinates,
        const IndexType SourceLocalSystemIndex,
        const IndexType SourceRank) const override
    {
        return Kratos::make_shared<BeamSplineMapperInterfaceInfo>(
            rCoordinates,
            SourceLocalSystemIndex,
            SourceRank,
            mLocalCoordTol);
    }

    InterfaceObject::ConstructionType GetInterfaceObjectType() const override
    {
        return InterfaceObject::ConstructionType::Geometry_Center;
    }

    void ProcessSearchResult(const InterfaceObject& rInterfaceObject) override;

    void ProcessSearchResultForApproximation(const InterfaceObject& rInterfaceObject) override;

    void GetValue(std::vector<int>& rValue, const InfoType ValueType) const override
    {
        rValue = mNodeIds;
    }

    void GetValue(int& rValue, const InfoType ValueType) const override
    {
        rValue = static_cast<int>(mPairingIndex);
    }

    void GetValue(GeometryType& rValue, const InfoType ValueType) const override
    {
        rValue = *(mpInterfaceObject->pGetBaseGeometry());
    }

    void GetValue(
        MatrixType& rRotationMatrixLocalToGlobal,
        VectorType& rProjectionPointValue,
        VectorType& rLinearShapeValues,
        GeometryType& rGeometryValue) const
    {
        rRotationMatrixLocalToGlobal = mRotationMatrixLocalToGlobal;
        rProjectionPointValue(0) = mProjectionOfPoint[0];
        rProjectionPointValue(1) = mProjectionOfPoint[1];
        rProjectionPointValue(2) = mProjectionOfPoint[2];
        rLinearShapeValues(0) = mLinearShapeFunctionValues[0];
        rLinearShapeValues(1) = mLinearShapeFunctionValues[1];
        rGeometryValue = *(mpInterfaceObject->pGetBaseGeometry());
    }

    void ComputeRotationMatrixInterfaceObject()
    {
        ComputeRotationMatrix();
    }

private:
    double mLocalCoordTol = 0.0;

    std::vector<int> mNodeIds;
    std::vector<double> mLinearShapeFunctionValues;
    Point mProjectionOfPoint;
    double mClosestProjectionDistance = std::numeric_limits<double>::max();
    ProjectionUtilities::PairingIndex mPairingIndex = ProjectionUtilities::PairingIndex::Unspecified;
    InterfaceObjectPointerType mpInterfaceObject;
    MatrixType mRotationMatrixLocalToGlobal;

    void SaveSearchResult(const InterfaceObject& rInterfaceObject, const bool ComputeApproximation);

    void ComputeRotationMatrix();

    friend class Serializer;

    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

class KRATOS_API(MAPPING_APPLICATION) BeamSplineMapperLocalSystem : public MapperLocalSystem
{
public:
    using MatrixType = Matrix;
    using VectorType = Vector;
    using GeometryType = Geometry<Node>;
    using NodePointerType = InterfaceObject::NodePointerType;

    struct ProjectionData
    {
        MatrixType RotationMatrixLocalToGlobal{3, 3};
        VectorType ProjectionPoint{3};
        VectorType LinearShapeValues{2};
        GeometryType BeamGeometry;
        NodePointerType pNode;
    };

    explicit BeamSplineMapperLocalSystem(NodePointerType pNode)
        : mpNode(pNode),
          mCachedEvaluationRotationVector(3, 0.0)
    {
    }

    ~BeamSplineMapperLocalSystem() override = default;

    void CalculateAll(
        MatrixType& rLocalMappingMatrix,
        EquationIdVectorType& rOriginIds,
        EquationIdVectorType& rDestinationIds,
        MapperLocalSystem::PairingStatus& rPairingStatus) const override
    {
        KRATOS_WARNING("BeamSplineMapperLocalSystem")
            << "CalculateAll() was called, but is not implemented for BeamSplineMapperLocalSystem." << std::endl;
    }

    CoordinatesArrayType& Coordinates() const override
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not initialized!" << std::endl;
        return mpNode->Coordinates();
    }

    ProjectionData CalculateProjectionData();

    void SaveEvaluationRotationVector(const VectorType& rEvaluationRotationVector)
    {
        KRATOS_ERROR_IF(rEvaluationRotationVector.size() != 3)
            << "Expected a 3D evaluation rotation vector." << std::endl;
        mCachedEvaluationRotationVector = rEvaluationRotationVector;
    }

    void GetEvaluationRotationVector(VectorType& rEvaluationRotationVector) const
    {
        rEvaluationRotationVector.resize(3, false);
        for (IndexType i = 0; i < 3; ++i) {
            rEvaluationRotationVector(i) = mCachedEvaluationRotationVector(i);
        }
    }

    MapperLocalSystemUniquePointer Create(NodePointerType pNode) const override
    {
        return Kratos::make_unique<BeamSplineMapperLocalSystem>(pNode);
    }

    void PairingInfo(std::ostream& rOStream, const int EchoLevel) const override;

private:
    NodePointerType mpNode;
    VectorType mCachedEvaluationRotationVector;
};

template<class TSparseSpace, class TDenseSpace>
class KRATOS_API(MAPPING_APPLICATION) BeamSplineMapper : public Mapper<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BeamSplineMapper);

    using BaseType = Mapper<TSparseSpace, TDenseSpace>;
    using InterfaceCommunicatorPointerType = Kratos::unique_ptr<InterfaceCommunicator>;
    using MapperInterfaceInfoUniquePointerType = typename InterfaceCommunicator::MapperInterfaceInfoUniquePointerType;
    using InterfaceVectorContainerType = InterfaceVectorContainer<TSparseSpace, TDenseSpace>;
    using InterfaceVectorContainerPointerType = Kratos::unique_ptr<InterfaceVectorContainerType>;
    using MapperLocalSystemPointer = Kratos::unique_ptr<MapperLocalSystem>;
    using MapperLocalSystemPointerVector = std::vector<MapperLocalSystemPointer>;
    using MapperUniquePointerType = typename BaseType::MapperUniquePointerType;
    using MatrixType = typename TDenseSpace::MatrixType;
    using VectorType = typename TDenseSpace::VectorType;
    using ComponentVariableType = Variable<double>;
    using GeometryType = Geometry<Node>;
    using NodePointerType = InterfaceObject::NodePointerType;

    BeamSplineMapper(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination)
        : mrModelPartOrigin(rModelPartOrigin),
          mrModelPartDestination(rModelPartDestination)
    {
    }

    BeamSplineMapper(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination, Parameters JsonParameters);

    void UpdateInterface(Kratos::Flags MappingOptions, double SearchRadius) override;

    void Map(
        const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
        const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
        const Variable<array_1d<double, 3>>& rDestinationVariableDisplacement,
        Kratos::Flags MappingOptions);

    void Map(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        KRATOS_ERROR << "This function is not supported for the BeamSplineMapper!" << std::endl;
    }

    void Map(
        const Variable<array_1d<double, 3>>& rOriginVariable,
        const Variable<array_1d<double, 3>>& rDestinationVariable,
        Kratos::Flags MappingOptions) override;

    void InverseMap(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        KRATOS_ERROR << "This function is not supported for the BeamSplineMapper!" << std::endl;
    }

    void InverseMap(
        const Variable<array_1d<double, 3>>& rOriginVariable,
        const Variable<array_1d<double, 3>>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        KRATOS_ERROR << "This function is not supported for the BeamSplineMapper!" << std::endl;
    }

    void InverseMap(
        const Variable<array_1d<double, 3>>& rOriginVariablesForces,
        const Variable<array_1d<double, 3>>& rOriginVariablesMoments,
        const Variable<array_1d<double, 3>>& rDestinationVariableForces,
        Kratos::Flags MappingOptions);

    MapperUniquePointerType Clone(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        Parameters JsonParameters) const override;

    std::string Info() const override;

    void PrintInfo(std::ostream& rOStream) const override;

    void PrintData(std::ostream& rOStream) const override;

private:
    using ComponentArrayType = std::array<std::vector<double>, 3>;

    struct BeamChainCacheData
    {
        std::string Key;
        std::vector<IndexType> SupportNodeIds;
        std::unordered_map<IndexType, IndexType> SupportNodeIdToLocalIndex;
        std::vector<double> LocalXCoordinates;
        double CoordinateHalfLength = 1.0;
        std::vector<MatrixType> SupportFramesLocalToGlobal;
        std::vector<array_1d<double, 3>> SupportReferenceCoordinates;
        MatrixType SplineSystemMatrix;
        MatrixType TransposedSplineSystemMatrix;
    };

    struct BeamChainSourceStateData
    {
        ComponentArrayType LocalDisplacements;
        ComponentArrayType LocalRotations;
        VectorType SplineCoefficientsY;
        VectorType SplineCoefficientsZ;
    };

    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;
    Parameters mMapperSettings;
    MapperLocalSystemPointerVector mMapperLocalSystems;
    InterfaceCommunicatorPointerType mpInterfaceCommunicator;
    double mLocalCoordTol = 0.0;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerOrigin;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerDestination;
    std::unordered_map<std::string, BeamChainCacheData> mBeamChainCache;
    std::unordered_map<IndexType, std::string> mNodeIdToBeamChainKey;

    void ValidateInput();

    void Initialize();

    void InitializeInterfaceCommunicator();

    void InitializeInterface(Kratos::Flags MappingOptions = Kratos::Flags());

    void BuildProblem(Kratos::Flags MappingOptions = Kratos::Flags());

    void CreateMapperLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems);

    void InitializeInformationBeams(
        const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
        const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
        const Variable<array_1d<double, 3>>& rDestinationVariableDisplacement);

    void InitializeInformationBeamsInverse(
        const Variable<array_1d<double, 3>>& rOriginVariablesForces,
        const Variable<array_1d<double, 3>>& rOriginVariablesMoments,
        const Variable<array_1d<double, 3>>& rDestinationVariableForces,
        const Kratos::Flags& rMappingOptions);

    void InitializeOriginForcesAndMoments(
        const Variable<array_1d<double, 3>>& rOriginVariablesForces,
        const Variable<array_1d<double, 3>>& rOriginVariablesMoments);

    double EvaluateKernelValue(const double Distance) const;

    double EvaluateKernelFirstDerivative(const double Distance) const;

    double EvaluateKernelSecondDerivative(const double Distance) const;

    void BuildSplineSystemMatrices(
        const std::vector<double>& rSourceCoordinates,
        MatrixType& rMss,
        MatrixType& rMssFirstDerivative,
        MatrixType& rMssSecondDerivative,
        MatrixType& rPs,
        MatrixType& rPsFirstDerivative) const;

    VectorType BuildEvaluationRow(
        const std::vector<double>& rSourceCoordinates,
        const double ProjectionCoordinate) const;

    VectorType BuildEvaluationDerivativeRow(
        const std::vector<double>& rSourceCoordinates,
        const double ProjectionCoordinate,
        const double CoordinateHalfLength) const;

    VectorType BuildRightHandSide(
        const std::vector<double>& rDisplacements,
        const std::vector<double>& rRotations,
        const double CoordinateHalfLength) const;

    void BuildLocalSourceData(
        const BeamChainCacheData& rBeamChainCacheData,
        const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
        const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
        ComponentArrayType& rLocalDisplacements,
        ComponentArrayType& rLocalRotations) const;

    const BeamChainCacheData& GetOrCreateBeamChainCacheData(
        const GeometryType& rBeamGeometry);

    void PrepareBeamChainCacheData();

    std::unordered_map<std::string, BeamChainSourceStateData> BuildAllBeamChainSourceStates(
        const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
        const Variable<array_1d<double, 3>>& rOriginVariablesRotations) const;

    std::string CreateBeamChainKey(
        const std::vector<IndexType>& rSupportNodeIds) const;

    MatrixType BuildSplineSystemMatrix(
        const std::vector<double>& rSourceCoordinates) const;

    MatrixType BuildTransposedMatrix(
        const MatrixType& rMatrix) const;

    void ComputeBeamChainSupport(
        const GeometryType& rBeamGeometry,
        std::vector<IndexType>& rSupportNodeIds,
        std::unordered_map<IndexType, IndexType>& rSupportNodeIdToLocalIndex) const;

    void ComputeSupportReferenceData(
        const std::vector<IndexType>& rSupportNodeIds,
        std::vector<double>& rLocalXCoordinates,
        double& rCoordinateHalfLength,
        std::vector<MatrixType>& rSupportFramesLocalToGlobal,
        std::vector<array_1d<double, 3>>& rSupportReferenceCoordinates) const;

    MatrixType BuildEvaluationFrameLocalToGlobal(
        const std::vector<IndexType>& rSupportNodeIds,
        const std::unordered_map<IndexType, IndexType>& rSupportNodeIdToLocalIndex,
        const std::vector<MatrixType>& rSupportFramesLocalToGlobal,
        const std::vector<array_1d<double, 3>>& rSupportReferenceCoordinates,
        const GeometryType& rBeamGeometry,
        const VectorType& rLinearShapeValues) const;

    VectorType SolveSplineCoefficients(
        const MatrixType& rSplineSystemMatrix,
        const VectorType& rRightHandSide) const;

    void MapDisplacements(
        const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
        const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
        const Variable<array_1d<double, 3>>& rDestinationVariableDisplacement,
        Kratos::Flags MappingOptions);

    VectorType TransformGlobalToLocal(
        const MatrixType& rRotationMatrixGlobalToLocal,
        const array_1d<double, 3>& rReferencePoint,
        const array_1d<double, 3>& rCoordinates) const;

    VectorType TransformVectorToLocal(
        const MatrixType& rRotationMatrixGlobalToLocal,
        const array_1d<double, 3>& rVector) const;

    VectorType TransformVectorToGlobal(
        const MatrixType& rRotationMatrixLocalToGlobal,
        const VectorType& rLocalVector) const;

    array_1d<double, 3> GetReferenceCoordinates(const Node& rNode) const;

    MatrixType ComputeInitialFrameLocalToGlobal(
        const array_1d<double, 3>& rTangent) const;

    array_1d<double, 3> RotateVectorAroundAxis(
        const array_1d<double, 3>& rVector,
        const array_1d<double, 3>& rAxis,
        const double Angle) const;

    VectorType EvaluatePointDisplacementLocal(
        const VectorType& rEvaluationRow,
        const VectorType& rEvaluationDerivativeRow,
        const VectorType& rSplineCoefficientsY,
        const VectorType& rSplineCoefficientsZ,
        const double AxialDisplacement,
        const double TorsionalRotation,
        const VectorType& rOffsetVectorLocal) const;

    MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const;

    Parameters GetMapperDefaultSettings() const;
};

}  // namespace Kratos
