//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
// Main authors:    Qinfei Ran - Initial development
// Contributor:

// This mapper uses the vector-valued rotational-recovery RBF formulation of
// Ahrem, Beckert and Wendland. The scalar kernel is used on the diagonal of
// the matrix-valued kernel Phi.

#pragma once

// System includes
#include <limits>
#include <array>
#include <string>
#include <unordered_map>

// External includes

// Project includes
#include "mappers/mapper.h"
#include "linear_solvers/linear_solver.h"
#include "custom_searching/interface_communicator.h"
#include "custom_utilities/interface_vector_container.h"
#include "custom_utilities/mapper_local_system.h"

#include "custom_utilities/projection_utilities.h"
#include "custom_utilities/beam_mapper_utilities.h"
#include "utilities/geometrical_projection_utilities.h"

namespace Kratos
{
class KRATOS_API(MAPPING_APPLICATION) BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo : public MapperInterfaceInfo
{
public:
    using InterfaceObjectPointerType = Kratos::shared_ptr<InterfaceObject>;
    using MatrixType = Matrix;
    using VectorType = Vector;
    using MapperInterfaceInfo::GetValue;

    explicit BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo(const double LocalCoordTol = 0.0)
        : mLocalCoordTol(LocalCoordTol)
    {
    }

    explicit BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo(
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
        return Kratos::make_shared<BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo>(mLocalCoordTol);
    }

    MapperInterfaceInfo::Pointer Create(
        const CoordinatesArrayType& rCoordinates,
        const IndexType SourceLocalSystemIndex,
        const IndexType SourceRank) const override
    {
        return Kratos::make_shared<BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo>(
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
        KRATOS_ERROR_IF_NOT(mpInterfaceObject)
            << "BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo does not hold an interface object." << std::endl;
        rValue = *(mpInterfaceObject->pGetBaseGeometry());
    }

    void GetValue(
        MatrixType& rRotationMatrix_L_G,
        VectorType& rProjectionPointValue,
        VectorType& rLinearShapeValues,
        GeometryType& rGeometryValue) const
    {
        KRATOS_ERROR_IF_NOT(mpInterfaceObject)
            << "BeamSplineMapperWithRecoveryOfRotationsInterfaceInfo does not hold an interface object." << std::endl;
        KRATOS_ERROR_IF(mLinearShapeFunctionValues.size() != 2)
            << "Expected two linear shape-function values for a 2-node beam, got "
            << mLinearShapeFunctionValues.size() << "." << std::endl;

        rRotationMatrix_L_G = mRotationMatrix_L_G;
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
    MatrixType mRotationMatrix_L_G;

    void SaveSearchResult(const InterfaceObject& rInterfaceObject, const bool ComputeApproximation);

    void ComputeRotationMatrix();

    friend class Serializer;

    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

class KRATOS_API(MAPPING_APPLICATION) BeamSplineMapperWithRecoveryOfRotationsLocalSystem : public MapperLocalSystem
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

    explicit BeamSplineMapperWithRecoveryOfRotationsLocalSystem(NodePointerType pNode)
        : mpNode(pNode),
          mCachedEvaluationRotationVector(3, 0.0)
    {
    }

    ~BeamSplineMapperWithRecoveryOfRotationsLocalSystem() override = default;

    void CalculateAll(
        MatrixType& rLocalMappingMatrix,
        EquationIdVectorType& rOriginIds,
        EquationIdVectorType& rDestinationIds,
        MapperLocalSystem::PairingStatus& rPairingStatus) const override
    {
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
        return Kratos::make_unique<BeamSplineMapperWithRecoveryOfRotationsLocalSystem>(pNode);
    }

    void PairingInfo(std::ostream& rOStream, const int EchoLevel) const override;

private:
    NodePointerType mpNode;
    VectorType mCachedEvaluationRotationVector;
};

template<class TSparseSpace, class TDenseSpace>
class KRATOS_API(MAPPING_APPLICATION) BeamSplineMapperWithRecoveryOfRotations : public Mapper<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BeamSplineMapperWithRecoveryOfRotations);

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
    using SparseMatrixType = typename TSparseSpace::MatrixType;
    using LinearSolverType = LinearSolver<TSparseSpace, TDenseSpace>;
    using LinearSolverSharedPointerType = Kratos::shared_ptr<LinearSolverType>;
    using ComponentVariableType = Variable<double>;
    using GeometryType = Geometry<Node>;
    using NodePointerType = InterfaceObject::NodePointerType;

    BeamSplineMapperWithRecoveryOfRotations(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination)
        : mrModelPartOrigin(rModelPartOrigin),
          mrModelPartDestination(rModelPartDestination)
    {
    }

    BeamSplineMapperWithRecoveryOfRotations(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination, Parameters JsonParameters);

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
        KRATOS_ERROR << "This function is not supported for the BeamSplineMapperWithRecoveryOfRotations!" << std::endl;
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
        KRATOS_ERROR << "This function is not supported for the BeamSplineMapperWithRecoveryOfRotations!" << std::endl;
    }

    void InverseMap(
        const Variable<array_1d<double, 3>>& rOriginVariable,
        const Variable<array_1d<double, 3>>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        KRATOS_ERROR << "This function is not supported for the BeamSplineMapperWithRecoveryOfRotations!" << std::endl;
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
    static constexpr IndexType Full3DPolynomialBasisSize = 4;
    static constexpr IndexType Full3DVectorPolynomialSize = 3 * Full3DPolynomialBasisSize;
    static constexpr IndexType LineAdaptedVectorPolynomialSize = 9;
    static constexpr IndexType EnrichedLineAdaptedVectorPolynomialSize = 12;
    static constexpr IndexType HighOrderLineAdaptedVectorPolynomialSize = 30;

    enum class FunctionalType
    {
        ValueX,
        ValueY,
        ValueZ,
        CurlX,
        CurlY,
        CurlZ
    };

    struct RecoveryOfRotationsCacheData
    {
        std::string Key;
        std::vector<IndexType> SupportNodeIds;
        std::unordered_map<IndexType, IndexType> SupportNodeIdToLocalIndex;
        std::vector<array_1d<double, 3>> SupportReferenceCoordinates;
        std::string PolynomialBasis;
        array_1d<double, 3> PolynomialReferencePoint = ZeroVector(3);
        array_1d<double, 3> PolynomialTangent = ZeroVector(3);
        MatrixType SystemMatrix;
        mutable SparseMatrixType SparseSystemMatrix;
        mutable SparseMatrixType TransposedSparseSystemMatrix;
    };

    struct RecoveryOfRotationsSourceStateData
    {
        VectorType Coefficients;
    };

    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;
    Parameters mMapperSettings;
    MapperLocalSystemPointerVector mMapperLocalSystems;
    InterfaceCommunicatorPointerType mpInterfaceCommunicator;
    double mLocalCoordTol = 0.0;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerOrigin;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerDestination;
    LinearSolverSharedPointerType mpLinearSolver = nullptr;
    std::unordered_map<std::string, RecoveryOfRotationsCacheData> mBeamChainCache;
    std::unordered_map<IndexType, std::string> mNodeIdToBeamChainKey;
    double mKernelRadius = 1.0;
    double mRegularization = 1.0e-8;
    std::string mKernelType = "cubic";
    std::string mPolynomialBasis = "auto";

    bool UsesLineAdaptedPolynomialBasis(const std::string& rPolynomialBasis) const;

    void ValidateInput();

    void Initialize();

    void InitializeInterfaceCommunicator();

    void CreateLinearSolver();

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

    double KernelValue(const array_1d<double, 3>& rDifference) const;

    array_1d<double, 3> KernelGradient(const array_1d<double, 3>& rDifference) const;

    MatrixType KernelHessian(const array_1d<double, 3>& rDifference) const;

    double ApplyFunctionalsToKernel(
        const FunctionalType RowFunctional,
        const FunctionalType ColumnFunctional,
        const array_1d<double, 3>& rRowCoordinates,
        const array_1d<double, 3>& rColumnCoordinates) const;

    void BuildAijBlock(
        const array_1d<double, 3>& rRowCoordinates,
        const array_1d<double, 3>& rColumnCoordinates,
        MatrixType& rBlock) const;

    void BuildPolynomialMatrix(
        const RecoveryOfRotationsCacheData& rCacheData,
        MatrixType& rPolynomialMatrix) const;

    void EvaluateLineAdaptedPolynomialModes(
        const RecoveryOfRotationsCacheData& rCacheData,
        const array_1d<double, 3>& rCoordinates,
        std::vector<array_1d<double, 3>>& rValues,
        std::vector<array_1d<double, 3>>& rCurls) const;

    std::string GetEffectivePolynomialBasis(const IndexType NumberOfSupportNodes) const;

    IndexType GetVectorPolynomialSize(const std::string& rPolynomialBasis) const;

    VectorType BuildRightHandSide(
        const RecoveryOfRotationsCacheData& rCacheData,
        const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
        const Variable<array_1d<double, 3>>& rOriginVariablesRotations) const;

    MatrixType BuildEvaluationMatrix(
        const RecoveryOfRotationsCacheData& rCacheData,
        const array_1d<double, 3>& rEvaluationCoordinates) const;

    VectorType EvaluateDisplacement(
        const MatrixType& rEvaluationMatrix,
        const VectorType& rCoefficients) const;

    const RecoveryOfRotationsCacheData& GetOrCreateBeamChainCacheData(
        const GeometryType& rBeamGeometry);

    void PrepareBeamChainCacheData();

    std::unordered_map<std::string, RecoveryOfRotationsSourceStateData> BuildAllBeamChainSourceStates(
        const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
        const Variable<array_1d<double, 3>>& rOriginVariablesRotations) const;

    std::string CreateBeamChainKey(
        const std::vector<IndexType>& rSupportNodeIds) const;

    MatrixType BuildSplineSystemMatrix(
        const RecoveryOfRotationsCacheData& rCacheData) const;

    SparseMatrixType BuildSparseMatrix(
        const MatrixType& rDenseMatrix) const;

    SparseMatrixType BuildTransposedSparseMatrix(
        const MatrixType& rDenseMatrix) const;

    void ComputeBeamChainSupport(
        const GeometryType& rBeamGeometry,
        std::vector<IndexType>& rSupportNodeIds,
        std::unordered_map<IndexType, IndexType>& rSupportNodeIdToLocalIndex) const;

    void ComputeSupportReferenceData(
        const std::vector<IndexType>& rSupportNodeIds,
        std::vector<array_1d<double, 3>>& rSupportReferenceCoordinates) const;

    void ComputeLineAdaptedPolynomialData(
        RecoveryOfRotationsCacheData& rCacheData) const;

    VectorType SolveSplineCoefficients(
        SparseMatrixType& rSystemMatrix,
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
        const MatrixType& rRotationMatrix_L_G,
        const VectorType& rLocalVector) const;

    array_1d<double, 3> GetReferenceCoordinates(const Node& rNode) const;

    MatrixType ComputeInitialFrameLocalToGlobal(
        const array_1d<double, 3>& rTangent) const;

    array_1d<double, 3> RotateVectorAroundAxis(
        const array_1d<double, 3>& rVector,
        const array_1d<double, 3>& rAxis,
        const double Angle) const;

    MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const;

    Parameters GetMapperDefaultSettings() const;
};

}  // namespace Kratos
