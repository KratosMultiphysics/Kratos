//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan I. Camarotti
//

#pragma once

// System includes
#include <variant>

// External includes

// Project includes
#include "mappers/mapper_define.h"
#include "interpolative_mapper_base.h"
#include "custom_utilities/mapper_backend.h"
#include "custom_utilities/mapper_local_system.h"
#include "custom_utilities/radial_basis_functions_utilities.h"
#include "utilities/rbf_shape_functions_utility.h"
#include "linear_solvers/linear_solver.h"
#include "factories/linear_solver_factory.h"
#include "custom_searching/interface_communicator.h"


namespace Kratos
{
///@name Kratos Classes
///@{

class RBFSupportAccumulator {
    public:

        ///@name Type Definitions
        ///@{
        using NodeType = Node;
        using GeometryType = Geometry<NodeType>;
        using CoordinatesArrayType = GeometryType::CoordinatesArrayType;

        struct Candidate {
            double mDistance{};
            IndexType mInterfaceEquationId;
            CoordinatesArrayType  mCoordinates = ZeroVector(3);

            friend class Serializer;
            void save(Serializer& rSerializer) const {
                rSerializer.save("Distance", mDistance);
                rSerializer.save("Coordinates", mCoordinates);
            }
            void load(Serializer& rSerializer) {
                rSerializer.load("Distance", mDistance);
                rSerializer.load("Coordinates", mCoordinates);
            }
        };

        explicit RBFSupportAccumulator(IndexType required_points)
                : mRequiredRBFSupportPoints(required_points) {}

        void AddNodeCandidate(IndexType eq_id, double distance, const CoordinatesArrayType& rCoordinates)
        {
            Candidate node_candidate;
            node_candidate.mDistance = distance;
            node_candidate.mInterfaceEquationId = eq_id;
            node_candidate.mCoordinates = rCoordinates;

            mCandidates.push_back(std::move(node_candidate));
            mNotOrdered = true;
        }

        void AddGeometryCandidate(IndexType interface_equation_id, double distance, const CoordinatesArrayType& rCoordinates)
        {
            Candidate geometry_candidate;
            geometry_candidate.mDistance = distance;
            geometry_candidate.mInterfaceEquationId = interface_equation_id;
            geometry_candidate.mCoordinates = rCoordinates;

            mCandidates.push_back(std::move(geometry_candidate)); 
            mNotOrdered = true;
        }

        void ReorderSupportPoints()
        {
            if (!mNotOrdered) return;

            std::stable_sort(mCandidates.begin(), mCandidates.end(),
                [](const Candidate& a, const Candidate& b){ return a.mDistance < b.mDistance; });

            mNumNeighbors = mCandidates.size();

            mNotOrdered = false;
        }

        IndexType NumNeighbors() const { return mNumNeighbors; }

        const std::vector<Candidate>& Candidates() const { return mCandidates; }

        IndexType GetRequiredRBFPoints() const {return mRequiredRBFSupportPoints;}

        bool HasEnough() const { return mNumNeighbors >= mRequiredRBFSupportPoints; }

        bool HasSome() const { return mNumNeighbors > 0; }

        // Get a matrix with the coordinates of the candidates
        Matrix GetCandidateCoordinatesMatrix() const
        {
            const IndexType n_points = mCandidates.size(); 
            Matrix coords(n_points, 3, 0.0);

            for (IndexType i = 0; i < n_points; ++i) {
                const CoordinatesArrayType& p = mCandidates[i].mCoordinates;
                coords(i, 0) = p[0];
                coords(i, 1) = p[1];
                coords(i, 2) = p[2];
            }

            return coords;
        }

    private:
        IndexType mRequiredRBFSupportPoints{0};
        std::vector<Candidate> mCandidates;
        IndexType           mNumNeighbors{0};
        bool                mNotOrdered{false};

        friend class Serializer;
        void save(Serializer& rSerializer) const {
            rSerializer.save("RequiredRBFSupportPoints", mRequiredRBFSupportPoints);
            rSerializer.save("Candidates", mCandidates);
            rSerializer.save("NeighborsNumber", mNumNeighbors);
            rSerializer.save("Dirty", mNotOrdered);
        }

        void load(Serializer& rSerializer) {
            rSerializer.load("RequiredRBFSupportPoints", mRequiredRBFSupportPoints);
            rSerializer.load("Candidates", mCandidates);
            rSerializer.load("NeighborsNumber", mNumNeighbors);
            rSerializer.load("Dirty", mNotOrdered);
        }
};   

class KRATOS_API(MAPPING_APPLICATION) RadialBasisFunctionsMapperInterfaceInfoIGA : public MapperInterfaceInfo
{
public:

    /// Constructor with support points
    RadialBasisFunctionsMapperInterfaceInfoIGA(const IndexType RequiredRBFSupportPoints)
        : mRequiredRBFSupportPoints(RequiredRBFSupportPoints),
        mRBFSupportAccumulator(RequiredRBFSupportPoints)
    {}

    explicit RadialBasisFunctionsMapperInterfaceInfoIGA(const CoordinatesArrayType& rCoordinates,
                                             const IndexType SourceLocalSystemIndex,
                                             const IndexType SourceRank,
                                             const IndexType RequiredRBFSupportPoints)
        : MapperInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank),  
        mRequiredRBFSupportPoints(RequiredRBFSupportPoints),
        mRBFSupportAccumulator(RequiredRBFSupportPoints) {}

    MapperInterfaceInfo::Pointer Create() const override
    {
        return Kratos::make_shared<RadialBasisFunctionsMapperInterfaceInfoIGA>(mRequiredRBFSupportPoints);
    }

    MapperInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                        const IndexType SourceLocalSystemIndex,
                                        const IndexType SourceRank) const override
    {
        return Kratos::make_shared<RadialBasisFunctionsMapperInterfaceInfoIGA>(
            rCoordinates,
            SourceLocalSystemIndex,
            SourceRank,
            mRequiredRBFSupportPoints);
    }

    InterfaceObject::ConstructionType GetInterfaceObjectType() const override
    {
        return InterfaceObject::ConstructionType::Geometry_Center;
    }

    void ProcessSearchResult(const InterfaceObject& rInterfaceObject) override;

    const RBFSupportAccumulator& GetRBFSupportAccumulator() const override
    {
        return mRBFSupportAccumulator;
    }

    void GetValue(IndexType& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mRequiredRBFSupportPoints;
    }

private:
    IndexType mRequiredRBFSupportPoints = 0;
    RBFSupportAccumulator mRBFSupportAccumulator;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.save("RequiredRBFSupportPoints", mRequiredRBFSupportPoints);
        rSerializer.save("RBFSupportAccumulator", mRBFSupportAccumulator);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.load("RequiredRBFSupportPoints", mRequiredRBFSupportPoints);
        rSerializer.load("RBFSupportAccumulator", mRBFSupportAccumulator);
    }
};

class KRATOS_API(MAPPING_APPLICATION) RadialBasisFunctionsMapperInterfaceInfoFEM : public MapperInterfaceInfo
{
public:

    RadialBasisFunctionsMapperInterfaceInfoFEM(const IndexType RequiredRBFSupportPoints)
        : mRequiredRBFSupportPoints(RequiredRBFSupportPoints),
          mRBFSupportAccumulator(RequiredRBFSupportPoints) {}

    explicit RadialBasisFunctionsMapperInterfaceInfoFEM(const CoordinatesArrayType& rCoordinates,
                                             const IndexType SourceLocalSystemIndex,
                                             const IndexType SourceRank,
                                             const IndexType RequiredRBFSupportPoints)
        : MapperInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank),  
          mRequiredRBFSupportPoints(RequiredRBFSupportPoints),
          mRBFSupportAccumulator(RequiredRBFSupportPoints) {}

    MapperInterfaceInfo::Pointer Create() const override
    {
        return Kratos::make_shared<RadialBasisFunctionsMapperInterfaceInfoFEM>(mRequiredRBFSupportPoints);
    }

    MapperInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                        const IndexType SourceLocalSystemIndex,
                                        const IndexType SourceRank) const override
    {
        return Kratos::make_shared<RadialBasisFunctionsMapperInterfaceInfoFEM>(
            rCoordinates,
            SourceLocalSystemIndex,
            SourceRank,
            mRequiredRBFSupportPoints);
    }

    InterfaceObject::ConstructionType GetInterfaceObjectType() const override
    {
        return InterfaceObject::ConstructionType::Node_Coords;
    }

    void ProcessSearchResult(const InterfaceObject& rInterfaceObject) override;

    const RBFSupportAccumulator& GetRBFSupportAccumulator() const override
    {
        return mRBFSupportAccumulator;
    }

    void GetValue(IndexType& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mRequiredRBFSupportPoints;
    }

private:
    IndexType mRequiredRBFSupportPoints = 0;
    RBFSupportAccumulator mRBFSupportAccumulator;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.save("RequiredRBFSupportPoints", mRequiredRBFSupportPoints);
        rSerializer.save("RBFSupportAccumulator", mRBFSupportAccumulator);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.load("RequiredRBFSupportPoints", mRequiredRBFSupportPoints);
        rSerializer.load("RBFSupportAccumulator", mRBFSupportAccumulator);
    }
};

class RadialBasisFunctionMapperLocalSystem : public MapperLocalSystem
{
public:
    // Default constructor
    RadialBasisFunctionMapperLocalSystem(): mpNode(nullptr), mpGeometry(nullptr), mRBFTypeString(""), mRBFTypeEnum(), mpPolynomialEquationIds(nullptr) {}

    // Construct from a node and a geometry (dummy constructor needed)
    explicit RadialBasisFunctionMapperLocalSystem(NodePointerType pNode, GeometryPointerType pGeometry, bool BuildOriginInterpolationMatrix, std::string RBFType, IndexType PolynomialDegree,
                                                  const std::vector<IndexType>* pPolynomialEquationIds)
            : mpNode(pNode), 
              mpGeometry(pGeometry),
              mBuildOriginInterpolationMatrix(BuildOriginInterpolationMatrix),
              mRBFTypeString(RBFType),
              mRBFTypeEnum(ParseRBFType(RBFType)),
              mPolynomialDegree(PolynomialDegree), 
              mNumberOfPolynomialTerms((PolynomialDegree + 3) * (PolynomialDegree + 2) * (PolynomialDegree + 1) / 6),
              mpPolynomialEquationIds(pPolynomialEquationIds){}

    // Construct from a node
    explicit RadialBasisFunctionMapperLocalSystem(NodePointerType pNode, bool BuildOriginInterpolationMatrix, std::string RBFType, IndexType PolynomialDegree, const std::vector<IndexType>* pPolynomialEquationIds): 
             mpNode(std::move(pNode)), 
             mpGeometry(nullptr),
             mBuildOriginInterpolationMatrix(BuildOriginInterpolationMatrix),
             mRBFTypeString(RBFType),
             mRBFTypeEnum(ParseRBFType(RBFType)),
             mPolynomialDegree(PolynomialDegree),
             mNumberOfPolynomialTerms((PolynomialDegree + 3) * (PolynomialDegree + 2) * (PolynomialDegree + 1) / 6),
             mpPolynomialEquationIds(pPolynomialEquationIds){}

    // Construct from a geometry (e.g. an integration-point “geometry” in IGA)
    explicit RadialBasisFunctionMapperLocalSystem(GeometryPointerType pGeometry, bool BuildOriginInterpolationMatrix, std::string RBFType, IndexType PolynomialDegree, const std::vector<IndexType>* pPolynomialEquationIds): 
             mpNode(nullptr), 
             mpGeometry(std::move(pGeometry)),
             mBuildOriginInterpolationMatrix(BuildOriginInterpolationMatrix),
             mRBFTypeString(RBFType),
             mRBFTypeEnum(ParseRBFType(RBFType)),
             mPolynomialDegree(PolynomialDegree),
             mNumberOfPolynomialTerms((PolynomialDegree + 3) * (PolynomialDegree + 2) * (PolynomialDegree + 1) / 6),
             mpPolynomialEquationIds(pPolynomialEquationIds){}

    void CalculateAll(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const override;

    CoordinatesArrayType& Coordinates() const override
    {
        KRATOS_DEBUG_ERROR_IF(!mpNode && !mpGeometry) << "LocalSystem not initialized" << std::endl;
        if (mpNode != nullptr){
            mCoordinates = mpNode->Coordinates();
            return mCoordinates;
        } else if (mpGeometry != nullptr) {
            mCoordinates = mpGeometry->Center();
            return mCoordinates;
        } else{
            KRATOS_ERROR << "RadialBasisFunctionMapperLocalSystem::Coordinates(): "
                 << "neither node nor geometry is set!" << std::endl;

            static CoordinatesArrayType dummy = ZeroVector(3);
            return dummy;
        }
    }

    MapperLocalSystemUniquePointer Create(NodePointerType pNode) const override
    {
       KRATOS_DEBUG_ERROR_IF(!pNode) << "Create(Node) received nullptr" << std::endl;
       return Kratos::make_unique<RadialBasisFunctionMapperLocalSystem>(pNode, mBuildOriginInterpolationMatrix, mRBFTypeString, mPolynomialDegree, mpPolynomialEquationIds);
    }

    MapperLocalSystemUniquePointer Create(GeometryPointerType pGeometry) const override
    {
        KRATOS_DEBUG_ERROR_IF(!pGeometry) << "Create(Geometry) received nullptr" << std::endl;
        return Kratos::make_unique<RadialBasisFunctionMapperLocalSystem>(pGeometry, mBuildOriginInterpolationMatrix, mRBFTypeString, mPolynomialDegree, mpPolynomialEquationIds);
    }

    void AddInterfaceInfo(MapperInterfaceInfoPointerType pInterfaceInfo) override 
    {
        mInterfaceInfos.clear();
        mInterfaceInfos.push_back(pInterfaceInfo);
    }

    /// Turn back information as a string.
    void PairingInfo(std::ostream& rOStream, const int EchoLevel) const override {KRATOS_ERROR << "Not implemented!"<<std::endl;}

private:
    NodePointerType     mpNode{nullptr};
    GeometryPointerType mpGeometry{nullptr};

    bool mBuildOriginInterpolationMatrix{};

    using RBFVariant = std::variant<
        RadialBasisFunctionsUtilities::InverseMultiquadric,
        RadialBasisFunctionsUtilities::Multiquadric,
        RadialBasisFunctionsUtilities::Gaussian,
        RadialBasisFunctionsUtilities::ThinPlateSpline,
        RadialBasisFunctionsUtilities::WendlandC2
    >;
    mutable RBFVariant mRBFKernel;
    mutable bool mIsRBFInitialized = false;
    mutable CoordinatesArrayType mCoordinates;

    std::string mRBFTypeString;
    RadialBasisFunctionsUtilities::RBFType mRBFTypeEnum;
    IndexType mPolynomialDegree = 0;
    IndexType mNumberOfPolynomialTerms;
    const std::vector<IndexType>* mpPolynomialEquationIds; // reference to the polynomial ids

    std::vector<double> EvaluatePolynomialBasis(const CoordinatesArrayType& rCoords, unsigned int degree) const;

    RadialBasisFunctionsUtilities::RBFType ParseRBFType(const std::string& rName) const
    {
        static const std::unordered_map<std::string, RadialBasisFunctionsUtilities::RBFType> rbf_type_map = {
            {"inverse_multiquadric", RadialBasisFunctionsUtilities::RBFType::InverseMultiquadric},
            {"multiquadric",         RadialBasisFunctionsUtilities::RBFType::Multiquadric},
            {"gaussian",             RadialBasisFunctionsUtilities::RBFType::Gaussian},
            {"thin_plate_spline",    RadialBasisFunctionsUtilities::RBFType::ThinPlateSpline},
            {"wendland_c2",          RadialBasisFunctionsUtilities::RBFType::WendlandC2}
        };

        const auto it = rbf_type_map.find(rName);
        if (it == rbf_type_map.end()) {
            KRATOS_ERROR << "Unrecognized RBF type: " << rName << std::endl;
        }
        return it->second;
    }

    void InitializeRBFKernel(const Matrix& coords) const
    {
        if (mIsRBFInitialized) {return;}

        switch (mRBFTypeEnum)
        {
            case RadialBasisFunctionsUtilities::RBFType::InverseMultiquadric:
            {
                const double h = RBFShapeFunctionsUtility::CalculateInverseMultiquadricShapeParameter(coords);
                mRBFKernel = RadialBasisFunctionsUtilities::InverseMultiquadric{h};
                break;
            }
            case RadialBasisFunctionsUtilities::RBFType::Multiquadric:
            {
                const double h = RBFShapeFunctionsUtility::CalculateInverseMultiquadricShapeParameter(coords);
                mRBFKernel = RadialBasisFunctionsUtilities::Multiquadric{h};
                break;
            }
            case RadialBasisFunctionsUtilities::RBFType::Gaussian:
            {
                const double h = RBFShapeFunctionsUtility::CalculateInverseMultiquadricShapeParameter(coords);
                mRBFKernel = RadialBasisFunctionsUtilities::Gaussian{h};
                break;
            }
            case RadialBasisFunctionsUtilities::RBFType::ThinPlateSpline:
            {
                mRBFKernel = RadialBasisFunctionsUtilities::ThinPlateSpline{};
                break;
            }
            case RadialBasisFunctionsUtilities::RBFType::WendlandC2:
            {
                const double h = RadialBasisFunctionsUtilities::CalculateWendlandC2SupportRadius(coords);
                mRBFKernel = RadialBasisFunctionsUtilities::WendlandC2{h};
                break;
            }
            default:
                KRATOS_ERROR << "Invalid RBF type. Available types are: inverse_multiquadric, multiquadric, gaussian, thin-plate spline and wendland C2" << std::endl;
        }
        mIsRBFInitialized = true;
    }
};

/// Radial Basis Functions (RBF) Mapper
/** This class implements the Radial Basis Functions mapping technique.
 * Each node on the destination side uses its N closest support points
 * on the origin side to build the RBF interpolation system,
 * optionally enriched with polynomial terms.
 *
 * The mapping matrix can be precomputed once and reused in every coupling step,
 * or a linear system can be solved at each step. For IGA-based mappings,
 * only the precomputed mapping matrix option is supported.
*/
template<class TSparseSpace, class TDenseSpace>
class RadialBasisFunctionMapper : public Mapper<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of RadialBasisFunctionMapper
    KRATOS_CLASS_POINTER_DEFINITION(RadialBasisFunctionMapper);

    using BaseType = Mapper<TSparseSpace, TDenseSpace>;

    using MapperLocalSystemPointer       = Kratos::unique_ptr<MapperLocalSystem>;
    using MapperLocalSystemPointerVector = std::vector<MapperLocalSystemPointer>;

    using InterfaceVectorContainerType        = InterfaceVectorContainer<TSparseSpace, TDenseSpace>;
    using InterfaceVectorContainerPointerType = Kratos::unique_ptr<InterfaceVectorContainerType>;

    using IndexType = std::size_t;

    using MapperUniquePointerType  = typename BaseType::MapperUniquePointerType;
    using MappingMatrixType        = typename BaseType::TMappingMatrixType;
    using MappingMatrixUniquePointerType = Kratos::unique_ptr<MappingMatrixType>;

    using LinearSolverType               = LinearSolver<TSparseSpace, TDenseSpace>;
    using LinearSolverSharedPointerType  = Kratos::shared_ptr<LinearSolverType>;

    using TSystemVectorType           = typename TSparseSpace::VectorType;
    using TSystemVectorUniquePointerType = Kratos::unique_ptr<TSystemVectorType>;

    using SparseMatrixType = typename TSparseSpace::MatrixType;
    using DenseMatrixType  = typename TDenseSpace::MatrixType;

    using InterfaceCommunicatorPointerType = Kratos::unique_ptr<InterfaceCommunicator>;
    using MapperInterfaceInfoUniquePointerType = typename InterfaceCommunicator::MapperInterfaceInfoUniquePointerType;

    using MappingSparseSpaceType = typename MapperDefinitions::SparseSpaceType;
    using DenseSpaceType = typename MapperDefinitions::DenseSpaceType;
    using MappingMatrixUtilitiesType = MappingMatrixUtilities<MappingSparseSpaceType, DenseSpaceType>;


    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    RadialBasisFunctionMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination)
                        : mrModelPartOrigin(rModelPartOrigin),
                          mrModelPartDestination(rModelPartDestination){}


    RadialBasisFunctionMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination,
                         Parameters JsonParameters);

    /// Destructor.
    ~RadialBasisFunctionMapper() override = default;

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(
        Kratos::Flags MappingOptions,
        double SearchRadius) override
    {
        KRATOS_ERROR << "UpdateInterface() is not implemented for this mapper" << std::endl;
    }

    void Map(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MappingOptions.Reset(MapperFlags::USE_TRANSPOSE);
            MappingOptions.Set(MapperFlags::INTERNAL_USE_TRANSPOSE, true);
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
        else if (MappingOptions.Is(MapperFlags::INTERNAL_USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
        }
    }

    void Map(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MappingOptions.Reset(MapperFlags::USE_TRANSPOSE);
            MappingOptions.Set(MapperFlags::INTERNAL_USE_TRANSPOSE, true);
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
        else if (MappingOptions.Is(MapperFlags::INTERNAL_USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            MapInternal(rOriginVariable, rDestinationVariable, MappingOptions);
        }
    }

    void InverseMap(
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
    }

    void InverseMap(
        const Variable< array_1d<double, 3> >& rOriginVariable,
        const Variable< array_1d<double, 3> >& rDestinationVariable,
        Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) {
            MapInternalTranspose(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else {
            GetInverseMapper()->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
    }

    ///@}
    ///@name Access
    ///@{

    MappingMatrixType& GetMappingMatrix() override
    {
        if (mMapperSettings["precompute_mapping_matrix"].GetBool()) return *(mpMappingMatrix.get());
        else KRATOS_ERROR << "'precompute_mapping_matrix' must be 'true' in your parameters to retrieve the computed mapping matrix!" << std::endl;
    }

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) const override
    {
        return Kratos::make_unique<RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>>(
            rModelPartOrigin,
            rModelPartDestination,
            JsonParameters);
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "RadialBasisFunctionMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RadialBasisFunctionMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

    // Get values
    // Always returns the true origin/destination regardless of 'destination_is_slave'
    ModelPart& GetInterfaceModelPartOrigin() override
    {
        return *mpCouplingInterfaceOrigin;
    }

    ModelPart& GetInterfaceModelPartDestination() override
    {
        return *mpCouplingInterfaceDestination;
    }

protected:
    ///@name Protected static Member Variables
    Parameters mLocalMapperSettings;

    /**
     * @brief This function origin model part
     * @return The origin model part
     */
    ModelPart& GetOriginModelPart()
    {
        return mrModelPartOrigin;
    }

    /**
     * @brief This function destination model part
     * @return The destination model part
     */
    ModelPart& GetDestinationModelPart()
    {
        return mrModelPartDestination;
    }

    /**
     * @brief This function return the interface vector container of the origin model part
     * @return The origin model part
    */
    InterfaceVectorContainerType* pGetInterfaceVectorContainerOrigin()
    {
        return mpInterfaceVectorContainerOrigin.get();
    }

    /**
     * @brief This function return the interface vector container of the destination model part
     * @return The destination model part
    */
    InterfaceVectorContainerType* pGetInterfaceVectorContainerDestination()
    {
        return mpInterfaceVectorContainerDestination.get();
    }

    void ValidateInput()
    {
        // backward compatibility
        if (mMapperSettings.Has("search_radius")) {
            KRATOS_WARNING("Mapper") << "DEPRECATION-WARNING: \"search_radius\" should be specified under \"search_settings\"!" << std::endl;
            const double search_radius = mMapperSettings["search_radius"].GetDouble();

            if (mMapperSettings.Has("search_settings")) {
                KRATOS_ERROR_IF(mMapperSettings["search_settings"].Has("search_radius")) << "\"search_radius\" specified twice, please only specify it in \"search_settings\"!" << std::endl;
            } else {
                mMapperSettings.AddValue("search_settings", Parameters());
            }

            mMapperSettings["search_settings"].AddEmptyValue("search_radius").SetDouble(search_radius);
            mMapperSettings.RemoveValue("search_radius");
        }

        if (mMapperSettings.Has("search_iterations")) {
            KRATOS_WARNING("Mapper") << "DEPRECATION-WARNING: \"search_iterations\" should be specified as \"max_num_search_iterations\" under \"search_settings\"!" << std::endl;
            const int search_iterations = mMapperSettings["search_iterations"].GetInt();

            if (mMapperSettings.Has("search_settings")) {
                KRATOS_ERROR_IF(mMapperSettings["search_settings"].Has("max_num_search_iterations")) << "\"search_iterations\" specified twice, please only specify it in \"search_settings\" (as \"max_num_search_iterations\")!" << std::endl;
            } else {
                mMapperSettings.AddValue("search_settings", Parameters());
            }

            mMapperSettings["search_settings"].AddEmptyValue("max_num_search_iterations").SetInt(search_iterations);
            mMapperSettings.RemoveValue("search_iterations");
        }

        MapperUtilities::CheckInterfaceModelParts(0);

        const Parameters mapper_default_settings(GetMapperDefaultSettings());

        mMapperSettings.ValidateAndAssignDefaults(mapper_default_settings);

        if (!mMapperSettings["search_settings"].Has("echo_level")) {
            // use the echo level of the mapper in case none was specified for the search
            mMapperSettings["search_settings"].AddEmptyValue("echo_level").SetInt(mMapperSettings["echo_level"].GetInt());
        }
    }

private:

    ///@name Private Operations
    ///@{
    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;
    bool mOriginIsIga = false;
    IndexType mRequiredRBFSupportPoints;

    ModelPart* mpCouplingMP = nullptr;
    ModelPart* mpCouplingInterfaceOrigin = nullptr;
    ModelPart* mpCouplingInterfaceDestination = nullptr;

    Parameters mMapperSettings;

    MapperLocalSystemPointerVector mMapperLocalSystemsOrigin;
    MapperLocalSystemPointerVector mMapperLocalSystemsDestination;

    InterfaceVectorContainerPointerType mpInterfaceVectorContainerOrigin;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerDestination;

    // Matrices for the mapping operation
    MappingMatrixUniquePointerType mpMappingMatrix;
    MappingMatrixUniquePointerType mpOriginInterpolationMatrix;
    MappingMatrixUniquePointerType mpDestinationEvaluationMatrix;

    MapperUniquePointerType mpInverseMapper = nullptr;

    LinearSolverSharedPointerType mpLinearSolver = nullptr;

    std::vector<IndexType> mPolynomialEquationIdsOrigin;
    std::vector<IndexType> mPolynomialEquationIdsDestination;

    std::string mRBFType;
    IndexType mPolynomialDegree;
    IndexType mNumberOfPolynomialTerms;

    MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo(const IndexType RequiredRBFSupportPoints) const 
    {
        if (mOriginIsIga){
            return Kratos::make_unique<RadialBasisFunctionsMapperInterfaceInfoIGA>(RequiredRBFSupportPoints);
        } else {
            return Kratos::make_unique<RadialBasisFunctionsMapperInterfaceInfoFEM>(RequiredRBFSupportPoints);
        }
    }

    void InitializeInterface(Kratos::Flags MappingOptions = Kratos::Flags());

    void AssignInterfaceEquationIds()
    {
        MapperUtilities::AssignInterfaceEquationIds(mpCouplingInterfaceDestination->GetCommunicator());
        MapperUtilities::AssignInterfaceEquationIds(mpCouplingInterfaceOrigin->GetCommunicator());
    }

    void AssignInterfaceEquationIdsIga()
    {
        MapperUtilities::AssignInterfaceEquationIds(mpCouplingInterfaceDestination->GetCommunicator());
        MapperUtilities::AssignInterfaceEquationIdsOnConditions(mpCouplingInterfaceOrigin->GetCommunicator());
    }

    void MapInternal(const Variable<double>& rOriginVariable,
                     const Variable<double>& rDestinationVariable,
                     Kratos::Flags MappingOptions);

    void MapInternalTranspose(const Variable<double>& rOriginVariable,
                              const Variable<double>& rDestinationVariable,
                              Kratos::Flags MappingOptions);

    void MapInternal(const Variable<array_1d<double, 3>>& rOriginVariable,
                     const Variable<array_1d<double, 3>>& rDestinationVariable,
                     Kratos::Flags MappingOptions);

    void MapInternalTranspose(const Variable<array_1d<double, 3>>& rOriginVariable,
                              const Variable<array_1d<double, 3>>& rDestinationVariable,
                              Kratos::Flags MappingOptions);

    void CreateLinearSolver();

    std::vector<IndexType> InitializePolynomialEquationIds(const ModelPart* pModelPart);

    std::vector<IndexType> InitializePolynomialEquationIdsIga(const ModelPart* pModelPart);

    void CalculateMappingMatrix();
    
    // For IGA, compute the mapping matrix mapping from origin control points to destination nodes
    std::unique_ptr<MappingMatrixType> ComputeMappingMatrixIgaOnControlPoints(const MappingMatrixType& rMappingMatrixGP, const ModelPart& rOriginModelPart) const;

    Parameters GetMapperDefaultSettings() const 
    {
        return Parameters(R"({
            "echo_level"                     : 0,
            "radial_basis_function_type"     : "thin_plate_spline",
            "additional_polynomial_degree"   : 0,
            "origin_is_iga"                  : false,
            "destination_is_iga"             : false,
            "precompute_mapping_matrix"      : true,
            "search_settings"                : {},
            "linear_solver_settings"         : {}
        })");
    }

    ///@}
    ///@name Private  Access
    ///@{

    MapperUniquePointerType& GetInverseMapper()
    {
        if (!mpInverseMapper) {
            InitializeInverseMapper();
        }
        return mpInverseMapper;
    }

    void InitializeInverseMapper()
    {
        KRATOS_ERROR << "Inverse Mapping is not supported yet!" << std::endl;
        mpInverseMapper = this->Clone(mrModelPartDestination,
                                      mrModelPartOrigin,
                                      mMapperSettings);
    }

    ///@}

}; // Class RadialBasisFunctionMapper

///@} addtogroup block
}  // namespace Kratos.