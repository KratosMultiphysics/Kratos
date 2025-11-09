//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
// Main authors:    E. G. Loera Villeda
// Contributor:     Juan I. Camarotti
//
// See PhD Thesis Tianyang Wang Chapter 5

#pragma once

// System includes

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
///@name Kratos Classes
///@{

class KRATOS_API(MAPPING_APPLICATION) BeamMapperInterfaceInfo : public MapperInterfaceInfo
{
public:
    using InterfaceObjectPointerType = Kratos::shared_ptr<InterfaceObject>;
    using MatrixType = Matrix;
    using VectorType = Vector;

    /// Default constructor.
    explicit BeamMapperInterfaceInfo(const double LocalCoordTol=0.0) : mLocalCoordTol(LocalCoordTol) {}

    explicit BeamMapperInterfaceInfo(const CoordinatesArrayType& rCoordinates,
                                 const IndexType SourceLocalSystemIndex,
                                 const IndexType SourceRank,
                                 const double LocalCoordTol=0.0)
        : MapperInterfaceInfo(rCoordinates, SourceLocalSystemIndex, SourceRank), mLocalCoordTol(LocalCoordTol) {}

    MapperInterfaceInfo::Pointer Create() const override
    {
        return Kratos::make_shared<BeamMapperInterfaceInfo>();
    }

    MapperInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                        const IndexType SourceLocalSystemIndex,
                                        const IndexType SourceRank) const override
    {
        return Kratos::make_shared<BeamMapperInterfaceInfo>(
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

    void GetValue(std::vector<int>& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mNodeIds;
    }

    void GetValue(double& rValue,
                  const InfoType ValueType) const override
    {
        rValue = mBeamMapperDistance;
    }

    void GetValue(std::vector<double>& rValue,
                  const InfoType ValueType) const override
    {
        std::vector<double> vector = mLinearShapeFunctionValues;
        vector.insert(vector.end(), mHermitianShapeFunctionValues.begin(), mHermitianShapeFunctionValues.end());
        vector.insert(vector.end(), mHermitianShapeFunctionValuesDerivatives.begin(), mHermitianShapeFunctionValuesDerivatives.end());

        rValue = vector;
    }

    void GetValue(int& rValue,
                  const InfoType ValueType) const override
    {   
        rValue = (int)mPairingIndex;
    }

    void GetValue(GeometryType& rValue,
                  const InfoType ValueType) const override
    {
        rValue = *(mpInterfaceObject->pGetBaseGeometry());
    }

    void GetValue(MatrixType& rRotMatrixValue, 
                  VectorType& rTransVectorValue,
                  VectorType& rLinearValue, 
                  VectorType& rHermitianValue, 
                  VectorType& rHermitianDerValue) const override 
    {
        rRotMatrixValue = mRotationMatrixOfBeam;

        rTransVectorValue(0) = mProjectionOfPoint[0];
        rTransVectorValue(1) = mProjectionOfPoint[1];
        rTransVectorValue(2) = mProjectionOfPoint[2];
        
        rLinearValue(0) = mLinearShapeFunctionValues[0];
        rLinearValue(1) = mLinearShapeFunctionValues[1];

        for (size_t i = 0; i < 4; i++){
            rHermitianValue(i) = mHermitianShapeFunctionValues[i];
            rHermitianDerValue(i) = mHermitianShapeFunctionValuesDerivatives[i];
        }
    }

    void ComputeRotationMatrixInterfaceObject()
    {
        ComputeRotationMatrix();
    }


private:
    double mLocalCoordTol; 
    
    std::vector<int> mNodeIds;
    std::vector<double> mLinearShapeFunctionValues;
    std::vector<double> mHermitianShapeFunctionValues;
    std::vector<double> mHermitianShapeFunctionValuesDerivatives;

    Point mProjectionOfPoint; // Point that results form the projection of the surface node on the beam

    double mClosestProjectionDistance = std::numeric_limits<double>::max(); // Distance between the surface node and the beam
    ProjectionUtilities::PairingIndex mPairingIndex = ProjectionUtilities::PairingIndex::Unspecified;

    InterfaceObjectPointerType mpInterfaceObject;

    MatrixType mRotationMatrixOfBeam;

    double mBeamMapperDistance = std::numeric_limits<double>::max();

    void SaveSearchResult(const InterfaceObject& rInterfaceObject,
                          const bool ComputeApproximation);

    // This computes the rotation matrix of the beam pointed by mpInterfaceObject
    void ComputeRotationMatrix(); 

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.save("NodeIds", mNodeIds);
        rSerializer.save("LSFValues", mLinearShapeFunctionValues);
        rSerializer.save("HSFValues", mHermitianShapeFunctionValues);
        rSerializer.save("HSFDValues", mHermitianShapeFunctionValuesDerivatives);
        rSerializer.save("ClosestProjectionDistance", mClosestProjectionDistance);
        rSerializer.save("PairingIndex", (int)mPairingIndex);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.load("NodeIds", mNodeIds);
        rSerializer.load("LSFValues", mLinearShapeFunctionValues);
        rSerializer.load("HSFValues", mHermitianShapeFunctionValues);
        rSerializer.load("HSFDValues", mHermitianShapeFunctionValuesDerivatives);
        rSerializer.load("ClosestProjectionDistance", mClosestProjectionDistance);
        int temp;
        rSerializer.load("PairingIndex", temp);
        mPairingIndex = (ProjectionUtilities::PairingIndex)temp;
    }
};

class KRATOS_API(MAPPING_APPLICATION) BeamMapperLocalSystem : public MapperLocalSystem
{
public:
    using VectorType = Vector;
    using BeamMapperInterfaceInfoPointerType = Kratos::shared_ptr<BeamMapperInterfaceInfo>;
    using GeometryType = Geometry<Node>;

    explicit BeamMapperLocalSystem(NodePointerType pNode) : mpNode(pNode) {
        VectorType zeroVector(3, 0.0);
        mRotationVectorOfSection = zeroVector;
    }

    void CalculateAll(MatrixType& rLocalMappingMatrix,
                      EquationIdVectorType& rOriginIds,
                      EquationIdVectorType& rDestinationIds,
                      MapperLocalSystem::PairingStatus& rPairingStatus) const override
    {
        KRATOS_WARNING("BeamMapperLocalSystem")
            << "CalculateAll() was called, but is not implemented for BeamMapperLocalSystem." << std::endl;
    }

    CoordinatesArrayType& Coordinates() const override
    {
        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;
        return mpNode->Coordinates();
    }

    void CalculateRotationMatrixInterfaceInfos(MatrixType& rRotationMatrix_G_B,
                                               VectorType& rTranslationVector_B_P,
                                               VectorType& rLinearShapeValues,
                                               VectorType& rHermitianShapeValues,
                                               VectorType& rHermitianDerShapeValues,
                                               GeometryType& rGeom,
                                               NodePointerType& pNode)
    {
        pNode = mpNode;

        for( auto& r_interface_info : mInterfaceInfos ){ 
            // Safely cast to BeamMapperInterfaceInfo
            auto beam_interface_info = std::dynamic_pointer_cast<BeamMapperInterfaceInfo>(r_interface_info);
            KRATOS_ERROR_IF_NOT(beam_interface_info) 
                << "Expected BeamMapperInterfaceInfo in mInterfaceInfos but got nullptr or wrong type." << std::endl;

            // Call beam-specific function
            beam_interface_info->ComputeRotationMatrixInterfaceObject();

            // Call base-class GetValue functions
            r_interface_info->GetValue(rRotationMatrix_G_B,
                                    rTranslationVector_B_P,
                                    rLinearShapeValues,
                                    rHermitianShapeValues,
                                    rHermitianDerShapeValues);

            r_interface_info->GetValue(rGeom, MapperInterfaceInfo::InfoType::Dummy);
            
            const Point point_to_proj(mpNode->Coordinates());
            Point projection_point;
            GeometricalProjectionUtilities::FastProjectOnLine(rGeom, point_to_proj, projection_point);
            rTranslationVector_B_P = projection_point;
        }
    }

    void SaveRotationVectorValue(const VectorType& rotationVector)
    {
        mRotationVectorOfSection = rotationVector;
    }

    void GetValue(VectorType& rotVectorValue)
    {
        rotVectorValue(0) = mRotationVectorOfSection(0);
        rotVectorValue(1) = mRotationVectorOfSection(1);
        rotVectorValue(2) = mRotationVectorOfSection(2);
    }

    MapperLocalSystemUniquePointer Create(NodePointerType pNode) const override
    {
        return Kratos::make_unique<BeamMapperLocalSystem>(pNode);
    }

    void PairingInfo(std::ostream& rOStream, const int EchoLevel) const override;

private:
    NodePointerType mpNode;
    VectorType mRotationVectorOfSection;
};

/// Beam Mapper
/** This class implements the Beam Mapping technique.
 * It is based on the work of T. Wang (PhD Thesis, Chapter 5)
 * and allows to map from a surface mesh to a beam mesh and vice versa.
 * The mapping is performed in a way that displacement and rotation
 * variables are mapped from the beam to the surface mesh and force
 * and moment variables are mapped from the surface mesh to the beam.
 * The mapping is consistent in a way that mechanical work is conserved.
 * This mapper can be used for example to couple a 3D continuum
 * structure with a 1D beam structure.
 * This mapper currently only works for linear beams (2-noded elements).
*/
template<class TSparseSpace, class TDenseSpace>
class KRATOS_API(MAPPING_APPLICATION) BeamMapper
    : public Mapper<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of BeamMapper
    KRATOS_CLASS_POINTER_DEFINITION(BeamMapper);

    using BaseType = Mapper<TSparseSpace, TDenseSpace>;

    using InterfaceCommunicatorPointerType = Kratos::unique_ptr<InterfaceCommunicator>;
    using MapperInterfaceInfoUniquePointerType = typename InterfaceCommunicator::MapperInterfaceInfoUniquePointerType;

    using InterfaceVectorContainerType = InterfaceVectorContainer<TSparseSpace, TDenseSpace>;
    using InterfaceVectorContainerPointerType = Kratos::unique_ptr<InterfaceVectorContainerType>;

    using MapperLocalSystemPointer = Kratos::unique_ptr<MapperLocalSystem>;
    using MapperLocalSystemPointerVector = std::vector<MapperLocalSystemPointer>;

    using MapperUniquePointerType = typename BaseType::MapperUniquePointerType;
    using TMappingMatrixType = typename BaseType::TMappingMatrixType;
    using TMappingMatrixUniquePointerType = Kratos::unique_ptr<TMappingMatrixType>;

    using MatrixType = typename TDenseSpace::MatrixType;
    using RotationMatrixVector = std::vector<MatrixType>;
    using VectorType = typename TDenseSpace::VectorType;

    using ComponentVariableType = Variable<double>;
    using GeometryType = Geometry<Node>;
    using GeometryPointerType = InterfaceObject::GeometryPointerType;

    using NodeType = InterfaceObject::NodeType;
    using NodePointerType = InterfaceObject::NodePointerType;
    ///@}
    ///@name Life Cycle
    ///@{

    // Default constructor, needed for registration
    BeamMapper(ModelPart& rModelPartOrigin,
               ModelPart& rModelPartDestination):
               mrModelPartOrigin(rModelPartOrigin),
               mrModelPartDestination(rModelPartDestination) {}

    BeamMapper(ModelPart& rModelPartOrigin,
               ModelPart& rModelPartDestination,
               Parameters JsonParameters):
               mrModelPartOrigin(rModelPartOrigin),
               mrModelPartDestination(rModelPartDestination),
               mMapperSettings(JsonParameters)
                          
    {
        KRATOS_TRY;

        mpInterfaceVectorContainerOrigin = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartOrigin);
        mpInterfaceVectorContainerDestination = Kratos::make_unique<InterfaceVectorContainerType>(rModelPartDestination);

        ValidateInput();
        
        mLocalCoordTol = JsonParameters["local_coord_tolerance"].GetDouble();
        KRATOS_ERROR_IF(mLocalCoordTol < 0.0) << "The local_coord_tolerance cannot be negative" << std::endl;
        
        //  In this function the search task is done, and the parameters for the local systems are stored
        // MapperInterfaceInfo has:
        // 1. A beam id to relate
        // 2. a t_B_P
        // 3. linear and hermitian interpolation values
        // 4. its own position in GCS space

        Initialize();

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(Kratos::Flags MappingOptions, double SearchRadius) override
    {
        // Set the Flags according to the type of remeshing
        if (MappingOptions.Is(MapperFlags::REMESHED)) {
            InitializeInterface(MappingOptions);
        }
        else {
            BuildProblem(MappingOptions); // in interpolative_mapper_base this is BuildMappingMatrix
        }
    }

    void Map(const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements, 
              const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
              const Variable<array_1d<double, 3>>& rDestinationVariable, 
              Kratos::Flags MappingOptions)
    {
        if (mMapperSettings["use_corotation"].GetBool() == false)
        {
            MapInternal( rOriginVariablesDisplacements, rOriginVariablesRotations, rDestinationVariable, MappingOptions );
        }
        else
        {
            MapInternalCorotation ( rOriginVariablesDisplacements, rOriginVariablesRotations, rDestinationVariable, MappingOptions );
        }        
    }

    void Map( const Variable<double>& rOriginVariable, const Variable<double>& rDestinationVariable,
              Kratos::Flags MappingOptions) override
    {
        KRATOS_ERROR << "This function is not supported for the Beam-Mapper!" << std::endl;
    }

    void Map(const Variable<array_1d<double, 3>>& rOriginVariable, const Variable<array_1d<double, 3>>& rDestinationVariable,
              Kratos::Flags MappingOptions) override
    {
        KRATOS_ERROR << "This function is not supported for the Beam-Mapper!" << std::endl;
    }

    void InverseMap(const Variable<double>& rOriginVariable, const Variable<double>& rDestinationVariable, 
                     Kratos::Flags MappingOptions) override
    {
        KRATOS_ERROR << "This function is not supported for the Beam-Mapper!" << std::endl;
    }

    void InverseMap(const Variable<array_1d<double, 3>>& rOriginVariable, const Variable<array_1d<double, 3>>& rDestinationVariable,
                     Kratos::Flags MappingOptions) override
    {
        KRATOS_ERROR << "This function is not supported for the Beam-Mapper!" << std::endl;
    }

    void InverseMap( const Variable<array_1d<double, 3>>& rOriginVariablesForces, 
                     const Variable<array_1d<double, 3>>& rOriginVariablesMoments,
                     const Variable<array_1d<double, 3>>& rDestinationVaribleForces,
                     Kratos::Flags MappingOptions)
    {
        InitializeOriginForcesAndMoments(rOriginVariablesForces, rOriginVariablesMoments);
        InitializeInformationBeamsInverse(rOriginVariablesForces, rOriginVariablesMoments, rDestinationVaribleForces, MappingOptions);
    }

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                                  ModelPart& rModelPartDestination,
                                  Parameters JsonParameters) const override
    {
        KRATOS_TRY;

        return Kratos::make_unique<BeamMapper<TSparseSpace, TDenseSpace>>(
            rModelPartOrigin,
            rModelPartDestination,
            JsonParameters);

        KRATOS_CATCH("");
    }

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
    std::string Info() const override
    {
        return "BeamMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BeamMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

    void ValidateInput(){
        Parameters mapper_default_settings(GetMapperDefaultSettings());
        mMapperSettings.ValidateAndAssignDefaults(mapper_default_settings);

        if (mMapperSettings["search_radius"].GetDouble() < 0.0) {
            const double search_radius = MapperUtilities::ComputeSearchRadius(
                                            mrModelPartOrigin,
                                            mrModelPartDestination,
                                            mMapperSettings["echo_level"].GetInt());
            mMapperSettings["search_radius"].SetDouble(search_radius);
        }
    };

    void Initialize(){
        InitializeInterfaceCommunicator();
        InitializeInterface();
    };

private:
    ///@name Member Variables
    ///@{
    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;

    Parameters mMapperSettings;

    MapperLocalSystemPointerVector mMapperLocalSystems;

    InterfaceCommunicatorPointerType mpInterfaceCommunicator;

    RotationMatrixVector mRotationMatrixOfBeams;

    double mLocalCoordTol;

    InterfaceVectorContainerPointerType mpInterfaceVectorContainerOrigin;
    InterfaceVectorContainerPointerType mpInterfaceVectorContainerDestination;

    ///@name Private Operations
    ///@{
    void InitializeInterfaceCommunicator()
    {
        mpInterfaceCommunicator = Kratos::make_unique<InterfaceCommunicator>(
            mrModelPartOrigin, mMapperLocalSystems, mMapperSettings["search_settings"]
        );
    }

    void InitializeInterface(Kratos::Flags MappingOptions = Kratos::Flags()){
        CreateMapperLocalSystems(mrModelPartDestination.GetCommunicator(), mMapperLocalSystems);
        BuildProblem(MappingOptions);
    };

    void InitializeInformationBeams(const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
                                    const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
                                    const Variable<array_1d<double, 3>>& rDestinationVariableDisplacement);
    
    void InitializeInformationBeamsCorotation(const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
                                    const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
                                    const Variable<array_1d<double, 3>>& rDestinationVariableDisplacement);

    void InitializeInformationBeamsInverse(const Variable< array_1d<double, 3> >& rOriginVariablesForces,
                                    const Variable<array_1d<double, 3>>& rOriginVariablesMoments,
                                    const Variable<array_1d<double, 3>>& rDestinationVariableForces,
                                    const Kratos::Flags& rMappingOptions);

    void InitializeOriginForcesAndMoments(const Variable<array_1d<double, 3>>& rOriginVariablesForces,
                                    const Variable<array_1d<double, 3>>& rOriginVariablesMoments);

    void BuildProblem(Kratos::Flags MappingOptions = Kratos::Flags()){
        MapperUtilities::AssignInterfaceEquationIds(mrModelPartOrigin.GetCommunicator());
        MapperUtilities::AssignInterfaceEquationIds(mrModelPartDestination.GetCommunicator());

        KRATOS_ERROR_IF_NOT(mpInterfaceCommunicator) << "mpInterfaceCommunicator is a nullptr" << std::endl;
        const MapperInterfaceInfoUniquePointerType p_ref_interface_info = GetMapperInterfaceInfo();
        mpInterfaceCommunicator->ExchangeInterfaceData(mrModelPartDestination.GetCommunicator(),
                                                    p_ref_interface_info);
    };

    void CreateMapperLocalSystems(
        const Communicator& rModelPartCommunicator,
        std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rLocalSystems)
    {
       MapperUtilities::CreateMapperLocalSystemsFromNodes(
            BeamMapperLocalSystem(nullptr),
            rModelPartCommunicator,
            rLocalSystems);
    }

    void MapInternal(const Variable< array_1d<double, 3> >& rOriginVariablesDisplacements,
                     const Variable< array_1d<double, 3> >& rOriginVariablesRotations,
                     const Variable< array_1d<double, 3> >& rDestinationVariableDisplacement,
                     Kratos::Flags MappingOptions)
    {   
        InitializeInformationBeams(rOriginVariablesDisplacements, rOriginVariablesRotations, rDestinationVariableDisplacement);
    }

    void MapInternalCorotation(const Variable< array_1d<double, 3> >& rOriginVariablesDisplacements,
                               const Variable< array_1d<double, 3> >& rOriginVariablesRotations,
                               const Variable< array_1d<double, 3> >& rDestinationVariableDisplacement,
                               Kratos::Flags MappingOptions)
    {
        InitializeInformationBeamsCorotation(rOriginVariablesDisplacements, rOriginVariablesRotations, rDestinationVariableDisplacement);
    }
    
    void CalculateRotationMatrixWithAngle(VectorType& rAxis, double& rAngle , MatrixType& rRotationMatrix);

    void GetRotationVector(const MatrixType& rRotationMatrix, VectorType& rRotationVector);

    MapperInterfaceInfoUniquePointerType GetMapperInterfaceInfo() const 
    {
        return Kratos::make_unique<BeamMapperInterfaceInfo>();
    }

    Parameters GetMapperDefaultSettings() const 
    {
        return Parameters( R"({
            "search_settings"              : {},
            "search_radius"            : -1.0,
            "search_iterations"        : 3,
            "local_coord_tolerance"    : 0.25,
            "echo_level"               : 0,
            "use_corotation"           : true
        })");
    }

    ///@}

}; // Class BeamMapper

///@} addtogroup block
}  // namespace Kratos.