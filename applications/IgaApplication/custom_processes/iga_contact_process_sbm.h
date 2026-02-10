//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Andrea Gorgi

#pragma once

// System includes

// External includes

// Project includes
#include "containers/pointer_vector.h"
#include "containers/model.h"
#include "geometries/geometry.h"

#include "processes/process.h"

#include "geometries/nurbs_coupling_geometry_2d.h"

#include "custom_conditions/sbm_contact_2D_condition.h"

#include "utilities/entities_utilities.h"

#include "spatial_containers/bins_dynamic.h"

#include <cmath>
#include <vector>
#include <string>

#include "spaces/ublas_space.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/* @class IgaContactProcessSbm
 * @ingroup IgaApplication
 **/
class KRATOS_API(IGA_APPLICATION) IgaContactProcessSbm
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IgaContactProcessSbm
    KRATOS_CLASS_POINTER_DEFINITION(IgaContactProcessSbm);

    typedef Node                                             NodeType;
    typedef Geometry<NodeType>                                  GeometryType;
    typedef GeometryType::Pointer                               GeometryPointerType;
    typedef typename GeometryType::GeometriesArrayType          GeometriesArrayType;
    typedef typename GeometryType::CoordinatesArrayType         CoordinatesArrayType;
    typedef typename GeometryType::IntegrationPointsArrayType   IntegrationPointsArrayType;

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    typedef typename Properties::Pointer PropertiesPointerType;

    typedef typename ModelPart::ElementsContainerType ElementsContainerType;
    typedef typename ModelPart::ConditionsContainerType ConditionsContainerType;

    // MODIFIED --------------------------------------
    typedef PointerVector<Node> ContainerNodeType;
    typedef PointerVector<Point> ContainerEmbeddedNodeType;

    typedef BrepCurveOnSurface<ContainerNodeType, true, ContainerEmbeddedNodeType> BrepCurveOnSurfaceType;

    typedef DenseVector<typename BrepCurveOnSurfaceType::Pointer> BrepCurveOnSurfaceArrayType;
    //--------------------------------------------------------------------------------------

    using PointType = Node;

    using ConditionsArrayType = ModelPart::ConditionsContainerType;
    using NodesArrayType = ModelPart::NodesContainerType;

    using NodePointerVector = GlobalPointersVector<NodeType>;

    using PointTypePointer = Node::Pointer;
    using PointVector = std::vector<PointType::Pointer>;
    using PointIterator = std::vector<PointType::Pointer>::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = std::vector<double>::iterator;
    using DynamicBins = BinsDynamic<3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;
    using PointerType = DynamicBins::PointerType;

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using SparseMatrixType = SparseSpaceType::MatrixType;
    
    using NurbsCurveGeometryType = NurbsCurveGeometry<2, PointerVector<Node>>;


    typedef typename Kratos::shared_ptr<Kratos::NurbsCouplingGeometry2D<PointType, PointerVector<NodeType>, false>> NurbsCouplingGeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    IgaContactProcessSbm(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~IgaContactProcessSbm() = default;

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    // void ExecuteInitialize() override {
    //     Execute();
    // };

    void ExecuteInitializeSolutionStep() override {
        Execute();
    };

    // void ExecuteBeforeSolutionLoop() override {
    //     Execute();
    // };

    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "main_model_part_name" : "ModelPart",
            "nurbs_volume_name" : "NurbsVolume",
            "embedded_model_part_name" : "IgaModelPart"
        })" );

        return default_parameters;
    }


        /// Creates conditions from geometries
    void CreateConditions(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rDestinationModelPart,
        std::string& rConditionName,
        SizeType& rIdCounter,
        PropertiesPointerType pPropMaster,
        PropertiesPointerType pPropSlave) const;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "IgaContactProcessSbm";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IgaContactProcessSbm";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:
    ///@name Iga functionalities
    ///@{

    Model* mpModel = nullptr;
    Parameters mParameters;
    SizeType mEchoLevel;
    bool mIntegrateOnTrueBoundary = false;

    ModelPart* mrSlaveModelPart = nullptr; 
    ModelPart* mrMasterModelPart = nullptr; 
    ModelPart* mrSlaveSkinModelPart = nullptr; 
    ModelPart* mrMasterSkinModelPart = nullptr; 
    
    ModelPart* mrContactModelPart = nullptr; 

    SparseMatrixType mSparseBrepMatrixSlave;
    SparseMatrixType mSparseBrepMatrixMaster;


    Properties::Pointer mpPropMaster;
    Properties::Pointer mpPropSlave;


    bool GetProjection(CoordinatesArrayType& slavePoint, GeometryType &slave_geometry, GeometryType &master_geometry, double& localProjection);

    ///@}
    ///@name Iga functionalities
    ///@{

    ///@}
    ///@name CAD functionalities
    ///@{

    /// Gets list of geometries from rModelPart
    void GetCadGeometryList(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters) const;

    ///@}

    ///@}
    ///@name Utility
    ///@{

    Parameters ReadParamatersFile(
        const std::string& rDataFileName) const;

    
    static void GetValuesVector(
        const Geometry<PointType>& rGeometry,
        Vector& rValues);

    
    static void GetDeformedPosition(
        const CoordinatesArrayType& rPointGlobalCoordinates, 
        const ModelPart& rModelPart,
        const SparseMatrixType& rSparseBrepMatrix,
        CoordinatesArrayType& rPointDeformedCoordinates);

    static void GetDeformedGradient(
        const CoordinatesArrayType& rPointGlobalCoordinates, 
        const ModelPart& rModelPart,
        const SparseMatrixType& SparseBrepMatrix,
        Matrix& rPointGradientDeformation,
        Matrix& rPointHessianDeformation);

    
    bool NewtonRaphsonCurveOnDeformed(
        CoordinatesArrayType& rProjectedPointLocalCoordinates,
        const CoordinatesArrayType& rPointGlobalCoordinates, 
        CoordinatesArrayType& rProjectedPointGlobalCoordinates,
        const NurbsCurveGeometryType& rPairedGeometry, 
        double& distance,
        const int rNumberOfInitialGuesses,
        const int MaxIterations = 20,
        const double Accuracy = 1e-6
        );
    
    static double ComputeTaylorTerm(double derivative, double dx, int n_k, double dy, int k);

    static unsigned long long Factorial(int n); 

    static std::vector<IndexType> FindClosestBrepId(
        const CoordinatesArrayType& rPointGlobalCoordinates, 
        const ModelPart& rModelPart,
        const SparseMatrixType& SparseBrepMatrix);



    
    bool ProjectPointViaRayTracingToMasterCurveDeformed(
        const CoordinatesArrayType& point_slave_deformed,
        const Vector& normal_slave_deformed,
        const NurbsCurveGeometryType& master_curve,
        const ModelPart& master_model_part,
        const SparseMatrixType& sparse_matrix_master,
        CoordinatesArrayType& rBestProjectedPoint,
        CoordinatesArrayType& rBestProjectedPointLocalCoords,
        double& rBestDistance,
        double tolerance = 1e-8,
        int max_iter = 100,
        int n_initial_guesses = 5);

    
    static bool ProjectToSkinBoundary(
        const ModelPart* ModelPart,
        std::string &rLayerName,
        const CoordinatesArrayType& rPoint,
        CoordinatesArrayType& rProjectedPoint,
        std::vector<array_1d<double, 3>>& rCurveDerivatives,
        int n_initial_guesses = 10);

    
    static bool ProjectBackToSurrogateBoundary(
        const ModelPart &rSurrogateModelPart,
        const CoordinatesArrayType& rSkinPoint,
        const SparseMatrixType& rSparseBrepMatrix,
        const Vector & rSkinNormal,
        IndexType& rBrepId,
        double& rSurrogateProjection
    );

    struct MasterSegmentData {
        BrepCurveOnSurfaceType::Pointer p_brep;
        double local_begin = 0.0;
        double local_end = 0.0;
        NodeType::Pointer p_node_begin = nullptr;
        NodeType::Pointer p_node_end = nullptr;
    };

    NodeType::Pointer FindExistingNode(
        const std::vector<NodeType::Pointer>& rNodes,
        const CoordinatesArrayType& rCoordinates,
        double tolerance) const;

    NodeType::Pointer CreateOrRetrieveContactNode(
        const CoordinatesArrayType& rCoordinates,
        std::vector<NodeType::Pointer>& rNodes,
        IndexType& rNextNodeId,
        double tolerance);

    static CoordinatesArrayType ConvertVectorToCoordinates(const Vector& rVector);

    static bool IsPointOnAxisAlignedSegment(
        const CoordinatesArrayType& rPoint,
        const CoordinatesArrayType& rSegmentBegin,
        const CoordinatesArrayType& rSegmentEnd,
        double tolerance);

    void CreateQuadratureGeometries(
        BrepCurveOnSurfaceType& rBrep,
        GeometriesArrayType& rGeometries,
        const std::vector<double>* pCustomSpans = nullptr);

    void CollectUniqueSlaveVertices(
        std::vector<CoordinatesArrayType>& rVertices,
        double tolerance) const;

    void SplitMasterSegmentsWithSlaveVertices(
        const std::vector<CoordinatesArrayType>& rSlaveVertices,
        std::vector<MasterSegmentData>& rMasterSegments,
        std::vector<NodeType::Pointer>& rMasterNodes,
        IndexType& rNextNodeId,
        double coordinate_tolerance,
        const std::string& master_layer_name,
        const std::string& slave_layer_name);

    void SplitMasterSegment(
        IndexType brep_id,
        double split_parameter,
        NodeType::Pointer p_split_node,
        std::vector<MasterSegmentData>& rMasterSegments,
        double tolerance);

    static void AppendUniqueVertex(
        std::vector<CoordinatesArrayType>& rVertices,
        const CoordinatesArrayType& rVertex,
        double tolerance);

    static bool ContainsVertex(
        const std::vector<CoordinatesArrayType>& rVertices,
        const CoordinatesArrayType& rVertex,
        double tolerance);

    void SetSurrogateNeighbourNodes(
        GeometryPointerType p_geometry,
        BrepCurveOnSurfaceType::Pointer p_brep_geometry,
        const std::string& rSlaveLayerName);
    
    void PrepareIntegrationOnTrueBoundary(
        ModelPart& rMasterSkinModelPart,
        ModelPart& rSlaveSkinModelPart) const;

    ///@}
    ///@name Input and output
    ///@{

    ///@}

}; // Class IgaContactProcessSbm

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  IgaContactProcessSbm& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const IgaContactProcessSbm& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
