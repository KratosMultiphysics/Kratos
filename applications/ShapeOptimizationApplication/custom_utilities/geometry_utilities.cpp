// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <functional>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "geometry_utilities.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/key_hash.h"
#include "shape_optimization_application.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "includes/global_pointer_variables.h"
#include "utilities/element_size_calculator.h"
#include "utilities/atomic_utilities.h"
#include "geometries/geometry_data.h"
#include "processes/calculate_nodal_area_process.h"

#include "spatial_containers/spatial_containers.h"

// ==============================================================================

namespace Kratos
{

void GeometryUtilities::ComputeUnitSurfaceNormals()
{
    KRATOS_TRY;

    const unsigned int domain_size = mrModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE);
    KRATOS_ERROR_IF(mrModelPart.NumberOfConditions() == 0) <<
        "> Normal calculation requires surface or line conditions to be defined!" << std::endl;
    KRATOS_ERROR_IF((domain_size == 3 && mrModelPart.ConditionsBegin()->GetGeometry().size() == 2)) <<
        "> Normal calculation of 2-noded conditions in 3D domains is not possible!" << std::endl;
    CalculateAreaNormalsFromConditions();
    CalculateUnitNormals();

    KRATOS_CATCH("");
}

void GeometryUtilities::ProjectNodalVariableOnUnitSurfaceNormals( const Variable<array_3d> &rNodalVariable )
{
    KRATOS_TRY;

    ProjectNodalVariableOnDirection(rNodalVariable, NORMALIZED_SURFACE_NORMAL);

    KRATOS_CATCH("");
}

void GeometryUtilities::ProjectNodalVariableOnDirection( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rDirectionVariable)
{
    KRATOS_TRY;

    for (ModelPart::NodeIterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); ++node_i)
    {
        array_3d &nodal_variable = node_i->FastGetSolutionStepValue(rNodalVariable);
        array_3d &node_normal = node_i->FastGetSolutionStepValue(rDirectionVariable);

        const double magnitude = inner_prod(nodal_variable, node_normal);
        noalias(nodal_variable) = magnitude * node_normal;
    }

    KRATOS_CATCH("");
}

void GeometryUtilities::ProjectNodalVariableOnTangentPlane( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rPlaneNormalVariable)
{
    KRATOS_TRY;

    for (ModelPart::NodeIterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); ++node_i)
    {
        array_3d &nodal_variable = node_i->FastGetSolutionStepValue(rNodalVariable);
        array_3d &node_normal = node_i->FastGetSolutionStepValue(rPlaneNormalVariable);

        const double magnitude = inner_prod(nodal_variable, node_normal);
        nodal_variable -= magnitude * node_normal;
    }

    KRATOS_CATCH("");
}

void GeometryUtilities::ExtractBoundaryNodes( std::string const& rBoundarySubModelPartName )
{
    KRATOS_TRY;

    ModelPart& r_boundary_model_part = mrModelPart.GetSubModelPart(rBoundarySubModelPartName);

    KRATOS_ERROR_IF(r_boundary_model_part.Nodes().size() != 0) << "ExtractBoundaryNodes: The boundary model part already has nodes!" << std::endl;

    // Some type-definitions
    typedef std::unordered_map<vector<unsigned int>, unsigned int, KeyHasherRange<vector<unsigned int>>, KeyComparorRange<vector<unsigned int>> > hashmap;

    // Create map to ask for number of faces for the given set of node ids representing one face in the model part
    hashmap n_boundaries_map;

    unsigned int domain_size = static_cast<unsigned int>(mrModelPart.GetProcessInfo()[DOMAIN_SIZE]);

    // Fill map that counts number of faces for given set of nodes
    for (auto& elem_i : mrModelPart.Elements())
    {
        KRATOS_ERROR_IF(elem_i.GetGeometry().WorkingSpaceDimension() < domain_size) << "ExtractBoundaryNodes: This function does only work"
            <<" for solid elements in 3D and surface elements in 2D!" << std::endl;

        Element::GeometryType::GeometriesArrayType boundaries = elem_i.GetGeometry().GenerateBoundariesEntities();

        for(unsigned int boundary=0; boundary<boundaries.size(); boundary++)
        {
            // Create vector that stores all node is of current face
            DenseVector<unsigned int> ids(boundaries[boundary].size());

            // Store node ids
            for(unsigned int i=0; i<boundaries[boundary].size(); i++)
                ids[i] = boundaries[boundary][i].Id();

            //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
            std::sort(ids.begin(), ids.end());

            // Fill the map
            n_boundaries_map[ids] += 1;
        }
    }

    // Vector to store all nodes on surface. Node ids may be listed several times
    std::vector<std::size_t> temp_boundary_node_ids;

    // Add surface nodes to sub-model part
    for(auto it=n_boundaries_map.begin(); it!=n_boundaries_map.end(); it++)
    {
        // If given node set represents face that is not overlapping with a face of another element, add it as skin element
        if(it->second == 1)
        {
            for(unsigned int i=0; i<it->first.size(); i++)
                temp_boundary_node_ids.push_back(it->first[i]);
        }
    }

    // Add nodes and remove double entries
    r_boundary_model_part.AddNodes(temp_boundary_node_ids);

    KRATOS_CATCH("");
}

void GeometryUtilities::ExtractEdgeNodes( std::string const& rEdgeSubModelPartName ) {
    KRATOS_TRY;

    KRATOS_ERROR_IF(mrModelPart.Conditions().size() == 0) << "ExtractEdgeNodes: No elements defined. Automatic edge detection will not find any edge nodes!" << std::endl;

    ModelPart& r_edge_model_part = mrModelPart.GetSubModelPart(rEdgeSubModelPartName);

    KRATOS_ERROR_IF(r_edge_model_part.Nodes().size() != 0) << "ExtractEdgeNodes: The edge model part already has nodes!" << std::endl;

    for (auto& r_node_i : mrModelPart.Nodes()) {
        auto& r_node_i_neighbours = r_node_i.GetValue(NEIGHBOUR_NODES);  // does not work with const
        for(const auto& r_node_j : r_node_i_neighbours) {
            bool node_i_and_j_share_edge = false;
            auto& r_element_neighbours = r_node_i.GetValue(NEIGHBOUR_CONDITIONS);  // does not work with const
            for(const auto& r_elem_k : r_element_neighbours) {
                Condition::GeometryType::GeometriesArrayType edges = r_elem_k.GetGeometry().GenerateEdges();
                for (const auto& edge : edges) {
                    bool edge_has_node_i = false;
                    bool edge_has_node_j = false;
                    for(const auto& r_node_m : edge) {
                        if (r_node_m.Id() == r_node_i.Id()) {
                            edge_has_node_i = true;
                        }
                        if (r_node_m.Id() == r_node_j.Id()) {
                            edge_has_node_j = true;
                        }
                    }
                    if (edge_has_node_i && edge_has_node_j) {
                        node_i_and_j_share_edge = true;
                        break;
                    }
                }
                if (node_i_and_j_share_edge) break;
            }
            if (!node_i_and_j_share_edge) continue;

            int count = 0;
            for(const auto& r_elem_k : r_element_neighbours) {
                const auto& r_element_geometry = r_elem_k.GetGeometry();
                for(const auto& r_node_l : r_element_geometry) {
                    if (r_node_l.Id() == r_node_j.Id()) count ++;
                }
            }
            if (count < 2){
                r_edge_model_part.AddNode(&r_node_i);
                break;
            }
        }
    }
    KRATOS_CATCH("");
}

std::tuple<std::vector<double>, std::vector<double>> GeometryUtilities::ComputeDistancesToBoundingModelPart(ModelPart& rBoundingModelPart)
{
    KRATOS_TRY;

    typedef Node NodeType;
    typedef NodeType::Pointer NodeTypePointer;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double>::iterator DoubleVectorIterator;
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    KRATOS_ERROR_IF(rBoundingModelPart.NumberOfElements() != 0) <<
        "ComputeDistancesToBoundingModelPart: Model part must only contain conditions!" << std::endl;

    NodeVector all_bounding_nodes;
    all_bounding_nodes.reserve(rBoundingModelPart.Nodes().size());
    for (ModelPart::NodesContainerType::iterator node_it = rBoundingModelPart.NodesBegin(); node_it != rBoundingModelPart.NodesEnd(); ++node_it)
    {
        all_bounding_nodes.push_back(*(node_it.base()));
    }
    const size_t bucket_size = 100;
    KDTree search_tree(all_bounding_nodes.begin(), all_bounding_nodes.end(), bucket_size);

    GeometryUtilities(rBoundingModelPart).ComputeUnitSurfaceNormals();

    std::tuple<std::vector<double>, std::vector<double>> distances_and_directions;
    std::vector<double>& r_signed_distances = std::get<0>(distances_and_directions);
    std::vector<double>& r_directions = std::get<1>(distances_and_directions);
    r_signed_distances.reserve(mrModelPart.NumberOfNodes());
    r_directions.reserve(mrModelPart.NumberOfNodes()*3);

    for (auto& r_node : mrModelPart.Nodes()){

        double distance;
        NodeTypePointer p_neighbor = search_tree.SearchNearestPoint(r_node, distance);

        const array_3d delta = r_node.Coordinates() - p_neighbor->Coordinates();
        const array_3d& bounding_normal = p_neighbor->FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
        const double projected_length = inner_prod(delta, bounding_normal);

        r_signed_distances.push_back(projected_length);

        r_directions.push_back(bounding_normal[0]);
        r_directions.push_back(bounding_normal[1]);
        r_directions.push_back(bounding_normal[2]);
    }

    return distances_and_directions;

    KRATOS_CATCH("");
}

double GeometryUtilities::ComputeVolume()
{
    KRATOS_TRY

    const double volume = block_for_each<SumReduction<double>>(
        mrModelPart.Elements(), [&](const ModelPart::ElementType& rElement) {
            return rElement.GetGeometry().DomainSize();
        });

    return mrModelPart.GetCommunicator().GetDataCommunicator().SumAll(volume);

    KRATOS_CATCH("");
}

void GeometryUtilities::ComputeVolumeShapeDerivatives(
    const Variable<array_3d>& rDerivativeVariable)
{
    KRATOS_TRY

    using CUInt = const unsigned int;
    using VolumeDerivativeMethodType = std::function<double(CUInt, CUInt, const GeometryType&)>;

    KRATOS_ERROR_IF_NOT(mrModelPart.HasNodalSolutionStepVariable(rDerivativeVariable))
        << rDerivativeVariable.Name() << " not found in solution step variables list in "
        << mrModelPart.FullName() << ".\n";

    VariableUtils().SetHistoricalVariableToZero(rDerivativeVariable, mrModelPart.Nodes());

    block_for_each(mrModelPart.Elements(), VolumeDerivativeMethodType(), [&](ModelPart::ElementType& rElement, VolumeDerivativeMethodType& rVolumeDerivativeMethodType){
        auto& r_geometry = rElement.GetGeometry();
        const auto& geometry_type = r_geometry.GetGeometryType();
        const SizeType dimension = r_geometry.WorkingSpaceDimension();

        switch (geometry_type) {
            case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                rVolumeDerivativeMethodType = [](CUInt NodeIndex, CUInt DirectionIndex,  const GeometryType& rGeometry) {
                    return 2.0 * ElementSizeCalculator<2, 3>::AverageElementSize(rGeometry) * ElementSizeCalculator<2, 3>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
                };
                break;
            case GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4:
                rVolumeDerivativeMethodType = [](CUInt NodeIndex, CUInt DirectionIndex,  const GeometryType& rGeometry) {
                    return 2.0 * ElementSizeCalculator<2, 4>::AverageElementSize(rGeometry) * ElementSizeCalculator<2, 4>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
                };
                break;
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                rVolumeDerivativeMethodType = [](CUInt NodeIndex, CUInt DirectionIndex,  const GeometryType& rGeometry) {
                    return 3.0 * std::pow(ElementSizeCalculator<3, 4>::AverageElementSize(rGeometry), 2) * ElementSizeCalculator<3, 4>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
                };
                break;
            case GeometryData::KratosGeometryType::Kratos_Prism3D6:
                rVolumeDerivativeMethodType = [](CUInt NodeIndex, CUInt DirectionIndex,  const GeometryType& rGeometry) {
                    return 3.0 * std::pow(ElementSizeCalculator<3, 6>::AverageElementSize(rGeometry), 2) * ElementSizeCalculator<3, 6>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
                };
                break;
            case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8:
                rVolumeDerivativeMethodType = [](CUInt NodeIndex, CUInt DirectionIndex,  const GeometryType& rGeometry) {
                    return 3.0 * std::pow(ElementSizeCalculator<3, 8>::AverageElementSize(rGeometry), 2) * ElementSizeCalculator<3, 8>::AverageElementSizeDerivative(NodeIndex, DirectionIndex, rGeometry);
                };
                break;
            default:
                KRATOS_ERROR << "Non supported geometry type." << std::endl;
        }

        for (SizeType c = 0; c < r_geometry.PointsNumber(); ++c) {
            auto& r_derivative_value = r_geometry[c].FastGetSolutionStepValue(rDerivativeVariable);

            for (SizeType k = 0; k < dimension; ++k) {
                const double derivative_value = rVolumeDerivativeMethodType(c, k, r_geometry);
                AtomicAdd(r_derivative_value[k], derivative_value);
            }
        }
    });

    mrModelPart.GetCommunicator().AssembleCurrentData(rDerivativeVariable);

    KRATOS_CATCH("");
}

void GeometryUtilities::CalculateNodalAreasFromConditions()
{
    KRATOS_TRY

    this->CalculateAreaNormalsFromConditions();

    //resetting the nodal area
    VariableUtils().SetHistoricalVariableToZero(NODAL_AREA, mrModelPart.Nodes());

    //calculating the normals and summing up at nodes
    block_for_each(mrModelPart.Nodes(), [&](Node& rNode)
    {
        const array_1d<double,3> normal = rNode.FastGetSolutionStepValue(NORMAL);
        rNode.FastGetSolutionStepValue(NODAL_AREA) = MathUtils<double>::Norm3(normal);

    });

    KRATOS_CATCH("")
}

void GeometryUtilities::CalculateAreaNormalsFromConditions()
{
    KRATOS_TRY

    //resetting the normals
    VariableUtils().SetHistoricalVariableToZero(NORMAL, mrModelPart.Nodes());

    //calculating the normals and summing up at nodes
    const array_1d<double,3> local_coords = ZeroVector(3);
    block_for_each(mrModelPart.Conditions(), [&](Condition& rCond)
    {
        auto& r_geometry = rCond.GetGeometry();

        const array_1d<double,3> normal = r_geometry.Normal(local_coords);
        const double coeff = 1.0/r_geometry.size();

        for(auto& node_i : r_geometry)
        {
            node_i.SetLock();
            noalias(node_i.FastGetSolutionStepValue(NORMAL)) += coeff * normal;
            node_i.UnSetLock();
        }
    });

    KRATOS_CATCH("")
}

void GeometryUtilities::CalculateUnitNormals()
{
    for (auto& node_i : mrModelPart.Nodes())
    {
        const array_1d<double,3>& area_normal = node_i.FastGetSolutionStepValue(NORMAL);
        array_3d& normalized_normal = node_i.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);

        const double norm2 = norm_2(area_normal);
        KRATOS_ERROR_IF(norm2<1e-10) << "CalculateUnitNormals: Norm2 of normal for node "
            << node_i.Id() << " is < 1e-10!" << std::endl;

        noalias(normalized_normal) = area_normal/norm2;
    }
}

double GeometryUtilities::CalculateAverageElementSize()
{

    double average_size = 0.0;

    for (auto& r_condition : mrModelPart.Conditions())
    {
        auto& r_geometry = r_condition.GetGeometry();

        average_size += r_geometry.Length();
        // const auto& geometry_type = r_geometry.GetGeometryType();


        // switch (geometry_type) {
        //     case GeometryData::KratosGeometryType::Kratos_Triangle3D3:
        //         average_size += ElementSizeCalculator<3, 3>::AverageElementSize(r_geometry);
        //     case GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4:
        //         average_size += ElementSizeCalculator<3, 4>::AverageElementSize(r_geometry);
        //     default:
        //         KRATOS_ERROR << "Non supported geometry type." << std::endl;
        // }
    }

    average_size /= mrModelPart.Conditions().size();

    return average_size;
}

void GeometryUtilities::CalculateGaussianCurvature() {

    KRATOS_TRY

    // check if surface normals are necessary (only for Taubin method)
    bool normals_necessary = false;
    for (const auto& rNode : mrModelPart.Nodes()) {
        std::string technique = this->GetCurvatureTechnique(rNode);
        if (technique == "Taubin") {
            normals_necessary = true;
            break;
        }
    }
    if (normals_necessary) this->ComputeUnitSurfaceNormals();

    // check if edge modelpart is necessary (only for Meyer method)
    bool edge_necessary = false;
    for (const auto& rNode : mrModelPart.Nodes()) {
        std::string technique = this->GetCurvatureTechnique(rNode);
        if (technique == "Meyer") {
            edge_necessary = true;
            break;
        }
    }
    if (edge_necessary) {
        ModelPart& r_edge_model_part = mrModelPart.HasSubModelPart(mrModelPart.Name() + "_edges") ? mrModelPart.GetSubModelPart(mrModelPart.Name() + "_edges") : mrModelPart.CreateSubModelPart(mrModelPart.Name() + "_edges");
        if (r_edge_model_part.Nodes().size() == 0) {
            this->ExtractEdgeNodes(mrModelPart.Name() + "_edges");
        }
    }

    block_for_each(mrModelPart.Nodes(), [&](NodeType &rNode) {

        double& r_gaussian_curvature = rNode.FastGetSolutionStepValue(GAUSSIAN_CURVATURE);

        std::string technique = this->GetCurvatureTechnique(rNode);

        // For quadratic elements (e.g. 6NTriangle)
        if (technique == "curvature_tensor") {
            r_gaussian_curvature = this->GaussianCurvatureForNodeFromTensor(rNode);
        }
        // For 3NTriangles by Meyer
        else if (technique == "Meyer") {
            r_gaussian_curvature = this->GaussianCurvatureForNodeMeyer(rNode);
        }
        // For 4NQuads by Taubin, variable naming according to Camprubi Estebo
        else if (technique == "Taubin") {
            r_gaussian_curvature = this->GaussianCurvatureForNodeTaubin(rNode);
        }
    });

    // Average gaussian curvature for edge nodes from neighbours (in case of 3NTriangle, Meyer)
    block_for_each(mrModelPart.Nodes(), [&](NodeType &rNode) {
        std::string technique = this->GetCurvatureTechnique(rNode);
        if (technique == "Meyer") {
            ModelPart& r_edge_model_part = mrModelPart.GetSubModelPart(mrModelPart.Name() + "_edges");
            const auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES).GetContainer();
            double& r_gaussian_curvature = rNode.FastGetSolutionStepValue(GAUSSIAN_CURVATURE);
            int counter = 0;
            if (r_edge_model_part.HasNode(rNode.Id())) {
                for (const auto& r_neighbour : r_neighbours) {
                    if (!r_edge_model_part.HasNode(r_neighbour->Id())) {
                        r_gaussian_curvature += r_neighbour->FastGetSolutionStepValue(GAUSSIAN_CURVATURE);
                        ++counter;
                    }
                }
                r_gaussian_curvature /= counter;
            }
        }
    });

    KRATOS_CATCH("");
}

double GeometryUtilities::GaussianCurvatureForNodeFromTensor(const NodeType &rNode) {

    const auto& r_condition_neighbours = rNode.GetValue(NEIGHBOUR_CONDITIONS).GetContainer();

    bool local_cartesian_basis_defined = false;
    Vector e_1 = ZeroVector(3);
    Vector e_2 = ZeroVector(3);
    Matrix b_total = ZeroMatrix(2,2);
    int counter = 0;
    for (const auto& r_condition_neighbour : r_condition_neighbours) {
        if (GeometryUtilities::CheckIfElementIsQuadratic(r_condition_neighbour)) {
            Matrix b = GeometryUtilities::CurvatureTensor(rNode, r_condition_neighbour);
            if (!local_cartesian_basis_defined) {
                GeometryUtilities::CartesianBaseVectors(rNode, r_condition_neighbour, e_1, e_2);
                local_cartesian_basis_defined = true;
            }
            Vector g_1 = ZeroVector(3);
            Vector g_2 = ZeroVector(3);
            GeometryUtilities::BaseVectors(rNode, r_condition_neighbour, g_1, g_2);
            Matrix b_transformed = ZeroMatrix(2,2);
            GeometryUtilities::TransformTensorCoefficients(b, b_transformed, g_1, g_2, e_1, e_2);
            b_total += b_transformed;
            ++counter;
        }
    }
    b_total /= counter;
    double gaussian_curvature = MathUtils<double>::Det2(b_total);

    return gaussian_curvature;
}

double GeometryUtilities::GaussianCurvatureForNodeMeyer(const NodeType &rNode) {

    const auto& r_condition_neighbours = rNode.GetValue(NEIGHBOUR_CONDITIONS).GetContainer();

    ModelPart& r_edge_model_part = mrModelPart.GetSubModelPart(mrModelPart.Name() + "_edges");
    if (r_edge_model_part.HasNode(rNode.Id())) {
        return 0;
    } else {
        double gaussian_curvature = 2*Globals::Pi;
        double A_mixed = 0.0;
        for (const auto& r_condition_neighbour : r_condition_neighbours) {
            double element_angle = 0;
            double element_a_mixed = 0;
            GeometryUtilities::InnerAngleAndMixedAreaOf3D3NTriangletAtNode(rNode, r_condition_neighbour, element_angle, element_a_mixed);

            gaussian_curvature -= element_angle;
            A_mixed += element_a_mixed;
        }
        gaussian_curvature /= A_mixed;
        return gaussian_curvature;
    }
}

double GeometryUtilities::GaussianCurvatureForNodeTaubin(const NodeType &rNode) {

    const auto& r_condition_neighbours = rNode.GetValue(NEIGHBOUR_CONDITIONS).GetContainer();

    Vector A3_p = rNode.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
    Vector A1_p = ZeroVector(3);
    Vector A2_p = ZeroVector(3);
    MathUtils<double>::OrthonormalBasis(A3_p, A1_p, A2_p);
    Vector curvature_tensor = ZeroVector(3);
    double area_sum = 0;
    for (const auto& r_condition_neighbour : r_condition_neighbours) {
        if (r_condition_neighbour->GetGeometry().GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4 ) continue;
        Vector kappa_Qi(3);
        Vector a_i(3);
        Vector b_i(3);
        Matrix A(3,3);
        int i = 0;
        for (const auto& r_element_node : r_condition_neighbour->GetGeometry()) {
            if (r_element_node.Id() == rNode.Id()) continue;
            Vector r_i = r_element_node - rNode.Coordinates();
            Vector t_i = r_i - A3_p * (MathUtils<double>::Dot3(A3_p, r_i));
            t_i /= MathUtils<double>::Norm3(t_i);
            kappa_Qi(i) = 2 * MathUtils<double>::Dot3(A3_p, r_i) / (MathUtils<double>::Norm3(r_i)*MathUtils<double>::Norm3(r_i));
            a_i(i) = MathUtils<double>::Dot3(t_i, A1_p);
            b_i(i) = MathUtils<double>::Dot3(t_i, A2_p);
            A(i, 0) = a_i(i) * a_i(i);
            A(i, 1) = 2 * a_i(i) * b_i(i);
            A(i, 2) = b_i(i) * b_i(i);
            ++i;
        }
        Vector curvature_tensor_element = ZeroVector(3);
        MathUtils<double>::Solve(A, curvature_tensor_element, kappa_Qi);
        area_sum += r_condition_neighbour->GetGeometry().Area();
        curvature_tensor += r_condition_neighbour->GetGeometry().Area() * curvature_tensor_element;
    }

    double gaussian_curvature = 0;
    if (area_sum > 0) {
        curvature_tensor /= area_sum;
        gaussian_curvature = curvature_tensor(0) * curvature_tensor(2) - curvature_tensor(1) * curvature_tensor(1);
    }
    return gaussian_curvature;
}

std::string GeometryUtilities::GetCurvatureTechnique(const NodeType &rNode) {

    const auto& r_condition_neighbours = rNode.GetValue(NEIGHBOUR_CONDITIONS).GetContainer();

    // Precedence: Curvature tensor > Taubin > Meyer
    if (this->CheckIfNodesHasQuadraticNeigbourElement(rNode)) {
        return "curvature_tensor";
    } else {
        for (const auto& r_condition_neighbour : r_condition_neighbours) {
            if (r_condition_neighbour->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4) {
                return "Taubin";
            }
        }
        return "Meyer";
    }
}

bool GeometryUtilities::CheckIfElementIsQuadratic(const Kratos::GlobalPointer<Kratos::Condition> pElement) {
    if (pElement->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D6) {
        return true;
    } else if (pElement->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8) {
        return true;
    } else if (pElement->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9) {
        return true;
    } else {
        return false;
    }
}

bool GeometryUtilities::CheckIfNodesHasQuadraticNeigbourElement(const NodeType &rNode) {
    const auto& r_condition_neighbours = rNode.GetValue(NEIGHBOUR_CONDITIONS).GetContainer();
    for (const auto& r_condition_neighbour : r_condition_neighbours) {
        if (this->CheckIfElementIsQuadratic(r_condition_neighbour)) {
            return true;
        }
    }
    return false;
}

void GeometryUtilities::LocalPointInElement(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement, Kratos::Point::CoordinatesArrayType& rLocalPoint) {
    Matrix local_points;
    pElement->GetGeometry().PointsLocalCoordinates(local_points);
    for (SizeType i = 0; i < pElement->GetGeometry().PointsNumber(); ++i) {
        if (pElement->GetGeometry()[i].Id() == rNode.Id()) {
            rLocalPoint[0] = local_points(i, 0);
            rLocalPoint[1] = local_points(i, 1);
            break;
        }
    }
}

void GeometryUtilities::BaseVectors(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement, Vector& rG1, Vector& rG2) {

    Kratos::Point::CoordinatesArrayType LocalPoint(2);
    GeometryUtilities::LocalPointInElement(rNode, pElement, LocalPoint);

    matrix<double> first_derivative;
    pElement->GetGeometry().ShapeFunctionsLocalGradients(first_derivative, LocalPoint);
    Vector g_1 = ZeroVector(3);
    Vector g_2 = ZeroVector(3);
    for (SizeType i = 0; i < pElement->GetGeometry().PointsNumber(); ++i) {
        g_1 += first_derivative(i, 0) * pElement->GetGeometry()[i].Coordinates();
        g_2 += first_derivative(i, 1) * pElement->GetGeometry()[i].Coordinates();
    }

    rG1 = g_1;
    rG2 = g_2;
}

void GeometryUtilities::CartesianBaseVectors(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement, Vector& rE1, Vector& rE2) {

    Vector g_1 = ZeroVector(3);
    Vector g_2 = ZeroVector(3);
    GeometryUtilities::BaseVectors(rNode, pElement, g_1, g_2);

    // local cartesian basis e_i
    Vector e_1 = g_1 / MathUtils<double>::Norm3(g_1);
    Vector e_2 = g_2 - MathUtils<double>::Dot3(g_2, e_1) * e_1;
    e_2 /= MathUtils<double>::Norm3(e_2);

    rE1 = e_1;
    rE2 = e_2;
}

Matrix GeometryUtilities::CurvatureTensor(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement) {
    Kratos::Point::CoordinatesArrayType LocalPoint(2);
    GeometryUtilities::LocalPointInElement(rNode, pElement, LocalPoint);

    Vector g_1 = ZeroVector(3);
    Vector g_2 = ZeroVector(3);
    GeometryUtilities::BaseVectors(rNode, pElement, g_1, g_2);

    GeometryData::ShapeFunctionsSecondDerivativesType second_derivative;
    pElement->GetGeometry().ShapeFunctionsSecondDerivatives(second_derivative, LocalPoint);

    Vector dg_1_du = ZeroVector(3);
    Vector dg_1_dv = ZeroVector(3);
    Vector dg_2_du = ZeroVector(3);
    Vector dg_2_dv = ZeroVector(3);

    for (SizeType i = 0; i < pElement->GetGeometry().PointsNumber(); ++i) {
        dg_1_du += second_derivative[i](0,0) * pElement->GetGeometry()[i].Coordinates();
        dg_1_dv += second_derivative[i](0,1) * pElement->GetGeometry()[i].Coordinates();
        dg_2_du += second_derivative[i](1,0) * pElement->GetGeometry()[i].Coordinates();
        dg_2_dv += second_derivative[i](1,1) * pElement->GetGeometry()[i].Coordinates();
    }

    Vector g_3 = ZeroVector(3);
    g_3 = MathUtils<double>::CrossProduct(g_1, g_2);
    g_3 *=  1 / MathUtils<double>::Norm3(g_3);

    // curvature tensor coefficients
    Matrix b(2,2);
    b(0, 0) = MathUtils<double>::Dot3(dg_1_du, g_3);
    b(1, 0) = MathUtils<double>::Dot3(dg_2_du, g_3);
    b(0, 1) = MathUtils<double>::Dot3(dg_1_dv, g_3);
    b(1, 1) = MathUtils<double>::Dot3(dg_2_dv, g_3);

    return b;
}

void GeometryUtilities::TransformTensorCoefficients(Matrix& rTensor, Matrix& rResultTensor, Vector&rG1, Vector&rG2, Vector&rE1, Vector&rE2) {

    Vector contra_G1 = ZeroVector(3);
    Vector contra_G2 = ZeroVector(3);
    Matrix covariant_metric(2,2);
    covariant_metric(0,0) = MathUtils<double>::Dot3(rG1, rG1);
    covariant_metric(1,0) = MathUtils<double>::Dot3(rG2, rG1);
    covariant_metric(0,1) = MathUtils<double>::Dot3(rG1, rG2);
    covariant_metric(1,1) = MathUtils<double>::Dot3(rG2, rG2);
    Matrix contravariant_metric(2,2);
    double det;
    MathUtils<double>::InvertMatrix2(covariant_metric, contravariant_metric, det);
    contra_G1 = contravariant_metric(0,0) * rG1 + contravariant_metric(1,0) * rG2;
    contra_G2 = contravariant_metric(0,1) * rG1 + contravariant_metric(1,1) * rG2;

    rResultTensor(0,0) = rTensor(0, 0) * MathUtils<double>::Dot3(rE1, contra_G1) * MathUtils<double>::Dot3(contra_G1, rE1);
    rResultTensor(0,0) += rTensor(1, 0) * MathUtils<double>::Dot3(rE1, contra_G2) * MathUtils<double>::Dot3(contra_G1, rE1);
    rResultTensor(0,0) += rTensor(0, 1) * MathUtils<double>::Dot3(rE1, contra_G1) * MathUtils<double>::Dot3(contra_G2, rE1);
    rResultTensor(0,0) += rTensor(1, 1) * MathUtils<double>::Dot3(rE1, contra_G2) * MathUtils<double>::Dot3(contra_G2, rE1);

    rResultTensor(1,0) = rTensor(0, 0) * MathUtils<double>::Dot3(rE2, contra_G1) * MathUtils<double>::Dot3(contra_G1, rE1);
    rResultTensor(1,0) += rTensor(1, 0) * MathUtils<double>::Dot3(rE2, contra_G2) * MathUtils<double>::Dot3(contra_G1, rE1);
    rResultTensor(1,0) += rTensor(0, 1) * MathUtils<double>::Dot3(rE2, contra_G1) * MathUtils<double>::Dot3(contra_G2, rE1);
    rResultTensor(1,0) += rTensor(1, 1) * MathUtils<double>::Dot3(rE2, contra_G2) * MathUtils<double>::Dot3(contra_G2, rE1);

    rResultTensor(0,1) = rTensor(0, 0) * MathUtils<double>::Dot3(rE1, contra_G1) * MathUtils<double>::Dot3(contra_G1, rE2);
    rResultTensor(0,1) += rTensor(1, 0) * MathUtils<double>::Dot3(rE1, contra_G2) * MathUtils<double>::Dot3(contra_G1, rE2);
    rResultTensor(0,1) += rTensor(0, 1) * MathUtils<double>::Dot3(rE1, contra_G1) * MathUtils<double>::Dot3(contra_G2, rE2);
    rResultTensor(0,1) += rTensor(1, 1) * MathUtils<double>::Dot3(rE1, contra_G2) * MathUtils<double>::Dot3(contra_G2, rE2);

    rResultTensor(1,1) = rTensor(0, 0) * MathUtils<double>::Dot3(rE2, contra_G1) * MathUtils<double>::Dot3(contra_G1, rE2);
    rResultTensor(1,1) += rTensor(1, 0) * MathUtils<double>::Dot3(rE2, contra_G2) * MathUtils<double>::Dot3(contra_G1, rE2);
    rResultTensor(1,1) += rTensor(0, 1) * MathUtils<double>::Dot3(rE2, contra_G1) * MathUtils<double>::Dot3(contra_G2, rE2);
    rResultTensor(1,1) += rTensor(1, 1) * MathUtils<double>::Dot3(rE2, contra_G2) * MathUtils<double>::Dot3(contra_G2, rE2);
}

void GeometryUtilities::InnerAngleAndMixedAreaOf3D3NTriangletAtNode(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement, double& rInnerAngle, double& rMixedArea) {

    const auto& node_i = rNode.Coordinates();

    Kratos::Point::CoordinatesArrayType node_j;
    Kratos::Point::CoordinatesArrayType node_k;
    if (pElement->GetGeometry()[0].Id() == rNode.Id()) {
        node_j = pElement->GetGeometry()[1].Coordinates();
        node_k = pElement->GetGeometry()[2].Coordinates();
    } else if (pElement->GetGeometry()[1].Id() == rNode.Id()) {
        node_j = pElement->GetGeometry()[2].Coordinates();
        node_k = pElement->GetGeometry()[0].Coordinates();
    } else if (pElement->GetGeometry()[2].Id() == rNode.Id()) {
        node_j = pElement->GetGeometry()[0].Coordinates();
        node_k = pElement->GetGeometry()[1].Coordinates();
    }

    const Kratos::Point::CoordinatesArrayType edge_ij = node_j - node_i;
    const Kratos::Point::CoordinatesArrayType edge_ik = node_k - node_i;;
    const double cos_theta_i = inner_prod(edge_ij, edge_ik) / (norm_2(edge_ij) * norm_2(edge_ik));
    const double theta_i = acos(cos_theta_i);
    rInnerAngle = theta_i;

    const Kratos::Point::CoordinatesArrayType edge_jk = node_k - node_j;
    const double theta_j = acos(inner_prod(-edge_ij, edge_jk) / (norm_2(-edge_ij) * norm_2(edge_jk)));
    const double theta_k = acos(inner_prod(-edge_ik, -edge_jk) / (norm_2(-edge_ik) * norm_2(-edge_jk)));

    // triangle is obtuse => 1/2 or 1/4 * triangle area
    if (theta_i > Globals::Pi / 2 || theta_j > Globals::Pi / 2 || theta_k > Globals::Pi / 2) {
        const double a = norm_2(edge_ij);
        const double b = norm_2(edge_ik);
        const double c = norm_2(edge_jk);
        const double s = (a+b+c) / 2.0;
        const double area = std::sqrt(s*(s-a)*(s-b)*(s-c));
        if (theta_i > Globals::Pi / 2) {
            rMixedArea += 0.5 * area;
        } else {
            rMixedArea += 0.25 * area;
        }
    }
    // triangle is non-obtuse => Voronoi area
    else {
        rMixedArea += 0.125 * (norm_2_square(edge_ij) * (cos(theta_k) / sin(theta_k)) + norm_2_square(edge_ik) * (cos(theta_j) / sin(theta_j)));
    }
}

}  // namespace Kratos.
