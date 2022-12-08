// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    David Schmölz
//                   Suneth Warnakulasuriya
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <string>
#include <unordered_map>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/global_pointer_variables.h"
#include "includes/model_part.h"
#include "processes/find_global_nodal_neighbours_for_entities_process.h"
#include "processes/find_global_nodal_entity_neighbours_process.h"
#include "shape_optimization_application.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"
#include "utilities/pointer_communicator.h"
#include "utilities/dense_householder_qr_decomposition.h"

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include "external_libraries/igl/include/igl/gaussian_curvature.h"
#include "external_libraries/igl/include/igl/massmatrix.h"
#include "external_libraries/igl/include/igl/invert_diag.h"
#include "external_libraries/igl/include/igl/writeOFF.h"

// ------------------------------------------------------------------------------
// Base mapper includes
// ------------------------------------------------------------------------------
#include "mapper_vertex_morphing.h"
#include "mapper_vertex_morphing_improved_integration.h"
#include "mapper_vertex_morphing_matrix_free.h"
#include "mapper_vertex_morphing_symmetric.h"

// ------------------------------------------------------------------------------
// Base class include
// ------------------------------------------------------------------------------
#include "mapper_vertex_morphing_adaptive_radius.h"

// ==============================================================================

namespace Kratos
{

template <class TBaseVertexMorphingMapper>
MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::MapperVertexMorphingAdaptiveRadius(
    ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart,
    Parameters MapperSettings)
    : BaseType(rOriginModelPart, rDestinationModelPart, MapperSettings),
        mrOriginModelPart(rOriginModelPart),
        mrDestinationModelPart(rDestinationModelPart),
        // TODO: delete filter_radius_factor
        mFilterRadiusFactor(MapperSettings["adaptive_filter_settings"]["filter_radius_factor"].GetDouble()),
        mMinimumFilterRadius(MapperSettings["adaptive_filter_settings"]["minimum_filter_radius"].GetDouble()),
        mNumberOfSmoothingIterations(MapperSettings["adaptive_filter_settings"]["filter_radius_smoothing_iterations"].GetInt()),
        mRadiusFunctionType(MapperSettings["adaptive_filter_settings"]["radius_function"].GetString()),
        mRadiusFunctionParameter(MapperSettings["adaptive_filter_settings"]["radius_function_parameter"].GetDouble()),
        mCurvatureLimit(MapperSettings["adaptive_filter_settings"]["curvature_limit"].GetDouble()),
        mMaxNumberOfNeighbors(MapperSettings["max_nodes_in_filter_radius"].GetInt())
{
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::Initialize()
{
    BaseType::Initialize();
    CalculateAdaptiveVertexMorphingRadius();

    KRATOS_INFO("ShapeOpt") << "filter_radius_factor:  " << mFilterRadiusFactor  << std::endl;
    KRATOS_INFO("ShapeOpt") << "minimum_filter_radius:  " << mMinimumFilterRadius  << std::endl;
    KRATOS_INFO("ShapeOpt") << "filter_radius_smoothing_iterations:  " << mNumberOfSmoothingIterations  << std::endl;

}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::Update()
{
    BaseType::Update();
    CalculateAdaptiveVertexMorphingRadius();
    // for(auto& node_i : mrDestinationModelPart.Nodes()) {
    //     CalculateCurvatureAtNode(node_i);
    // }
}

template <class TBaseVertexMorphingMapper>
std::string MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::Info() const
{
    return BaseType::Info() + "AdaptiveRadius";
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::PrintInfo(std::ostream &rOStream) const
{
    rOStream << BaseType::Info() << "AdaptiveRadius";
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::PrintData(std::ostream &rOStream) const
{
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateCurvatureIGL()
{
    // count number of quad elements
    // TODO: deal with quads 4, 8, 9 nodes
    // TODO: deal with triangles 6 nodes
    // int numQuadElms = 0;
    // for (int i = 0; i < mrDestinationModelPart.NumberOfElements; ++i) {
    //     ConditionType::Pointer condition_i = mrDestinationModelPart.mpConditions[i];
    //     if (condition_i -> NumberOfNodes() == 4) numQuadElms++;
    // }

    // Load a mesh in OFF format
    Eigen::MatrixXd V(mrDestinationModelPart.NumberOfNodes(), 3);
    // Eigen::MatrixXi F(mrDestinationModelPart->mNumElems + numQuadElms, 3);
    unsigned int i = 0;
    for (auto& node_i : mrDestinationModelPart.Nodes()) {
        V(i,0) = node_i.X();
        V(i,1) = node_i.Y();
        V(i,2) = node_i.Z();
        ++i;
    }
    KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "Size of V : " << V.size() << "..." << std::endl;
    // for (int i = 0; i < mrDestinationModelPart.NumberOfNodes(); ++i) {
    //     NodeType::Pointer node_i = mrDestinationModelPart.GetNode()[i];
    //     V(i,0) = node_i->X();
    //     V(i,1) = node_i->Y();
    //     V(i,2) = node_i->Z();
    // }

    // int counter = 0;
    Eigen::MatrixXi F(mrDestinationModelPart.NumberOfConditions(), 3);
    unsigned int j = 0;
    for (auto& condition_j : mrDestinationModelPart.Conditions()) {
        NodeType& r_node_1 = condition_j.GetGeometry()[0];
        NodeType& r_node_2 = condition_j.GetGeometry()[1];
        NodeType& r_node_3 = condition_j.GetGeometry()[2];
        // KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "node_1 : " << r_node_1 << "..." << std::endl;
        // KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "node_2 : " << r_node_2 << "..." << std::endl;
        // KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "node_3 : " << r_node_3 << "..." << std::endl;
        // KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "node_1_index : " << r_node_1.GetValue(MAPPING_ID) << "..." << std::endl;
        // KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "node_2_index : " << r_node_2.GetValue(MAPPING_ID) << "..." << std::endl;
        // KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "node_3_index : " << r_node_3.GetValue(MAPPING_ID) << "..." << std::endl;
        F(j,0) = r_node_1.GetValue(MAPPING_ID);
        F(j,1) = r_node_2.GetValue(MAPPING_ID);
        F(j,2) = r_node_3.GetValue(MAPPING_ID);
        ++j;
        // counter ++;
    }
    KRATOS_INFO("ShapeOpt:::ADAPTIVE VM") << "Size of F : " << F.size() << "..." << std::endl;
    Eigen::VectorXd K;
    igl::gaussian_curvature(V,F,K);
    Eigen::SparseMatrix<double> M,Minv;
    // SparseMatrixType M,Minv;
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
    igl::invert_diag(M,Minv);
    K = (Minv*K).eval();

    igl::writeOFF("test.off",V,F);

    i = 0;
    for (auto& node_i : mrDestinationModelPart.Nodes()) {
        double& gaussian_curvature_i = node_i.FastGetSolutionStepValue(GAUSSIAN_CURVATURE);
        gaussian_curvature_i = K[i];
        ++i;
    }
    // // if quad elements in mesh, add a second triangle
    // for (int j = 0; j < mpPart->mNumElems; ++j) {
    //     Element::ConstPointer element_j = mpPart->mElements[j];
    //     if (element_j -> NumberOfNodes() == 4) {
    //         F(counter,0) = mpPart->GetNodeIndex(element_j -> mNodes[0]->mId);
    //         F(counter,1) = mpPart->GetNodeIndex(element_j -> mNodes[2]->mId);
    //         F(counter,2) = mpPart->GetNodeIndex(element_j -> mNodes[3]->mId);
    //         counter ++;
    //     }
    // }
}


template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateCurvatureAtNode(NodeType &rNode)
{
    KRATOS_TRY

    const auto& r_node_i = rNode.Coordinates();
    const auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES).GetContainer();
    const auto& r_condition_neighbours = rNode.GetValue(NEIGHBOUR_CONDITIONS).GetContainer();

    double& r_gaussian_curvature = rNode.FastGetSolutionStepValue(GAUSSIAN_CURVATURE);

    // TODO: ACTIVATE AGAIN / FIND A WAY TO DEAL WITH EDGE NODES!
    // check if node lies on edge
    // TODO: for all element types
    // for (const auto& r_neighbour : r_neighbours) {
    //     int count = 0;
    //     for (const auto& r_condition_neighbour : r_condition_neighbours) {
    //         for (int i = 0; i < r_condition_neighbour->GetGeometry().PointsNumber(); ++i) {
    //             if (r_condition_neighbour->GetGeometry()[i].Id() == r_neighbour->Id()) count++;
    //         }
    //     }
    //     if (count < 2) {
    //         r_gaussian_curvature = 0;
    //         return;
    //     }
    // }

    // r_gaussian_curvature = 2*Globals::Pi;
    r_gaussian_curvature = 0;
    double A_mixed = 0.0;

    for (const auto& r_condition_neighbour : r_condition_neighbours) {
        // TODO: find way of dealing with quadratic/bi-quadratic elements
        // KRATOS_INFO("ShapeOpt") << "KratosGeometryType = " << r_condition_neighbour->GetGeometry().GetGeometryType() << std::endl;
        if (r_condition_neighbour->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D6) {
            r_gaussian_curvature += this->CalculateCurvatureOf6Node3DTriangletAtNode(rNode, r_condition_neighbour);// * r_condition_neighbour->GetGeometry().Area();
            // A_mixed += r_condition_neighbour->GetGeometry().Area();
            ++A_mixed;
            // KRATOS_INFO("ShapeOpt") << "r_gaussian_curvature at node id = " << rNode.Id() << " : " << r_gaussian_curvature << std::endl;
        } else {
            double element_angle = 0;
            double element_a_mixed = 0;
            this->CalculateInnerAngleAndMixedAreaOfElementAtNode(rNode, r_condition_neighbour, element_angle, element_a_mixed);

            r_gaussian_curvature -= element_angle;
            A_mixed += element_a_mixed;
        }
    }

    // KRATOS_INFO("ShapeOpt") << "Final 2pi-inner angle at node id = " << rNode.Id() << " : " << r_gaussian_curvature << std::endl;
    r_gaussian_curvature /= A_mixed;
    // KRATOS_INFO("ShapeOpt") << "Final mixed area at node id = " << rNode.Id() << " : " << A_mixed << std::endl;

    KRATOS_CATCH("");
}

template <class TBaseVertexMorphingMapper>
double MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::Cotan(const double& rAngle) {
    return cos(rAngle) / sin(rAngle);
}

template <class TBaseVertexMorphingMapper>
double MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CurvatureFunction(const double& rCurvature, const double& rElementSize) {

    // Using the function delta_K = a / (kappa + b)
    // two equations for parameter a and b:
    // (I) a / b = delta_K(r_min)
    // (II) a / (kappa_limit + b) = mRadiusFunctionParameter
    // This is based on the "analytic" relation (found out by numerical experiment) between curvature change and filter radius: delta_K ~ 4 / r**4
    // => r = pow(4 / delta_K, 0.25)
    if (mRadiusFunctionType == "analytic")
    {
        double delta_K_max = 4 / pow(mMinimumFilterRadius/rElementSize, 4);

        double b = (mRadiusFunctionParameter * mCurvatureLimit) / (delta_K_max - mRadiusFunctionParameter);
        double a = b * delta_K_max;
        double delta_K = a / (abs(rCurvature) + b);

        double filter_radius = (delta_K > 0) ? rElementSize * pow(4/delta_K, 0.25) : mMinimumFilterRadius;
        return filter_radius;
    }
    // Using linear function r = r_min + a * h * kappa
    else if (mRadiusFunctionType == "linear")
    {
        return mMinimumFilterRadius + mRadiusFunctionParameter * rElementSize * abs(rCurvature);
    }
    // Using square root function r = r_min + a * h * sqrt(kappa)
    else if (mRadiusFunctionType == "square_root") {
        return mMinimumFilterRadius + mRadiusFunctionParameter * rElementSize * sqrt(abs(rCurvature));
    }
    // Using fourth root function r = r_min + a * h * kappa**0.25
    else if (mRadiusFunctionType == "fourth_root") {
        return mMinimumFilterRadius + mRadiusFunctionParameter * rElementSize * pow(abs(rCurvature), 0.25);
    } else {
        KRATOS_ERROR << "ShapeOpt Adaptive Filter: Curvature function type " << mRadiusFunctionType << " not supported for adaptive filter vertex morphing method." << std::endl;
    }
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateInnerAngleAndMixedAreaOfElementAtNode(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement, double& rInnerAngle, double& rMixedArea) {

    if (pElement->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3) {
        this->CalculateInnerAngleAndMixedAreaOf3Node3DTriangletAtNode(rNode, pElement, rInnerAngle, rMixedArea);
    } else if (pElement->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D6) {
        this->CalculateInnerAngleAndMixedAreaOf6Node3DTriangletAtNode(rNode, pElement, rInnerAngle, rMixedArea);
    } else if (pElement->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4) {
        this->CalculateInnerAngleAndMixedAreaOf4Node3DQuadrilateraltAtNode(rNode, pElement, rInnerAngle, rMixedArea);
    } else if (pElement->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8) {
        this->CalculateInnerAngleAndMixedAreaOf8Node3DQuadrilateraltAtNode(rNode, pElement, rInnerAngle, rMixedArea);
    } else {
        KRATOS_ERROR << "ShapeOpt Adaptive Filter: Geometry type not supported for adaptive filter vertex morphing method." << std::endl;
    }
}

/**
 *      v
 *      ^
 *      |
 *      k
 *      |`\
 *      |  `\
 *      |    `\
 *      |      `\
 *      |        `\
 *      i----------j --> u
 */
template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateInnerAngleAndMixedAreaOfTriangleAtNodeI(const Kratos::Point::CoordinatesArrayType& rNodeI, const Kratos::Point::CoordinatesArrayType& rNodeJ, const Kratos::Point::CoordinatesArrayType& rNodeK, double& rInnerAngle, double& rMixedArea) {

    const Kratos::Point::CoordinatesArrayType edge_ij = rNodeJ - rNodeI;
    const Kratos::Point::CoordinatesArrayType edge_ik = rNodeK - rNodeI;;
    const double cos_theta_i = inner_prod(edge_ij, edge_ik) / (norm_2(edge_ij) * norm_2(edge_ik));
    const double theta_i = acos(cos_theta_i);
    rInnerAngle = theta_i;

    const Kratos::Point::CoordinatesArrayType edge_jk = rNodeK - rNodeJ;
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
        rMixedArea += 0.125 * (norm_2_square(edge_ij) * this->Cotan(theta_k) + norm_2_square(edge_ik) * this->Cotan(theta_j));
    }
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateInnerAngleAndMixedAreaOf3Node3DTriangletAtNode(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement, double& rInnerAngle, double& rMixedArea) {

    const auto& r_node_i = rNode.Coordinates();

    Kratos::Point::CoordinatesArrayType r_node_j;
    Kratos::Point::CoordinatesArrayType r_node_k;

    if (pElement->GetGeometry()[0].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[1].Coordinates();
        r_node_k = pElement->GetGeometry()[2].Coordinates();
    } else if (pElement->GetGeometry()[1].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[2].Coordinates();
        r_node_k = pElement->GetGeometry()[0].Coordinates();
    } else if (pElement->GetGeometry()[2].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[0].Coordinates();
        r_node_k = pElement->GetGeometry()[1].Coordinates();
    }
    this->CalculateInnerAngleAndMixedAreaOfTriangleAtNodeI(r_node_i, r_node_j, r_node_k, rInnerAngle, rMixedArea);
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateInnerAngleAndMixedAreaOf6Node3DTriangletAtNode(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement, double& rInnerAngle, double& rMixedArea) {

    const auto& r_node_i = rNode.Coordinates();

    Kratos::Point::CoordinatesArrayType r_node_j;
    Kratos::Point::CoordinatesArrayType r_node_k;
    Kratos::Point::CoordinatesArrayType r_node_corner_j;
    Kratos::Point::CoordinatesArrayType r_node_corner_k;

    bool corner_node = true;
    if (pElement->GetGeometry()[0].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[3].Coordinates();
        r_node_k = pElement->GetGeometry()[5].Coordinates();
    } else if (pElement->GetGeometry()[1].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[4].Coordinates();
        r_node_k = pElement->GetGeometry()[3].Coordinates();
    } else if (pElement->GetGeometry()[2].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[5].Coordinates();
        r_node_k = pElement->GetGeometry()[4].Coordinates();
    } else if (pElement->GetGeometry()[3].Id() == rNode.Id()) {
        corner_node = false;
        r_node_j = pElement->GetGeometry()[4].Coordinates();
        r_node_corner_j = pElement->GetGeometry()[1].Coordinates();
        r_node_k = pElement->GetGeometry()[5].Coordinates();
        r_node_corner_k = pElement->GetGeometry()[0].Coordinates();
    } else if (pElement->GetGeometry()[4].Id() == rNode.Id()) {
        corner_node = false;
        r_node_j = pElement->GetGeometry()[5].Coordinates();
        r_node_corner_j = pElement->GetGeometry()[2].Coordinates();
        r_node_k = pElement->GetGeometry()[3].Coordinates();
        r_node_corner_k = pElement->GetGeometry()[1].Coordinates();
    } else if (pElement->GetGeometry()[5].Id() == rNode.Id()) {
        corner_node = false;
        r_node_j = pElement->GetGeometry()[3].Coordinates();
        r_node_corner_j = pElement->GetGeometry()[0].Coordinates();
        r_node_k = pElement->GetGeometry()[4].Coordinates();
        r_node_corner_k = pElement->GetGeometry()[2].Coordinates();
    }

    if (corner_node) {
        this->CalculateInnerAngleAndMixedAreaOfTriangleAtNodeI(r_node_i, r_node_j, r_node_k, rInnerAngle, rMixedArea);
    } else {
        double inner_angle = 0.0;
        this->CalculateInnerAngleAndMixedAreaOfTriangleAtNodeI(r_node_i, r_node_j, r_node_k, inner_angle, rMixedArea);
        rInnerAngle += inner_angle;
        this->CalculateInnerAngleAndMixedAreaOfTriangleAtNodeI(r_node_i, r_node_corner_j, r_node_j, inner_angle, rMixedArea);
        rInnerAngle += inner_angle;
        this->CalculateInnerAngleAndMixedAreaOfTriangleAtNodeI(r_node_i, r_node_k, r_node_corner_k, inner_angle, rMixedArea);
        rInnerAngle += inner_angle;
    }
}

/**
 * The node ordering corresponds with:
 *          2
 *          |`\
 *          |  `\
 *          5    `4
 *          |      `\
 *          |        `\
 *          0-----3----1
 */
template <class TBaseVertexMorphingMapper>
double MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateCurvatureOf6Node3DTriangletAtNode(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement) {

    Kratos::Point::CoordinatesArrayType rLocalPoint;

    if (pElement->GetGeometry()[0].Id() == rNode.Id()) {
        rLocalPoint[0] = 0.0;
        rLocalPoint[1] = 0.0;
    } else if (pElement->GetGeometry()[1].Id() == rNode.Id()) {
        rLocalPoint[0] = 1.0;
        rLocalPoint[1] = 0.0;
    } else if (pElement->GetGeometry()[2].Id() == rNode.Id()) {
        rLocalPoint[0] = 0.0;
        rLocalPoint[1] = 1.0;
    } else if (pElement->GetGeometry()[3].Id() == rNode.Id()) {
        rLocalPoint[0] = 0.5;
        rLocalPoint[1] = 0.0;
    } else if (pElement->GetGeometry()[4].Id() == rNode.Id()) {
        rLocalPoint[0] = 0.5;
        rLocalPoint[1] = 0.5;
    } else if (pElement->GetGeometry()[5].Id() == rNode.Id()) {
        rLocalPoint[0] = 0.0;
        rLocalPoint[1] = 0.5;
    }

    matrix<double> first_derivative;
    pElement->GetGeometry().ShapeFunctionsLocalGradients(first_derivative, rLocalPoint);
    GeometryData::ShapeFunctionsSecondDerivativesType second_derivative;
    pElement->GetGeometry().ShapeFunctionsSecondDerivatives(second_derivative, rLocalPoint);
    matrix<double> df;
    Vector g_1(3);
    Vector g_2(3);
    Vector n(3);
    Vector n_tilde(3);
    double n_hat;
    Vector dn_tilde_du(3);
    Vector dn_tilde_dv(3);

    Vector dg_1_du(3);
    Vector dg_1_dv(3);
    Vector dg_2_du(3);
    Vector dg_2_dv(3);
    df.resize(3, 2, true);
    for (int i = 0; i < pElement->GetGeometry().PointsNumber(); ++i) {
        g_1 += first_derivative(i, 0) * pElement->GetGeometry()[i].Coordinates();
        g_2 += first_derivative(i, 1) * pElement->GetGeometry()[i].Coordinates();
        for (int dim = 0; dim < 3; ++dim) {
            df(dim, 0) += first_derivative(i, 0) * pElement->GetGeometry()[i].Coordinates()[dim];
            df(dim, 1) += first_derivative(i, 1) * pElement->GetGeometry()[i].Coordinates()[dim];
        }
        dg_1_du += second_derivative[i](0,0) * pElement->GetGeometry()[i].Coordinates();
        dg_1_dv += second_derivative[i](0,1) * pElement->GetGeometry()[i].Coordinates();
        dg_2_du += second_derivative[i](1,0) * pElement->GetGeometry()[i].Coordinates();
        dg_2_dv += second_derivative[i](1,1) * pElement->GetGeometry()[i].Coordinates();
    }
    n_tilde = MathUtils<double>::CrossProduct(g_1, g_2);
    n_hat = MathUtils<double>::Norm3(n_tilde);
    n = n_tilde / n_hat;
    // g_1 /= MathUtils<double>::Norm3(g_1);
    // g_2 /= MathUtils<double>::Norm3(g_2);

    Vector dn_du = MathUtils<double>::CrossProduct(dg_1_du, g_2) + MathUtils<double>::CrossProduct(g_1, dg_2_du);
    Vector dn_dv = MathUtils<double>::CrossProduct(dg_1_dv, g_2) + MathUtils<double>::CrossProduct(g_1, dg_2_dv);

    Matrix dn(3, 2);
    for (int i = 0; i < 3; ++i) {
        dn(i, 0) = dn_du(i);
        dn(i, 1) = dn_dv(i);
    }
    // solve least square problem df * S = dn for shape operator S [2x2] matrix
    // with direct QR solver
    DenseHouseholderQRDecomposition<UblasSpace<double, Matrix, Vector>> qr_decomposition;
    qr_decomposition.Compute(df);
    Matrix S(2, 2);
    qr_decomposition.Solve(dn, S);
    // eigenvalues of S => principal curvatures kappa_1 and kappa_2
    double m = (S(0,0) + S(1,1)) / 2;
    double p = S(0,0) * S(1,1) - S(0,1) * S(1,0);
    double lambda_1 = m + sqrt(m*m - p);
    double lambda_2 = m - sqrt(m*m - p);

    // curvature tensor coefficients (see Klingbeil)
    Matrix b(2,2);
    b(0, 0) = MathUtils<double>::Dot3(dg_1_du, n);
    b(1, 0) = MathUtils<double>::Dot3(dg_2_du, n);
    b(0, 1) = MathUtils<double>::Dot3(dg_1_dv, n);
    b(1, 1) = MathUtils<double>::Dot3(dg_2_dv, n);

    // transform curvature tensor coefficients to standardbasis
    Vector contra_g_1(3);
    Vector contra_g_2(3);
    Matrix covariant_metric(2,2);
    covariant_metric(0,0) = MathUtils<double>::Dot3(g_1, g_1);
    covariant_metric(1,0) = MathUtils<double>::Dot3(g_2, g_1);
    covariant_metric(0,1) = MathUtils<double>::Dot3(g_1, g_2);
    covariant_metric(1,1) = MathUtils<double>::Dot3(g_2, g_2);
    Matrix contravariant_metric = MathUtils<double>::InvertMatrix2(covariant_metric);
    contra_g_1 = contravariant_metric(0,0) * g_1 + contravariant_metric(0,1) * g_2;
    contra_g_2 = contravariant_metric(1,0) * g_1 + contravariant_metric(1,1) * g_2;
    Vector e_1(1, 0, 0);
    Vector e_2(0, 1, 0);
    Matrix b_e(2,2);
    b_e(0,0) = b(0, 0) * MathUtils<double>::Dot3(e_1, contra_g_1) * MathUtils<double>::Dot3(contra_g_1, e_1);
    b_e(0,0) += b(1, 0) * MathUtils<double>::Dot3(e_1, contra_g_2) * MathUtils<double>::Dot3(contra_g_1, e_1);
    b_e(0,0) += b(0, 1) * MathUtils<double>::Dot3(e_1, contra_g_1) * MathUtils<double>::Dot3(contra_g_2, e_1);
    b_e(0,0) += b(1, 1) * MathUtils<double>::Dot3(e_1, contra_g_2) * MathUtils<double>::Dot3(contra_g_2, e_1);

    b_e(1,0) = b(0, 0) * MathUtils<double>::Dot3(e_2, contra_g_1) * MathUtils<double>::Dot3(contra_g_1, e_1);
    b_e(1,0) += b(1, 0) * MathUtils<double>::Dot3(e_2, contra_g_2) * MathUtils<double>::Dot3(contra_g_1, e_1);
    b_e(1,0) += b(0, 1) * MathUtils<double>::Dot3(e_2, contra_g_1) * MathUtils<double>::Dot3(contra_g_2, e_1);
    b_e(1,0) += b(1, 1) * MathUtils<double>::Dot3(e_2, contra_g_2) * MathUtils<double>::Dot3(contra_g_2, e_1);

    b_e(0,1) = b(0, 0) * MathUtils<double>::Dot3(e_1, contra_g_1) * MathUtils<double>::Dot3(contra_g_1, e_2);
    b_e(0,1) += b(1, 0) * MathUtils<double>::Dot3(e_1, contra_g_2) * MathUtils<double>::Dot3(contra_g_1, e_2);
    b_e(0,1) += b(0, 1) * MathUtils<double>::Dot3(e_1, contra_g_1) * MathUtils<double>::Dot3(contra_g_2, e_2);
    b_e(0,1) += b(1, 1) * MathUtils<double>::Dot3(e_1, contra_g_2) * MathUtils<double>::Dot3(contra_g_2, e_2);

    b_e(1,1) = b(0, 0) * MathUtils<double>::Dot3(e_2, contra_g_1) * MathUtils<double>::Dot3(contra_g_1, e_2);
    b_e(1,1) += b(1, 0) * MathUtils<double>::Dot3(e_2, contra_g_2) * MathUtils<double>::Dot3(contra_g_1, e_2);
    b_e(1,1) += b(0, 1) * MathUtils<double>::Dot3(e_2, contra_g_1) * MathUtils<double>::Dot3(contra_g_2, e_2);
    b_e(1,1) += b(1, 1) * MathUtils<double>::Dot3(e_2, contra_g_2) * MathUtils<double>::Dot3(contra_g_2, e_2);

    double gaussian_curvature = MathUtils<double>::Det2(b_e);

    if (rNode.Id() == 1) {
        KRATOS_INFO("ShapeOpt") << "Node id = " << rNode.Id() << " || Element id = " << pElement->GetId() << " : " << std::endl;
        KRATOS_INFO("") << "g_1 = " << g_1 << " || g_2 = " << g_2 << std::endl;
        // KRATOS_INFO("") << "df = " << df <<  std::endl;
        // KRATOS_INFO("") << "dn = " << dn <<  std::endl;
        // KRATOS_INFO("") << "S = " << S <<  std::endl;
        // KRATOS_INFO("") << "m = " << m << " || p = " << p << std::endl;
        // KRATOS_INFO("") << "lambda_1 = " << lambda_1 << " || lambda_2 = " << lambda_2 << std::endl;
        KRATOS_INFO("") << "normal n = " << b << std::endl;
        KRATOS_INFO("") << "curvature tensor b = " << b << std::endl;
        KRATOS_INFO("") << "gaussian curvature = " << gaussian_curvature << std::endl;
    }

    return gaussian_curvature;
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateInnerAngleAndMixedAreaOf4Node3DQuadrilateraltAtNode(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement, double& rInnerAngle, double& rMixedArea) {

    // TODO: check angle at corner and decide how to triangulate based on that
    const auto& r_node_i = rNode.Coordinates();

    Kratos::Point::CoordinatesArrayType r_node_j;
    Kratos::Point::CoordinatesArrayType r_node_k;
    Kratos::Point::CoordinatesArrayType r_node_l;

    bool one_triangle = false;
    if (pElement->GetGeometry()[0].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[1].Coordinates();
        r_node_k = pElement->GetGeometry()[2].Coordinates();
        r_node_l = pElement->GetGeometry()[3].Coordinates();
    } else if (pElement->GetGeometry()[2].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[3].Coordinates();
        r_node_k = pElement->GetGeometry()[0].Coordinates();
        r_node_l = pElement->GetGeometry()[1].Coordinates();
    } else if (pElement->GetGeometry()[1].Id() == rNode.Id()) {
        one_triangle = true;
        r_node_j = pElement->GetGeometry()[2].Coordinates();
        r_node_k = pElement->GetGeometry()[0].Coordinates();
    } else if (pElement->GetGeometry()[3].Id() == rNode.Id()) {
        one_triangle = true;
        r_node_j = pElement->GetGeometry()[0].Coordinates();
        r_node_k = pElement->GetGeometry()[2].Coordinates();
    }

    if (one_triangle) {
        this->CalculateInnerAngleAndMixedAreaOfTriangleAtNodeI(r_node_i, r_node_j, r_node_k, rInnerAngle, rMixedArea);
    } else {
        double inner_angle = 0.0;
        this->CalculateInnerAngleAndMixedAreaOfTriangleAtNodeI(r_node_i, r_node_j, r_node_k, inner_angle, rMixedArea);
        rInnerAngle += inner_angle;
        this->CalculateInnerAngleAndMixedAreaOfTriangleAtNodeI(r_node_i, r_node_k, r_node_l, inner_angle, rMixedArea);
        rInnerAngle += inner_angle;
    }
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateInnerAngleAndMixedAreaOf8Node3DQuadrilateraltAtNode(const NodeType& rNode, const Kratos::GlobalPointer<Kratos::Condition> pElement, double& rInnerAngle, double& rMixedArea) {

    const auto& r_node_i = rNode.Coordinates();

    Kratos::Point::CoordinatesArrayType r_node_j;
    Kratos::Point::CoordinatesArrayType r_node_k;
    Kratos::Point::CoordinatesArrayType r_node_l;
    Kratos::Point::CoordinatesArrayType r_node_m;
    Kratos::Point::CoordinatesArrayType r_node_n;
    Kratos::Point::CoordinatesArrayType r_center_point = Kratos::ZeroVector(3);
    for (int i = 0; i < pElement->GetGeometry().PointsNumber(); ++i) {
        r_center_point += pElement->GetGeometry().ShapeFunctionValue(i, Kratos::ZeroVector(2)) * pElement->GetGeometry()[i].Coordinates();
    }

    bool corner_node = true;
    bool four_triangles = false;
    if (pElement->GetGeometry()[0].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[4].Coordinates();
        r_node_k = r_center_point;
        r_node_l = pElement->GetGeometry()[7].Coordinates();
    } else if (pElement->GetGeometry()[1].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[5].Coordinates();
        r_node_k = r_center_point;
        r_node_l = pElement->GetGeometry()[4].Coordinates();
    } else if (pElement->GetGeometry()[2].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[6].Coordinates();
        r_node_k = r_center_point;
        r_node_l = pElement->GetGeometry()[5].Coordinates();
    } else if (pElement->GetGeometry()[3].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[7].Coordinates();
        r_node_k = r_center_point;
        r_node_l = pElement->GetGeometry()[6].Coordinates();
    } else if (pElement->GetGeometry()[4].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[1].Coordinates();
        r_node_k = r_center_point;
        r_node_l = pElement->GetGeometry()[0].Coordinates();
    } else if (pElement->GetGeometry()[6].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[3].Coordinates();
        r_node_k = r_center_point;
        r_node_l = pElement->GetGeometry()[2].Coordinates();
    } else if (pElement->GetGeometry()[5].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[2].Coordinates();
        r_node_k = r_center_point;
        r_node_l = pElement->GetGeometry()[1].Coordinates();
    } else if (pElement->GetGeometry()[7].Id() == rNode.Id()) {
        r_node_j = pElement->GetGeometry()[0].Coordinates();
        r_node_k = r_center_point;
        r_node_l = pElement->GetGeometry()[3].Coordinates();
    }

    if (corner_node) {
        double inner_angle = 0.0;
        this->CalculateInnerAngleAndMixedAreaOfTriangleAtNodeI(r_node_i, r_node_j, r_node_k, inner_angle, rMixedArea);
        rInnerAngle += inner_angle;
        this->CalculateInnerAngleAndMixedAreaOfTriangleAtNodeI(r_node_i, r_node_k, r_node_l, inner_angle, rMixedArea);
        rInnerAngle += inner_angle;
    } else {
        double inner_angle = 0.0;
        this->CalculateInnerAngleAndMixedAreaOfTriangleAtNodeI(r_node_i, r_node_j, r_node_k, inner_angle, rMixedArea);
        rInnerAngle += inner_angle;
        this->CalculateInnerAngleAndMixedAreaOfTriangleAtNodeI(r_node_i, r_node_k, r_node_l, inner_angle, rMixedArea);
        rInnerAngle += inner_angle;
        this->CalculateInnerAngleAndMixedAreaOfTriangleAtNodeI(r_node_i, r_node_l, r_node_m, inner_angle, rMixedArea);
        rInnerAngle += inner_angle;
        this->CalculateInnerAngleAndMixedAreaOfTriangleAtNodeI(r_node_i, r_node_m, r_node_n, inner_angle, rMixedArea);
        rInnerAngle += inner_angle;
    }
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateCurvatureBasedFilterRadius()
{
    KRATOS_TRY

    class GlobalPointerAdder
    {
    public:
        typedef GlobalPointersVector<NodeType> value_type;
        typedef GlobalPointersVector<NodeType> return_type;

        return_type gp_vector;
        return_type GetValue()
        {
            gp_vector.Unique();
            return gp_vector;
        }

        void LocalReduce(const value_type &rGPVector)
        {
            for (auto &r_gp : rGPVector.GetContainer())
            {
                this->gp_vector.push_back(r_gp);
            }
        }
        void ThreadSafeReduce(GlobalPointerAdder &rOther)
        {
#pragma omp critical
            {
                for (auto &r_gp : rOther.gp_vector.GetContainer())
                {
                    this->gp_vector.push_back(r_gp);
                }
            }
        }
    };

    const auto &r_data_communicator = mrDestinationModelPart.GetCommunicator().GetDataCommunicator();

    FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType> find_neighbours_process(
        mrDestinationModelPart,
        NEIGHBOUR_NODES);
    find_neighbours_process.Execute();

    FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType> find_condition_neighbours_process(
        mrDestinationModelPart,
        NEIGHBOUR_CONDITIONS);
    find_condition_neighbours_process.Execute();

    GlobalPointersVector<NodeType> all_global_pointers =
        block_for_each<GlobalPointerAdder>(mrDestinationModelPart.Nodes(), [&](NodeType &rNode)
                                            { return rNode.GetValue(NEIGHBOUR_NODES); });

    GlobalPointerCommunicator<NodeType> pointer_comm(r_data_communicator, all_global_pointers);

    auto coordinates_proxy = pointer_comm.Apply(
        [](const GlobalPointer<NodeType> &rGP)
        { return rGP->Coordinates(); });

    double step_size = 1;

    block_for_each(mrDestinationModelPart.Nodes(), [&](NodeType &rNode) {
        double max_distance = -1.0;
        const auto& r_coordinates_origin_node = rNode.Coordinates();
        const auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES).GetContainer();
        for (const auto& r_neighbour : r_neighbours) {
            const auto& r_coordinates_neighbour_node = coordinates_proxy.Get(r_neighbour);
            const double distance = norm_2(r_coordinates_origin_node - r_coordinates_neighbour_node);
            max_distance = std::max(max_distance, distance);
        }
        CalculateCurvatureAtNode(rNode);
        double& r_gaussian_curvature = rNode.FastGetSolutionStepValue(GAUSSIAN_CURVATURE);

        double& vm_radius = rNode.FastGetSolutionStepValue(VERTEX_MORPHING_RADIUS);
        vm_radius = this->CurvatureFunction(r_gaussian_curvature, max_distance);
        double& vm_radius_raw = rNode.FastGetSolutionStepValue(VERTEX_MORPHING_RADIUS_RAW);
        vm_radius_raw = this->CurvatureFunction(r_gaussian_curvature, max_distance);
    });

    KRATOS_CATCH("");
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateNeighbourBasedFilterRadius()
{
    KRATOS_TRY

    class GlobalPointerAdder
    {
    public:
        typedef GlobalPointersVector<NodeType> value_type;
        typedef GlobalPointersVector<NodeType> return_type;

        return_type gp_vector;
        return_type GetValue()
        {
            gp_vector.Unique();
            return gp_vector;
        }

        void LocalReduce(const value_type &rGPVector)
        {
            for (auto &r_gp : rGPVector.GetContainer())
            {
                this->gp_vector.push_back(r_gp);
            }
        }
        void ThreadSafeReduce(GlobalPointerAdder &rOther)
        {
#pragma omp critical
            {
                for (auto &r_gp : rOther.gp_vector.GetContainer())
                {
                    this->gp_vector.push_back(r_gp);
                }
            }
        }
    };

    const auto &r_data_communicator = mrDestinationModelPart.GetCommunicator().GetDataCommunicator();

    FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType> find_neighbours_process(
        mrDestinationModelPart,
        NEIGHBOUR_NODES);
    find_neighbours_process.Execute();

    FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType> find_condition_neighbours_process(
        mrDestinationModelPart,
        NEIGHBOUR_CONDITIONS);
    find_condition_neighbours_process.Execute();

    GlobalPointersVector<NodeType> all_global_pointers =
        block_for_each<GlobalPointerAdder>(mrDestinationModelPart.Nodes(), [&](NodeType &rNode)
                                            { return rNode.GetValue(NEIGHBOUR_NODES); });

    GlobalPointerCommunicator<NodeType> pointer_comm(r_data_communicator, all_global_pointers);

    auto coordinates_proxy = pointer_comm.Apply(
        [](const GlobalPointer<NodeType> &rGP)
        { return rGP->Coordinates(); });

    block_for_each(mrDestinationModelPart.Nodes(), [&](NodeType &rNode) {
        double max_distance = -1.0;
        const auto& r_coordinates_origin_node = rNode.Coordinates();
        const auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES).GetContainer();
        for (const auto& r_neighbour : r_neighbours) {
            const auto& r_coordinates_neighbour_node = coordinates_proxy.Get(r_neighbour);
            const double distance = norm_2(r_coordinates_origin_node - r_coordinates_neighbour_node);
            max_distance = std::max(max_distance, distance);
        }

        // KRATOS_INFO("ShapeOpt") << "max distance of node id = " << rNode.Id() << " : " << max_distance << std::endl;
        double& vm_radius = rNode.FastGetSolutionStepValue(VERTEX_MORPHING_RADIUS);
        vm_radius = max_distance * mFilterRadiusFactor;
        double& vm_radius_raw = rNode.FastGetSolutionStepValue(VERTEX_MORPHING_RADIUS_RAW);
        vm_radius_raw = max_distance * mFilterRadiusFactor;

        // rNode.SetValue(VERTEX_MORPHING_RADIUS, max_distance * mFilterRadiusFactor);
        // rNode.SetValue(VERTEX_MORPHING_RADIUS_RAW, max_distance * mFilterRadiusFactor);
        // KRATOS_INFO("ShapeOpt") << "VERTEX_MORPHING_RADIUS on node id = " << rNode.Id() << " : " << max_distance * mFilterRadiusFactor << std::endl;
    });

    KRATOS_CATCH("");
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::SmoothenNeighbourBasedFilterRadius()
{
    KRATOS_TRY

    const IndexType number_of_nodes = mrDestinationModelPart.NumberOfNodes();
    Vector raw_radius(number_of_nodes);
    Vector temp_radius(number_of_nodes);
    IndexPartition<IndexType>(number_of_nodes).for_each([&](const IndexType iNode) {
        raw_radius[iNode] = (mrDestinationModelPart.NodesBegin() + iNode)->FastGetSolutionStepValue(VERTEX_MORPHING_RADIUS_RAW);
    });

    for (IndexType iter = 0; iter < mNumberOfSmoothingIterations; ++iter) {
        IndexPartition<IndexType>(number_of_nodes).for_each([&](const IndexType iNode) {
            const auto& r_node = *(mrDestinationModelPart.NodesBegin() + iNode);
            double current_radius = r_node.FastGetSolutionStepValue(VERTEX_MORPHING_RADIUS);
            const double current_raw_radius = raw_radius[iNode];
            if (current_raw_radius > current_radius) {
                temp_radius[iNode] = current_raw_radius;
            } else {
                temp_radius[iNode] = current_radius;
            }
        });

        IndexPartition<IndexType>(number_of_nodes).for_each([&](const IndexType iNode) {
            auto& r_node_i = *(mrDestinationModelPart.NodesBegin() + iNode);
            double& radius = r_node_i.FastGetSolutionStepValue(VERTEX_MORPHING_RADIUS);

            NodeVector neighbor_nodes(mMaxNumberOfNeighbors);
            std::vector<double> resulting_squared_distances(mMaxNumberOfNeighbors);

            const IndexType number_of_neighbors = mpSearchTree->SearchInRadius(
                                                    r_node_i,
                                                    radius,
                                                    neighbor_nodes.begin(),
                                                    resulting_squared_distances.begin(),
                                                    mMaxNumberOfNeighbors);

            std::vector<double> list_of_weights(number_of_neighbors, 0.0);
            double sum_of_weights = 0.0;
            this->ComputeWeightForAllNeighbors(r_node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
            radius = 0.0;

            for(IndexType neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; ++neighbor_itr) {
                const NodeType& r_node_j = *neighbor_nodes[neighbor_itr];
                const double radius_j = temp_radius[r_node_j.GetValue(MAPPING_ID)];
                radius += radius_j* list_of_weights[neighbor_itr] / sum_of_weights;
            }
        });
    }

    KRATOS_CATCH("");
}

template <class TBaseVertexMorphingMapper>
double MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::GetVertexMorphingRadius(const NodeType &rNode) const
{
    // KRATOS_INFO("ShapeOpt") << "VERTEX_MORPHING_RADIUS on node id = " << rNode.Id() << " : " << rNode.GetValue(VERTEX_MORPHING_RADIUS) << std::endl;
    return std::max(rNode.FastGetSolutionStepValue(VERTEX_MORPHING_RADIUS), mMinimumFilterRadius);
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateAdaptiveVertexMorphingRadius()
{
    BuiltinTimer timer;
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ShapeOpt") << "Starting calculation of adaptive vertext morphing radius for  " << mrDestinationModelPart.FullName() << "..." << std::endl;

    AssignMappingIds();
    CreateListOfNodesInOriginModelPart();
    CreateSearchTreeWithAllNodesInOriginModelPart();
    // CalculateNeighbourBasedFilterRadius();
    CalculateCurvatureBasedFilterRadius();
    SmoothenNeighbourBasedFilterRadius();

    KRATOS_INFO("ShapeOpt") << "Finished calculation of adaptive vertext morphing radius in " << timer.ElapsedSeconds() << " s." << std::endl;
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CreateSearchTreeWithAllNodesInOriginModelPart()
{
    BuiltinTimer timer;
    KRATOS_INFO("ShapeOpt") << "Creating search tree to perform mapping..." << std::endl;
    mpSearchTree = Kratos::make_unique<KDTree>(mListOfNodesInOriginModelPart.begin(), mListOfNodesInOriginModelPart.end(), mBucketSize);
    KRATOS_INFO("ShapeOpt") << "Search tree created in: " << timer.ElapsedSeconds() << " s" << std::endl;
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CreateListOfNodesInOriginModelPart()
{
    mListOfNodesInOriginModelPart.resize(mrOriginModelPart.Nodes().size());
    int counter = 0;
    for (ModelPart::NodesContainerType::iterator node_it = mrOriginModelPart.NodesBegin(); node_it != mrOriginModelPart.NodesEnd(); ++node_it)
    {
        NodeTypePointer pnode = *(node_it.base());
        mListOfNodesInOriginModelPart[counter++] = pnode;
    }
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::ComputeWeightForAllNeighbors(
    const ModelPart::NodeType& destination_node,
    const NodeVector& neighbor_nodes,
    const unsigned int number_of_neighbors,
    std::vector<double>& list_of_weights,
    double& sum_of_weights)
{
    // KRATOS_INFO("ShapeOpt") << "VERTEX_MORPHING_RADIUS on node id = " << destination_node.Id() << " : " << destination_node.GetValue(VERTEX_MORPHING_RADIUS) << std::endl;

    for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++) {
        const NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
        const double weight = this->mpFilterFunction->ComputeWeight( destination_node.Coordinates(), neighbor_node.Coordinates(), GetVertexMorphingRadius(destination_node) );

        list_of_weights[neighbor_itr] = weight;
        sum_of_weights += weight;
    }
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::AssignMappingIds()
{
    unsigned int i = 0;

    // Note: loop in the same order as in AllocateMatrix(), to avoid reallocations of the matrix.
    for(auto& node_i : mrOriginModelPart.Nodes())
        node_i.SetValue(MAPPING_ID,i++);

    i = 0;
    for(auto& node_i : mrDestinationModelPart.Nodes())
        node_i.SetValue(MAPPING_ID,i++);
}

// template instantiations
template class MapperVertexMorphingAdaptiveRadius<MapperVertexMorphing>;
template class MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingMatrixFree>;
template class MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingImprovedIntegration>;
template class MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingSymmetric>;

} // namespace Kratos.
