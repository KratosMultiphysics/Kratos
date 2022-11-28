//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "geometries/oriented_bounding_box.h"
#include "utilities/math_utils.h"
#include "utilities/mortar_utilities.h"
#include "utilities/intersection_utilities.h"

namespace Kratos
{
template<std::size_t TDim>
OrientedBoundingBox<TDim>::OrientedBoundingBox(
    const array_1d<double, 3>& rCenterCoords,
    const array_1d<array_1d<double, 3>, TDim>& rOrientationVectors,
    const array_1d<double, TDim>& rHalfLength
    ) : mPointCenter(rCenterCoords),
        mOrientationVectors(rOrientationVectors),
        mHalfLength(rHalfLength)
{}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
OrientedBoundingBox<TDim>::OrientedBoundingBox(
    const array_1d<double, 3>& rCenterCoords,
    const array_1d<array_1d<double, 3>, TDim>& rAxisCoordinates
    ) : mPointCenter(rCenterCoords)
{
    // Iterate the axis coordinates
    for (std::size_t i = 0; i < TDim; ++i) {
        mOrientationVectors[i] = rAxisCoordinates[i] - rCenterCoords;
        mHalfLength[i] = norm_2(mOrientationVectors[i]);
        mOrientationVectors[i] /= mHalfLength[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
OrientedBoundingBox<2>::OrientedBoundingBox(
    const GeometryType& rGeometry,
    const double BoundingBoxFactor,
    const bool BuildFromBoundingBox
    )
{
    // Check the intersection of each edge of the object bounding box against the intersecting object bounding box
    NodeType geometry_low_node, geometry_high_node;
    rGeometry.BoundingBox(geometry_low_node, geometry_high_node);

    // Creating OBB
    noalias(mOrientationVectors[0]) = geometry_high_node - geometry_low_node;
    const double norm_orientation_vectors = norm_2(mOrientationVectors[0]);
    if (norm_orientation_vectors > ZeroTolerance)
        mOrientationVectors[0] /= norm_orientation_vectors;
    else
        KRATOS_ERROR << "Zero norm on OrientedBoundingBox direction" << std::endl;
    mOrientationVectors[1][0] = mOrientationVectors[0][1];
    mOrientationVectors[1][1] = - mOrientationVectors[0][0];
    mOrientationVectors[1][2] = 0.0;
    noalias(mPointCenter) = 0.5 * (geometry_low_node.Coordinates() + geometry_high_node.Coordinates());
    mHalfLength[0] = 0.5 * norm_orientation_vectors + BoundingBoxFactor;
    mHalfLength[1] = BoundingBoxFactor;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
OrientedBoundingBox<3>::OrientedBoundingBox(
    const GeometryType& rGeometry,
    const double BoundingBoxFactor,
    const bool BuildFromBoundingBox
    )
{
    if (BuildFromBoundingBox || rGeometry.WorkingSpaceDimension() == rGeometry.LocalSpaceDimension()) {
        // Check the intersection of each face of the object bounding box against the intersecting object bounding box
        NodeType geometry_low_node, geometry_high_node;
        rGeometry.BoundingBox(geometry_low_node, geometry_high_node);

        // Creating OBB
        noalias(mOrientationVectors[0]) = geometry_high_node - geometry_low_node;
        const double norm_mOrientationVectors = norm_2(mOrientationVectors[0]);
        if (norm_mOrientationVectors > ZeroTolerance)
            mOrientationVectors[0] /= norm_mOrientationVectors;
        else
            KRATOS_ERROR << "Zero norm in OrientedBoundingBox direction" << std::endl;
        MathUtils<double>::OrthonormalBasis(mOrientationVectors[0], mOrientationVectors[1], mOrientationVectors[2]);

        // Compute center
        noalias(mPointCenter) = rGeometry.Center().Coordinates();// 0.5 * (geometry_low_node.Coordinates() + geometry_high_node.Coordinates());
    } else { // Build from orthonormal base (only for surfaces)
        // Compute the center, normal and base vectors
        noalias(mPointCenter) = rGeometry.Center().Coordinates();
        GeometryType::CoordinatesArrayType aux_coords;
        rGeometry.PointLocalCoordinates(aux_coords, mPointCenter);
        noalias(mOrientationVectors[0]) = rGeometry.UnitNormal(aux_coords);
        MathUtils<double>::OrthonormalBasis(mOrientationVectors[0], mOrientationVectors[1], mOrientationVectors[2]);

        // Getting the farthest node
        double distance = 0.0;
        const Point center = rGeometry.Center();
        Point aux_point;
        IndexType aux_i_node = 0;
        for (IndexType i_node = 0; i_node < rGeometry.size(); ++i_node) {
            noalias(aux_point.Coordinates()) = rGeometry[i_node].Coordinates();
            MortarUtilities::RotatePoint(aux_point, center, mOrientationVectors[1], mOrientationVectors[2], false);
            const double aux_distance = norm_2(aux_point.Coordinates()- mPointCenter);
            if (distance < aux_distance) {
                aux_i_node = i_node;
                distance = aux_distance;
            }
        }

        const array_1d<double, 3> aux_vector = mPointCenter - rGeometry[aux_i_node].Coordinates();
        mOrientationVectors[1] = inner_prod(aux_vector, mOrientationVectors[1]) * mOrientationVectors[1] + inner_prod(aux_vector, mOrientationVectors[2]) * mOrientationVectors[2];
        mOrientationVectors[1] /= norm_2(mOrientationVectors[1]);

        // Construct  mOrientationVectors[2]  using a cross  product
        MathUtils<double>::UnitCrossProduct(mOrientationVectors[2], mOrientationVectors[1] , mOrientationVectors[0]);
    }

    double distance_0 = 0.0;
    double distance_1 = 0.0;
    double distance_2 = 0.0;
    for (auto& r_node : rGeometry) {
        const array_1d<double, 3> vector_points = r_node.Coordinates() - mPointCenter;
        distance_0 = std::max(distance_0, std::abs(inner_prod(vector_points, mOrientationVectors[0])));
        distance_1 = std::max(distance_1, std::abs(inner_prod(vector_points, mOrientationVectors[1])));
        distance_2 = std::max(distance_2, std::abs(inner_prod(vector_points, mOrientationVectors[2])));
    }
    mHalfLength[0] = distance_0 + BoundingBoxFactor;
    mHalfLength[1] = distance_1 + BoundingBoxFactor;
    mHalfLength[2] = distance_2 + BoundingBoxFactor;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
const array_1d<double, 3>& OrientedBoundingBox<TDim>::GetCenter() const
{
    return mPointCenter;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void OrientedBoundingBox<TDim>::SetCenter(const array_1d<double, 3>& rCenterCoords)
{
    noalias(mPointCenter) = rCenterCoords;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
const array_1d<array_1d<double, 3>, TDim>& OrientedBoundingBox<TDim>::GetOrientationVectors() const
{
    return mOrientationVectors;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void OrientedBoundingBox<TDim>::SetOrientationVectors(const array_1d<array_1d<double, 3>, TDim>& rOrientationVectors)
{
    for (std::size_t i_dim = 0; i_dim < TDim; ++i_dim) {
        noalias(mOrientationVectors[i_dim]) = rOrientationVectors[i_dim];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
const array_1d<double, TDim>& OrientedBoundingBox<TDim>::GetHalfLength() const
{
    return mHalfLength;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void OrientedBoundingBox<TDim>::SetHalfLength(const array_1d<double, TDim>& rHalfLength)
{
    mHalfLength = rHalfLength;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bool OrientedBoundingBox<2>::IsInside(const OrientedBoundingBox<2>& rOtherOrientedBoundingBox)  const
{
    // Signs
    constexpr static std::array<double, 4> sign_components_X2D = {{-1.0, 1.0, 1.0, -1.0}};
    constexpr static std::array<double, 4> sign_components_Y2D = {{-1.0, -1.0, 1.0, 1.0}};

    // Getting nodes from second
    const auto& r_second_obb_center = rOtherOrientedBoundingBox.GetCenter();
    const auto& r_second_obb_half_length = rOtherOrientedBoundingBox.GetHalfLength();
    const auto& r_second_obb_orientation_vectors = rOtherOrientedBoundingBox.GetOrientationVectors();

    // Checking each point
    for (std::size_t i_point = 0; i_point < 4; ++i_point) {
        array_1d<double, 3> second_point = r_second_obb_center + sign_components_X2D[i_point] *  r_second_obb_orientation_vectors[0] * r_second_obb_half_length[0] + sign_components_Y2D[i_point] * r_second_obb_orientation_vectors[1] * r_second_obb_half_length[1];

        // Check if inside
        if (CheckIsInside2D(second_point)) {
            return true;
        }
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bool OrientedBoundingBox<3>::IsInside(const OrientedBoundingBox<3>& rOtherOrientedBoundingBox) const
{
    // Getting inverted rotation matrix
    BoundedMatrix<double, 4, 4> rotation_matrix, inverted_rotation_matrix;
    for (std::size_t i = 0; i < 3; ++i) {
        rotation_matrix(i, 0) = mOrientationVectors[0][i];
        rotation_matrix(i, 1) = mOrientationVectors[1][i];
        rotation_matrix(i, 2) = mOrientationVectors[2][i];
        rotation_matrix(i, 3) = 0.0;
        rotation_matrix(3, i) = 0.0;
    }
    rotation_matrix(3, 3) = 1.0;

    double det_rotation_matrix;
    MathUtils<double>::InvertMatrix(rotation_matrix, inverted_rotation_matrix, det_rotation_matrix);

    // Signs
    constexpr static std::array<double, 8> sign_components_X3D = {{-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0}};
    constexpr static std::array<double, 8> sign_components_Y3D = {{-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0}};
    constexpr static std::array<double, 8> sign_components_Z3D = {{-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0}};

    // Getting nodes from second
    const auto& r_second_obb_center = rOtherOrientedBoundingBox.GetCenter();
    const auto& r_second_obb_half_length = rOtherOrientedBoundingBox.GetHalfLength();
    const auto& r_second_obb_orientation_vectors = rOtherOrientedBoundingBox.GetOrientationVectors();

    // Checking each point
    for (std::size_t i_point = 0; i_point < 8; ++i_point) {
        array_1d<double, 3> second_point = r_second_obb_center + sign_components_X3D[i_point] *  r_second_obb_orientation_vectors[0] * r_second_obb_half_length[0] + sign_components_Y3D[i_point] * r_second_obb_orientation_vectors[1] * r_second_obb_half_length[1] + sign_components_Z3D[i_point] * r_second_obb_orientation_vectors[2] * r_second_obb_half_length[2];

        // Check if inside
        if (CheckIsInside3D(second_point, inverted_rotation_matrix)) {
            return true;
        }
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool OrientedBoundingBox<TDim>::HasIntersection(
    const OrientedBoundingBox<TDim>& rOtherOrientedBoundingBox,
    const OBBHasIntersectionType OBBType
    ) const
{
    if (OBBType == OBBHasIntersectionType::Direct) {
        return DirectHasIntersection(rOtherOrientedBoundingBox);
    } else if (OBBType == OBBHasIntersectionType::SeparatingAxisTheorem) {
        return SeparatingAxisTheoremHasIntersection(rOtherOrientedBoundingBox);
    } else {
        KRATOS_ERROR << "OBBType not well defined: " << static_cast<int>(OBBType) << std::endl;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool OrientedBoundingBox<TDim>::DirectHasIntersection(const OrientedBoundingBox<TDim>& rOtherOrientedBoundingBox) const
{
    /* We check if points are inside */
    // Checking one combination
    if (this->IsInside(rOtherOrientedBoundingBox))
        return true;

    // Checking the other
    if (rOtherOrientedBoundingBox.IsInside(*this))
        return true;

    /* We check for intersection of edges/faces */
    auto geom1 = this->GetEquivalentGeometry();
    auto geom2 = rOtherOrientedBoundingBox.GetEquivalentGeometry();

    // Id 2D we check edges
    if constexpr (TDim == 2) {
        const auto r_edges_1 = geom1.GenerateEdges();
        const auto r_edges_2 = geom2.GenerateEdges();
        Point int_pt(0.0,0.0,0.0);
        for (auto& r_edge_1 : r_edges_1) {
            for (auto& r_edge_2 : r_edges_2) {
                const int int_id = IntersectionUtilities::ComputeLineLineIntersection(
                    r_edge_1,
                    r_edge_2[0].Coordinates(),
                    r_edge_2[1].Coordinates(),
                    int_pt.Coordinates());

                if (int_id != 0){
                    return true;
                }
            }
        }
    } else { // If 3D we check faces
        const auto faces_1 = geom1.GenerateFaces();
        const auto faces_2 = geom2.GenerateFaces();
        for (auto& r_face_1 : faces_1) {
            for (auto& r_face_2 : faces_2) {
                if (r_face_1.HasIntersection(r_face_2)){
                    return true;
                }
            }
        }
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bool OrientedBoundingBox<2>::SeparatingAxisTheoremHasIntersection(const OrientedBoundingBox<2>& rOtherOrientedBoundingBox) const
{
    // Auxiliary values
    const auto& r_orientation_vectors_2 = rOtherOrientedBoundingBox.GetOrientationVectors();
    const array_1d<double, 3> relative_position = rOtherOrientedBoundingBox.GetCenter() - mPointCenter;

    array_1d<double, 3> cross_product_1, cross_product_2, cross_product_3, cross_product_4;

    MathUtils<double>::CrossProduct(cross_product_1, mOrientationVectors[0], r_orientation_vectors_2[0]);
    MathUtils<double>::CrossProduct(cross_product_2, mOrientationVectors[0], r_orientation_vectors_2[1]);
    MathUtils<double>::CrossProduct(cross_product_3, mOrientationVectors[1], r_orientation_vectors_2[0]);
    MathUtils<double>::CrossProduct(cross_product_4, mOrientationVectors[1], r_orientation_vectors_2[1]);

    return !(GetSeparatingPlane2D(relative_position, mOrientationVectors[0], rOtherOrientedBoundingBox) ||
        GetSeparatingPlane2D(relative_position, mOrientationVectors[1], rOtherOrientedBoundingBox) ||
        GetSeparatingPlane2D(relative_position, r_orientation_vectors_2[0], rOtherOrientedBoundingBox) ||
        GetSeparatingPlane2D(relative_position, r_orientation_vectors_2[1], rOtherOrientedBoundingBox) ||
        GetSeparatingPlane2D(relative_position, cross_product_1, rOtherOrientedBoundingBox) ||
        GetSeparatingPlane2D(relative_position, cross_product_2, rOtherOrientedBoundingBox) ||
        GetSeparatingPlane2D(relative_position, cross_product_3, rOtherOrientedBoundingBox) ||
        GetSeparatingPlane2D(relative_position, cross_product_4, rOtherOrientedBoundingBox)
        );
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bool OrientedBoundingBox<3>::SeparatingAxisTheoremHasIntersection(const OrientedBoundingBox<3>& rOtherOrientedBoundingBox) const
{
    // Auxiliary values
    const auto& r_orientation_vectors_2 = rOtherOrientedBoundingBox.GetOrientationVectors();
    const array_1d<double, 3> relative_position = rOtherOrientedBoundingBox.GetCenter() - mPointCenter;

    array_1d<double, 3> cross_product_1, cross_product_2, cross_product_3, cross_product_4, cross_product_5, cross_product_6, cross_product_7, cross_product_8, cross_product_9;

    MathUtils<double>::CrossProduct(cross_product_1, mOrientationVectors[0], r_orientation_vectors_2[0]);
    MathUtils<double>::CrossProduct(cross_product_2, mOrientationVectors[0], r_orientation_vectors_2[1]);
    MathUtils<double>::CrossProduct(cross_product_3, mOrientationVectors[0], r_orientation_vectors_2[2]);
    MathUtils<double>::CrossProduct(cross_product_4, mOrientationVectors[1], r_orientation_vectors_2[0]);
    MathUtils<double>::CrossProduct(cross_product_5, mOrientationVectors[1], r_orientation_vectors_2[1]);
    MathUtils<double>::CrossProduct(cross_product_6, mOrientationVectors[1], r_orientation_vectors_2[2]);
    MathUtils<double>::CrossProduct(cross_product_7, mOrientationVectors[2], r_orientation_vectors_2[0]);
    MathUtils<double>::CrossProduct(cross_product_8, mOrientationVectors[2], r_orientation_vectors_2[1]);
    MathUtils<double>::CrossProduct(cross_product_9, mOrientationVectors[2], r_orientation_vectors_2[2]);

    return !(GetSeparatingPlane3D(relative_position, mOrientationVectors[0], rOtherOrientedBoundingBox) ||
        GetSeparatingPlane3D(relative_position, mOrientationVectors[1], rOtherOrientedBoundingBox) ||
        GetSeparatingPlane3D(relative_position, mOrientationVectors[2], rOtherOrientedBoundingBox) ||
        GetSeparatingPlane3D(relative_position, r_orientation_vectors_2[0], rOtherOrientedBoundingBox) ||
        GetSeparatingPlane3D(relative_position, r_orientation_vectors_2[1], rOtherOrientedBoundingBox) ||
        GetSeparatingPlane3D(relative_position, r_orientation_vectors_2[2], rOtherOrientedBoundingBox) ||
        GetSeparatingPlane3D(relative_position, cross_product_1, rOtherOrientedBoundingBox) ||
        GetSeparatingPlane3D(relative_position, cross_product_2, rOtherOrientedBoundingBox) ||
        GetSeparatingPlane3D(relative_position, cross_product_3, rOtherOrientedBoundingBox) ||
        GetSeparatingPlane3D(relative_position, cross_product_4, rOtherOrientedBoundingBox) ||
        GetSeparatingPlane3D(relative_position, cross_product_5, rOtherOrientedBoundingBox) ||
        GetSeparatingPlane3D(relative_position, cross_product_6, rOtherOrientedBoundingBox) ||
        GetSeparatingPlane3D(relative_position, cross_product_7, rOtherOrientedBoundingBox) ||
        GetSeparatingPlane3D(relative_position, cross_product_8, rOtherOrientedBoundingBox) ||
        GetSeparatingPlane3D(relative_position, cross_product_9, rOtherOrientedBoundingBox)
        );
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool OrientedBoundingBox<TDim>::GetSeparatingPlane(
    const array_1d<double, 3>& rRelativePosition,
    const array_1d<double, 3>& rPlane,
    const OrientedBoundingBox& rOtherOrientedBoundingBox
    ) const
{
    // 2D/3D computation
    if constexpr (TDim == 2) {
        return GetSeparatingPlane2D(rRelativePosition, rPlane, rOtherOrientedBoundingBox);
    } else {
        return GetSeparatingPlane3D(rRelativePosition, rPlane, rOtherOrientedBoundingBox);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool OrientedBoundingBox<TDim>::GetSeparatingPlane2D(
    const array_1d<double, 3>& rRelativePosition,
    const array_1d<double, 3>& rPlane,
    const OrientedBoundingBox& rOtherOrientedBoundingBox
    ) const
{
    const auto& r_half_length_2 = rOtherOrientedBoundingBox.GetHalfLength();
    const auto& r_orientation_vectors_2 = rOtherOrientedBoundingBox.GetOrientationVectors();
    return (std::abs(inner_prod(rRelativePosition, rPlane)) >
            (std::abs(inner_prod((mOrientationVectors[0] * mHalfLength[0]), rPlane)) +
            std::abs(inner_prod((mOrientationVectors[1] * mHalfLength[1]), rPlane)) +
            std::abs(inner_prod((r_orientation_vectors_2[0] * r_half_length_2[0]), rPlane)) +
            std::abs(inner_prod((r_orientation_vectors_2[1] * r_half_length_2[1]), rPlane))));
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool OrientedBoundingBox<TDim>::GetSeparatingPlane3D(
    const array_1d<double, 3>& rRelativePosition,
    const array_1d<double, 3>& rPlane,
    const OrientedBoundingBox& rOtherOrientedBoundingBox
    ) const
{
    const auto& r_half_length_2 = rOtherOrientedBoundingBox.GetHalfLength();
    const auto& r_orientation_vectors_2 = rOtherOrientedBoundingBox.GetOrientationVectors();
    return (std::abs(inner_prod(rRelativePosition, rPlane)) >
            (std::abs(inner_prod((mOrientationVectors[0] * mHalfLength[0]), rPlane)) +
            std::abs(inner_prod((mOrientationVectors[1] * mHalfLength[1]), rPlane)) +
            std::abs(inner_prod((mOrientationVectors[2] * mHalfLength[2]), rPlane)) +
            std::abs(inner_prod((r_orientation_vectors_2[0] * r_half_length_2[0]), rPlane)) +
            std::abs(inner_prod((r_orientation_vectors_2[1] * r_half_length_2[1]), rPlane)) +
            std::abs(inner_prod((r_orientation_vectors_2[2] * r_half_length_2[2]), rPlane))));
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Quadrilateral2D4<Point> OrientedBoundingBox<2>::GetEquivalentGeometry() const
{
    // Signs
    constexpr static std::array<double, 4> sign_components_X2D = {{-1.0, 1.0, 1.0, -1.0}};
    constexpr static std::array<double, 4> sign_components_Y2D = {{-1.0, -1.0, 1.0, 1.0}};

    // Create a quad points
    std::vector<Point::Pointer> points(4);
    array_1d<double, 3> auxiliar_coords;
    for (std::size_t i_point = 0; i_point < 4; ++i_point) {
        noalias(auxiliar_coords) = mPointCenter + sign_components_X2D[i_point] *  mOrientationVectors[0] * mHalfLength[0] + sign_components_Y2D[i_point] * mOrientationVectors[1] * mHalfLength[1];
        points[i_point] = Kratos::make_shared<Point>(auxiliar_coords);
    }

    Quadrilateral2D4<Point> quadrilateral(PointerVector<Point>{points});

    return quadrilateral;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Hexahedra3D8<Point> OrientedBoundingBox<3>::GetEquivalentGeometry() const
{
    // Signs
    constexpr static std::array<double, 8> sign_components_X3D = {{-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0}};
    constexpr static std::array<double, 8> sign_components_Y3D = {{-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0}};
    constexpr static std::array<double, 8> sign_components_Z3D = {{-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0}};

    // Create a hexa points
    std::vector<Point::Pointer> points(8);
    array_1d<double, 3> auxiliar_coords;
    for (std::size_t i_point = 0; i_point < 8; ++i_point) {
        noalias(auxiliar_coords) = mPointCenter + sign_components_X3D[i_point] *  mOrientationVectors[0] * mHalfLength[0] + sign_components_Y3D[i_point] * mOrientationVectors[1] * mHalfLength[1] + sign_components_Z3D[i_point] * mOrientationVectors[2] * mHalfLength[2];
        points[i_point] = Kratos::make_shared<Point>(auxiliar_coords);
    }

    Hexahedra3D8<Point> hexahedra(PointerVector<Point>{points});

    return hexahedra;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void OrientedBoundingBox<2>::GetEquivalentRotatedGeometry(OutputType& rGeometry)
{
    for (std::size_t i_point = 0; i_point < 4; ++i_point) {
        array_1d<double, 3>& r_coordinates = rGeometry[i_point].Coordinates();
        RotateNode2D(r_coordinates);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void OrientedBoundingBox<3>::GetEquivalentRotatedGeometry(OutputType& rGeometry)
{
    // Getting inverted rotation matrix
    BoundedMatrix<double, 4, 4> rotation_matrix, inverted_rotation_matrix;
    for (std::size_t i = 0; i < 3; ++i) {
        rotation_matrix(i, 0) = mOrientationVectors[0][i];
        rotation_matrix(i, 1) = mOrientationVectors[1][i];
        rotation_matrix(i, 2) = mOrientationVectors[2][i];
        rotation_matrix(i, 3) = 0.0;
        rotation_matrix(3, i) = 0.0;
    }
    rotation_matrix(3, 3) = 1.0;

    double det_rotation_matrix;
    MathUtils<double>::InvertMatrix(rotation_matrix, inverted_rotation_matrix, det_rotation_matrix);

    for (std::size_t i_point = 0; i_point < 8; ++i_point) {
        array_1d<double, 3>& r_coordinates = rGeometry[i_point].Coordinates();
        RotateNode3D(r_coordinates, inverted_rotation_matrix);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void OrientedBoundingBox<TDim>::RotateNode2D(array_1d<double, 3>& rCoords) const
{
    // Compute angle
    const double angle = - std::atan2(mOrientationVectors[0][1], mOrientationVectors[0][0]);

    // Avoid if no rotation
    if (std::abs(angle) < ZeroTolerance) {
        return void();
    }

    array_1d<double, 2> old_coords;
    old_coords[0] = rCoords[0];
    old_coords[1] = rCoords[1];

    // Rotate
    old_coords[0] -= mPointCenter[0];
    old_coords[1] -= mPointCenter[1];
    rCoords[0] = std::cos(angle) * old_coords[0] - std::sin(angle) * old_coords[1] + mPointCenter[0];
    rCoords[1] = std::cos(angle) * old_coords[1] + std::sin(angle) * old_coords[0] + mPointCenter[1];
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void OrientedBoundingBox<TDim>::RotateNode3D(
    array_1d<double, 3>& rCoords,
    const BoundedMatrix<double, 4, 4>& rInvertedRotationMatrix
    ) const
{
    array_1d<double, 4> old_coords;
    old_coords[0] = rCoords[0] - mPointCenter[0];
    old_coords[1] = rCoords[1] - mPointCenter[1];
    old_coords[2] = rCoords[2] - mPointCenter[2];
    old_coords[3] = 1.0;

    const array_1d<double, 4> new_coords = prod(rInvertedRotationMatrix, old_coords);

    rCoords[0] = new_coords[0];
    rCoords[1] = new_coords[1];
    rCoords[2] = new_coords[2];

    // Restore movement
    noalias(rCoords) += mPointCenter;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool OrientedBoundingBox<TDim>::CheckIsInside2D(array_1d<double, 3>& rCoords) const
{
    // We move to X-Y alignment
    RotateNode2D(rCoords);

    return (std::abs(rCoords[0] - mPointCenter[0]) <= mHalfLength[0] + ZeroTolerance) && (std::abs(rCoords[1] - mPointCenter[1]) <= mHalfLength[1] + ZeroTolerance);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool OrientedBoundingBox<TDim>::CheckIsInside3D(
    array_1d<double, 3>& rCoords,
    BoundedMatrix<double, 4, 4> rInvertedRotationMatrix
    ) const
{
    // We move to X-Y-Z alignment
    RotateNode3D(rCoords, rInvertedRotationMatrix);

    return (std::abs(rCoords[0] - mPointCenter[0]) <= mHalfLength[0] + ZeroTolerance) && (std::abs(rCoords[1] - mPointCenter[1]) <= mHalfLength[1] + ZeroTolerance) && (std::abs(rCoords[2] - mPointCenter[2]) <= mHalfLength[2] + ZeroTolerance);
}

/***********************************************************************************/
/***********************************************************************************/

template class OrientedBoundingBox<2>;
template class OrientedBoundingBox<3>;

}  // namespace Kratos.
