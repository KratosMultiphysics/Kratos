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
#include "includes/obb.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"

namespace Kratos
{
template<std::size_t TDim>
OBB<TDim>::OBB(
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
const array_1d<double, 3>& OBB<TDim>::GetCenter() const
{
    return mPointCenter;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void OBB<TDim>::SetCenter(const array_1d<double, 3>& rCenterCoords)
{
    noalias(mPointCenter) = rCenterCoords;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
const array_1d<array_1d<double, 3>, TDim>& OBB<TDim>::GetOrientationVectors() const
{
    return mOrientationVectors;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void OBB<TDim>::SetOrientationVectors(const array_1d<array_1d<double, 3>, TDim>& rOrientationVectors)
{
    for (std::size_t i_dim = 0; i_dim < TDim; ++i_dim) {
        noalias(mOrientationVectors[i_dim]) = rOrientationVectors[i_dim];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
const array_1d<double, TDim>& OBB<TDim>::GetHalfLength() const
{
    return mHalfLength;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void OBB<TDim>::SetHalfLength(const array_1d<double, TDim>& rHalfLength)
{
    mHalfLength = rHalfLength;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bool OBB<2>::HasIntersection(const OBB<2>& rOtherOBB)
{
    // Signs
    constexpr static std::array<double, 4> SignComponentsX2D = {1.0, -1.0, -1.0, 1.0};
    constexpr static std::array<double, 4> SignComponentsY2D = {1.0, 1.0, -1.0, -1.0};

    // Getting nodes from second
    const auto& r_second_obb_center = rOtherOBB.GetCenter();
    const auto& r_second_obb_half_length = rOtherOBB.GetHalfLength();
    const auto& r_second_obb_orientation_vectors = rOtherOBB.GetOrientationVectors();

    // Checking each point
    for (std::size_t i_point = 0; i_point < 4; ++i_point) {
        array_1d<double, 3> second_point = r_second_obb_center + SignComponentsX2D[i_point] *  r_second_obb_orientation_vectors[0] * r_second_obb_half_length[0] + SignComponentsY2D[i_point] * r_second_obb_orientation_vectors[1] * r_second_obb_half_length[1];

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
bool OBB<3>::HasIntersection(const OBB<3>& rOtherOBB)
{
    // Signs
    constexpr static std::array<double, 8> SignComponentsX3D = {-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0};
    constexpr static std::array<double, 8> SignComponentsY3D = {-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0};
    constexpr static std::array<double, 8> SignComponentsZ3D = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0};

    // Getting nodes from second
    const auto& r_second_obb_center = rOtherOBB.GetCenter();
    const auto& r_second_obb_half_length = rOtherOBB.GetHalfLength();
    const auto& r_second_obb_orientation_vectors = rOtherOBB.GetOrientationVectors();

    // Checking each point
    for (std::size_t i_point = 0; i_point < 8; ++i_point) {
        array_1d<double, 3> second_point = r_second_obb_center + SignComponentsX3D[i_point] *  r_second_obb_orientation_vectors[0] * r_second_obb_half_length[0] + SignComponentsY3D[i_point] * r_second_obb_orientation_vectors[1] * r_second_obb_half_length[1] + SignComponentsZ3D[i_point] * r_second_obb_orientation_vectors[2] * r_second_obb_half_length[2];

        // Check if inside
        if (CheckIsInside3D(second_point)) {
            return true;
        }
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Geometry<Point> OBB<2>::GetEquiavelentGeometry()
{
    // Signs
    constexpr static std::array<double, 4> SignComponentsX2D = {1.0, -1.0, -1.0, 1.0};
    constexpr static std::array<double, 4> SignComponentsY2D = {1.0, 1.0, -1.0, -1.0};

    // Create a quad points
    PointerVector<Point> points(4);
    for (std::size_t i_point = 0; i_point < 4; ++i_point) {
        points[i_point] = Point(mPointCenter + SignComponentsX2D[i_point] *  mOrientationVectors[0] * mHalfLength[0] + SignComponentsY2D[i_point] * mOrientationVectors[1] * mHalfLength[1]);
    }

    Quadrilateral2D4<Point> quadrilateral(points);

    return quadrilateral;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Geometry<Point> OBB<3>::GetEquiavelentGeometry()
{
    // Signs
    constexpr static std::array<double, 8> SignComponentsX3D = {-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0};
    constexpr static std::array<double, 8> SignComponentsY3D = {-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0};
    constexpr static std::array<double, 8> SignComponentsZ3D = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0};

    // Create a hexa points
    PointerVector<Point> points(8);
    for (std::size_t i_point = 0; i_point < 8; ++i_point) {
        points[i_point] = Point(mPointCenter + SignComponentsX3D[i_point] *  mOrientationVectors[0] * mHalfLength[0] + SignComponentsY3D[i_point] * mOrientationVectors[1] * mHalfLength[1] + SignComponentsZ3D[i_point] * mOrientationVectors[2] * mHalfLength[2]);
    }

    Hexahedra3D8<Point> hexahedra(points);

    return hexahedra;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void OBB<TDim>::RotateNode2D(array_1d<double, 3>& rCoords)
{
    array_1d<double, 2> old_coords;
    old_coords[0] = rCoords[0];
    old_coords[1] = rCoords[1];

    // Compute angle
    const double angle = std::atan2(mOrientationVectors[0][1], mOrientationVectors[0][0]);

    // Rotate
    rCoords[0] = std::cos(angle) * (old_coords[0] - mPointCenter[0]) - std::sin(angle) * (old_coords[1] - mPointCenter[1]);
    rCoords[1] = std::sin(angle) * (old_coords[1] - mPointCenter[1]) + std::cos(angle) * (old_coords[0] - mPointCenter[0]);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool OBB<TDim>::CheckIsInside2D(array_1d<double, 3>& rCoords)
{
    // We move to X-Y alignment
    RotateNode2D(rCoords);

    return (rCoords[1] > mPointCenter[1] - mHalfLength[1]) && (rCoords[1] < mPointCenter[1] + mHalfLength[1]) && (rCoords[0] > mPointCenter[0] - mHalfLength[0]) && (rCoords[0] < mPointCenter[0] + mHalfLength[0]);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool OBB<TDim>::CheckIsInside3D(array_1d<double, 3>& rCoords)
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template class OBB<2>;
template class OBB<3>;

}  // namespace Kratos.
