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
#include "utilities/math_utils.h"

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
bool OBB<2>::IsInside(const OBB<2>& rOtherOBB)  const
{
    // Signs
    constexpr static std::array<double, 4> sign_components_X2D = {-1.0, 1.0, 1.0, -1.0};
    constexpr static std::array<double, 4> sign_components_Y2D = {-1.0, -1.0, 1.0, 1.0};

    // Getting nodes from second
    const auto& r_second_obb_center = rOtherOBB.GetCenter();
    const auto& r_second_obb_half_length = rOtherOBB.GetHalfLength();
    const auto& r_second_obb_orientation_vectors = rOtherOBB.GetOrientationVectors();

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
bool OBB<3>::IsInside(const OBB<3>& rOtherOBB) const
{
    // Signs
    constexpr static std::array<double, 8> sign_components_X3D = {-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0};
    constexpr static std::array<double, 8> sign_components_Y3D = {-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0};
    constexpr static std::array<double, 8> sign_components_Z3D = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0};

    // Getting nodes from second
    const auto& r_second_obb_center = rOtherOBB.GetCenter();
    const auto& r_second_obb_half_length = rOtherOBB.GetHalfLength();
    const auto& r_second_obb_orientation_vectors = rOtherOBB.GetOrientationVectors();

    // Checking each point
    for (std::size_t i_point = 0; i_point < 8; ++i_point) {
        array_1d<double, 3> second_point = r_second_obb_center + sign_components_X3D[i_point] *  r_second_obb_orientation_vectors[0] * r_second_obb_half_length[0] + sign_components_Y3D[i_point] * r_second_obb_orientation_vectors[1] * r_second_obb_half_length[1] + sign_components_Z3D[i_point] * r_second_obb_orientation_vectors[2] * r_second_obb_half_length[2];

        // Check if inside
        if (CheckIsInside3D(second_point)) {
            return true;
        }
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool OBB<TDim>::HasIntersection(const OBB<TDim>& rOtherOBB) const
{
    // Checking one combination
    if (this->IsInside(rOtherOBB))
        return true;

    // Checking the other
    if (rOtherOBB.IsInside(*this))
        return true;

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
Quadrilateral2D4<Point> OBB<2>::GetEquivalentGeometry()
{
    // Signs
    constexpr static std::array<double, 4> sign_components_X2D = {-1.0, 1.0, 1.0, -1.0};
    constexpr static std::array<double, 4> sign_components_Y2D = {-1.0, -1.0, 1.0, 1.0};

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
Hexahedra3D8<Point> OBB<3>::GetEquivalentGeometry()
{
    // Signs
    constexpr static std::array<double, 8> sign_components_X3D = {-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0};
    constexpr static std::array<double, 8> sign_components_Y3D = {-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0};
    constexpr static std::array<double, 8> sign_components_Z3D = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0};

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
void OBB<2>::GetEquivalentRotatedGeometry(OutpuType& rGeometry)
{
    for (std::size_t i_point = 0; i_point < 4; ++i_point) {
        array_1d<double, 3>& r_coordinates = rGeometry[i_point].Coordinates();
        RotateNode2D(r_coordinates);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void OBB<3>::GetEquivalentRotatedGeometry(OutpuType& rGeometry)
{
    for (std::size_t i_point = 0; i_point < 8; ++i_point) {
        array_1d<double, 3>& r_coordinates = rGeometry[i_point].Coordinates();
        RotateNode3D(r_coordinates);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void OBB<TDim>::RotateNode2D(array_1d<double, 3>& rCoords) const
{
    // Compute angle
    const double angle = - std::atan2(mOrientationVectors[0][1], mOrientationVectors[0][0]);

    // Avoid if no rotation
    if (std::abs(angle) < std::numeric_limits<double>::epsilon()) {
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
void OBB<TDim>::RotateNode3D(array_1d<double, 3>& rCoords) const
{
    array_1d<double, 4> old_coords;
    old_coords[0] = rCoords[0] - mPointCenter[0];
    old_coords[1] = rCoords[1] - mPointCenter[1];
    old_coords[2] = rCoords[2] - mPointCenter[2];
    old_coords[3] = 1.0;

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

    const array_1d<double, 4> new_coords = prod(inverted_rotation_matrix, old_coords);

    rCoords[0] = new_coords[0];
    rCoords[1] = new_coords[1];
    rCoords[2] = new_coords[2];

    // Restore movement
    noalias(rCoords) += mPointCenter;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool OBB<TDim>::CheckIsInside2D(array_1d<double, 3>& rCoords) const
{
    // We move to X-Y alignment
    RotateNode2D(rCoords);

    return (std::abs(rCoords[0] - mPointCenter[0]) <= mHalfLength[0] + ZeroTolerance) && (std::abs(rCoords[1] - mPointCenter[1]) <= mHalfLength[1] + ZeroTolerance);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool OBB<TDim>::CheckIsInside3D(array_1d<double, 3>& rCoords) const
{
    // We move to X-Y-Z alignment
    RotateNode3D(rCoords);

    return (std::abs(rCoords[0] - mPointCenter[0]) <= mHalfLength[0] + ZeroTolerance) && (std::abs(rCoords[1] - mPointCenter[1]) <= mHalfLength[1] + ZeroTolerance) && (std::abs(rCoords[2] - mPointCenter[2]) <= mHalfLength[2] + ZeroTolerance);
}

/***********************************************************************************/
/***********************************************************************************/

template class OBB<2>;
template class OBB<3>;

}  // namespace Kratos.
