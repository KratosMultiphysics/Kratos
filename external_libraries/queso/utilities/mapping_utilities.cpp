//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

// STL includes
#include <cstdlib>

// Project includes
#include "queso/utilities/mapping_utilities.h"

namespace queso {

// Static member operations for Mapping
PointType Mapping::PointFromGlobalToParam( const PointType& rGlobalCoord, const BoundingBoxType& rBoundXYZ, const BoundingBoxType& rBoundUVW){
    const auto delta_xyz = Math::Subtract( rBoundXYZ.second, rBoundXYZ.first );
    const auto delta_uvw = Math::Subtract( rBoundUVW.second, rBoundUVW.first );

    return PointType{ ( (rGlobalCoord[0] - rBoundXYZ.first[0]) / std::abs(delta_xyz[0]) * std::abs(delta_uvw[0]) ) + rBoundUVW.first[0],
                      ( (rGlobalCoord[1] - rBoundXYZ.first[1]) / std::abs(delta_xyz[1]) * std::abs(delta_uvw[1]) ) + rBoundUVW.first[1],
                      ( (rGlobalCoord[2] - rBoundXYZ.first[2]) / std::abs(delta_xyz[2]) * std::abs(delta_uvw[2]) ) + rBoundUVW.first[2] };
}

PointType Mapping::PointFromParamToGlobal( const PointType& rLocalCoord, const BoundingBoxType& rBoundXYZ, const BoundingBoxType& rBoundUVW ) {
    const auto delta_xyz = Math::Subtract( rBoundXYZ.second,  rBoundXYZ.first );
    const auto delta_uvw = Math::Subtract( rBoundUVW.second,  rBoundUVW.first );

    return PointType{ ( (rLocalCoord[0] - rBoundUVW.first[0]) / std::abs(delta_uvw[0]) * std::abs(delta_xyz[0]) ) + rBoundXYZ.first[0],
                      ( (rLocalCoord[1] - rBoundUVW.first[1]) / std::abs(delta_uvw[1]) * std::abs(delta_xyz[1]) ) + rBoundXYZ.first[1],
                      ( (rLocalCoord[2] - rBoundUVW.first[2]) / std::abs(delta_uvw[2]) * std::abs(delta_xyz[2]) ) + rBoundXYZ.first[2] };
}

Vector3i Mapping::GetMatrixIndicesFromVectorIndex(const IndexType Index, const Vector3i& rNumberOfElements) {
    Vector3i result;
    const IndexType index_in_row_column_plane = Index % (rNumberOfElements[0]*rNumberOfElements[1]);
    result[0] = index_in_row_column_plane % rNumberOfElements[0]; // row
    result[1] = index_in_row_column_plane / rNumberOfElements[0]; // column
    result[2] = Index / (rNumberOfElements[0]*rNumberOfElements[1]);   // depth

    return result;
}

IndexType Mapping::GetVectorIndexFromMatrixIndices(const IndexType RowIndex, const IndexType ColumnIndex, const IndexType DepthIndex, const Vector3i& rNumberOfElements) {
    return DepthIndex * (rNumberOfElements[1]*rNumberOfElements[0]) + ColumnIndex * rNumberOfElements[0] + RowIndex;
}

IndexType Mapping::GetVectorIndexFromMatrixIndices(const Vector3i& rIndices, const Vector3i& rNumberOfElements) {
    return rIndices[2] * (rNumberOfElements[1]*rNumberOfElements[0]) + rIndices[1] * rNumberOfElements[0] + rIndices[0];
}
std::pair<PointType, PointType> Mapping::GetBoundingBoxFromIndex(IndexType Index, const PointType& rLowerBound, const PointType& rUpperBound, const Vector3i& rNumberOfElements)  {
    const auto indices = GetMatrixIndicesFromVectorIndex(Index, rNumberOfElements);
    return GetBoundingBoxFromIndex( indices[0], indices[1], indices[2], rLowerBound, rUpperBound, rNumberOfElements);
}

std::pair<PointType, PointType> Mapping::GetBoundingBoxFromIndex(const Vector3i& rIndices, const PointType& rLowerBound, const PointType& rUpperBound, const Vector3i& rNumberOfElements )  {
    return GetBoundingBoxFromIndex( rIndices[0], rIndices[1], rIndices[2], rLowerBound, rUpperBound, rNumberOfElements);
}

std::pair<PointType, PointType> Mapping::GetBoundingBoxFromIndex(IndexType i, IndexType j, IndexType k, const PointType& rLowerBound, const PointType& rUpperBound, const Vector3i& rNumberOfElements)  {
    const PointType indices_d{ static_cast<double>(i), static_cast<double>(j), static_cast<double>(k) };
    PointType delta;
    delta[0] = std::abs(rUpperBound[0] - rLowerBound[0]) / (rNumberOfElements[0]);
    delta[1] = std::abs(rUpperBound[1] - rLowerBound[1]) / (rNumberOfElements[1]);
    delta[2] = std::abs(rUpperBound[2] - rLowerBound[2]) / (rNumberOfElements[2]);
    return std::make_pair( Math::Add( rLowerBound, Math::MultElementWise(delta, indices_d)),
                           Math::Add( rLowerBound, Math::MultElementWise(delta, Math::Add({1.0, 1.0, 1.0}, indices_d))) );
}

// Member operations for Mapper
Vector3i Mapper::GetMatrixIndicesFromVectorIndex(const IndexType Index ) const {
    return Mapping::GetMatrixIndicesFromVectorIndex(Index, mNumberOfElements);
}

IndexType Mapper::GetVectorIndexFromMatrixIndices(const IndexType RowIndex, const IndexType ColumnIndex, const IndexType DepthIndex ) const {
    return Mapping::GetVectorIndexFromMatrixIndices(RowIndex, ColumnIndex, DepthIndex, mNumberOfElements );
}

IndexType Mapper::GetVectorIndexFromMatrixIndices(const Vector3i& rIndices ) const {
    return Mapping::GetVectorIndexFromMatrixIndices(rIndices, mNumberOfElements );
}

std::pair<PointType, PointType> Mapper::GetBoundingBoxXYZFromIndex(IndexType Index) const {
    return Mapping::GetBoundingBoxFromIndex(Index, mBoundXYZ.first, mBoundXYZ.second, mNumberOfElements);
}

std::pair<PointType, PointType> Mapper::GetBoundingBoxXYZFromIndex(const Vector3i& Indices) const {
    return Mapping::GetBoundingBoxFromIndex(Indices, mBoundXYZ.first, mBoundXYZ.second, mNumberOfElements);
}

std::pair<PointType, PointType> Mapper::GetBoundingBoxXYZFromIndex(IndexType i, IndexType j, IndexType k) const {
    return Mapping::GetBoundingBoxFromIndex(i, j, k, mBoundXYZ.first, mBoundXYZ.second, mNumberOfElements);
}

std::pair<PointType, PointType> Mapper::GetBoundingBoxUVWFromIndex(IndexType Index) const {
    if( mBSplineMesh ){
        return Mapping::GetBoundingBoxFromIndex(Index, mBoundUVW.first, mBoundUVW.second, mNumberOfElements);
    }
    return mBoundUVW;
}

std::pair<PointType, PointType> Mapper::GetBoundingBoxUVWFromIndex(const Vector3i& Indices) const {
    if( mBSplineMesh ){
        return Mapping::GetBoundingBoxFromIndex(Indices, mBoundUVW.first, mBoundUVW.second, mNumberOfElements);
    }
    return mBoundUVW;
}

std::pair<PointType, PointType> Mapper::GetBoundingBoxUVWFromIndex(IndexType i, IndexType j, IndexType k) const {
    if( mBSplineMesh ){
        return Mapping::GetBoundingBoxFromIndex(i, j, k, mBoundUVW.first, mBoundUVW.second, mNumberOfElements);
    }
    return mBoundUVW;
}

IndexType Mapper::NumberOfElements() const {
    return mNumberOfElements[0]*mNumberOfElements[1]*mNumberOfElements[2];
}

} // End namespace queso