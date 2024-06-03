// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef MAPPING_UTILITIES_INCLUDE_H
#define MAPPING_UTILITIES_INCLUDE_H

//// Project includes
#include "queso/includes/define.hpp"
#include "queso/includes/parameters.h"

namespace queso {

///@name QuESo Classes
///@{

///
/**
 * @class  Mapping
 * @author Manuel Messmer
 * @brief  Provides operations two map between spaces. Static interface.
 * @see Mapper for non-static interface.
*/
class Mapping {

public:

    ///@name Public Operations
    ///@{

    /// @brief Maps point from global to parametric space.
    /// @param rGlobalCoord Point to map.
    /// @param rBoundsXYZ physical bounds of background mesh.
    /// @param rBoundsUVW parametric bounds of background mesh.
    /// @return PointType
    static PointType PointFromGlobalToParam( const PointType& rGlobalCoord, const BoundingBoxType& rBoundsXYZ, const BoundingBoxType& rBoundsUVW);

    /// @brief Maps point from parametric to global space.
    /// @param rLocalCoord Point to map.
    /// @param rBoundsXYZ of background mesh.
    /// @param rBoundsUVW of background mesh.
    /// @return PointType.
    static PointType PointFromParamToGlobal( const PointType& rLocalCoord, const BoundingBoxType& rBoundsXYZ, const BoundingBoxType& rBoundsUVW);

    /// @brief Maps univariate index to trivariate indices i -> (i,j,k)
    /// @see GetVectorIndexFromMatrixIndices
    /// @param Index
    /// @param rNumberOfElements in background mesh.
    /// @return Vector3i
    static Vector3i GetMatrixIndicesFromVectorIndex(const IndexType Index, const Vector3i& rNumberOfElements);

    /// @brief Maps trivariate indices to univariate index (i,j,k) -> i.
    /// @see GetMatrixIndicesFromVectorIndex.
    /// @param RowIndex
    /// @param ColumnIndex
    /// @param DepthIndex
    /// @param rNumberOfElements in background mesh.
    /// @return IndexType.
    static IndexType GetVectorIndexFromMatrixIndices(const IndexType RowIndex, const IndexType ColumnIndex, const IndexType DepthIndex, const Vector3i& rNumberOfElements);

    /// @brief Maps trivariate indices to univariate index (i,j,k) -> i.
    /// @see GetMatrixIndicesFromVectorIndex.
    /// @param rIndices
    /// @param rNumberOfElements in background mesh.
    /// @return IndexType
    static IndexType GetVectorIndexFromMatrixIndices(const Vector3i& rIndices, const Vector3i& rNumberOfElements);

    /// @brief Creates bounding box from given index.
    /// @param Index
    /// @param rLowerBound of background mesh.
    /// @param rUpperBound of background mesh.
    /// @param rNumberOfElements in background mesh.
    /// @return BoundingBoxType
    static BoundingBoxType GetBoundingBoxFromIndex(IndexType Index, const PointType& rLowerBound, const PointType& rUpperBound, const Vector3i& rNumberOfElements);

    /// @brief Creates bounding box from given indices.
    /// @param Indices
    /// @param rLowerBound of background mesh.
    /// @param rUpperBound of background mesh.
    /// @param rNumberOfElements in background mesh.
    /// @return BoundingBoxType
    static BoundingBoxType GetBoundingBoxFromIndex(const Vector3i& Indices, const PointType& rLowerBound, const PointType& rUpperBound, const Vector3i& rNumberOfElements );

    /// @brief Creates bounding box from given indices.
    /// @param i
    /// @param j
    /// @param k
    /// @param rLowerBound of background mesh.
    /// @param rUpperBound of background mesh.
    /// @param rNumberOfElements in background mesh.
    /// @return BoundingBoxType
    static BoundingBoxType GetBoundingBoxFromIndex(IndexType i, IndexType j, IndexType k, const PointType& rLowerBound, const PointType& rUpperBound, const Vector3i& rNumberOfElements);

    ///@}

}; // End class Mapping

/**
 * @class  Mapper
 * @author Manuel Messmer
 * @brief  Provides operations two map between spaces. Non-Static interface.
 * @see Mapping for static interface.
*/
class Mapper {
public:
    ///@name Life Cycle
    ///@{

    /// @brief Constructor
    /// @param rParameters
    Mapper( const Parameters& rParameters ) :
        mBoundXYZ( std::make_pair(rParameters.LowerBoundXYZ(), rParameters.UpperBoundXYZ()) ),
        mBoundUVW( std::make_pair(rParameters.LowerBoundUVW(), rParameters.UpperBoundUVW()) ),
        mNumberOfElements(rParameters.NumberOfElements()), mBSplineMesh(rParameters.Get<bool>("b_spline_mesh"))
    {
    }

    ///@}
    ///@name Public Operations
    ///@{


    /// @brief Maps univariate index to trivariate indices i -> (i,j,k).
    /// @see GetVectorIndexFromMatrixIndices.
    /// @param Index
    /// @return Vector3i.
    Vector3i GetMatrixIndicesFromVectorIndex(const IndexType Index) const;

    /// @brief Maps trivariate indices to univariate index (i,j,k) -> i.
    /// @see GetMatrixIndicesFromVectorIndex.
    /// @param RowIndex
    /// @param ColumnIndex
    /// @param DepthIndex
    /// @return IndexType.
    IndexType GetVectorIndexFromMatrixIndices(const IndexType RowIndex, const IndexType ColumnIndex, const IndexType DepthIndex ) const;

    /// @brief Maps trivariate indices to univariate index (i,j,k) -> i.
    /// @see GetMatrixIndicesFromVectorIndex.
    /// @param rIndices
    /// @return IndexType.
    IndexType GetVectorIndexFromMatrixIndices(const Vector3i& rIndices ) const;

    /// @brief Creates bounding box in physical space from given index.
    /// @param Index
    /// @return BoundingBoxType.
    BoundingBoxType GetBoundingBoxXYZFromIndex(IndexType Index) const;

    /// @brief Creates bounding box in physical space from given indices.
    /// @param Indices
    /// @return BoundingBoxType.
    BoundingBoxType GetBoundingBoxXYZFromIndex(const Vector3i& Indices) const;

    /// @brief Creates bounding box in physical space from given indices.
    /// @param i
    /// @param j
    /// @param k
    /// @return BoundingBoxType
    BoundingBoxType GetBoundingBoxXYZFromIndex(IndexType i, IndexType j, IndexType k) const;

    /// @brief Creates bounding box in parametric space from given index.
    /// @param Index
    /// @return BoundingBoxType.
    BoundingBoxType GetBoundingBoxUVWFromIndex(IndexType Index) const;

    /// @brief Creates bounding box in parametric space from given indices.
    /// @param Indices
    /// @return BoundingBoxType.
    BoundingBoxType GetBoundingBoxUVWFromIndex(const Vector3i& Indices) const;

    /// @brief Creates bounding box in parametric space from given indices.
    /// @param i
    /// @param j
    /// @param k
    /// @return BoundingBoxType
    BoundingBoxType GetBoundingBoxUVWFromIndex(IndexType i, IndexType j, IndexType k) const;

    /// @brief Returns global number of elements (including inactive elements).
    /// @return IndexType.
    IndexType NumberOfElements() const;

    ///@}
private:
    ///@name Private Members
    ///@{
    const BoundingBoxType mBoundXYZ;
    const BoundingBoxType mBoundUVW;
    const Vector3i mNumberOfElements;
    const bool mBSplineMesh;
    ///@}
}; // End class Mapper.
///@} End queso classes.

} // End namespace queso

#endif