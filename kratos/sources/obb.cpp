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

namespace Kratos
{
template<std::size_t TDim>
OBB<TDim>::OBB(
    const array_1d<double, 3>& rCenterCoords,
    const array_1d<double, 3>& rOrientationVector,
    const double HalfDiagonal
    ) : mPointCenter(rCenterCoords),
        mOrientationVector(rOrientationVector),
        mHalfDiagonal(HalfDiagonal)
{}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
array_1d<double, 3> OBB<TDim>::GetCenter()
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
array_1d<double, 3>& OBB<TDim>::GetOrientationVector()
{
    return mOrientationVector;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void OBB<TDim>::SetOrientationVector(const array_1d<double, 3>& rOrientationVector)
{
    noalias(mOrientationVector) = rOrientationVector;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
double& OBB<TDim>::GetHalfDiagonal()
{
    return mHalfDiagonal;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void OBB<TDim>::SetHalfDiagonal(const double HalfDiagonal)
{
    mHalfDiagonal = HalfDiagonal;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bool OBB<2>::HasIntersection(const OBB<2>& rOtherOBB)
{
    // TODO Finish
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
bool OBB<3>::HasIntersection(const OBB<3>& rOtherOBB)
{
    // TODO Finish
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template class OBB<2>;
template class OBB<3>;

}  // namespace Kratos.
