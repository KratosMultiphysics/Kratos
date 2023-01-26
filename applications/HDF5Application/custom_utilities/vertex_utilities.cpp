//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

// Internal includes
#include "vertex_utilities.h"


namespace Kratos
{
namespace HDF5
{


BruteForcePointLocatorAdaptor::BruteForcePointLocatorAdaptor(ModelPart& rModelPart,
                                                             const Globals::Configuration configuration,
                                                             const double tolerance)
    : mLocator(rModelPart),
      mrModelPart(rModelPart),
      mConfiguration(configuration),
      mTolerance(tolerance)
{
}


const Element::WeakPointer BruteForcePointLocatorAdaptor::FindElement(const Point& rPoint) const
{
    KRATOS_TRY

    Kratos::Vector shape_function_values;

    const int element_id =  mLocator.FindElement(
        rPoint,
        shape_function_values,
        mConfiguration,
        mTolerance);

    // Awkward syntax due to the intrusive pointers of Element

    return -1 < element_id ? mrModelPart.pGetElement(element_id).get() : nullptr;

    KRATOS_CATCH("");
}


} // namespace HDF5
} // namespace Kratos