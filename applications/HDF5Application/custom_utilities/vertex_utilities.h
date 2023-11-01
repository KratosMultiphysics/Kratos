//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

#pragma once

// Project includes
#include "includes/node.h"
#include "includes/variables.h"

// Core includes
#include "includes/element.h"
#include "utilities/brute_force_point_locator.h"
#include "includes/model_part.h"

namespace Kratos
{
namespace HDF5
{

/** Interface for point locators
 *  @details narrow down the scope of a locator's tasks (eg.: find element containing a point)
 */
struct KRATOS_API(HDF5_APPLICATION) PointLocatorAdaptor
{
    KRATOS_CLASS_POINTER_DEFINITION(PointLocatorAdaptor);

    virtual ~PointLocatorAdaptor() {}

    virtual const Element::WeakPointer FindElement(const Point& rPoint) const = 0;
}; // struct PointLocatorAdaptor


/// BruteForcePointLocator with configuration and tolerance persistence
class KRATOS_API(HDF5_APPLICATION) BruteForcePointLocatorAdaptor final : public PointLocatorAdaptor
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BruteForcePointLocatorAdaptor);

    BruteForcePointLocatorAdaptor(ModelPart& rModelPart,
                                  const Globals::Configuration configuration,
                                  const double tolerance);

    const Element::WeakPointer FindElement(const Point& rPoint) const override;

private:
    BruteForcePointLocator mLocator;

    const ModelPart& mrModelPart;

    const Globals::Configuration mConfiguration;

    const double mTolerance;
}; // class BruteForcePointLocatorAdaptor


} // namespace HDF5
} // namespace Kratos