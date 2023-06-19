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

#ifndef KRATOS_HDF5APPLICATION_VERTEX_UTILITIES_H
#define KRATOS_HDF5APPLICATION_VERTEX_UTILITIES_H

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


#define KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(TValue)   \
    virtual TValue GetValue(const Node& rNode, const Variable<TValue>& rVariable) const = 0


#define KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(TValue)                                          \
    TValue GetValue(const Node& rNode, const Variable<TValue>& rVariable) const override final \
    {                                                                                             \
        return rNode.GetSolutionStepValue(rVariable);                                             \
    }


#define KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(TValue)                                      \
    TValue GetValue(const Node& rNode, const Variable<TValue>& rVariable) const override final \
    {                                                                                             \
        return rNode.GetValue(rVariable);                                                         \
    }


/** Base class for getting historical/non-historical variables from nodes
 *  @details the goal is to aggregate @ref{Node::GetValue} and @ref{Node::GetSolutionStepValue}
 *  to GetValue at compile time. Derived classes implement GetValue.
 *  @note overloads for every variable type are necessary since virtual templates are illegal.
 *  As a result, variables with new types will have to be added here maually.
 */
struct NodalVariableGetter
{
    KRATOS_CLASS_POINTER_DEFINITION(NodalVariableGetter);

    using Array3 = array_1d<double,3>;
    using Array4 = array_1d<double,4>;
    using Array6 = array_1d<double,6>;
    using Array9 = array_1d<double,9>;

    virtual ~NodalVariableGetter() {}

    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(bool);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(int);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(double);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(Kratos::Vector);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(Array3);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(Array4);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(Array6);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(Array9);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(Kratos::Matrix);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(DenseVector<int>);
}; // struct NodalVariableGetter


struct HistoricalVariableGetter : public NodalVariableGetter
{
    KRATOS_CLASS_POINTER_DEFINITION(HistoricalVariableGetter);
    using NodalVariableGetter::Array3;
    using NodalVariableGetter::Array4;
    using NodalVariableGetter::Array6;
    using NodalVariableGetter::Array9;

    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(bool);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(int);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(double);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(Kratos::Vector);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(Array3);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(Array4);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(Array6);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(Array9);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(Kratos::Matrix);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(DenseVector<int>);
}; // struct HistoricalVariableGetter


struct NonHistoricalVariableGetter : public NodalVariableGetter
{
    KRATOS_CLASS_POINTER_DEFINITION(NonHistoricalVariableGetter);
    using NodalVariableGetter::Array3;
    using NodalVariableGetter::Array4;
    using NodalVariableGetter::Array6;
    using NodalVariableGetter::Array9;

    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(bool);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(int);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(double);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(Kratos::Vector);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(Array3);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(Array4);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(Array6);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(Array9);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(Kratos::Matrix);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(DenseVector<int>);
}; // struct NonHistoricalVariableGetter


#undef KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER
#undef KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER
#undef KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER


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

#endif // KRATOS_HDF5APPLICATION_NODAL_VARIABLE_GETTER_H