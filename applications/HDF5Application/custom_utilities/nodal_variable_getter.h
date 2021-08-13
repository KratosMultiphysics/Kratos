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

#ifndef KRATOS_HDF5Application_NODAL_VARIABLE_GETTER_H
#define KRATOS_HDF5Application_NODAL_VARIABLE_GETTER_H

// Project includes
#include "includes/node.h"
#include "includes/variables.h"


namespace Kratos
{
namespace HDF5
{
namespace Detail
{


#define KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(TValue)   \
    virtual TValue GetValue(const Node<3>& rNode, const Variable<TValue>& rVariable) const = 0


#define KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(TValue)                                          \
    TValue GetValue(const Node<3>& rNode, const Variable<TValue>& rVariable) const override final \
    {                                                                                             \
        return rNode.GetSolutionStepValue(rVariable);                                             \
    }


#define KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(TValue)                                      \
    TValue GetValue(const Node<3>& rNode, const Variable<TValue>& rVariable) const override final \
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

    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(bool);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(int);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(double);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(Kratos::Vector);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(Array3);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(AdjointExtensions::Pointer);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(Kratos::Matrix);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(std::string);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(DenseVector<int>);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(Quaternion<double>);
    KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER(TableStreamUtility::Pointer);
};


struct HistoricalVariableGetter : public NodalVariableGetter
{
    KRATOS_CLASS_POINTER_DEFINITION(HistoricalVariableGetter);
    using NodalVariableGetter::Array3;

    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(bool);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(int);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(double);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(Kratos::Vector);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(Array3);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(AdjointExtensions::Pointer);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(Kratos::Matrix);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(std::string);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(DenseVector<int>);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(Quaternion<double>);
    KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER(TableStreamUtility::Pointer);
};


struct NonHistoricalVariableGetter : public NodalVariableGetter
{
    KRATOS_CLASS_POINTER_DEFINITION(NonHistoricalVariableGetter);
    using NodalVariableGetter::Array3;

    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(bool);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(int);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(double);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(Kratos::Vector);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(Array3);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(AdjointExtensions::Pointer);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(Kratos::Matrix);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(std::string);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(DenseVector<int>);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(Quaternion<double>);
    KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER(TableStreamUtility::Pointer);
};


#undef KRATOS_DEFINE_NON_HISTORICAL_VARIABLE_GETTER
#undef KRATOS_DEFINE_HISTORICAL_VARIABLE_GETTER
#undef KRATOS_DECLARE_VIRTUAL_VARIABLE_GETTER


} // namespace Detail
} // namespace HDF5
} // namespace Kratos

#endif // KRATOS_HDF5Application_NODAL_VARIABLE_GETTER_H