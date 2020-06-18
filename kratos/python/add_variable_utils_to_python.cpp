//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes

// External includes

// Project includes
#include "python/add_variable_utils_to_python.h"
#include "includes/define_python.h"
#include "processes/process.h"

// Variable utilities
#include "utilities/variable_utils.h"

namespace Kratos {
namespace Python {

void VariableUtilsUpdateCurrentPosition(
    VariableUtils &rVariableUtils,
    const ModelPart::NodesContainerType &rNodes
    )
{
    rVariableUtils.UpdateCurrentPosition(rNodes);
}

void VariableUtilsUpdateCurrentPositionWithVariable(
    VariableUtils &rVariableUtils,
    const ModelPart::NodesContainerType &rNodes,
    const VariableUtils::ArrayVarType &rUpdateVariable
    )
{
    rVariableUtils.UpdateCurrentPosition(rNodes, rUpdateVariable);
}

void VariableUtilsUpdateCurrentPositionWithVariableAndPosition(
    VariableUtils &rVariableUtils,
    const ModelPart::NodesContainerType &rNodes,
    const VariableUtils::ArrayVarType &rUpdateVariable,
    const IndexType BufferPosition
    )
{
    rVariableUtils.UpdateCurrentPosition(rNodes, rUpdateVariable, BufferPosition);
}

template<class TVarType>
void VariableUtilsCopyModelPartNodalVar(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart,
    const unsigned int BuffStep = 0)
{
    rVariableUtils.CopyModelPartNodalVar(rVariable, rOriginModelPart, rDestinationModelPart, BuffStep);
}

template<class TVarType>
void VariableUtilsCopyModelPartNodalVarWithDestination(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TVarType &rDestinationVariable,
    const ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart,
    const unsigned int BuffStep = 0)
{
    rVariableUtils.CopyModelPartNodalVar(rVariable, rDestinationVariable, rOriginModelPart, rDestinationModelPart, BuffStep);
}

template<class TVarType>
void CopyModelPartNodalVarToNonHistoricalVar(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart,
    const unsigned int BuffStep = 0)
{
    rVariableUtils.CopyModelPartNodalVarToNonHistoricalVar(rVariable, rOriginModelPart, rDestinationModelPart, BuffStep);
}

template<class TVarType>
void CopyModelPartNodalVarToNonHistoricalVarWithDestination(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TVarType &rDestinationVariable,
    const ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart,
    const unsigned int BuffStep = 0)
{
    rVariableUtils.CopyModelPartNodalVarToNonHistoricalVar(rVariable, rDestinationVariable, rOriginModelPart, rDestinationModelPart, BuffStep);
}

template <class TDataType, class TContainerType, class TVarType = Variable<TDataType>>
void VariableUtilsSetNonHistoricalVariable(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TDataType &rValue,
    TContainerType &rContainer)
{
    rVariableUtils.SetNonHistoricalVariable(rVariable, rValue, rContainer);
}

template <class TDataType, class TContainerType, class TVarType = Variable<TDataType>>
void VariableUtilsSetNonHistoricalVariableForFlag(
    VariableUtils &rVariableUtils,
    const TVarType &rVariable,
    const TDataType &rValue,
    TContainerType &rContainer,
    const Flags Flag,
    const bool CheckValue = true)
{
    rVariableUtils.SetNonHistoricalVariable(rVariable, rValue, rContainer, Flag, CheckValue);
}

void AddVariableUtilsToPython(pybind11::module &m)
{
    namespace py = pybind11;

    py::class_<VariableUtils>(m, "VariableUtils")
        .def(py::init<>())
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<bool>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<double>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<array_1d<double, 3>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<array_1d<double, 4>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<array_1d<double, 6>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<array_1d<double, 9>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<Quaternion<double>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<Vector>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVar<Variable<Matrix>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<bool>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<double>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<array_1d<double, 3>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<array_1d<double, 4>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<array_1d<double, 6>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<array_1d<double, 9>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<Quaternion<double>>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<Vector>>)
        .def("CopyModelPartNodalVar", VariableUtilsCopyModelPartNodalVarWithDestination<Variable<Matrix>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<bool>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<double>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<array_1d<double, 3>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<array_1d<double, 4>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<array_1d<double, 6>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<array_1d<double, 9>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<Quaternion<double>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<Vector>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVar<Variable<Matrix>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<bool>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<double>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<array_1d<double, 3>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<array_1d<double, 4>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<array_1d<double, 6>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<array_1d<double, 9>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<Quaternion<double>>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<Vector>>)
        .def("CopyModelPartNodalVarToNonHistoricalVar", CopyModelPartNodalVarToNonHistoricalVarWithDestination<Variable<Matrix>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<bool>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<double>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double, 3>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double, 4>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double, 6>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<array_1d<double, 9>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<Quaternion<double>>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<Vector>>)
        .def("CopyModelPartElementalVar", &VariableUtils::CopyModelPartElementalVar<Variable<Matrix>>)
        .def("SetVectorVar", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 3>> &rVariable, const array_1d<double, 3> &rValue, ModelPart::NodesContainerType &rNodes) {rVariableUtils.SetVariable(rVariable, rValue, rNodes);})
        .def("SetVectorVar", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 3>> &rVariable, array_1d<double, 3> &rValue, ModelPart::NodesContainerType &rNodes, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, BuffStep);})
        .def("SetVectorVar", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 3>> &rVariable, const array_1d<double, 3> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag);})
        .def("SetVectorVar", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 3>> &rVariable, array_1d<double, 3> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue);})
        .def("SetVectorVar", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 3>> &rVariable, const array_1d<double, 3> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue, BuffStep);})
        .def("SetScalarVar", [](VariableUtils &rVariableUtils, const Variable<double> &rVariable, const double &rValue, ModelPart::NodesContainerType &rNodes) {rVariableUtils.SetVariable(rVariable, rValue, rNodes);})
        .def("SetScalarVar", [](VariableUtils &rVariableUtils, const Variable<double> &rVariable, double &rValue, ModelPart::NodesContainerType &rNodes, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, BuffStep);})
        .def("SetScalarVar", [](VariableUtils &rVariableUtils, const Variable<double> &rVariable, const double &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag);})
        .def("SetScalarVar", [](VariableUtils &rVariableUtils, const Variable<double> &rVariable, double &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue);})
        .def("SetScalarVar", [](VariableUtils &rVariableUtils, const Variable<double> &rVariable, const double &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue, BuffStep);})
        .def("SetNonHistoricalScalarVar", VariableUtilsSetNonHistoricalVariable<int, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalScalarVar", VariableUtilsSetNonHistoricalVariable<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVectorVar", VariableUtilsSetNonHistoricalVariable<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<int> &rVariable, const int &rValue, ModelPart::NodesContainerType &rNodes) {rVariableUtils.SetVariable(rVariable, rValue, rNodes);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<int> &rVariable, int &rValue, ModelPart::NodesContainerType &rNodes, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<bool> &rVariable, const bool &rValue, ModelPart::NodesContainerType &rNodes) {rVariableUtils.SetVariable(rVariable, rValue, rNodes);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<bool> &rVariable, bool &rValue, ModelPart::NodesContainerType &rNodes, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<double> &rVariable, const double &rValue, ModelPart::NodesContainerType &rNodes) {rVariableUtils.SetVariable(rVariable, rValue, rNodes);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<double> &rVariable, double &rValue, ModelPart::NodesContainerType &rNodes, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 3>> &rVariable, const array_1d<double, 3> &rValue, ModelPart::NodesContainerType &rNodes) {rVariableUtils.SetVariable(rVariable, rValue, rNodes);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 3>> &rVariable, array_1d<double, 3> &rValue, ModelPart::NodesContainerType &rNodes, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 4>> &rVariable, const array_1d<double, 4> &rValue, ModelPart::NodesContainerType &rNodes) {rVariableUtils.SetVariable(rVariable, rValue, rNodes);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 4>> &rVariable, array_1d<double, 4> &rValue, ModelPart::NodesContainerType &rNodes, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 6>> &rVariable, const array_1d<double, 6> &rValue, ModelPart::NodesContainerType &rNodes) {rVariableUtils.SetVariable(rVariable, rValue, rNodes);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 6>> &rVariable, array_1d<double, 6> &rValue, ModelPart::NodesContainerType &rNodes, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 9>> &rVariable, const array_1d<double, 9> &rValue, ModelPart::NodesContainerType &rNodes) {rVariableUtils.SetVariable(rVariable, rValue, rNodes);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 9>> &rVariable, array_1d<double, 9> &rValue, ModelPart::NodesContainerType &rNodes, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<Vector> &rVariable, const Vector &rValue, ModelPart::NodesContainerType &rNodes) {rVariableUtils.SetVariable(rVariable, rValue, rNodes);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<Vector> &rVariable, Vector &rValue, ModelPart::NodesContainerType &rNodes, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<Matrix> &rVariable, const Matrix &rValue, ModelPart::NodesContainerType &rNodes) {rVariableUtils.SetVariable(rVariable, rValue, rNodes);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<Matrix> &rVariable, Matrix &rValue, ModelPart::NodesContainerType &rNodes, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<Quaternion<double>> &rVariable, const Quaternion<double> &rValue, ModelPart::NodesContainerType &rNodes) {rVariableUtils.SetVariable(rVariable, rValue, rNodes);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<Quaternion<double>> &rVariable, Quaternion<double> &rValue, ModelPart::NodesContainerType &rNodes, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<int> &rVariable, const int &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<int> &rVariable, int &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<int> &rVariable, const int &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<bool> &rVariable, const bool &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<bool> &rVariable, bool &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<bool> &rVariable, const bool &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<double> &rVariable, const double &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<double> &rVariable, double &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<double> &rVariable, const double &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 3>> &rVariable, const array_1d<double, 3> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 3>> &rVariable, array_1d<double, 3> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 3>> &rVariable, const array_1d<double, 3> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 4>> &rVariable, const array_1d<double, 4> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 4>> &rVariable, array_1d<double, 4> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 4>> &rVariable, const array_1d<double, 4> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 6>> &rVariable, const array_1d<double, 6> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 6>> &rVariable, array_1d<double, 6> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 6>> &rVariable, const array_1d<double, 6> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 9>> &rVariable, const array_1d<double, 9> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 9>> &rVariable, array_1d<double, 9> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<array_1d<double, 9>> &rVariable, const array_1d<double, 9> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<Vector> &rVariable, const Vector &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<Vector> &rVariable, Vector &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<Vector> &rVariable, const Vector &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<Matrix> &rVariable, const Matrix &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<Matrix> &rVariable, Matrix &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<Matrix> &rVariable, const Matrix &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue, BuffStep);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<Quaternion<double>> &rVariable, const Quaternion<double> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<Quaternion<double>> &rVariable, const Quaternion<double> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue);})
        .def("SetVariable", [](VariableUtils &rVariableUtils, const Variable<Quaternion<double>> &rVariable, const Quaternion<double> &rValue, ModelPart::NodesContainerType &rNodes, const Flags Flag, const bool CheckValue, const IndexType BuffStep) {rVariableUtils.SetVariable(rVariable, rValue, rNodes, Flag, CheckValue, BuffStep);})
        .def("SetHistoricalVariableToZero", &VariableUtils::SetHistoricalVariableToZero<int>)
        .def("SetHistoricalVariableToZero", &VariableUtils::SetHistoricalVariableToZero<double>)
        .def("SetHistoricalVariableToZero", &VariableUtils::SetHistoricalVariableToZero<array_1d<double, 3>>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<int, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<int, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<int, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<double, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 3>, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 4>, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 6>, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariableToZero", &VariableUtils::SetNonHistoricalVariableToZero<array_1d<double, 9>, ModelPart::MasterSlaveConstraintContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<int, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<bool, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 4>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 6>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 9>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<Quaternion<double>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<Vector, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<Matrix, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<bool, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<double, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 3>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 4>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 6>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 9>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<Quaternion<double>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<Vector, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<Matrix, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<bool, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<double, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 3>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 4>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 6>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<array_1d<double, 9>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<Quaternion<double>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<Vector, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariable<Matrix, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<bool, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<double, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 4>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 6>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 9>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<Quaternion<double>, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<Vector, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<Matrix, ModelPart::NodesContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<bool, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<double, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 3>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 4>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 6>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 9>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<Quaternion<double>, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<Vector, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<Matrix, ModelPart::ConditionsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<bool, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<double, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 3>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 4>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 6>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<array_1d<double, 9>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<Quaternion<double>, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<Vector, ModelPart::ElementsContainerType>)
        .def("SetNonHistoricalVariable", VariableUtilsSetNonHistoricalVariableForFlag<Matrix, ModelPart::ElementsContainerType>)
        .def("ClearNonHistoricalData", &VariableUtils::ClearNonHistoricalData<ModelPart::NodesContainerType>)
        .def("ClearNonHistoricalData", &VariableUtils::ClearNonHistoricalData<ModelPart::ConditionsContainerType>)
        .def("ClearNonHistoricalData", &VariableUtils::ClearNonHistoricalData<ModelPart::ElementsContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::NodesContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::ConditionsContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::ElementsContainerType>)
        .def("SetFlag", &VariableUtils::SetFlag<ModelPart::MasterSlaveConstraintContainerType>)
        .def("ResetFlag", &VariableUtils::ResetFlag<ModelPart::NodesContainerType>)
        .def("ResetFlag", &VariableUtils::ResetFlag<ModelPart::ConditionsContainerType>)
        .def("ResetFlag", &VariableUtils::ResetFlag<ModelPart::ElementsContainerType>)
        .def("ResetFlag", &VariableUtils::ResetFlag<ModelPart::MasterSlaveConstraintContainerType>)
        .def("FlipFlag", &VariableUtils::FlipFlag<ModelPart::NodesContainerType>)
        .def("FlipFlag", &VariableUtils::FlipFlag<ModelPart::ConditionsContainerType>)
        .def("FlipFlag", &VariableUtils::FlipFlag<ModelPart::ElementsContainerType>)
        .def("FlipFlag", &VariableUtils::FlipFlag<ModelPart::MasterSlaveConstraintContainerType>)
        .def("SaveScalarVar", &VariableUtils::SaveVariable<double>)              // To be removed
        .def("SaveVectorVar", &VariableUtils::SaveVariable<array_1d<double, 3>>) // To be removed
        .def("SaveVariable", &VariableUtils::SaveVariable<bool>)
        .def("SaveVariable", &VariableUtils::SaveVariable<double>)
        .def("SaveVariable", &VariableUtils::SaveVariable<array_1d<double, 3>>)
        .def("SaveVariable", &VariableUtils::SaveVariable<array_1d<double, 4>>)
        .def("SaveVariable", &VariableUtils::SaveVariable<array_1d<double, 6>>)
        .def("SaveVariable", &VariableUtils::SaveVariable<array_1d<double, 9>>)
        .def("SaveVariable", &VariableUtils::SaveVariable<Vector>)
        .def("SaveVariable", &VariableUtils::SaveVariable<Matrix>)
        .def("SaveScalarNonHistoricalVar", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::NodesContainerType>)              // To be removed
        .def("SaveVectorNonHistoricalVar", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 3>, ModelPart::NodesContainerType>) // To be removed
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<bool, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 3>, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 4>, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 6>, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 9>, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<Vector, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<Matrix, ModelPart::NodesContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<bool, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 3>, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 4>, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 6>, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 9>, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<Vector, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<Matrix, ModelPart::ConditionsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<bool, ModelPart::ElementsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<double, ModelPart::ElementsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 3>, ModelPart::ElementsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 4>, ModelPart::ElementsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 6>, ModelPart::ElementsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<array_1d<double, 9>, ModelPart::ElementsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<Vector, ModelPart::ElementsContainerType>)
        .def("SaveNonHistoricalVariable", &VariableUtils::SaveNonHistoricalVariable<Matrix, ModelPart::ElementsContainerType>)
        .def("SelectNodeList", &VariableUtils::SelectNodeList)
        .def("CopyScalarVar", &VariableUtils::CopyVariable<double>)                                                                    // To be removed
        .def("CopyVectorVar", &VariableUtils::CopyVariable<array_1d<double, 3>>)                                                       // To be removed
        .def("CopyVariable", &VariableUtils::CopyVariable<bool>)
        .def("CopyVariable", &VariableUtils::CopyVariable<double>)
        .def("CopyVariable", &VariableUtils::CopyVariable<array_1d<double, 3>>)
        .def("CopyVariable", &VariableUtils::CopyVariable<array_1d<double, 4>>)
        .def("CopyVariable", &VariableUtils::CopyVariable<array_1d<double, 6>>)
        .def("CopyVariable", &VariableUtils::CopyVariable<array_1d<double, 9>>)
        .def("CopyVariable", &VariableUtils::CopyVariable<Vector>)
        .def("CopyVariable", &VariableUtils::CopyVariable<Matrix>)
        .def("ApplyFixity", &VariableUtils::ApplyFixity<Variable<double>>)
        .def("ApplyVector", &VariableUtils::ApplyVector<Variable<double>>)
        .def("SumHistoricalNodeScalarVariable", &VariableUtils::SumHistoricalVariable<double>)
        .def("SumHistoricalNodeVectorVariable", &VariableUtils::SumHistoricalVariable<array_1d<double, 3>>)
        .def("SumNonHistoricalNodeScalarVariable", &VariableUtils::SumNonHistoricalNodeScalarVariable<Variable<double>>)
        .def("SumNonHistoricalNodeVectorVariable", &VariableUtils::SumNonHistoricalNodeVectorVariable)
        .def("SumConditionScalarVariable", &VariableUtils::SumConditionScalarVariable<Variable<double>>)
        .def("SumConditionVectorVariable", &VariableUtils::SumConditionVectorVariable)
        .def("SumElementScalarVariable", &VariableUtils::SumElementScalarVariable<Variable<double>>)
        .def("SumElementVectorVariable", &VariableUtils::SumElementVectorVariable)
        .def("AddDof", &VariableUtils::AddDof<Variable<double>>)
        .def("AddDof", &VariableUtils::AddDofWithReaction<Variable<double>>)
        .def("CheckVariableKeys", &VariableUtils::CheckVariableKeys)
        .def("CheckDofs", &VariableUtils::CheckDofs)
        .def("UpdateCurrentToInitialConfiguration", &VariableUtils::UpdateCurrentToInitialConfiguration)
        .def("UpdateInitialToCurrentConfiguration", &VariableUtils::UpdateInitialToCurrentConfiguration)
        .def("UpdateCurrentPosition", VariableUtilsUpdateCurrentPosition)
        .def("UpdateCurrentPosition", VariableUtilsUpdateCurrentPositionWithVariable)
        .def("UpdateCurrentPosition", VariableUtilsUpdateCurrentPositionWithVariableAndPosition);
}

} // namespace Python.
} // Namespace Kratos
