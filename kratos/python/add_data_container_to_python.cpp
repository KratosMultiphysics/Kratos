//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes

// External includes
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Project includes
#include "add_data_container_to_python.h"
#include "containers/array_1d.h"
#include "containers/data_container/data_container.h"
#include "containers/data_container/sparse_data_value_policy.h"
#include "includes/define_python.h"

namespace Kratos::Python
{

namespace py = pybind11;

namespace
{

/**
 * @brief Python-side view over the std::string data of a DataContainer.
 * @details std::string values cannot be exposed through the NumPy buffer protocol; this
 * proxy re-fetches the span from the container on every access (making it robust against
 * sparse storage updates) and keeps the container Python object alive.
 */
class StringDataSpan
{
public:
    StringDataSpan(py::object Container, DataAccessor<std::string> Accessor)
        : mContainer(Container), mAccessor(Accessor) {}

    std::size_t size() const
    {
        return GetSpan().size();
    }

    std::string GetItem(py::ssize_t Index) const
    {
        auto span = GetSpan();
        return span[CheckedIndex(Index, span.size())];
    }

    void SetItem(py::ssize_t Index, const std::string& rValue)
    {
        auto span = GetSpan();
        span[CheckedIndex(Index, span.size())] = rValue;
    }

    std::vector<std::string> ToList() const
    {
        auto span = GetSpan();
        return std::vector<std::string>(span.begin(), span.end());
    }

private:
    Kratos::span<std::string> GetSpan() const
    {
        return mContainer.cast<DataContainer&>().GetDataSpan(mAccessor);
    }

    static std::size_t CheckedIndex(py::ssize_t Index, std::size_t Size)
    {
        if (Index < 0) Index += static_cast<py::ssize_t>(Size);
        if (Index < 0 || static_cast<std::size_t>(Index) >= Size) throw py::index_error();
        return static_cast<std::size_t>(Index);
    }

    py::object mContainer;                 /// The DataContainer Python object (kept alive by the proxy)
    DataAccessor<std::string> mAccessor;   /// The accessor of the string data
};

using DataContainerBinderType = py::class_<DataContainer, DataContainer::Pointer>;

/**
 * @brief Bind the DataAccessor / DataValuePolicy / SparseDataValuePolicy triple and the
 * typed DataContainer methods for one value type.
 * @details Class names follow the Variable naming precedent (DoubleVariable ->
 * DoubleDataValuePolicy etc.); for array_1d<double, 3> the dimension goes into the
 * suffix (Array1DDataValuePolicy3, matching Array1DVariable3).
 */
template<class TValueType>
void AddDataContainerTypeInterface(
    py::module& m,
    DataContainerBinderType& rContainerBinder,
    const std::string& rPrefix,
    const std::string& rSuffix = "")
{
    // Accessor: produced by DataContainer.Add / DataContainer.GetAccessor, no Python constructor
    py::class_<DataAccessor<TValueType>>(m, (rPrefix + "DataAccessor" + rSuffix).c_str(),
        "Lightweight handle to the data of one variable in a DataContainer. Obtained from DataContainer.Add or DataContainer.GetAccessor.")
        .def("GetIndex", &DataAccessor<TValueType>::GetIndex,
             "Index of the data chunk inside its container.")
        .def("GetStepAccessor", &DataAccessor<TValueType>::GetStepAccessor,
             py::arg("category"), py::arg("step_before_current") = 0,
             "Accessor to the same data re-parameterized to another step.")
        .def("GetStepCategory", &DataAccessor<TValueType>::GetStepCategory,
             "Step category this accessor refers to.")
        .def("GetStepBeforeCurrent", &DataAccessor<TValueType>::GetStepBeforeCurrent,
             "How many steps before the current one this accessor refers to.")
        .def(py::self == py::self);

    // Dense value policy
    py::class_<DataValuePolicy<TValueType>, DataValuePolicyBase>(m, (rPrefix + "DataValuePolicy" + rSuffix).c_str(),
        "Value policy storing one dense value per entity. The optional zero argument customizes the value used to initialize new entries.")
        .def(py::init<>())
        .def(py::init<const TValueType&>(), py::arg("zero"))
        .def("Zero", &DataValuePolicy<TValueType>::Zero, py::return_value_policy::copy,
             "Value used to zero-initialize entries governed by this policy.");

    // Sparse value policy
    py::class_<SparseDataValuePolicy<TValueType>, DataValuePolicy<TValueType>>(m, (rPrefix + "SparseDataValuePolicy" + rSuffix).c_str(),
        "Value policy for sparse storage driven by an integer index variable: entities with a non-negative index own a value, the others do not. Chunks start empty and grow via DataContainer.UpdateSparseStorage / DataContainer.AddToSparseStorage.")
        .def(py::init<DataAccessor<int>>(), py::arg("index_accessor"))
        .def(py::init<DataAccessor<int>, const TValueType&>(), py::arg("index_accessor"), py::arg("zero"));

    // Typed container methods (overloads dispatch on the bound Variable / policy types)
    rContainerBinder.def("Add",
        [](DataContainer& rSelf, const Variable<TValueType>& rVariable, const DataValuePolicy<TValueType>& rValuePolicy, const DataHistoryPolicyBase& rHistoryPolicy) {
            return rSelf.Add(rVariable, rValuePolicy, rHistoryPolicy);
        },
        py::arg("variable"), py::arg("value_policy"), py::arg("history_policy"),
        "Add a data chunk for the given registered variable and policies; returns an accessor to its data.");

    rContainerBinder.def("Add",
        [](DataContainer& rSelf, const Variable<TValueType>& rVariable, const DataValuePolicy<TValueType>& rValuePolicy) {
            return rSelf.Add(rVariable, rValuePolicy);
        },
        py::arg("variable"), py::arg("value_policy"),
        "Add a non-historical data chunk for the given registered variable; returns an accessor to its data.");

    rContainerBinder.def("GetAccessor",
        [](DataContainer& rSelf, const Variable<TValueType>& rVariable, const DataValuePolicy<TValueType>& rValuePolicy, const DataHistoryPolicyBase& rHistoryPolicy) {
            return rSelf.GetAccessor(rVariable, rValuePolicy, rHistoryPolicy);
        },
        py::arg("variable"), py::arg("value_policy"), py::arg("history_policy"),
        "Get an accessor for existing data with matching variable and policy types.");

    rContainerBinder.def("GetAccessor",
        [](DataContainer& rSelf, const Variable<TValueType>& rVariable, const DataValuePolicy<TValueType>& rValuePolicy) {
            return rSelf.GetAccessor(rVariable, rValuePolicy);
        },
        py::arg("variable"), py::arg("value_policy"),
        "Get an accessor for existing non-historical data with a matching variable and value policy type.");

    rContainerBinder.def("Has",
        [](const DataContainer& rSelf, const Variable<TValueType>& rVariable) {
            return rSelf.Has(rVariable);
        },
        py::arg("variable"),
        "Check if the variable exists in the container.");

    rContainerBinder.def("Has",
        [](const DataContainer& rSelf, const Variable<TValueType>& rVariable, const DataValuePolicyBase& rValuePolicy) {
            return rSelf.Has(rVariable, rValuePolicy);
        },
        py::arg("variable"), py::arg("value_policy"),
        "Check if the variable exists in the container with the given value policy type.");

    rContainerBinder.def("Has",
        [](const DataContainer& rSelf, const Variable<TValueType>& rVariable, const DataValuePolicyBase& rValuePolicy, const DataHistoryPolicyBase& rHistoryPolicy) {
            return rSelf.Has(rVariable, rValuePolicy, rHistoryPolicy);
        },
        py::arg("variable"), py::arg("value_policy"), py::arg("history_policy"),
        "Check if the variable exists in the container with the given value and history policy types.");

    // GetDataSpan: NumPy view for numeric types, proxy object for strings.
    // The container Python object is passed as the array base, so NumPy neither frees the
    // memory nor lets the container die while the array is referenced. The view aliases
    // the chunk storage directly: UpdateSparseStorage / AddToSparseStorage / Initialize
    // reallocate chunks, invalidating previously obtained arrays - re-fetch after them.
    if constexpr (std::is_same_v<TValueType, array_1d<double, 3>>) {
        static_assert(sizeof(array_1d<double, 3>) == 3 * sizeof(double),
                      "array_1d<double, 3> must be laid out as 3 contiguous doubles to be exposed as a NumPy view.");
        rContainerBinder.def("GetDataSpan",
            [](py::object Self, const DataAccessor<TValueType>& rAccessor) {
                auto& r_self = Self.cast<DataContainer&>();
                auto data_span = r_self.GetDataSpan(rAccessor);
                return py::array_t<double>(
                    {static_cast<py::ssize_t>(data_span.size()), static_cast<py::ssize_t>(3)},
                    {static_cast<py::ssize_t>(sizeof(array_1d<double, 3>)), static_cast<py::ssize_t>(sizeof(double))},
                    reinterpret_cast<double*>(data_span.data()),
                    Self);
            },
            py::arg("accessor"),
            "NumPy view of shape (n_entities, 3) over the data of the given accessor's step. The view aliases the chunk memory; re-fetch it after sparse storage updates.");
    } else if constexpr (std::is_same_v<TValueType, std::string>) {
        rContainerBinder.def("GetDataSpan",
            [](py::object Self, const DataAccessor<TValueType>& rAccessor) {
                return StringDataSpan(Self, rAccessor);
            },
            py::arg("accessor"),
            "StringDataSpan view over the string data of the given accessor's step.");
    } else {
        rContainerBinder.def("GetDataSpan",
            [](py::object Self, const DataAccessor<TValueType>& rAccessor) {
                auto& r_self = Self.cast<DataContainer&>();
                auto data_span = r_self.GetDataSpan(rAccessor);
                return py::array_t<TValueType>(
                    {static_cast<py::ssize_t>(data_span.size())},
                    {static_cast<py::ssize_t>(sizeof(TValueType))},
                    data_span.data(),
                    Self);
            },
            py::arg("accessor"),
            "NumPy view over the data of the given accessor's step (one value per entity). The view aliases the chunk memory; re-fetch it after sparse storage updates.");
    }
}

} // anonymous namespace

void AddDataContainerToPython(pybind11::module& m)
{
    // Step categories used by the history policies and accessors
    py::enum_<StepCategory>(m, "StepCategory")
        .value("AnyStep", StepCategory::AnyStep)
        .value("TimeStep", StepCategory::TimeStep)
        .value("IterationStep", StepCategory::IterationStep)
        .value("SubStep", StepCategory::SubStep);

    // History policies
    py::class_<DataHistoryPolicyBase>(m, "DataHistoryPolicyBase",
        "Base class of the history policies governing how many step slots a chunk stores.")
        .def("GetTotalNumberOfSteps", &DataHistoryPolicyBase::GetTotalNumberOfSteps,
             "Number of stored steps including the current one.")
        .def("GetLatestStepIndex", &DataHistoryPolicyBase::GetLatestStepIndex,
             "Index of the current step slot.")
        .def("GetStepCategory", &DataHistoryPolicyBase::GetStepCategory,
             "Step category tracked by this policy.");

    py::class_<NonHistoricalDataPolicy, DataHistoryPolicyBase>(m, "NonHistoricalDataPolicy",
        "History policy storing a single value (no history) per entity.")
        .def(py::init<>());

    py::class_<HistoricalDataPolicy, DataHistoryPolicyBase>(m, "HistoricalDataPolicy",
        "History policy storing a ring buffer of steps of one category (e.g. time steps).")
        .def(py::init<StepCategory, std::size_t>(), py::arg("category"), py::arg("buffer_size"));

    // Value policy base (concrete policies are bound per value type)
    py::class_<DataValuePolicyBase>(m, "DataValuePolicyBase",
        "Base class of the value policies governing how chunk values are stored.")
        .def("IsSparse", &DataValuePolicyBase::IsSparse,
             "Whether the policy stores data sparsely.");

    // Python-side view for std::string spans
    py::class_<StringDataSpan>(m, "StringDataSpan",
        "View over the string data of a DataContainer; re-fetches the underlying span on every access.")
        .def("__len__", &StringDataSpan::size)
        .def("__getitem__", &StringDataSpan::GetItem)
        .def("__setitem__", &StringDataSpan::SetItem)
        .def("ToList", &StringDataSpan::ToList,
             "Copy the values into a Python list.");

    // The container itself
    DataContainerBinderType container_binder(m, "DataContainer",
        "Chunked, type-erased, variable-keyed storage of entity data. Variables must be registered Kratos variables.");
    container_binder.def(py::init<>());
    container_binder.def(py::init<std::size_t>(), py::arg("chunk_size"));
    container_binder.def("Initialize", &DataContainer::Initialize, py::arg("other"), py::arg("chunk_size"),
        "Mirror the chunk structure of another container with fresh zero-initialized chunks of the given chunk size.");
    container_binder.def("UpdateSparseStorage", &DataContainer::UpdateSparseStorage, py::arg("index_accessor"),
        "Rebuild all sparse chunks keyed on the given index accessor to the number of active (non-negative) indices. Discards their previous values.");
    container_binder.def("AddToSparseStorage", &DataContainer::AddToSparseStorage, py::arg("index_accessor"), py::arg("entity_indices"),
        "Assign sparse indices to the given (currently inactive) entities and grow the sparse chunks keyed on the index accessor, preserving existing values.");
    container_binder.def("CloneStepData", &DataContainer::CloneStepData, py::arg("step_category"),
        "Clone the current step onto the next step slot for all chunks of the given category.");
    container_binder.def("__str__", PrintObject<DataContainer>);

    // Typed interfaces, following the Variable naming precedent
    AddDataContainerTypeInterface<double>(m, container_binder, "Double");
    AddDataContainerTypeInterface<int>(m, container_binder, "Integer");
    AddDataContainerTypeInterface<array_1d<double, 3>>(m, container_binder, "Array1D", "3");
    AddDataContainerTypeInterface<std::string>(m, container_binder, "String");
}

}  // namespace Kratos::Python.
