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
#include <pybind11/pybind11.h>

// Project includes
#include "future/python/add_entity_data_to_python.h"
#include "future/containers/entity_data_container.h"
#include "future/containers/model_part_data_container.h"
#include "includes/define_python.h"

namespace Kratos::Future::Python
{

namespace py = pybind11;

namespace
{

using EntityDataContainerBinderType = py::class_<EntityDataContainer, EntityDataContainer::Pointer>;

/**
 * @brief Bind the typed EntityDataContainer methods for one value type.
 * @details Mirrors the Phase I AddDataContainerTypeInterface naming (Double/Integer/Array1D...3/String suffixes are implicit through the bound Variable types: pybind dispatches the overloads on the Variable<T> argument).
 */
template<class TValueType>
void AddEntityDataTypeInterface(EntityDataContainerBinderType& rBinder)
{
    rBinder.def("AddVariable",
        [](EntityDataContainer& rSelf, const Variable<TValueType>& rVariable) {
            return rSelf.AddVariable(rVariable);
        },
        py::arg("variable"),
        "Add a non-historical dense variable (one value per registered entity), zero-initialized with the variable zero.");

    rBinder.def("AddHistoricalVariable",
        [](EntityDataContainer& rSelf, const Variable<TValueType>& rVariable) {
            return rSelf.AddHistoricalVariable(rVariable);
        },
        py::arg("variable"),
        "Add a TimeStep-buffered dense variable using the container's buffer size. Advance the buffer with CloneStepData(StepCategory.TimeStep).");

    rBinder.def("Add",
        [](EntityDataContainer& rSelf, const Variable<TValueType>& rVariable, const DataValuePolicy<TValueType>& rValuePolicy, const DataHistoryPolicyBase& rHistoryPolicy) {
            return rSelf.Add(rVariable, rValuePolicy, rHistoryPolicy);
        },
        py::arg("variable"), py::arg("value_policy"), py::arg("history_policy"),
        "Add a variable with full control over the value and history policies (Phase I semantics; the caller is responsible for a meaningful policy zero).");

    rBinder.def("GetAccessor",
        [](const EntityDataContainer& rSelf, const Variable<TValueType>& rVariable) {
            return rSelf.GetAccessor(rVariable);
        },
        py::arg("variable"),
        "Accessor for a previously added variable; use with GetDataContainer().GetDataSpan for bulk NumPy access, indexing rows with Index(entity).");

    rBinder.def("Has",
        [](const EntityDataContainer& rSelf, const Variable<TValueType>& rVariable) {
            return rSelf.Has(rVariable);
        },
        py::arg("variable"),
        "Check if the variable has been added to this container.");

    rBinder.def("GetValue",
        [](EntityDataContainer& rSelf, EntityDataContainer::IndexType EntityId, const Variable<TValueType>& rVariable, int StepBeforeCurrent) {
            return rSelf.GetValue(EntityId, rVariable, StepBeforeCurrent);
        },
        py::arg("entity_id"), py::arg("variable"), py::arg("step_before_current") = 0,
        "Value of the variable for the entity with the given id (step_before_current selects the buffer step of historical variables).");

    rBinder.def("SetValue",
        [](EntityDataContainer& rSelf, EntityDataContainer::IndexType EntityId, const Variable<TValueType>& rVariable, const TValueType& rValue, int StepBeforeCurrent) {
            rSelf.SetValue(EntityId, rVariable, rValue, StepBeforeCurrent);
        },
        py::arg("entity_id"), py::arg("variable"), py::arg("value"), py::arg("step_before_current") = 0,
        "Assign the value of the variable for the entity with the given id.");
}

/**
 * @brief Bind the entity-typed registration/access sugar for one entity kind.
 */
template<class TEntityType>
void AddEntityOverloads(EntityDataContainerBinderType& rBinder)
{
    rBinder.def("RegisterEntity",
        [](EntityDataContainer& rSelf, const TEntityType& rEntity) {
            return rSelf.RegisterEntity(rEntity);
        },
        py::arg("entity"),
        "Register the entity (by its Id), returning its dense slot. Idempotent.");

    rBinder.def("Index",
        [](const EntityDataContainer& rSelf, const TEntityType& rEntity) {
            return rSelf.Index(rEntity);
        },
        py::arg("entity"),
        "Dense slot of the registered entity.");
}

} // anonymous namespace

void AddEntityDataToPython(pybind11::module& m)
{
    EntityDataContainerBinderType entity_data_binder(m, "EntityDataContainer",
        "Experimental side-car storage mapping Id-addressed Kratos entities onto a DataContainer. Parallel to (and fully independent from) the entities' existing data storage.");
    entity_data_binder.def(py::init<>());
    entity_data_binder.def(py::init<std::size_t>(), py::arg("buffer_size"));
    entity_data_binder.def(py::init<std::size_t, std::size_t>(), py::arg("buffer_size"), py::arg("chunk_size"));
    entity_data_binder.def("RegisterEntityId", &EntityDataContainer::RegisterEntityId, py::arg("entity_id"),
        "Register an entity id, returning its dense slot. Idempotent. Growth invalidates previously obtained NumPy views.");
    entity_data_binder.def("UnregisterEntityId", &EntityDataContainer::UnregisterEntityId, py::arg("entity_id"),
        "Unregister an entity id. Phase II limitation: the slot becomes a hole (not reclaimed); re-registering assigns a new slot.");
    entity_data_binder.def("HasEntity", &EntityDataContainer::HasEntity, py::arg("entity_id"),
        "Check if the entity id is registered.");
    entity_data_binder.def("Index",
        [](const EntityDataContainer& rSelf, EntityDataContainer::IndexType EntityId) {
            return rSelf.Index(EntityId);
        },
        py::arg("entity_id"),
        "Dense slot of the registered entity id.");
    entity_data_binder.def("NumberOfEntities", &EntityDataContainer::NumberOfEntities,
        "Number of registered entities.");
    entity_data_binder.def("Capacity", &EntityDataContainer::Capacity,
        "Number of slots per step of the dense chunks (the span length).");
    entity_data_binder.def("CloneStepData", &EntityDataContainer::CloneStepData, py::arg("step_category"),
        "Clone the current step onto the next step slot for all chunks of the given category.");
    entity_data_binder.def("GetBufferSize", &EntityDataContainer::GetBufferSize,
        "Number of steps stored by historical variables (fixed at construction).");
    entity_data_binder.def("GetDataContainer", &EntityDataContainer::GetDataContainer, py::return_value_policy::reference_internal,
        "The underlying DataContainer (use its GetDataSpan for NumPy access; index rows with Index(entity)). Views are invalidated by registration growth.");
    entity_data_binder.def("__str__", PrintObject<EntityDataContainer>);

    AddEntityDataTypeInterface<double>(entity_data_binder);
    AddEntityDataTypeInterface<int>(entity_data_binder);
    AddEntityDataTypeInterface<array_1d<double, 3>>(entity_data_binder);
    AddEntityDataTypeInterface<std::string>(entity_data_binder);

    AddEntityOverloads<Node>(entity_data_binder);
    AddEntityOverloads<Element>(entity_data_binder);
    AddEntityOverloads<Condition>(entity_data_binder);
    AddEntityOverloads<ModelPart::GeometryType>(entity_data_binder);
    AddEntityOverloads<MasterSlaveConstraint>(entity_data_binder);

    py::class_<ModelPartDataContainer, ModelPartDataContainer::Pointer>(m, "ModelPartDataContainer",
        "Experimental DataContainer-backed storage for the entities of a ModelPart (snapshot at construction; call Update after adding entities). The Model must outlive this object.")
        .def(py::init<ModelPart&>(), py::arg("model_part"), py::keep_alive<1, 2>())
        .def(py::init<ModelPart&, std::size_t>(), py::arg("model_part"), py::arg("chunk_size"), py::keep_alive<1, 2>())
        .def("Update", &ModelPartDataContainer::Update,
            "Register the entities added to the ModelPart since construction (removals are not tracked).")
        .def("CloneStepData", &ModelPartDataContainer::CloneStepData, py::arg("step_category"),
            "Clone the current step for all entity kinds; call right after ModelPart.CloneTimeStep.")
        .def("Nodes", &ModelPartDataContainer::Nodes, py::return_value_policy::reference_internal,
            "The nodes' EntityDataContainer.")
        .def("Elements", &ModelPartDataContainer::Elements, py::return_value_policy::reference_internal,
            "The elements' EntityDataContainer (keyed by element Id; independent per element, unlike the legacy shared-geometry path).")
        .def("Conditions", &ModelPartDataContainer::Conditions, py::return_value_policy::reference_internal,
            "The conditions' EntityDataContainer (keyed by condition Id).")
        .def("Geometries", &ModelPartDataContainer::Geometries, py::return_value_policy::reference_internal,
            "The geometries' EntityDataContainer.")
        .def("MasterSlaveConstraints", &ModelPartDataContainer::MasterSlaveConstraints, py::return_value_policy::reference_internal,
            "The master-slave constraints' EntityDataContainer.")
        .def("GetModelPart", &ModelPartDataContainer::GetModelPart, py::return_value_policy::reference_internal,
            "The managed model part.")
        .def("__str__", PrintObject<ModelPartDataContainer>);
}

}  // namespace Kratos::Future::Python.
