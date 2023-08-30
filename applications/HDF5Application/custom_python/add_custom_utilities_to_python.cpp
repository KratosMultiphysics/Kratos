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
//                   Suneth Warnakulasuriya
//

// System includes

// External includes
#include "pybind11/stl.h"

// Project includes
#include "includes/define_python.h"

// Application includes
#include "custom_utilities/vertex.h"
#include "custom_utilities/vertex_utilities.h"
#include "custom_utilities/container_io_utils.h"
#include "custom_utilities/mesh_location_container.h"

// Include base h
#include "add_custom_utilities_to_python.h"


namespace Kratos
{
namespace Python
{

class PointLocatorAdaptorTrampoline : public HDF5::PointLocatorAdaptor
{
public:
    using HDF5::PointLocatorAdaptor::PointLocatorAdaptor;

    const Element::WeakPointer FindElement(const Point& rPoint) const override
    {
        using ReturnType = const Element::WeakPointer;
        using BaseType = HDF5::PointLocatorAdaptor;

        PYBIND11_OVERRIDE_PURE(
            ReturnType,
            BaseType,
            FindElement,
            rPoint);
    }
}; // class PointLocatorAdaptorTrampoline


void AddCustomUtilitiesToPython(pybind11::module& rModule)
{
    pybind11::class_<HDF5::PointLocatorAdaptor, HDF5::PointLocatorAdaptor::Pointer, PointLocatorAdaptorTrampoline>(rModule, "PointLocatorAdaptor")
        .def(pybind11::init<>())
        .def("FindElement", &HDF5::PointLocatorAdaptor::FindElement)
        ;

    pybind11::class_<HDF5::BruteForcePointLocatorAdaptor, HDF5::BruteForcePointLocatorAdaptor::Pointer, HDF5::PointLocatorAdaptor>(rModule, "BruteForcePointLocatorAdaptor")
        .def(pybind11::init<ModelPart&, const Globals::Configuration, const double>())
        ;

    #define KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(...)                                                                                     \
        .def("GetValue", [](const HDF5::Detail::Vertex& rVertex, const Variable<__VA_ARGS__>& rVariable) {                                          \
                return rVertex.GetValue(rVariable, HDF5::Internals::NonHistoricalIO{});                                                             \
            },                                                                                                                                      \
            pybind11::arg("variable"))                                                                                                              \
        .def("GetSolutionStepValue", [](const HDF5::Detail::Vertex& rVertex, const Variable<__VA_ARGS__>& rVariable, const IndexType StepIndex) {   \
                return rVertex.GetValue(rVariable, HDF5::Internals::HistoricalIO(StepIndex));                                                       \
            },                                                                                                                                      \
            pybind11::arg("variable"),                                                                                                              \
            pybind11::arg("step_index") = 0)                                                                                                        \


    pybind11::class_<HDF5::Detail::Vertex, HDF5::Detail::Vertex::Pointer, Point>(rModule, "Vertex")
        .def(pybind11::init<const array_1d<double,3>&, const HDF5::PointLocatorAdaptor&, std::size_t>(), pybind11::arg("position"), pybind11::arg("locator"), pybind11::arg("vertex_id"))
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(int)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(double)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(array_1d<double, 3>)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(array_1d<double, 4>)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(array_1d<double, 6>)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(array_1d<double, 9>)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(Kratos::Vector)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(Kratos::Matrix)
        .def_static("MakeShared", &HDF5::Detail::Vertex::MakeShared, pybind11::arg("position"), pybind11::arg("locator"), pybind11::arg("vertex_id"))
        .def("IsLocated", &HDF5::Detail::Vertex::IsLocated)
        .def("GetID", &HDF5::Detail::Vertex::GetID)
        ;

    pybind11::class_<HDF5::Detail::VertexContainerType, HDF5::Detail::VertexContainerType::Pointer>(rModule, "VertexContainer")
        .def(pybind11::init<>())
        .def("push_back", pybind11::overload_cast<const HDF5::Detail::VertexContainerType::pointer&>(&HDF5::Detail::VertexContainerType::push_back))
        ;

    #undef KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING

    pybind11::class_<HDF5::MeshLocationContainer, HDF5::MeshLocationContainer::Pointer>(rModule, "MeshLocationContainer")
        .def(pybind11::init<>())
        .def("Set", &HDF5::MeshLocationContainer::Set, pybind11::arg("hdf5_rank_id"), pybind11::arg("hdf5_process_id"), pybind11::arg("hdf5_mesh_location"))
        .def("Has", &HDF5::MeshLocationContainer::Has, pybind11::arg("hdf5_rank_id"), pybind11::arg("hdf5_process_id"))
        .def("Get", &HDF5::MeshLocationContainer::Get, pybind11::arg("hdf5_rank_id"), pybind11::arg("hdf5_process_id"))
        .def("__str__", PrintObject<HDF5::MeshLocationContainer>)
        ;

    rModule.def("HasProcessId", &HDF5::HasProcessId, pybind11::arg("parameters"));
    rModule.def("GetProcessId", &HDF5::GetProcessId, pybind11::arg("parameters"));
    rModule.def("AddProcessId", &HDF5::AddProcessId, pybind11::arg("parameters"), pybind11::arg("hdf5_rank_id"), pybind11::arg("hdf5_process_id"));
    rModule.def("HasMeshLocationContainer", &HDF5::HasMeshLocationContainer, pybind11::arg("model_part"));
    rModule.def("GetMeshLocationContainer", pybind11::overload_cast<ModelPart&>(&HDF5::GetMeshLocationContainer), pybind11::arg("model_part"));
    rModule.def("SetMeshLocationContainer", &HDF5::SetMeshLocationContainer, pybind11::arg("model_part"), pybind11::arg("mesh_location_container"));
}


} // namespace Python
} // namespace Kratos