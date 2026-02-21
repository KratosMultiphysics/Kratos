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

// External includes
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/stl/filesystem.h>

// Internal includes
#include "add_custom_utilities_to_python.h"

// Project includes
#include "custom_utilities/vertex.h"
#include "custom_utilities/vertex_utilities.h"
#include "custom_utilities/container_io_utils.h"
#include "custom_utilities/data_type_utilities.h"
#include "custom_utilities/journal.hpp"


namespace Kratos::Python {


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

    #undef KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING

    pybind11::class_<HDF5::Detail::VertexContainerType, HDF5::Detail::VertexContainerType::Pointer>(rModule, "VertexContainer")
        .def(pybind11::init<>())
        .def("push_back", pybind11::overload_cast<const HDF5::Detail::VertexContainerType::pointer&>(&HDF5::Detail::VertexContainerType::push_back))
        ;

    auto sub_module = rModule.def_submodule("Utilities");
    sub_module.def("GetListOfAvailableVariables", &HDF5::Internals::GetListOfAvailableVariables<ModelPart::NodesContainerType>, pybind11::arg("container"), pybind11::arg("data_communicator"));
    sub_module.def("GetListOfAvailableVariables", &HDF5::Internals::GetListOfAvailableVariables<ModelPart::ConditionsContainerType>, pybind11::arg("container"), pybind11::arg("data_communicator"));
    sub_module.def("GetListOfAvailableVariables", &HDF5::Internals::GetListOfAvailableVariables<ModelPart::ElementsContainerType>, pybind11::arg("container"), pybind11::arg("data_communicator"));
    sub_module.def("GetListOfAvailableVariables", &HDF5::Internals::GetListOfAvailableVariables<ModelPart::PropertiesContainerType>, pybind11::arg("container"), pybind11::arg("data_communicator"));

    pybind11::class_<JournalBase, JournalBase::Pointer>(rModule, "JournalBase")
        .def(pybind11::init<const std::filesystem::path&>())
        .def(pybind11::init<const std::filesystem::path&,const JournalBase::Extractor&>())
        .def("GetFilePath",
             &JournalBase::GetFilePath,
             "Get the path to the underlying file.")
        .def("SetExtractor",
             pybind11::overload_cast<const JournalBase::Extractor&>(&JournalBase::SetExtractor),
             "Set the extractor handling the conversion from Model to string.")
        .def("Push",
             pybind11::overload_cast<const Model&>(&JournalBase::Push),
             "Insert a new entry at the end, extracted from the input model.")
        .def("EraseIf",
             &JournalBase::EraseIf,
             "Erase all lines from the associated file matching the provided predicate.")
        .def("Clear",
             &JournalBase::Clear,
             "Delete the journal file.")
        .def("__len__",
             &JournalBase::size)
        .def("__iter__",
             [](const JournalBase& rJournal){return pybind11::make_iterator(rJournal.begin(), rJournal.end());},
             pybind11::keep_alive<0, 1>())
        ;

    pybind11::class_<Journal, Journal::Pointer>(rModule, "Journal")
        .def(pybind11::init<const std::filesystem::path&>())
        .def(pybind11::init<const std::filesystem::path&,const Journal::Extractor&>())
        .def("GetFilePath",
             &Journal::GetFilePath,
             "Get the path to the underlying file.")
        .def("SetExtractor",
             pybind11::overload_cast<const Journal::Extractor&>(&Journal::SetExtractor),
             "Set the extractor handling the conversion from Model to Parameters.")
        .def("Push",
             &Journal::Push,
             "Insert a new entry at the end, extracted from the input model.")
        .def("EraseIf",
             &Journal::EraseIf,
             "Erase all lines from the associated file matching the provided predicate.")
        .def("Clear",
             &Journal::Clear,
             "Delete the journal file.")
        .def("__len__",
             &Journal::size)
        .def("__iter__",
             [](const Journal& rJournal){return pybind11::make_iterator(rJournal.begin(), rJournal.end());},
             pybind11::keep_alive<0, 1>())
        ;

    #ifdef KRATOS_BUILD_TESTING // <== defined through CMake if cpp test sources are built
    auto testing_module = rModule.def_submodule("HDF5Testing", "Module for testing Python bindings.");
    testing_module.def("TestJournal", [](const Model& rModel, Journal& rJournal) -> void {
            constexpr const std::size_t iteration_count = 1e2;
            IndexPartition<>(iteration_count).for_each([&rModel,&rJournal](std::size_t index){
                rJournal.Push(rModel);
            });
        })
        ;
    #endif
}


} // namespace Kratos::Python
