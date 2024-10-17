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

// Internal includes
#include "add_custom_utilities_to_python.h"

// Project includes
#include "custom_utilities/vertex.h"
#include "custom_utilities/vertex_utilities.h"


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

    #define KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(TValue) \
        .def("GetValue", [](const HDF5::Detail::Vertex& rVertex, const Variable<TValue>& rVariable) {return rVertex.GetValue(rVariable);})

    using Array3 = array_1d<double,3>;
    using Array4 = array_1d<double,4>;
    using Array6 = array_1d<double,6>;
    using Array9 = array_1d<double,9>;

    pybind11::class_<HDF5::Detail::Vertex, HDF5::Detail::Vertex::Pointer, Point>(rModule, "Vertex")
        .def(pybind11::init<const array_1d<double,3>&, const HDF5::PointLocatorAdaptor&, bool>())
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(bool)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(int)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(double)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(Kratos::Vector)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(Array3)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(Array4)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(Array6)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(Array9)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(Kratos::Matrix)
        KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING(DenseVector<int>)
        .def_static("MakeShared", &HDF5::Detail::Vertex::MakeShared)
        .def("IsLocated", &HDF5::Detail::Vertex::IsLocated)
        .def("GetID", &HDF5::Detail::Vertex::GetID)
        ;

    pybind11::class_<HDF5::Detail::VertexContainerType, HDF5::Detail::VertexContainerType::Pointer>(rModule, "VertexContainer")
        .def(pybind11::init<>())
        .def("push_back", pybind11::overload_cast<const HDF5::Detail::VertexContainerType::pointer&>(&HDF5::Detail::VertexContainerType::push_back))
        ;

    #undef KRATOS_DEFINE_VERTEX_GETVALUE_OVERLOAD_BINDING
}


} // namespace Python
} // namespace Kratos