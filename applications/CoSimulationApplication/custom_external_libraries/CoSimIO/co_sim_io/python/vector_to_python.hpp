//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef CO_SIM_IO_VECTOR_TO_PYHON_INCLUDED
#define CO_SIM_IO_VECTOR_TO_PYHON_INCLUDED

// Exposure of the CoSimIO to Python

// System includes
#include <vector>
#include <sstream>

// pybind includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

//
#include "includes/define.hpp"

namespace CoSimIO {

// small wrapper around std::vector for using it in Python
// Note: deliberately not using PYBIND11_MAKE_OPAQUE as this
// can be problematic if integrated in other libraries
template<typename TDataType>
class VectorWrapper
{
public:
    VectorWrapper() = default;
    explicit VectorWrapper(const std::size_t I_Size) { mVector.resize(I_Size); }
    explicit VectorWrapper(const std::vector<TDataType>& I_Vector)
        : mVector(I_Vector) {};

    VectorWrapper(const VectorWrapper& Other) : mVector(Other.mVector) {}
    VectorWrapper& operator=(const VectorWrapper&) = delete;

    std::vector<TDataType>& Vector() {return mVector;}
    const std::vector<TDataType>& Vector() const {return mVector;}

private:
    std::vector<TDataType> mVector;
};

} // namespace CoSimIO

namespace { // anonymous namespace

template<typename TDataType>
void AddVectorWithTypeToPython(pybind11::module& m, const std::string& Name)
{
    namespace py = pybind11;

    using VectorType = CoSimIO::VectorWrapper<TDataType>;

    const std::string full_name = Name+"Vector";

    py::class_<VectorType>(m, full_name.c_str())
        .def(py::init<>())
        .def(py::init<const VectorType&>())
        .def(py::init( [](const py::list& l){
            VectorType vec(l.size());
            for(std::size_t i=0; i<l.size(); ++i) {
                vec.Vector()[i] = py::cast<TDataType>(l[i]);
            }
            return vec;
        }))
        .def(py::init( [Name](py::buffer b){
            py::buffer_info info = b.request();
            CO_SIM_IO_ERROR_IF(info.format != py::format_descriptor<TDataType>::value) << "Expected a " << Name << " array!";
            CO_SIM_IO_ERROR_IF(info.ndim != 1) << "Buffer dimension of 1 is required, got: " << info.ndim << std::endl;
            VectorType vec(info.shape[0]);
            for (int i=0; i<info.shape[0]; ++i) {
                vec.Vector()[i] = static_cast<TDataType *>(info.ptr)[i];
            }
            return vec;
        }))

        .def("__len__", [](VectorType& v)
            { return v.Vector().size(); } )
        .def("size", [](VectorType& v)
            { return v.Vector().size(); } )
        .def("resize", [](VectorType& v, const std::size_t& size)
            { v.Vector().resize(size); } )
        .def("append", [](VectorType& v, const TDataType& val)
            { v.Vector().push_back(val); } )
        .def("__setitem__", [](VectorType& v, const std::size_t I_Index, const TDataType& val)
            { v.Vector()[I_Index] = val; } )
        .def("__getitem__", [](VectorType& v, const std::size_t I_Index)
            { return v.Vector()[I_Index]; } )
        .def("__iter__",     [](VectorType& v)
            {return py::make_iterator(v.Vector().begin(), v.Vector().end());},  py::keep_alive<0,1>())

        .def("__str__",   [](const VectorType& v)
            {
                std::stringstream ss;
                const std::size_t size = v.Vector().size();

                ss << "[";
                if(size>0) ss << v.Vector()[0];
                if(size>1) {
                    for(std::size_t i=1; i<size; ++i)
                        ss<<", "<<v.Vector()[i];
                }
                ss << "]";

                return ss.str(); } )
        ;
}

} // anonymous namespace

namespace CoSimIO {

void AddCoSimIOVectorToPython(pybind11::module& m)
{
    AddVectorWithTypeToPython<int>(m, "Int");
    AddVectorWithTypeToPython<double>(m, "Double");
}

} // namespace CoSimIO

#endif // CO_SIM_IO_VECTOR_TO_PYHON_INCLUDED
