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
#include "../impl/define.hpp"

// has to be done outside of any namespace!
PYBIND11_MAKE_OPAQUE(std::vector<int>)
PYBIND11_MAKE_OPAQUE(std::vector<double>)

namespace {

template<typename TDataType>
void AddVectorWithTypeToPython(pybind11::module& m, const std::string& Name)
{
    namespace py = pybind11;

    using VectorType = std::vector<TDataType>;

    const std::string full_name = Name+"Vector";

    py::class_<VectorType>(m, full_name.c_str())
        .def(py::init<>())
        .def(py::init<const VectorType&>())
        .def(py::init( [](const py::list& l){
            VectorType vec(l.size());
            for(std::size_t i=0; i<l.size(); ++i) {
                vec[i] = py::cast<TDataType>(l[i]);
            }
            return vec;
        }))
        .def(py::init( [Name](py::buffer b){
            py::buffer_info info = b.request();
            CO_SIM_IO_ERROR_IF(info.format != py::format_descriptor<TDataType>::value) << "Expected a " << Name << " array!";
            CO_SIM_IO_ERROR_IF(info.ndim != 1) << "Buffer dimension of 1 is required, got: " << info.ndim << std::endl;
            VectorType vec(info.shape[0]);
            for (int i=0; i<info.shape[0]; ++i) {
                vec[i] = static_cast<TDataType *>(info.ptr)[i];
            }
            return vec;
        }))

        .def("__len__", [](VectorType& v)
            { return v.size(); } )
        .def("size", [](VectorType& v)
            { return v.size(); } )
        .def("resize", [](VectorType& v, const std::size_t& size)
            { v.resize(size); } )
        .def("append", [](VectorType& v, const TDataType& val)
            { v.push_back(val); } )
        .def("__setitem__", [](VectorType& v, const std::size_t I_Index, const TDataType& val)
            { v[I_Index] = val; } )
        .def("__getitem__", [](VectorType& v, const std::size_t I_Index)
            { return v[I_Index]; } )
        .def("__iter__",     [](VectorType& v)
            {return py::make_iterator(v.begin(), v.end());},  py::keep_alive<0,1>())

        .def("__str__",   [](const VectorType& v)
            {
                std::stringstream ss;
                const std::size_t size = v.size();

                ss << "[";
                if(size>0) ss << v[0];
                if(size>1) {
                    for(std::size_t i=1; i<size; ++i)
                        ss<<", "<<v[i];
                }
                ss << "]";

                return ss.str(); } )
        ;
}

}

void AddCoSimIOVectorToPython(pybind11::module& m)
{
    AddVectorWithTypeToPython<int>(m, "Int");
    AddVectorWithTypeToPython<double>(m, "Double");
}

#endif // CO_SIM_IO_VECTOR_TO_PYHON_INCLUDED
