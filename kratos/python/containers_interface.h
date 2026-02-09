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

#pragma once

// System includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// External includes

// Project includes
#include "includes/define_python.h"
#include "utilities/parallel_utilities.h"

namespace Kratos::Python {
namespace py = pybind11;

template <class TContainerType> class ContainerInterface {
public:
  ContainerInterface() {};
  virtual ~ContainerInterface() {};

  void CreateInterface(pybind11::module &m, std::string ContainerName) {
    py::class_<TContainerType, typename TContainerType::Pointer>(
        m, ContainerName.c_str())
        .def(py::init<>())
        .def("__len__", [](TContainerType &self) { return self.size(); })
        .def(
            "__setitem__",
            [](TContainerType &self, unsigned int i,
               typename TContainerType::value_type &value) { self[i] = value; })
        .def("__getitem__",
             [](TContainerType &self, unsigned int i) { return self(i); })
        .def(
            "__iter__",
            [](TContainerType &self) {
              return py::make_iterator(self.begin(), self.end());
            },
            py::keep_alive<0, 1>())
        .def("append", [](TContainerType &self,
                          typename TContainerType::value_type &value) {
          self.push_back(value);
        });
  }
};

template <class TContainerType> class MapInterface {
public:
  MapInterface() {};
  virtual ~MapInterface() {};

  void CreateInterface(pybind11::module &m, std::string ContainerName) {
    py::class_<TContainerType, typename TContainerType::Pointer>(
        m, ContainerName.c_str())
        .def(py::init<>())
        .def("__len__", [](TContainerType &self) { return self.size(); })
        .def("__contains__",
             [](TContainerType &self,
                const typename TContainerType::value_type &value) {
               return (self.find(value.Id()) != self.end());
             })
        .def("__contains__",
             [](TContainerType &self, unsigned int i) {
               return (self.find(i) != self.end());
             })
        .def(
            "__setitem__",
            [](TContainerType &self, unsigned int i,
               typename TContainerType::value_type &value) { self[i] = value; })
        .def("__getitem__",
             [](TContainerType &self, unsigned int i) { return self(i); })
        .def(
            "__iter__",
            [](TContainerType &self) {
              return py::make_iterator(self.begin(), self.end());
            },
            py::keep_alive<0, 1>())
        .def("clear", [](TContainerType &self) { self.clear(); });
  }
};

template <class TContainerType> class PointerVectorPythonInterface {
public:
  PointerVectorPythonInterface() {};
  virtual ~PointerVectorPythonInterface() {};

  void CreateInterface(pybind11::module &m, std::string ContainerName) {
    py::class_<TContainerType, typename TContainerType::Pointer>(
        m, ContainerName.c_str())
        .def(py::init<>())
        .def("__len__", [](TContainerType &self) { return self.size(); })
        .def(
            "__setitem__",
            [](TContainerType &self, unsigned int i,
               typename TContainerType::value_type &value) { self[i] = value; })
        .def("__getitem__",
             [](TContainerType &self, unsigned int i) { return self(i); })
        .def(
            "__iter__",
            [](TContainerType &self) {
              return py::make_iterator(self.begin(), self.end());
            },
            py::keep_alive<0, 1>())
        .def("append", [](TContainerType &self,
                          typename TContainerType::value_type::Pointer pValue) {
          self.push_back(pValue);
        });
  }
};

template <class TContainerType> class PointerVectorSetPythonInterface {
public:
  PointerVectorSetPythonInterface() {};
  virtual ~PointerVectorSetPythonInterface() {};

  void CreateInterface(pybind11::module &m, std::string ContainerName) {
    py::class_<TContainerType, typename TContainerType::Pointer>(
        m, ContainerName.c_str())
        .def(py::init<>())
        .def("__len__", [](TContainerType &self) { return self.size(); })
        .def("__contains__",
             [](TContainerType &self,
                const typename TContainerType::value_type &value) {
               return (self.find(value.Id()) != self.end());
             })
        .def("__contains__",
             [](TContainerType &self, unsigned int i) {
               return (self.find(i) != self.end());
             })
        .def("__setitem__",
             [](TContainerType &self,
                typename TContainerType::value_type &value) {
               self[value.Id()] = value;
             })
        .def(
            "__setitem__",
            [](TContainerType &self, typename TContainerType::pointer &pvalue) {
              self(pvalue->Id()) = pvalue;
            })
        .def("__getitem__",
             [](TContainerType &self, unsigned int i) { return self(i); })
        .def(
            "__iter__",
            [](TContainerType &self) {
              return py::make_iterator(self.ptr_begin(), self.ptr_end());
            },
            py::keep_alive<0, 1>()) // TODO: decide if here we should use
                                    // ptr_iterators or iterators
        .def("append",
             [](TContainerType &self, typename TContainerType::pointer value) {
               self.push_back(value);
             })
        .def("clear", [](TContainerType &self) { self.clear(); })
        .def("back",
             [](TContainerType &rThis) {
               KRATOS_ERROR_IF(rThis.empty())
                   << "attempting to call \"back()\" on an empty container.";
               return *(rThis.ptr_end() - 1);
             })
        .def("IdsList", [](TContainerType &self) {
          const int n = self.size();
          py::array_t<int> ids(n);
          auto ids_ptr = static_cast<int *>(ids.request().ptr);
          IndexPartition<int>(n).for_each(
              [&](int i) { ids_ptr[i] = (self.ptr_begin() + i)->get()->Id(); });
          return ids;
        });
  }
};

} // namespace Kratos::Python.
