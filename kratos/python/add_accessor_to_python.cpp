//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "python/add_accessor_to_python.h"
#include "includes/properties.h"
#include "includes/accessor.h"

namespace Kratos::Python
{

void AddAccessorToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<Accessor, Accessor::UniquePointer>(m, "Accessor")
        .def(py::init<>())
        // .def("GetValue", [&](Accessor &self, Variable<double> &rVariable, Properties &rProperties, Accessor::GeometryType &rGeometry, Vector &rShapeFunctionVector, ProcessInfo &rProcessInfo)
        //     { 
        //         return self.GetValue(rVariable, rProperties, rGeometry, rShapeFunctionVector, rProcessInfo); 
        //     })
        // .def("GetValue", [&](Accessor &self, Variable<int> &rVariable, Properties &rProperties, Accessor::GeometryType &rGeometry, Vector &rShapeFunctionVector, ProcessInfo &rProcessInfo)
        //     { 
        //         return self.GetValue(rVariable, rProperties, rGeometry, rShapeFunctionVector, rProcessInfo); 
        //     })
        // .def("GetValue", [&](Accessor &self, Variable<bool> &rVariable, Properties &rProperties, Accessor::GeometryType &rGeometry, Vector &rShapeFunctionVector, ProcessInfo &rProcessInfo)
        //     { 
        //         return self.GetValue(rVariable, rProperties, rGeometry, rShapeFunctionVector, rProcessInfo); 
        //     })
        // .def("GetValue", [&](Accessor &self, Variable<Vector> &rVariable, Properties &rProperties, Accessor::GeometryType &rGeometry, Vector &rShapeFunctionVector, ProcessInfo &rProcessInfo)
        //     { 
        //         return self.GetValue(rVariable, rProperties, rGeometry, rShapeFunctionVector, rProcessInfo); 
        //     })
        // .def("GetValue", [&](Accessor &self, Variable<Matrix> &rVariable, Properties &rProperties, Accessor::GeometryType &rGeometry, Vector &rShapeFunctionVector, ProcessInfo &rProcessInfo)
        //     { 
        //         return self.GetValue(rVariable, rProperties, rGeometry, rShapeFunctionVector, rProcessInfo); 
        //     })
        // .def("GetValue", [&](Accessor &self, Variable<std::string> &rVariable, Properties &rProperties, Accessor::GeometryType &rGeometry, Vector &rShapeFunctionVector, ProcessInfo &rProcessInfo)
        //     { 
        //         return self.GetValue(rVariable, rProperties, rGeometry, rShapeFunctionVector, rProcessInfo); 
        //     })
        // .def("GetValue", [&](Accessor &self, Variable<array_1d<double, 3>> &rVariable, Properties &rProperties, Accessor::GeometryType &rGeometry, Vector &rShapeFunctionVector, ProcessInfo &rProcessInfo)
        //     { 
        //         return self.GetValue(rVariable, rProperties, rGeometry, rShapeFunctionVector, rProcessInfo); 
        //     })
        // .def("GetValue", [&](Accessor &self, Variable<array_1d<double, 4>> &rVariable, Properties &rProperties, Accessor::GeometryType &rGeometry, Vector &rShapeFunctionVector, ProcessInfo &rProcessInfo)
        //     { 
        //         return self.GetValue(rVariable, rProperties, rGeometry, rShapeFunctionVector, rProcessInfo); 
        //     })
        // .def("GetValue", [&](Accessor &self, Variable<array_1d<double, 6>> &rVariable, Properties &rProperties, Accessor::GeometryType &rGeometry, Vector &rShapeFunctionVector, ProcessInfo &rProcessInfo)
        //     { 
        //         return self.GetValue(rVariable, rProperties, rGeometry, rShapeFunctionVector, rProcessInfo); 
        //     })
        // .def("GetValue", [&](Accessor &self, Variable<array_1d<double, 9>> &rVariable, Properties &rProperties, Accessor::GeometryType &rGeometry, Vector &rShapeFunctionVector, ProcessInfo &rProcessInfo)
        //     { 
        //         return self.GetValue(rVariable, rProperties, rGeometry, rShapeFunctionVector, rProcessInfo); 
        //     })
        ;
}

}  // namespace Kratos::Python.

