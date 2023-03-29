//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "add_controllers_to_python.h"
#include "controllers/controller.h"

namespace Kratos::Python
{

class ControllerTrampoline : public Controller
{
public:
    //Inherit the constructors
    using Controller::Controller;

    bool Evaluate() override
    {
        using ReturnType = bool;
        using BaseType = Controller;
        PYBIND11_OVERRIDE_PURE(
            ReturnType,
            BaseType,
            Evaluate);
    }

    Controller::Pointer Create(
        Model& rModel,
        Parameters ThisParameters) const override
    {
        using ReturnType = Controller::Pointer;
        using BaseType = Controller;
        PYBIND11_OVERRIDE_PURE(
            ReturnType,
            BaseType,
            Create,
            rModel,
            ThisParameters);
    }
}; // class ControllerTrampoline

void AddControllersToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<Controller, Controller::Pointer, ControllerTrampoline>(m,"Controller")
        .def(py::init<>())
        .def("Create", &Controller::Create)
        .def("Evaluate", &Controller::Evaluate)
        .def("GetDefaultParameters", &Controller::GetDefaultParameters)
        .def("__str__", PrintObject<Controller>)
        ;
}

} // Namespace Kratos::Python
