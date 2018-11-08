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
//                   Ruben Zorrilla
//


// System includes

// External includes


// Project includes
#include "includes/define_python.h"
#include "includes/constitutive_law.h"

//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

// 3D constitutive laws
#include "custom_constitutive/euler_3d_law.h"
#include "custom_constitutive/bingham_3d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"
#include "custom_constitutive/herschel_bulkley_3d_law.h"
#include "custom_constitutive/newtonian_two_fluid_3d_law.h"

// 2D constitutive laws
#include "custom_constitutive/euler_2d_law.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_two_fluid_2d_law.h"

namespace Kratos
{

namespace Python
{

void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_< Euler2DLaw, Euler2DLaw::Pointer, ConstitutiveLaw >(m,"Euler2DLaw")
    .def(  py::init<>() );

    py::class_< Euler3DLaw, Euler3DLaw::Pointer, ConstitutiveLaw >(m,"Euler3DLaw")
    .def( py::init<>() );

    py::class_< Bingham3DLaw, Bingham3DLaw::Pointer, ConstitutiveLaw >(m,"Bingham3DLaw")
    .def( py::init<>() );

    py::class_< Newtonian2DLaw, Newtonian2DLaw::Pointer, ConstitutiveLaw >(m,"Newtonian2DLaw")
    .def( py::init<>() );

    py::class_< Newtonian3DLaw, Newtonian3DLaw::Pointer, ConstitutiveLaw >(m,"Newtonian3DLaw")
    .def( py::init<>() );

    py::class_< HerschelBulkley3DLaw, HerschelBulkley3DLaw::Pointer, ConstitutiveLaw >(m,"HerschelBulkley3DLaw")
    .def( py::init<>() );

    py::class_< NewtonianTwoFluid3DLaw, NewtonianTwoFluid3DLaw::Pointer, ConstitutiveLaw >(m,"NewtonianTwoFluid3DLaw")
    .def( py::init<>() );

    py::class_< NewtonianTwoFluid2DLaw, NewtonianTwoFluid2DLaw::Pointer, ConstitutiveLaw >(m,"NewtonianTwoFluid2DLaw")
    .def( py::init<>() );
}

}  // namespace Python.
}  // namespace Kratos.
