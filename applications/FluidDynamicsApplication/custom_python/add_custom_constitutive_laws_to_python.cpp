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

namespace Kratos
{

namespace Python
{

using namespace pybind11;


void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{

    class_< Euler2DLaw, Euler2DLaw::Pointer, ConstitutiveLaw >(m,"Euler2DLaw")
    .def(  init<>() );

    class_< Euler3DLaw, Euler3DLaw::Pointer, ConstitutiveLaw >(m,"Euler3DLaw")
    .def( init<>() );
    
    class_< Bingham3DLaw, Bingham3DLaw::Pointer, ConstitutiveLaw >(m,"Bingham3DLaw")
    .def( init<>() );

    class_< Newtonian2DLaw, Newtonian2DLaw::Pointer, ConstitutiveLaw >(m,"Newtonian2DLaw")
    .def( init<>() );
     
    class_< Newtonian3DLaw, Newtonian3DLaw::Pointer, ConstitutiveLaw >(m,"Newtonian3DLaw")
    .def( init<>() );
  
    class_< HerschelBulkley3DLaw, HerschelBulkley3DLaw::Pointer, ConstitutiveLaw >(m,"HerschelBulkley3DLaw")
    .def( init<>() );

    class_< NewtonianTwoFluid3DLaw, NewtonianTwoFluid3DLaw::Pointer, ConstitutiveLaw >(m,"NewtonianTwoFluid3DLaw")
    .def( init<>() );
    
}

}  // namespace Python.
}  // namespace Kratos.
