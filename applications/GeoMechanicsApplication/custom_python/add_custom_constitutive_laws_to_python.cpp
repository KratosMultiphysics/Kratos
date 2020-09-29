// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//


// System includes

// Project includes
#include "includes/constitutive_law.h"

//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

//constitutive laws
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"
#include "custom_constitutive/bilinear_cohesive_2D_law.hpp"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    class_< BilinearCohesive3DLaw, BilinearCohesive3DLaw::Pointer, ConstitutiveLaw >
    (m, "BilinearCohesive3DLaw")
    .def( init<>() );
    class_< BilinearCohesive2DLaw, BilinearCohesive2DLaw::Pointer, ConstitutiveLaw >
    (m, "BilinearCohesive2DLaw")
    .def( init<>() ) ;

}

}  // namespace Python.
}  // namespace Kratos.
