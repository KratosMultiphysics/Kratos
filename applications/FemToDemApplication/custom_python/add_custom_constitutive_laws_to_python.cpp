//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

// System includes
#include <pybind11/pybind11.h>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"

#include "python/variable_indexing_python.h"
#include "python/add_mesh_to_python.h"


//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

//constitutive laws
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_constitutive/linear_plane_stress.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_plane_strain_2d.h"

namespace Kratos
{
namespace Python
{
    void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
    {
        py::class_<LinearPlaneStress, typename LinearPlaneStress::Pointer, ConstitutiveLaw >
            (m, "LinearPlaneStress").def(py::init<>() )
            ;

        py::class_<LinearPlaneStrain, typename LinearPlaneStrain::Pointer, ConstitutiveLaw >
            (m, "LinearPlaneStrain").def(py::init<>() )
            ;

        py::class_<ElasticIsotropic3D, typename ElasticIsotropic3D::Pointer, ConstitutiveLaw >
            (m, "ElasticIsotropic3D").def(py::init<>() )
            ;

        py::class_<HyperElasticIsotropicNeoHookean3D, typename HyperElasticIsotropicNeoHookean3D::Pointer, ConstitutiveLaw >
            (m, "HyperElasticIsotropicNeoHookean3D").def(py::init<>() )
            ;

        py::class_<HyperElasticIsotropicNeoHookeanPlaneStrain2D, typename HyperElasticIsotropicNeoHookeanPlaneStrain2D::Pointer, ConstitutiveLaw >
            (m, "HyperElasticIsotropicNeoHookeanPlaneStrain2D").def(py::init<>() )
            ;
    }
}  // namespace Python.
}  // namespace Kratos.