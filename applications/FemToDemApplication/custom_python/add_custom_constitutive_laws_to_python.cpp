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
        py::class_<LinearPlaneStressFEMDEM, typename LinearPlaneStressFEMDEM::Pointer, ConstitutiveLaw >
            (m, "LinearPlaneStressFEMDEM").def(py::init<>() )
            ;

        py::class_<LinearPlaneStrainFEMDEM, typename LinearPlaneStrainFEMDEM::Pointer, ConstitutiveLaw >
            (m, "LinearPlaneStrainFEMDEM").def(py::init<>() )
            ;

        py::class_<ElasticIsotropic3DFEMDEM, typename ElasticIsotropic3DFEMDEM::Pointer, ConstitutiveLaw >
            (m, "ElasticIsotropic3DFEMDEM").def(py::init<>() )
            ;

        py::class_<HyperElasticIsotropicNeoHookean3DFEMDEM, typename HyperElasticIsotropicNeoHookean3DFEMDEM::Pointer, ConstitutiveLaw >
            (m, "HyperElasticIsotropicNeoHookean3DFEMDEM").def(py::init<>() )
            ;

        py::class_<HyperElasticIsotropicNeoHookeanPlaneStrain2DFEMDEM, typename HyperElasticIsotropicNeoHookeanPlaneStrain2DFEMDEM::Pointer, ConstitutiveLaw >
            (m, "HyperElasticIsotropicNeoHookeanPlaneStrain2DFEMDEM").def(py::init<>() )
            ;
    }
}  // namespace Python.
}  // namespace Kratos.