/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hurga $
//   Date:                $Date: 2009-02-02 14:03:23 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_ADD_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED


// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "custom_python/add_constitutive_laws_to_python.h"
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "constitutive_laws/dummy_constitutive_law.h"
#include "constitutive_laws/tutorial_damage_model.h"
#include "constitutive_laws/isotropic_2d.h"
#include "constitutive_laws/isotropic_3d.h"
#include "constitutive_laws/neo_hookean_3d.h"
#include "constitutive_laws/hyperelastic_3d.h"
#include "constitutive_laws/hyperelastic_2d.h"
// #include "constitutive_laws/viscoelastic_2d.h" // new VISCOELASTICITY
// #include "constitutive_laws/viscofibers_2d.h" // new VISCOELASTIC Fibers
// #include "constitutive_laws/viscofibers_hypermatrix_2d.h" // new VISCOELASTIC Fibers and Hyperelastic Matrix
#include "constitutive_laws/von_mises_3d.h"
#include "constitutive_laws/hypoelastic_2d.h"
#include "constitutive_laws/plane_strain.h"
#include "constitutive_laws/plane_stress.h"
#include "constitutive_laws/fluid_2d.h"
#include "constitutive_laws/external_isotropic_3d.h"
#include "constitutive_laws/drucker_prager.h"
#include "constitutive_laws/cam_clay_3d.h"
//#include "constitutive_laws/isotropic_elastic_large_strain.h"
#include "constitutive_laws/hooks_law.h"
#include "constitutive_laws/isotropic_planestress_wrinkling.h"
#include "constitutive_laws/isotropic_damage_2d.h"
#include "constitutive_laws/isotropic_rankine_damage_2d.h"
#include "constitutive_laws/isotropic_rankine_damage_3d.h"
#include "constitutive_laws/isotropic_damage_3d.h"
#include "constitutive_laws/isotropic_damage_implex.h"
#include "constitutive_laws/plasticity_2d.h"
#include "constitutive_laws/plane_stress_J2.h"
#include "constitutive_laws/brittle_material_2d.h"
#include "constitutive_laws/orthotropic_3d.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"
#include "python/add_mesh_to_python.h"
#include "python/variable_indexing_python.h"
#include "fluency_criteria/fluency_criteria.h"
#include "constitutive_laws/mohr_coulomb_plane_strain.h"

#include "constitutive_laws/von_mises_3d.h"


namespace Kratos
{

namespace Python
{

typedef ConstitutiveLaw ConstitutiveLawBaseType;
typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;
typedef FluencyCriteria::Pointer FluencyCriteriaPointer;
typedef SofteningHardeningCriteria::Pointer SofteningHardeningCriteriaPointer;
typedef Properties::Pointer PropertiesPointer;

typedef std::vector<ConstitutiveLaw::Pointer> MaterialsContainer;
typedef ConstitutiveLaw::Pointer  ConstitutiveLawPointer;

void Push_Back_Constitutive_Laws( MaterialsContainer& ThisMaterialsContainer,
                                  ConstitutiveLawPointer ThisConstitutiveLaw )
{
    ThisMaterialsContainer.push_back( ThisConstitutiveLaw );
}

void  AddConstitutiveLawsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_< MaterialsContainer >(m, "MaterialsContainer")
    .def(py::init<>() )
    .def( "PushBack", Push_Back_Constitutive_Laws )
    ;

    py::class_< DummyConstitutiveLaw, ConstitutiveLawBaseType >
    (m, "DummyConstitutiveLaw" )
    .def( py::init<>() )
    ;

    py::class_< TutorialDamageModel, ConstitutiveLawBaseType >
    (m, "TutorialDamageModel" )
    .def( py::init<>() )
    ;

    py::class_< Isotropic2D, ConstitutiveLawBaseType >
    (m, "Isotropic2D" )
    .def( py::init<>() )
    // .def("Clone",              &Isotropic2D::Clone)
    ;

    py::class_< PlaneStrain, ConstitutiveLawBaseType >
    (m, "PlaneStrain" )
    .def( py::init<>() )
    ;

    py::class_< PlaneStress, ConstitutiveLawBaseType >
    (m, "PlaneStress" )
    .def( py::init<>() )
    ;

    py::class_< MohrCoulombPlaneStrain, ConstitutiveLawBaseType >
    (m, "MohrCoulombPlaneStrain")
    .def( py::init<>() )
    ;

    py::class_< Isotropic3D, ConstitutiveLawBaseType >
    (m, "Isotropic3D" )
    .def( py::init<>() )
    ;

    py::class_< DruckerPrager, ConstitutiveLawBaseType >
    (m, "DruckerPrager" )
    .def( py::init<>() )
    ;

    py::class_< Orthotropic3D, ConstitutiveLawBaseType >
    (m, "Orthotropic3D" )
    .def( py::init<>() )
    ;

    py::class_< Isotropic_Damage_2D, ConstitutiveLawBaseType >
    (m, "IsotropicDamage2D" )
    .def( py::init<>() )
    .def( py::init<FluencyCriteriaPointer, SofteningHardeningCriteriaPointer, PropertiesPointer>() )
    ;

    py::class_< Isotropic_Damage_3D, ConstitutiveLawBaseType >
    (m, "IsotropicDamage3D" )
    .def( py::init<>() )
    .def( py::init<FluencyCriteriaPointer, SofteningHardeningCriteriaPointer, PropertiesPointer>() )
    ;

    py::class_< IsotropicDamageIMPLEX, ConstitutiveLawBaseType >
    (m, "IsotropicDamageIMPLEX" )
    .def( py::init<>() )
    ;


    py::class_<Plasticity2D, ConstitutiveLawBaseType >
    (m, "Plasticity2D" )
    .def( py::init<>() )
    .def( py::init<FluencyCriteriaPointer, PropertiesPointer>() )
    ;

    py::class_<PlaneStressJ2, ConstitutiveLawBaseType >
    (m, "PlaneStressJ2" )
    .def( py::init<>() )
    ;

//             py::class_<Plasticity3D, bases< ConstitutiveLawBaseType >
//             ("Plasticity3D",
//             py::init<>() )
//                 .def(py::init<FluencyCriteriaPointer,SofteningHardeningCriteriaPointer, PropertiesPointer>())
//                  ;





    py::class_<BrittleMaterial2D, ConstitutiveLawBaseType >
    (m, "BrittleMaterial2D" )
    .def( py::init<>() )
    .def( py::init<FluencyCriteriaPointer, PropertiesPointer>() )
    ;


    py::class_<IsotropicRankineDamage2D, ConstitutiveLawBaseType >
    (m, "IsotropicRankineDamage2D" )
    .def( py::init<>() )
    ;

    py::class_<IsotropicRankineDamage3D, ConstitutiveLawBaseType >
    (m, "IsotropicRankineDamage3D" )
    .def( py::init<>() )
    ;

    py::class_< VonMises3D, ConstitutiveLawBaseType >
    (m, "VonMises3D" )
    .def( py::init<>() )
    ;

    py::class_< Hypoelastic2D, ConstitutiveLawBaseType >
    (m, "Hypoelastic2D" )
    .def( py::init<>() )
    ;

    py::class_< Fluid2D, ConstitutiveLawBaseType >
    (m, "Fluid2D" )
    .def( py::init<>() )
    ;

    py::class_< ExternalIsotropic3D, ConstitutiveLawBaseType >
    (m, "ExternalIsotropic3D" )
    .def( py::init<>() )
    ;

    py::class_< HooksLaw, ConstitutiveLawBaseType >
    (m, "HooksLaw" )
    .def( py::init<>() )
    ;

    py::class_< IsotropicPlaneStressWrinkling, ConstitutiveLawBaseType >
    (m, "IsotropicPlaneStressWrinkling" )
    .def( py::init<>() )
    ;

    py::class_< Hyperelastic3D, ConstitutiveLawBaseType >
    (m, "Hyperelastic3D" )
    .def( py::init<>() )
    ;


    py::class_< Hyperelastic2D, ConstitutiveLawBaseType >
    (m, "Hyperelastic2D" )
    .def( py::init<>() )
    ;

    py::class_< CamClay3D, ConstitutiveLawBaseType >
    (m, "CamClay3D" )
    .def( py::init<>() )
    ;

    py::class_< NeoHookean3D, ConstitutiveLawBaseType >
    (m, "NeoHookean3D" )
    .def( py::init<>() )
    ;

    /*
           py::class_< Viscofibers2D, bases< ConstitutiveLawBaseType >
                   ("Viscofibers2D",
                    py::init<>() )
                   ;


           py::class_< Viscoelastic2D, bases< ConstitutiveLawBaseType >
                   ("Viscoelastic2D",
                    py::init<>() )
                   ;


    py::class_< Viscofibers_Hypermatrix2D, bases< ConstitutiveLawBaseType >
                   ("Viscofibers_Hypermatrix2D",
                    py::init<>() )
                   ;
     */
//    py::class_<Plane_Stress_Damage_Orthotropic_2D  , bases< ConstitutiveLawBaseType >
//    ("PlaneStressDamageOrthotropic2D",
//    py::init<>() )
//    //.def(py::init<FluencyCriteriaType const&>())
//                         .def(py::init<FluencyCriteriaPointer>())
//    ;
    /*
       py::class_<ComposeMaterial , bases< ConstitutiveLawBaseType >
       ("ComposeMaterial",
       py::init<>() )
                            .def(py::init<MaterialsContainer>())
       ;*/


}
}  // namespace Python.
}  // namespace Kratos.
#endif // KRATOS_ADD_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED defined
