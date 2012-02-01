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
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "custom_python/add_constitutive_laws_to_python.h"
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "constitutive_laws/isotropic_2d.h"
#include "constitutive_laws/isotropic_3d.h"
#include "constitutive_laws/hyperelastic_3d.h"
#include "constitutive_laws/hyperelastic_2d.h"
// #include "constitutive_laws/viscoelastic_2d.h" // new VISCOELASTICITY
// #include "constitutive_laws/viscofibers_2d.h" // new VISCOELASTIC Fibers
// #include "constitutive_laws/viscofibers_hypermatrix_2d.h" // new VISCOELASTIC Fibers and Hyperelastic Matrix
#include "constitutive_laws/von_mises_3d.h"
#include "constitutive_laws/hypoelastic_2d.h"
#include "constitutive_laws/plane_strain.h"
#include "constitutive_laws/fluid_2d.h"
#include "constitutive_laws/external_isotropic_3d.h"
#include "constitutive_laws/drucker_prager.h"
//#include "constitutive_laws/isotropic_elastic_large_strain.h"
#include "constitutive_laws/hooks_law.h"
#include "constitutive_laws/drucker_prager_law.h"
#include "constitutive_laws/isotropic_planestress_wrinkling.h"
#include "constitutive_laws/isotropic_damage_2d.h"
#include "constitutive_laws/isotropic_rankine_damage_2d.h"
#include "constitutive_laws/isotropic_damage_3d.h"
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
#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"
#include "fluency_criteria/fluency_criteria.h"
#include "constitutive_laws/mohr_coulomb_plane_strain.h"

#include "constitutive_laws/von_mises_3d.h"


namespace Kratos
{
    namespace Python
    {
        using namespace boost::python;
        
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
            ThisMaterialsContainer.push_back(ThisConstitutiveLaw);
        }
        
        void  AddConstitutiveLawsToPython()
        {
            class_< MaterialsContainer >("MaterialsContainer", init<>() )
                    .def("PushBack", Push_Back_Constitutive_Laws)
                    ;
            
            class_< Isotropic2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("Isotropic2D",
                     init<>() )
                    // .def("Clone",              &Isotropic2D::Clone)
                    ;
            
            class_< PlaneStrain, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("PlaneStrain",
                     init<>() )
                    ;
            
            class_< MohrCoulombPlaneStrain, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("MohrCoulombPlaneStrain",
                     init<>() )
                    ;
            
            class_< Isotropic3D, bases< ConstitutiveLawBaseType >,  boost::noncopyable >
                    ("Isotropic3D", 
                     init<>() ) 
                    ;
            
	    class_< Orthotropic3D, bases< ConstitutiveLawBaseType >,  boost::noncopyable >
                    ("Orthotropic3D", 
                     init<>() ) 
                    ;
		    
            class_< Isotropic_Damage_2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("IsotropicDamage2D",
                     init<>() )
                    .def(init<FluencyCriteriaPointer,SofteningHardeningCriteriaPointer, PropertiesPointer>())
                    ;
            
            class_<Plasticity2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("Plasticity2D",
                     init<>() )
                    .def(init<FluencyCriteriaPointer, PropertiesPointer>())
                    ;
            
            class_<PlaneStressJ2, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("PlaneStressJ2",
                     init<>() )
                    ;
            
//             class_<Plasticity3D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
//             ("Plasticity3D",
//             init<>() )
//                 .def(init<FluencyCriteriaPointer,SofteningHardeningCriteriaPointer, PropertiesPointer>())
//                  ;

            

            

            class_<BrittleMaterial2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("BrittleMaterial2D",
                     init<>() )
                     .def(init<FluencyCriteriaPointer, PropertiesPointer>())
                    ;
            
            
	    class_<IsotropicRankineDamage2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("IsotropicRankineDamage2D",
                     init<>() )
                     .def(init<>())
                    ;
		    

            class_< VonMises3D, bases< ConstitutiveLawBaseType >,  boost::noncopyable >
                    ("VonMises3D", 
                     init<>() ) 
                    ;
            
            class_< Hypoelastic2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("Hypoelastic2D",
                     init<>() )
                    ;
            
            class_< Fluid2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("Fluid2D",
                     init<>() )
                    ;
            
            class_< ExternalIsotropic3D, bases< ConstitutiveLawBaseType >,  boost::noncopyable >
                    ("ExternalIsotropic3D", 
                     init<>() ) 
                    ;
            
            class_< HooksLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("HooksLaw",
                     init<>() )
                    ;
            
            class_< DruckerPragerLaw, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("DruckerPrager",
                     init<>() )
                    ;
            
            class_< IsotropicPlaneStressWrinkling, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("IsotropicPlaneStressWrinkling",
                     init<>() )
                    ;
            
	    class_< Hyperelastic3D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
			("Hyperelastic3D",
			init<>() )
			;           

             
            class_< Hyperelastic2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("Hyperelastic2D",
                     init<>() )
                    ;

	    /*    
            class_< Viscofibers2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("Viscofibers2D",
                     init<>() )
                    ;


            class_< Viscoelastic2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("Viscoelastic2D",
                     init<>() )
                    ;

    
	    class_< Viscofibers_Hypermatrix2D, bases< ConstitutiveLawBaseType >, boost::noncopyable >
                    ("Viscofibers_Hypermatrix2D",
                     init<>() )
                    ;
		    */
// 			class_<Plane_Stress_Damage_Orthotropic_2D  , bases< ConstitutiveLawBaseType >, boost::noncopyable >
// 			("PlaneStressDamageOrthotropic2D",
// 			init<>() )
// 			//.def(init<FluencyCriteriaType const&>())
//                         .def(init<FluencyCriteriaPointer>())
// 			;
/*
			class_<ComposeMaterial , bases< ConstitutiveLawBaseType >, boost::noncopyable >
			("ComposeMaterial",
			init<>() )
                        .def(init<MaterialsContainer>())
			;*/

        
        }
    }  // namespace Python.
}  // namespace Kratos.
#endif // KRATOS_ADD_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED defined
