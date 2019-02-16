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
//   Last Modified by:    $Author: Nelson $
//   Date:                $Date: 2009-02-02 14:03:23 $
//   Revision:            $Revision: 1.5 $
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "python/variable_indexing_python.h"
#include "custom_python/add_fluency_criteria_to_python.h"
#include "fluency_criteria/fluency_criteria.h"
#include "fluency_criteria/energy_yield_function.h"
#include "fluency_criteria/isotropic_rankine_yield_function.h"
//#include "fluency_criteria/tresca_yield_function.h"
#include "fluency_criteria/von_mises_yield_function.h"
#include "fluency_criteria/modified_morh_coulomb_yield_function.h"
#include "fluency_criteria/morh_coulomb_yield_function.h"
#include "fluency_criteria/standard_morh_coulomb_yield_function.h"

//#include "fluency_criteria/rankine_yield_function.h"
//#include "fluency_criteria/drucker_prager_yield_function.h"

#include "soft_hard_behavior/softening_hardening_criteria.h"
#include "soft_hard_behavior/exponencial_softening.h"
#include "soft_hard_behavior/linear_softening.h"
#include "soft_hard_behavior/cohesion_softening.h"
#include "soft_hard_behavior/friction_softening.h"
#include "soft_hard_behavior/dilatancy_softening.h"


#include "spaces/ublas_space.h"





namespace Kratos
{

namespace Python
{

typedef FluencyCriteria                   FluencyCriteriaType;
typedef SofteningHardeningCriteria        SofteningHardeningType;
typedef Morh_Coulomb_Yield_Function       MorhCoulombType;
typedef Isotropic_Rankine_Yield_Function  RankineType;


typedef FluencyCriteria::Pointer            FluencyPointerType;
typedef SofteningHardeningCriteria::Pointer SofteningHardeningPointerType;
typedef MorhCoulombType::Pointer            MorhCoulombPointerType;
typedef RankineType::Pointer                RankinePointerType;

void  AddFluencyCriteriaToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_< FluencyCriteriaType >
    (m, "FluencyCriteriaType")
    .def( py::init<>() )
    ;

    py::class_< SofteningHardeningType >
    (m, "SofteningHardeningType")
    .def( py::init<>() )
    ;

    py::enum_<myState>(m, "State")
    .value("Plane_Stress", Plane_Stress)
    .value("Plane_Strain", Plane_Strain)
    .value("Tri_D", Tri_D)
    ;

    py::enum_<myPotencialPlastic>(m, "PotencialPlastic")
    .value("Not_Associated", Not_Associated)
    .value("Associated", Associated)
    ;


    /*
    py::class_<Rankine_Yield_Function, bases< FluencyCriteriaType > >
    ("RankineYieldFunction",
    py::init<myState> () )
    ;
    */

    py::class_<Isotropic_Rankine_Yield_Function, FluencyCriteriaType >
    (m, "IsotropicRankineYieldFunction")
    .def( py::init<SofteningHardeningPointerType, myState> () )
    ;

    /*
    py::class_<Tresca_Yield_Function, bases< FluencyCriteriaType > >
    ("TrescaYieldFunction",
    py::init<myState> () )
    ;
    */

    py::class_<Von_Mises_Yield_Function, FluencyCriteriaType >
    (m, "VonMisesYieldFunction")
    .def( py::init<const myState&, const SofteningHardeningPointerType&> () )
    //py::init<myState, myPotencialPlastic> () )
    ;

    py::class_<Modified_Morh_Coulomb_Yield_Function, FluencyCriteriaType >
    (m, "ModifiedMorhCoulombYieldFunction")
    .def( py::init<myState, MorhCoulombPointerType, RankinePointerType> () )
    /* py::init<const SofteningHardeningPointerType&,
         const SofteningHardeningPointerType&,
         const SofteningHardeningPointerType&,
         const SofteningHardeningPointerType&,
         const myState,
         const myPotencialPlastic> () )
         */
    ;


    py::class_<Morh_Coulomb_Yield_Function, FluencyCriteriaType >
    (m, "MorhCoulombYieldFunction")
    .def(py::init<const SofteningHardeningPointerType&,
     const SofteningHardeningPointerType&,
     const SofteningHardeningPointerType&,
     const myState,
     const myPotencialPlastic> () )
    ;

    py::class_<Standard_Morh_Coulomb_Yield_Function, FluencyCriteriaType >
    (m, "StandardMorhCoulombYieldFunction")
    .def(py::init<const SofteningHardeningPointerType&,
     const myState> () )
    ;

    /*
    py::class_<Drucker_Prager_Yield_Function, bases< FluencyCriteriaType > >
    ("DruckerPragerYieldFunction",
    py::init<myState> () )
    ;
                 */


    py::class_<Energy_Yield_Function, FluencyCriteriaType >
    (m, "EnergyYieldFunction")
    .def( py::init<myState> () )
    ;


    py::class_<Exponential_Softening, SofteningHardeningType >
    (m, "ExponentialSoftening")
    .def( py::init<> () )
    ;

    py::class_<Linear_Softening, SofteningHardeningType >
    (m, "LinearSoftening")
    .def( py::init<> () )
    ;

    py::class_<Cohesion_Softening, SofteningHardeningType >
    (m, "CohesionSoftening")
    .def( py::init<> () )
    ;

    py::class_<Friction_Softening, SofteningHardeningType >
    (m, "FrictionSoftening")
    .def( py::init<> () )
    ;

    py::class_<Dilatancy_Softening, SofteningHardeningType >
    (m, "DilatancySoftening")
    .def( py::init<> () )
    ;


}


}  // namespace Python.

}  // namespace Kratos.

