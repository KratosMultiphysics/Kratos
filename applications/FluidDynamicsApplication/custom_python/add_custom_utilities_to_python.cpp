/*
==============================================================================
KratosFluiDynamicsApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


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
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

#include "custom_utilities/dynamic_smagorinsky_utilities.h"
#include "custom_utilities/periodic_condition_utilities.h"
#include "custom_utilities/fractional_step_settings.h"
#include "custom_utilities/fractional_step_settings_periodic.h"
#include "utilities/split_tetrahedra.h"

namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    class_<DynamicSmagorinskyUtils>("DynamicSmagorinskyUtils", init<ModelPart&,unsigned int>())
    .def("StoreCoarseMesh",&DynamicSmagorinskyUtils::StoreCoarseMesh)
    .def("CalculateC",&DynamicSmagorinskyUtils::CalculateC)
    .def("CorrectFlagValues",&DynamicSmagorinskyUtils::CorrectFlagValues)
    ;

    typedef void (PeriodicConditionUtilities::*AddDoubleVariableType)(Properties&,Variable<double>&);
    typedef void (PeriodicConditionUtilities::*AddVariableComponentType)(Properties&,VariableComponent< VectorComponentAdaptor< array_1d<double, 3> > >&);

    AddDoubleVariableType AddDoubleVariable = &PeriodicConditionUtilities::AddPeriodicVariable;
    AddVariableComponentType AddVariableComponent = &PeriodicConditionUtilities::AddPeriodicVariable;

    class_<PeriodicConditionUtilities>("PeriodicConditionUtilities", init<ModelPart&,unsigned int>())
            .def("SetUpSearchStructure",&PeriodicConditionUtilities::SetUpSearchStructure)
            .def("DefinePeriodicBoundary",&PeriodicConditionUtilities::DefinePeriodicBoundary)
            .def("AddPeriodicVariable",AddDoubleVariable)
            .def("AddPeriodicVariable",AddVariableComponent)
    ;

    typedef SolverSettings<SparseSpaceType,LocalSpaceType,LinearSolverType> BaseSettingsType;

    class_ < BaseSettingsType, boost::noncopyable >
    ( "BaseSettingsType",no_init );

    enum_<FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::StrategyLabel>("StrategyLabel")
    .value("Velocity",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::Velocity)
    .value("Pressure",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::Pressure)
    //.value("EddyViscosity",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::EddyViscosity)
    ;

    enum_<FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::TurbulenceModelLabel>("TurbulenceModelLabel")
    .value("SpalartAllmaras",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SpalartAllmaras)
    ;

    typedef void (FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*SetStrategyByParamsType)(FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::StrategyLabel const&,LinearSolverType::Pointer,const double,const unsigned int);
    typedef void (FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*BuildTurbModelType)(FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::TurbulenceModelLabel const&, LinearSolverType::Pointer, const double, const unsigned int);
    typedef void (FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*PassTurbModelType)(Process::Pointer);
    SetStrategyByParamsType ThisSetStrategyOverload = &FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetStrategy;
    BuildTurbModelType SetTurbModel_Build = &FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetTurbulenceModel;
    PassTurbModelType SetTurbModel_Pass = &FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetTurbulenceModel;

    class_< FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>,bases<BaseSettingsType>, boost::noncopyable>
            ("FractionalStepSettings",init<ModelPart&,unsigned int,unsigned int,bool,bool,bool>())
    .def("SetStrategy",ThisSetStrategyOverload)
    .def("SetTurbulenceModel",SetTurbModel_Build)
    .def("SetTurbulenceModel",SetTurbModel_Pass)
    .def("GetStrategy",&FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::pGetStrategy)
    .def("SetEchoLevel",&FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetEchoLevel)
    ;
    
    class_< FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>,bases<BaseSettingsType>, boost::noncopyable>
            ("FractionalStepSettingsPeriodic",init<ModelPart&,unsigned int,unsigned int,bool,bool,bool,const Kratos::Variable<int>&>())
    .def("SetStrategy",&FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetStrategy)
    .def("GetStrategy",&FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::pGetStrategy)
    .def("SetEchoLevel",&FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetEchoLevel)
    ;

}





}  // namespace Python.

} // Namespace Kratos

