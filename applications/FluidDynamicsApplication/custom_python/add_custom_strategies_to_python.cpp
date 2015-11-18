/*
==============================================================================
KratosFluidDynamicsApplication
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
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


// System includes

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_utilities/solver_settings.h"

#include "spaces/ublas_space.h"

// builder_and_solvers
#include "custom_strategies/builder_and_solvers/residualbased_block_builder_and_solver_periodic.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/fs_strategy.h"
#include "custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme_turbulent.h"
#include "custom_strategies/strategies/residualbased_predictorcorrector_velocity_bdf_scheme_turbulent.h"
#include "custom_strategies/strategies/gear_scheme.h"

//linear solvers
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

namespace Python
{
using namespace boost::python;

void  AddCustomStrategiesToPython()
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;

    //********************************************************************
    //********************************************************************

    class_< ResidualBasedBlockBuilderAndSolverPeriodic< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            bases< ResidualBasedBlockBuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > >,
            boost::noncopyable >
            ("ResidualBasedBlockBuilderAndSolverPeriodic", init<LinearSolverType::Pointer, const Variable<int>& >());

    class_< FSStrategy< SparseSpaceType,LocalSpaceType, LinearSolverType >, bases<BaseSolvingStrategyType>, boost::noncopyable >
            ("FSStrategy",init<ModelPart&,LinearSolverType::Pointer,LinearSolverType::Pointer,bool,bool,double,double,int,int,unsigned int,unsigned int,bool>())
            .def(init< ModelPart&, SolverSettings< SparseSpaceType,LocalSpaceType, LinearSolverType >&, bool >() )
            .def(init< ModelPart&, SolverSettings< SparseSpaceType,LocalSpaceType, LinearSolverType >&, bool, const Kratos::Variable<int>& >() )
            .def("CalculateReactions",&FSStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::CalculateReactions)
            .def("AddIterationStep",&FSStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::AddIterationStep)
            .def("ClearExtraIterationSteps",&FSStrategy<SparseSpaceType,LocalSpaceType,LinearSolverType>::ClearExtraIterationSteps)
            ;


    class_< ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent< SparseSpaceType, LocalSpaceType >,
            bases< BaseSchemeType >,  boost::noncopyable >
            ("ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent",init<double,double,unsigned int,Process::Pointer >() )
            .def(init<double,double,unsigned int >())// constructor without a turbulence model
            .def(init<double,double,unsigned int,const Kratos::Variable<int>&>())// constructor without a turbulence model for periodic boundary conditions
            .def(init<double,double,unsigned int,Kratos::Variable<double>&>())// constructor with a non-default flag for slip conditions
            ;


    class_< ResidualBasedPredictorCorrectorBDFSchemeTurbulent< SparseSpaceType, LocalSpaceType >,
            bases< BaseSchemeType >,  boost::noncopyable >
            ("ResidualBasedPredictorCorrectorBDFSchemeTurbulent",init<unsigned int,Process::Pointer >() )
            .def(init<unsigned int >())// constructor without a turbulence model
            .def(init<unsigned int,Kratos::Variable<double>&>())// constructor with a non-default flag for slip conditions
            ;
            
    class_< GearScheme< SparseSpaceType, LocalSpaceType >,
            bases< BaseSchemeType >,  boost::noncopyable >
            ("GearScheme",init<>()) // default constructor
            .def(init< Process::Pointer >()) // constructor passing a turbulence model
            ;

}

}  // namespace Python.

} // Namespace Kratos

