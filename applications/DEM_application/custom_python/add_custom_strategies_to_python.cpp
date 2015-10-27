/*
==============================================================================
KratosTestApplication 
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
#include "custom_python/add_custom_strategies_to_python.h"

#include "spaces/ublas_space.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/explicit_solver_strategy.h"
#include "custom_strategies/strategies/explicit_solver_continuum.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

//convergence criteria
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//schemes
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/forward_euler_scheme.h"
#include "custom_strategies/schemes/constant_average_acceleration_scheme.h"
#include "custom_strategies/schemes/mid_point_scheme.h"

//builder_and_solvers
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

//Search
#include "custom_utilities/omp_dem_search.h"
#include "custom_utilities/dem_fem_search.h"

// Create and Destroy
#include "custom_utilities/create_and_destroy.h"

namespace Kratos
{

    namespace Python
    {   
        using namespace boost::python;

        void  AddCustomStrategiesToPython()
        {
          typedef UblasSpace<double, CompressedMatrix, Vector>                          SparseSpaceType;
          typedef UblasSpace<double, Matrix, Vector>                                    LocalSpaceType;

          typedef LinearSolver<SparseSpaceType, LocalSpaceType >                        LinearSolverType;
          typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >  BaseSolvingStrategyType;
          //typedef Scheme< SparseSpaceType, LocalSpaceType >                             BaseSchemeType;
          
          typedef OMP_DEMSearch                                                         OmpDemSearchType;
          typedef DEMSearch<OmpDemSearchType >                                          DemSearchType;

          class_<IntegrationScheme, boost::noncopyable >
                    ("IntegrationScheme", init< >())
                  ;

          class_<ParticleCreatorDestructor, boost::noncopyable >
                    ("ParticleCreatorDestructor", init<>())
                  ;

          class_< ForwardEulerScheme, bases<IntegrationScheme>,  boost::noncopyable>
          (
                    "ForwardEulerScheme", init<>()
                  )
                  ;

          class_< MidPointScheme, bases<IntegrationScheme>,  boost::noncopyable>
          (
                    "MidPointScheme", init<>()
                  )
                  ;

          class_< ConstAverageAccelerationScheme, bases<IntegrationScheme>,  boost::noncopyable>
          (
                    "ConstAverageAccelerationScheme", init<>()
                  )
                  ;
                    
          class_<DemSearchType, bases<SpatialSearch>, boost::noncopyable>
                    ("OMP_DEMSearch", init<>())
                    ;
        
          typedef ExplicitSolverStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType > ExplicitSolverStrategyType;
          typedef ContinuumExplicitSolverStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType > ContinuumExplicitSolverStrategyType;
          
          class_< ExplicitSolverSettings, boost::noncopyable >
          (
            "ExplicitSolverSettings", init<>() )
          .def_readwrite("r_model_part",&ExplicitSolverSettings::r_model_part)
          .def_readwrite("contact_model_part",&ExplicitSolverSettings::contact_model_part)
          .def_readwrite("fem_model_part",&ExplicitSolverSettings::fem_model_part)
          .def_readwrite("cluster_model_part",&ExplicitSolverSettings::cluster_model_part)
          .def_readwrite("inlet_model_part",&ExplicitSolverSettings::inlet_model_part)
          ;
          
          class_< ExplicitSolverStrategyType, bases< BaseSolvingStrategyType >,  boost::noncopyable>
          (
          "ExplicitSolverStrategy", init< ExplicitSolverSettings&, double, double, double, bool, int, double, double, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, IntegrationScheme::Pointer, SpatialSearch::Pointer>())
                  .def("Initialize", &ExplicitSolverStrategyType::Initialize)
                  .def("InitialTimeStepCalculation", &ExplicitSolverStrategyType::InitialTimeStepCalculation)
                  .def("PrepareElementsForPrinting", &ContinuumExplicitSolverStrategyType::PrepareElementsForPrinting)
          ;
          
          class_< ContinuumExplicitSolverStrategyType, bases< ExplicitSolverStrategyType >,  boost::noncopyable>
          (
          "ContinuumExplicitSolverStrategy", init< ExplicitSolverSettings&, double, double, double, bool, int, double, double, ParticleCreatorDestructor::Pointer,DEM_FEM_Search::Pointer, IntegrationScheme::Pointer, SpatialSearch::Pointer>())
                  .def("Initialize", &ContinuumExplicitSolverStrategyType::Initialize)
                  .def("InitialTimeStepCalculation", &ContinuumExplicitSolverStrategyType::InitialTimeStepCalculation)
                  .def("PrepareContactElementsForPrinting", &ContinuumExplicitSolverStrategyType::PrepareContactElementsForPrinting)                  
          
          ;

        }

    }  // namespace Python.

} // Namespace Kratos
