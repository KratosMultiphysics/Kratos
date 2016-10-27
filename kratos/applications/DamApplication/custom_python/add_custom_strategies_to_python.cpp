//   
//   Project Name:        
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: $
//

// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 

// Project includes
#include "includes/define.h"
#include "containers/flags.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "spaces/ublas_space.h"
#include "includes/kratos_parameters.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"

//builders and solvers

//schemes
#include "custom_strategies/schemes/incrementalupdate_static_smoothing_scheme.hpp"
#include "custom_strategies/schemes/bossak_displacement_smoothing_scheme.hpp"

//strategies
#include "custom_strategies/strategies/dam_eulerian_convection_diffusion_strategy.hpp"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

void  AddCustomStrategiesToPython()
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;


    //custom scheme types
    typedef IncrementalUpdateStaticSmoothingScheme< SparseSpaceType, LocalSpaceType >  IncrementalUpdateStaticSmoothingSchemeType;     
    typedef BossakDisplacementSmoothingScheme< SparseSpaceType, LocalSpaceType >  BossakDisplacementSmoothingSchemeType;
    
    //custom strategies types
    typedef DamEulerianConvectionDiffusionStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > DamEulerianConvectionDiffusionStrategyType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    class_< IncrementalUpdateStaticSmoothingSchemeType, bases< BaseSchemeType >, boost::noncopyable >("IncrementalUpdateStaticSmoothingScheme",
        init< >());

    class_< BossakDisplacementSmoothingSchemeType, bases< BaseSchemeType >, boost::noncopyable >("BossakDisplacementSmoothingScheme",
        init< double >());
            
    class_< DamEulerianConvectionDiffusionStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >("DamEulerianConvectionDiffusionStrategy", 
        init < ModelPart&, LinearSolverType::Pointer, Parameters&, bool, int >());

}

}  // namespace Python.
} // Namespace Kratos

