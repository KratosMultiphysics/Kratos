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

//Trilinos includes
//~ #include "mpi.h"
//~ #include "Epetra_FECrsMatrix.h"
//~ #include "Epetra_FEVector.h"
//~ #include "trilinos_space.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"

//builders and solvers

//schemes
#include "custom_strategies/schemes/incrementalupdate_static_smoothing_scheme.hpp"
#include "custom_strategies/schemes/incrementalupdate_static_damped_smoothing_scheme.hpp"
#include "custom_strategies/schemes/bossak_displacement_smoothing_scheme.hpp"
#include "custom_strategies/schemes/dam_UP_scheme.hpp"
#include "custom_strategies/schemes/dam_P_scheme.hpp"
//~ #include "custom_strategies/schemes/trilinos_incrementalupdate_static_damped_scheme.hpp"
//~ #include "custom_strategies/schemes/trilinos_dam_UP_scheme.hpp"

//strategies


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

    //custom scheme types
    typedef IncrementalUpdateStaticSmoothingScheme< SparseSpaceType, LocalSpaceType >  IncrementalUpdateStaticSmoothingSchemeType; 
    typedef IncrementalUpdateStaticDampedSmoothingScheme< SparseSpaceType, LocalSpaceType >  IncrementalUpdateStaticDampedSmoothingSchemeType;
    typedef BossakDisplacementSmoothingScheme< SparseSpaceType, LocalSpaceType >  BossakDisplacementSmoothingSchemeType;
    typedef DamUPScheme< SparseSpaceType, LocalSpaceType >  DamUPSchemeType;
    typedef DamPScheme< SparseSpaceType, LocalSpaceType >  DamPSchemeType;
    
    //~ typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    //~ typedef Scheme< TrilinosSparseSpaceType, LocalSpaceType > TrilinosBaseSchemeType;
    //~ typedef TrilinosIncrementalUpdateStaticDampedScheme<TrilinosSparseSpaceType, LocalSpaceType> TrilinosIncrementalUpdateStaticDampedSchemeType;
    //~ typedef TrilinosDamUPScheme<TrilinosSparseSpaceType, LocalSpaceType> TrilinosDamUPSchemeType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Schemes
    class_< IncrementalUpdateStaticSmoothingSchemeType, bases< BaseSchemeType >, boost::noncopyable >("IncrementalUpdateStaticSmoothingScheme",
        init< >());

    class_< IncrementalUpdateStaticDampedSmoothingSchemeType, bases< BaseSchemeType >, boost::noncopyable >("IncrementalUpdateStaticDampedSmoothingScheme",
        init< double, double >());

    class_< BossakDisplacementSmoothingSchemeType, bases< BaseSchemeType >, boost::noncopyable >("BossakDisplacementSmoothingScheme",
        init< double, double, double >());

	class_< DamUPSchemeType, bases< BaseSchemeType >,  boost::noncopyable >("DamUPScheme",
        init< double, double, double, double >());
    
    class_< DamPSchemeType, bases< BaseSchemeType >,  boost::noncopyable >("DamPScheme",
        init< double, double >());

    //~ class_< TrilinosIncrementalUpdateStaticDampedSchemeType, bases<TrilinosBaseSchemeType>, boost::noncopyable >( "TrilinosIncrementalUpdateStaticDampedScheme", 
        //~ init< double >() );

	//~ class_< TrilinosDamUPSchemeType, bases< TrilinosBaseSchemeType >,  boost::noncopyable >("TrilinosDamUPScheme",
        //~ init< double, double, double, double >());
}

}  // namespace Python.
} // Namespace Kratos

