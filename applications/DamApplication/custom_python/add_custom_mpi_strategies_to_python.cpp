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
#include "custom_python/add_custom_mpi_strategies_to_python.h"
#include "spaces/ublas_space.h"
#include "includes/kratos_parameters.h"

//Trilinos includes
#include "mpi.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "trilinos_space.h"

//linear solvers

//strategies

//builders and solvers

//schemes
#include "custom_strategies/schemes/trilinos_incrementalupdate_static_damped_scheme.hpp"
#include "custom_strategies/schemes/trilinos_dam_UP_scheme.hpp"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

void  AddCustomMPIStrategiesToPython()
{
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;

    typedef Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosBaseSchemeType;
    
    typedef TrilinosIncrementalUpdateStaticDampedScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosIncrementalUpdateStaticDampedSchemeType;
    typedef TrilinosDamUPScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosDamUPSchemeType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Schemes
    class_< TrilinosIncrementalUpdateStaticDampedSchemeType, bases<TrilinosBaseSchemeType>, boost::noncopyable >( "TrilinosIncrementalUpdateStaticDampedScheme", 
        init< double >() );
	class_< TrilinosDamUPSchemeType, bases< TrilinosBaseSchemeType >,  boost::noncopyable >("TrilinosDamUPScheme",
        init< double, double, double, double >());
}

}  // namespace Python.
} // Namespace Kratos
