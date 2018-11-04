//   
//   Project Name:        
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: $
//

// External includes
#include "spaces/ublas_space.h"

// Project includes
#include "custom_python/add_custom_mpi_strategies_to_python.h"
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

using namespace pybind11;

void  AddCustomMPIStrategiesToPython(pybind11::module& m)
{
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;

    typedef Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosBaseSchemeType;
    
    typedef TrilinosIncrementalUpdateStaticDampedScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosIncrementalUpdateStaticDampedSchemeType;
    typedef TrilinosDamUPScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosDamUPSchemeType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Schemes
    class_< TrilinosIncrementalUpdateStaticDampedSchemeType, typename TrilinosIncrementalUpdateStaticDampedSchemeType::Pointer, TrilinosBaseSchemeType>
    (m,  "TrilinosIncrementalUpdateStaticDampedScheme")
    .def(init< double >());
	class_< TrilinosDamUPSchemeType, typename TrilinosDamUPSchemeType::Pointer, TrilinosBaseSchemeType >
    (m, "TrilinosDamUPScheme")
    .def(init< double, double, double, double >());
}

}  // namespace Python.
} // Namespace Kratos
