//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:              March 2017 $
//   Revision:            $Revision:                 1.0 $
//

// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 

// Project includes
#include "includes/define.h"
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
#include "custom_strategies/schemes/trilinos_newmark_quasistatic_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/trilinos_newmark_quasistatic_damped_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/trilinos_newmark_dynamic_U_Pw_scheme.hpp"


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
    
    typedef TrilinosNewmarkQuasistaticUPwScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosNewmarkQuasistaticUPwSchemeType;
    typedef TrilinosNewmarkQuasistaticDampedUPwScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosNewmarkQuasistaticDampedUPwSchemeType;
    typedef TrilinosNewmarkDynamicUPwScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosNewmarkDynamicUPwSchemeType;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Schemes
    class_< TrilinosNewmarkQuasistaticUPwSchemeType, bases<TrilinosBaseSchemeType>, boost::noncopyable >( "TrilinosNewmarkQuasistaticUPwScheme", 
        init< double, double, double >() );
    class_< TrilinosNewmarkQuasistaticDampedUPwSchemeType,bases< TrilinosBaseSchemeType >, boost::noncopyable >("TrilinosNewmarkQuasistaticDampedUPwScheme",
        init<  double, double, double, double, double >());
    class_< TrilinosNewmarkDynamicUPwSchemeType,bases< TrilinosBaseSchemeType >, boost::noncopyable >("TrilinosNewmarkDynamicUPwScheme",
        init<  double, double, double, double, double >());

}

}  // namespace Python.
} // Namespace Kratos
