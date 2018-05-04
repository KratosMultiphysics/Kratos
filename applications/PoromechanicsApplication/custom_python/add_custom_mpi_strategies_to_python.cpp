//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:              March 2017 $
//   Revision:            $Revision:                 1.0 $
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
#include "custom_strategies/schemes/trilinos_newmark_quasistatic_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/trilinos_newmark_quasistatic_damped_U_Pw_scheme.hpp"
#include "custom_strategies/schemes/trilinos_newmark_dynamic_U_Pw_scheme.hpp"


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
    
    typedef TrilinosNewmarkQuasistaticUPwScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosNewmarkQuasistaticUPwSchemeType;
    typedef TrilinosNewmarkQuasistaticDampedUPwScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosNewmarkQuasistaticDampedUPwSchemeType;
    typedef TrilinosNewmarkDynamicUPwScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosNewmarkDynamicUPwSchemeType;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Schemes
    class_< TrilinosNewmarkQuasistaticUPwSchemeType, typename TrilinosNewmarkQuasistaticUPwSchemeType::Pointer, TrilinosBaseSchemeType >
    (m, "TrilinosNewmarkQuasistaticUPwScheme")
    .def(init< double, double, double >() );
    class_< TrilinosNewmarkQuasistaticDampedUPwSchemeType, typename TrilinosNewmarkQuasistaticDampedUPwSchemeType::Pointer, TrilinosBaseSchemeType >
    (m, "TrilinosNewmarkQuasistaticDampedUPwScheme")
    .def(init<  double, double, double, double, double >());
    class_< TrilinosNewmarkDynamicUPwSchemeType, typename TrilinosNewmarkDynamicUPwSchemeType::Pointer, TrilinosBaseSchemeType >
    (m, "TrilinosNewmarkDynamicUPwScheme")
    .def(init<  double, double, double, double, double >());

}

}  // namespace Python.
} // Namespace Kratos
