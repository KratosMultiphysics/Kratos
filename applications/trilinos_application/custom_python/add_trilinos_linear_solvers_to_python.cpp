//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

#include "custom_python/add_trilinos_strategies_to_python.h"

//Trilinos includes
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"


// Project includes
#include "trilinos_application.h"
#include "trilinos_space.h"
#include "spaces/ublas_space.h"
#include "includes/model_part.h"


//linear solvers
#include "linear_solvers/linear_solver.h"

//teuchos parameter list
#include "Teuchos_ParameterList.hpp"

#include "external_includes/epetra_default_utility.h"
#include "external_includes/aztec_solver.h"
#include "external_includes/amesos_solver.h"
#include "external_includes/ml_solver.h"

#include "external_includes/amgcl_mpi_solver.h"
#include "external_includes/amgcl_mpi_schur_complement_solver.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;

typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

void Solve(TrilinosLinearSolverType& solver, 
            TrilinosSparseSpaceType::MatrixType& rA, 
            TrilinosSparseSpaceType::VectorType& rX, 
            TrilinosSparseSpaceType::VectorType& rB
            
            )
{
    solver.Solve(rA,rX,rB);
}

void  AddLinearSolvers(pybind11::module& m)
{
    
    class_<TrilinosLinearSolverType, TrilinosLinearSolverType::Pointer > (m,"TrilinosLinearSolver")
    .def(init<>())
    .def("Solve", Solve);

    class_<EpetraDefaultSetter>(m,"EpetraDefaultSetter").def( init<>())
    .def("SetDefaults", &EpetraDefaultSetter::SetDefaults);

    typedef AztecSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AztecSolverType;
    class_<AztecSolverType, typename AztecSolverType::Pointer, TrilinosLinearSolverType >
    (m,"AztecSolver")
    .def( init< Teuchos::ParameterList&, std::string, Teuchos::ParameterList&, double, int, int >())
    .def(init<Parameters>())
    .def("SetScalingType", &AztecSolverType::SetScalingType)
    ;

    typedef AmesosSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AmesosSolverType;
    class_<AmesosSolverType, typename AmesosSolverType::Pointer, TrilinosLinearSolverType >
    (m,"AmesosSolver").def( init<const std::string&, Teuchos::ParameterList& >())
    .def(init<Parameters>())
    .def(init<Parameters>())
    .def_static("HasSolver", &AmesosSolverType::HasSolver)
    ;

    typedef MultiLevelSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > MLSolverType;
    class_<MLSolverType, typename MLSolverType::Pointer, TrilinosLinearSolverType >
    (m,"MultiLevelSolver").def( init<Teuchos::ParameterList&, Teuchos::ParameterList&, double, int >())
    .def(init<Parameters>())
    .def("SetScalingType", &MLSolverType::SetScalingType)
    .def("SetReformPrecAtEachStep", &MLSolverType::SetReformPrecAtEachStep)
    ;

    typedef AmgclMPISolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AmgclMPISolverType;
    class_<AmgclMPISolverType, typename AmgclMPISolverType::Pointer, TrilinosLinearSolverType >
    (m,"AmgclMPISolver").def( init<Parameters>()) //init<double, int,int,bool >())
    .def("SetDoubleParameter", &AmgclMPISolverType::SetDoubleParameter)
    .def("SetIntParameter", &AmgclMPISolverType::SetIntParameter)
    ;
    
    typedef AmgclMPISchurComplementSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AmgclMPISchurComplementSolverType;
    class_<AmgclMPISchurComplementSolverType, typename AmgclMPISchurComplementSolverType::Pointer, TrilinosLinearSolverType >
    (m,"AmgclMPISchurComplementSolver").def( init<Parameters>()) 
    ;
    
    enum_<AztecScalingType>(m,"AztecScalingType")
    .value("NoScaling", NoScaling)
    .value("LeftScaling", LeftScaling)
    .value("SymmetricScaling", SymmetricScaling)
    ;

    enum_<MLSolverType::ScalingType>(m,"MLSolverScalingType")
    .value("NoScaling", MLSolverType::NoScaling)
    .value("LeftScaling", MLSolverType::LeftScaling)
    ;
    
}


} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
