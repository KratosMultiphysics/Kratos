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

#include "custom_python/add_trilinos_convergence_criterias_to_python.h"

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
#include "includes/define.h"
#include "trilinos_application.h"
#include "trilinos_space.h"
#include "spaces/ublas_space.h"
#include "add_trilinos_linear_solvers_to_python.h"
#include "includes/model_part.h"

//convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#include "solving_strategies/convergencecriterias/or_criteria.h"
//#include "solving_strategies/convergencecriterias/displacement_criteria.h"
//
#include "custom_strategies/convergencecriterias/trilinos_displacement_criteria.h"
#include "custom_strategies/convergencecriterias/trilinos_residual_criteria.h"
#include "custom_strategies/convergencecriterias/trilinos_up_criteria.h"

//teuchos parameter list
#include "Teuchos_ParameterList.hpp"

namespace Kratos
{

namespace Python
{

using namespace pybind11;


void  AddConvergenceCriterias(pybind11::module& m)
{
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    //typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;

    //typedef Epetra_FECrsMatrix FECrsMatrix;


    //********************************************************************
    //********************************************************************
    //convergence criteria base class
    typedef ConvergenceCriteria< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosConvergenceCriteria;

    typedef typename ConvergenceCriteria< TrilinosSparseSpaceType, TrilinosLocalSpaceType > ::Pointer TrilinosConvergenceCriteriaPointer;

    class_< TrilinosConvergenceCriteria, TrilinosConvergenceCriteriaPointer > (m,"TrilinosConvergenceCriteria")
    .def(init<>())
    .def("SetActualizeRHSFlag", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::SetActualizeRHSFlag)
    .def("GetActualizeRHSflag", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::GetActualizeRHSflag)
    .def("PreCriteria", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::PreCriteria)
    .def("PostCriteria", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::PostCriteria)
    .def("Initialize", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Initialize)
    .def("InitializeSolutionStep", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::InitializeSolutionStep)
    .def("FinalizeSolutionStep", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::FinalizeSolutionStep)
    .def("Check", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Check)
    .def("SetEchoLevel", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::SetEchoLevel)
    ;

    class_< TrilinosDisplacementCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            typename TrilinosDisplacementCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Pointer,
            TrilinosConvergenceCriteria>(m,"TrilinosDisplacementCriteria")
            .def(init< double, double >());

    class_< TrilinosUPCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            typename TrilinosUPCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Pointer,
            TrilinosConvergenceCriteria >
            (m,"TrilinosUPCriteria")
            .def(init< double, double, double, double >());

    class_< TrilinosResidualCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            typename TrilinosResidualCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Pointer,
            TrilinosConvergenceCriteria >
            (m,"TrilinosResidualCriteria")
            .def(init< TrilinosSparseSpaceType::DataType, TrilinosSparseSpaceType::DataType >());

    class_<And_Criteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            typename And_Criteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Pointer,
            TrilinosConvergenceCriteria>
            (m,"TrilinosAndCriteria")
            .def(init<TrilinosConvergenceCriteriaPointer, TrilinosConvergenceCriteriaPointer > ());

    class_<Or_Criteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            typename Or_Criteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Pointer,
            TrilinosConvergenceCriteria>
            (m,"TrilinosOrCriteria")
            .def(init<TrilinosConvergenceCriteriaPointer, TrilinosConvergenceCriteriaPointer > ());
}


} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
