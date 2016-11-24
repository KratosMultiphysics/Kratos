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
#include <boost/python.hpp>

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
#include "custom_strategies/convergencecriterias/trilinos_up_criteria.h"

//teuchos parameter list
#include "Teuchos_ParameterList.hpp"

namespace Kratos
{

namespace Python
{

using namespace boost::python;


void  AddConvergenceCriterias()
{
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    //typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;

    typedef Epetra_FECrsMatrix FECrsMatrix;


    //********************************************************************
    //********************************************************************
    //convergence criteria base class
    typedef ConvergenceCriteria< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosConvergenceCriteria;
    typedef ConvergenceCriteria< TrilinosSparseSpaceType, TrilinosLocalSpaceType > ::Pointer TTrilinosConvergenceCriteriaPointer;
    class_< TrilinosConvergenceCriteria, boost::noncopyable > ("TrilinosConvergenceCriteria", init<>())
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
            bases< TrilinosConvergenceCriteria >,
            boost::noncopyable >
            ("TrilinosDisplacementCriteria", init< double, double, Epetra_MpiComm& >());

    class_< TrilinosUPCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            bases< TrilinosConvergenceCriteria >,
            boost::noncopyable >
            ("TrilinosUPCriteria", init< double, double, double, double, Epetra_MpiComm& >());

    class_< ResidualCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            bases< TrilinosConvergenceCriteria >,
            boost::noncopyable >
            ("TrilinosResidualCriteria", init< TrilinosSparseSpaceType::DataType, TrilinosSparseSpaceType::DataType >());
            
    class_<And_Criteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            bases< TrilinosConvergenceCriteria >,
            boost::noncopyable >
            ("TrilinosAndCriteria", init<TTrilinosConvergenceCriteriaPointer, TTrilinosConvergenceCriteriaPointer > ());

    class_<Or_Criteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            bases< TrilinosConvergenceCriteria >,
            boost::noncopyable >
            ("TrilinosOrCriteria", init<TTrilinosConvergenceCriteriaPointer, TTrilinosConvergenceCriteriaPointer > ());
}


} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
