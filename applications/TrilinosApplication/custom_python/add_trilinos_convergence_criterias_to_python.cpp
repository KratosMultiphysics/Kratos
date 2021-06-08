//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include "custom_python/add_trilinos_convergence_criterias_to_python.h"

//Trilinos includes
#include "Epetra_FEVector.h"

// Project includes
#include "trilinos_space.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#include "solving_strategies/convergencecriterias/or_criteria.h"

// Application includes
#include "custom_strategies/convergencecriterias/trilinos_displacement_criteria.h"
#include "custom_strategies/convergencecriterias/trilinos_residual_criteria.h"
#include "custom_strategies/convergencecriterias/trilinos_mixed_generic_criteria.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

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

    py::class_< TrilinosConvergenceCriteria, TrilinosConvergenceCriteriaPointer > (m,"TrilinosConvergenceCriteria")
    .def(py::init<>())
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

    py::class_< TrilinosDisplacementCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            typename TrilinosDisplacementCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Pointer,
            TrilinosConvergenceCriteria>(m,"TrilinosDisplacementCriteria")
            .def(py::init< double, double >());

    py::class_< TrilinosResidualCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            typename TrilinosResidualCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Pointer,
            TrilinosConvergenceCriteria >
            (m,"TrilinosResidualCriteria")
            .def(py::init< TrilinosSparseSpaceType::DataType, TrilinosSparseSpaceType::DataType >());

    py::class_<And_Criteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            typename And_Criteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Pointer,
            TrilinosConvergenceCriteria>
            (m,"TrilinosAndCriteria")
            .def(py::init<TrilinosConvergenceCriteriaPointer, TrilinosConvergenceCriteriaPointer > ());

    py::class_<Or_Criteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >,
            typename Or_Criteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::Pointer,
            TrilinosConvergenceCriteria>
            (m,"TrilinosOrCriteria")
            .def(py::init<TrilinosConvergenceCriteriaPointer, TrilinosConvergenceCriteriaPointer > ());

    typedef typename TrilinosMixedGenericCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType>::ConvergenceVariableListType ConvergenceVariableListType;
    py::class_<
        TrilinosMixedGenericCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename TrilinosMixedGenericCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        TrilinosConvergenceCriteria>(m, "TrilinosMixedGenericCriteria")
        .def(py::init< const ConvergenceVariableListType& >());
}


} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
