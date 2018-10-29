//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes
#include <complex>

// External includes

// Project includes
#include "python/add_linear_solvers_to_python.h"
#include "includes/define_python.h"
#include "includes/kratos_parameters.h"
#include "spaces/ublas_space.h"

#include "linear_solvers/direct_solver.h"
#include "linear_solvers/iterative_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "external_includes/superlu_solver.h"
#include "external_includes/superlu_iterative_solver.h"
#include "external_includes/gmres_solver.h"

#ifdef INCLUDE_FEAST
  #include "external_includes/feast_solver.h"
#endif

#ifdef INCLUDE_PASTIX
  #include "external_includes/pastix_solver.h"
  #include "external_includes/pastix_complex_solver.h"
#endif

#include "externalsolvers_application.h"
#include "includes/standard_linear_solver_factory.h"

namespace Kratos {
namespace Python {

template <class TDataType>
using TSpaceType = UblasSpace<TDataType, compressed_matrix<TDataType>, vector<TDataType>>;
template <class TDataType>
using TLocalSpaceType = UblasSpace<TDataType, matrix<TDataType>, vector<TDataType>>;
template <class TDataType>
using TLinearSolverType = LinearSolver<TSpaceType<TDataType>, TLocalSpaceType<TDataType>>;
template <class TDataType>
using TDirectSolverType = DirectSolver<TSpaceType<TDataType>, TLocalSpaceType<TDataType>>;

void  AddLinearSolversToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef TUblasSparseSpace<double> SpaceType;
    typedef TUblasDenseSpace<double> LocalSpaceType;
    typedef TUblasSparseSpace<std::complex<double>> ComplexSpaceType;
    typedef TUblasDenseSpace<std::complex<double>> ComplexLocalSpaceType;
    typedef TLinearSolverType<double> LinearSolverType;
    typedef TDirectSolverType<double> DirectSolverType;
    typedef SuperLUSolver<SpaceType,  LocalSpaceType> SuperLUSolverType;
    typedef SuperLUIterativeSolver<SpaceType,  LocalSpaceType> SuperLUIterativeSolverType;
    typedef IterativeSolver<SpaceType, LocalSpaceType> IterativeSolverType;
    typedef GMRESSolver<SpaceType, LocalSpaceType> GMRESSolverType;
    typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;

    //***************************************************************************
    // Linear solvers
    //***************************************************************************
#ifdef INCLUDE_FEAST
    typedef FEASTSolver<SpaceType, LocalSpaceType> FEASTSolverType;                          //SOME PROBLEM WITH THE SKYLINE_CUSTOM ... TO BE FIXED
    py::class_<FEASTSolverType, FEASTSolverType::Pointer, LinearSolverType >
        (m, "FEASTSolver")
        .def(py::init<Parameters>() )
        .def(py::init<Parameters, TLinearSolverType<std::complex<double>>::Pointer>())
        ;
#endif

    py::class_<SuperLUSolverType, typename SuperLUSolverType::Pointer,DirectSolverType>
    (m, "SuperLUSolver")
        .def(py::init<>() )
        .def(py::init<Parameters>());

    py::class_<SuperLUIterativeSolverType, typename SuperLUIterativeSolverType::Pointer,LinearSolverType>
    (m, "SuperLUIterativeSolver")
        .def(py::init<>() )
        .def(py::init<double,int,int,double,double,double>())
        .def(py::init<Parameters>())
        ;

#ifdef INCLUDE_PASTIX
    typedef PastixSolver<SpaceType,  LocalSpaceType> PastixSolverType;
    py::class_<PastixSolverType, typename PastixSolverType::Pointer, LinearSolverType>
    (m, "PastixSolver")
        .def(py::init<int,bool>() )
        .def(py::init<double,int,int,int,bool>())
        .def(py::init<Parameters>());
        ;

    typedef PastixComplexSolver<TSpaceType<std::complex<double>>, TLocalSpaceType<std::complex<double>>> PastixComplexSolverType;
    py::class_<PastixComplexSolverType, typename PastixComplexSolverType::Pointer, TDirectSolverType<std::complex<double>>>
    (m,"PastixComplexSolver")
        .def(py::init<Parameters&>())
        ;
#endif

    py::class_<GMRESSolverType,typename GMRESSolverType::Pointer, IterativeSolverType>
    (m, "GMRESSolver")
        .def(py::init<Parameters >())
        .def(py::init<Parameters,  PreconditionerType::Pointer >())
        .def(py::init<double>())
        .def(py::init<double, unsigned int>())
        .def(py::init<double, unsigned int,  PreconditionerType::Pointer>())
        .def("__str__", PrintObject<GMRESSolverType>)
        ;

//     ExternalSolversApplicationRegisterLinearSolvers();

}

}  // namespace Python.


//Must put this definition here to avoid a problem with multiply defined symbols when including the external C libraries
ExternalSolversApplicationRegisterLinearSolvers::ExternalSolversApplicationRegisterLinearSolvers()
{
    typedef TUblasSparseSpace<double> SpaceType;
    typedef TUblasDenseSpace<double> LocalSpaceType;
    //typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
    typedef SuperLUSolver<SpaceType,  LocalSpaceType> SuperLUSolverType;
    typedef SuperLUIterativeSolver<SpaceType,  LocalSpaceType> SuperLUIterativeSolverType;
    typedef GMRESSolver<SpaceType, LocalSpaceType> GMRESSolverType;

    //REGISTERING SOLVERS
    static auto GMRESSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,GMRESSolverType>();
    static auto SuperLUSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,SuperLUSolverType>();
    static auto SuperLUIterativeSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,SuperLUIterativeSolverType>();

    KRATOS_REGISTER_LINEAR_SOLVER("GMRESSolver", GMRESSolverFactory);
    KRATOS_REGISTER_LINEAR_SOLVER("Super_LU", SuperLUSolverFactory); // NOTE: This is duplicated by retrocompatibility
    KRATOS_REGISTER_LINEAR_SOLVER("SuperLUSolver", SuperLUSolverFactory);
    KRATOS_REGISTER_LINEAR_SOLVER("SuperLUIterativeSolver", SuperLUIterativeSolverFactory);

#ifdef INCLUDE_PASTIX
    typedef TUblasSparseSpace<std::complex<double>> ComplexSpaceType;
    typedef TUblasDenseSpace<std::complex<double>> ComplexLocalSpaceType;
    typedef PastixSolver<SpaceType,  LocalSpaceType> PastixSolverType;
    static auto PastixSolverFactory = StandardLinearSolverFactory<SpaceType,LocalSpaceType,PastixSolverType>();
    KRATOS_REGISTER_LINEAR_SOLVER("PastixSolver", PastixSolverFactory);
    typedef PastixComplexSolver<ComplexSpaceType, ComplexLocalSpaceType> PastixComplexSolverType;
    static auto PastixComplexSolverFactory = StandardLinearSolverFactory<ComplexSpaceType, ComplexLocalSpaceType, PastixComplexSolverType>();
    KRATOS_REGISTER_COMPLEX_LINEAR_SOLVER("PastixComplexSolver", PastixComplexSolverFactory);
#endif

#ifdef INCLUDE_FEAST
    typedef FEASTSolver<SpaceType, LocalSpaceType> FEASTSolverType;
    static auto FEASTSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,FEASTSolverType>();
    KRATOS_REGISTER_LINEAR_SOLVER("FEASTSolver", FEASTSolverFactory);
#endif
}

} // Namespace Kratos

