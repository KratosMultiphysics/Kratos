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

// External includes
#include <complex>

namespace Kratos
{

namespace Python
{
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
    typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
    typedef DirectSolver<SpaceType,  LocalSpaceType> DirectSolverType;
    typedef SuperLUSolver<SpaceType,  LocalSpaceType> SuperLUSolverType;
    typedef SuperLUIterativeSolver<SpaceType,  LocalSpaceType> SuperLUIterativeSolverType;
    typedef IterativeSolver<SpaceType, LocalSpaceType> IterativeSolverType;
    typedef GMRESSolver<SpaceType, LocalSpaceType> GMRESSolverType;
    typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;

    using namespace pybind11;


    //***************************************************************************
    //linear solvers
    //***************************************************************************
#ifdef INCLUDE_FEAST
    typedef FEASTSolver<SpaceType, LocalSpaceType> FEASTSolverType;                          //SOME PROBLEM WITH THE SKYLINE_CUSTOM ... TO BE FIXED
    class_<FEASTSolverType, FEASTSolverType::Pointer, LinearSolverType >
        (m, "FEASTSolver")
        .def(init<Parameters::Pointer>() )
        .def(init<Parameters::Pointer, TLinearSolverType<std::complex<double>>::Pointer>())
        ;
#endif


    class_<SuperLUSolverType, typename SuperLUSolverType::Pointer,DirectSolverType>
    (m, "SuperLUSolver")
      .def(init<>() )
      .def(init<Parameters>());

    class_<SuperLUIterativeSolverType, typename SuperLUIterativeSolverType::Pointer,LinearSolverType>
    (m, "SuperLUIterativeSolver")
    .def(init<>() )
    .def(init<double,int,int,double,double,double>())
    .def(init<Parameters>())
    ;

#ifdef INCLUDE_PASTIX
    typedef PastixSolver<SpaceType,  LocalSpaceType> PastixSolverType;
    class_<PastixSolverType, typename PastixSolverType::Pointer, LinearSolverType>
    (m, "PastixSolver")
    .def(init<int,bool>() )
    .def(init<double,int,int,int,bool>())
    .def(init<Parameters>());
    ;
    typedef PastixComplexSolver<TSpaceType<std::complex<double>>, TLocalSpaceType<std::complex<double>>> PastixComplexSolverType;
    class_<PastixComplexSolverType, typename PastixComplexSolverType::Pointer, TDirectSolverType<std::complex<double>>>
    (m,"PastixComplexSolver")
    .def(init<Parameters&>())
    ;
#endif

    class_<GMRESSolverType,typename GMRESSolverType::Pointer, IterativeSolverType>
    (m, "GMRESSolver")
    .def(init<Parameters >())
    .def(init<Parameters,  PreconditionerType::Pointer >())
    .def(init<double>())
    .def(init<double, unsigned int>())
    .def(init<double, unsigned int,  PreconditionerType::Pointer>())
    .def("__str__", KRATOS_DEF_PYTHON_STR(GMRESSolverType))
    ;
}

}  // namespace Python.

} // Namespace Kratos

