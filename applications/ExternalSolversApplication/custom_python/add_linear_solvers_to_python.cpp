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

// External includes
#include <boost/python.hpp>
#include <complex>


// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "python/add_equation_systems_to_python.h"
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
#include "includes/linear_solver_factory.h"  
#include "custom_utilities/register_linear_solvers.h"

namespace Kratos
{

namespace Python
{
template <class TDataType>
using TLinearSolverType = LinearSolver<TUblasSparseSpace<TDataType>, TUblasDenseSpace<TDataType>>;
template <class TDataType>
using TDirectSolverType = DirectSolver<TUblasSparseSpace<TDataType>, TUblasDenseSpace<TDataType>>;

void  AddLinearSolversToPython()
{
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

    using namespace boost::python;


    //***************************************************************************
    //linear solvers
    //***************************************************************************
#ifdef INCLUDE_FEAST
    typedef FEASTSolver<SpaceType, LocalSpaceType> FEASTSolverType;
    class_<FEASTSolverType, FEASTSolverType::Pointer, bases<LinearSolverType>, boost::noncopyable >
        ( "FEASTSolver", init<Parameters>() )
        .def(init<Parameters, TLinearSolverType<std::complex<double>>::Pointer>())
        ;
#endif    
          
    
    class_<SuperLUSolverType, bases<DirectSolverType>, boost::noncopyable >
    ( "SuperLUSolver",
      init<>() )
      .def(init<Parameters>());
      
    class_<SuperLUIterativeSolverType, bases<LinearSolverType>, boost::noncopyable >
    ( "SuperLUIterativeSolver",init<>() )
    .def(init<double,int,int,double,double,double>())
    .def(init<Parameters>())
    ;

#ifdef INCLUDE_PASTIX
    typedef PastixSolver<SpaceType,  LocalSpaceType> PastixSolverType;
    class_<PastixSolverType, bases<LinearSolverType>, boost::noncopyable >
    ( "PastixSolver",init<int,bool>() )
    .def(init<double,int,int,int,bool>())
    .def(init<Parameters>());
    ;
    typedef PastixComplexSolver<ComplexSpaceType, ComplexLocalSpaceType> PastixComplexSolverType;
    class_<PastixComplexSolverType, bases<TDirectSolverType<std::complex<double>>>, boost::noncopyable >
    ("PastixComplexSolver",init<Parameters>())
    ;
#endif
    
    class_<GMRESSolverType, bases<IterativeSolverType>, boost::noncopyable >
    ( "GMRESSolver")
    .def(init<Parameters >())
    .def(init<double>())
    .def(init<double, unsigned int>())
    .def(init<double, unsigned int,  PreconditionerType::Pointer>())
    .def(self_ns::str(self))
    ;
    
    
    
//     RegisterLinearSolvers();
    
}

    

}  // namespace Python.




//Must put this definition here to avoid a problem with multiply defined symbols when including the external C libraries
RegisterLinearSolvers::RegisterLinearSolvers()
    {
        typedef TUblasSparseSpace<double> SpaceType;
        typedef TUblasDenseSpace<double> LocalSpaceType;
        //typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
        typedef SuperLUSolver<SpaceType,  LocalSpaceType> SuperLUSolverType;
        typedef SuperLUIterativeSolver<SpaceType,  LocalSpaceType> SuperLUIterativeSolverType;
        typedef GMRESSolver<SpaceType, LocalSpaceType> GMRESSolverType;

        //REGISTERING SOLVERS
        typedef LinearSolverFactoryBase<SpaceType,  LocalSpaceType> LinearSolverFactoryBaseType;
        
        static auto GMRESSolverFactory= LinearSolverFactory<SpaceType,LocalSpaceType,GMRESSolverType>();
        static auto SuperLUSolverFactory= LinearSolverFactory<SpaceType,LocalSpaceType,SuperLUSolverType>();
        static auto SuperLUIterativeSolverFactory= LinearSolverFactory<SpaceType,LocalSpaceType,SuperLUIterativeSolverType>();
            
        KratosComponents<LinearSolverFactoryBaseType>::Add(std::string("GMRESSolver"), GMRESSolverFactory);
        KratosComponents<LinearSolverFactoryBaseType>::Add(std::string("SuperLUSolver"), SuperLUSolverFactory);
        KratosComponents<LinearSolverFactoryBaseType>::Add(std::string("SuperLUIterativeSolver"), SuperLUIterativeSolverFactory);
        
    #ifdef INCLUDE_PASTIX
        typedef TUblasSparseSpace<std::complex<double>> ComplexSpaceType;
        typedef TUblasDenseSpace<std::complex<double>> ComplexLocalSpaceType;
        typedef LinearSolverFactoryBase<ComplexSpaceType, ComplexLocalSpaceType> ComplexLinearSolverFactoryBaseType;
        typedef PastixSolver<SpaceType,  LocalSpaceType> PastixSolverType;
        static auto PastixSolverFactory = LinearSolverFactory<SpaceType,LocalSpaceType,PastixSolverType>();    
        KratosComponents<LinearSolverFactoryBaseType>::Add(std::string("PastixSolver"), PastixSolverFactory);
        typedef PastixComplexSolver<ComplexSpaceType, ComplexLocalSpaceType> PastixComplexSolverType;
        static auto PastixComplexSolverFactory = LinearSolverFactory<ComplexSpaceType, ComplexLocalSpaceType, PastixComplexSolverType>();
        KratosComponents<ComplexLinearSolverFactoryBaseType>::Add(std::string("PastixComplexSolver"), PastixComplexSolverFactory);
    #endif
        
    #ifdef INCLUDE_FEAST
        typedef FEASTSolver<SpaceType, LocalSpaceType> FEASTSolverType;
        static auto FEASTSolverFactory= LinearSolverFactory<SpaceType,LocalSpaceType,FEASTSolverType>();
        KratosComponents<LinearSolverFactoryBaseType>::Add(std::string("FEASTSolver"), FEASTSolverFactory);
    #endif
    }
    
} // Namespace Kratos

