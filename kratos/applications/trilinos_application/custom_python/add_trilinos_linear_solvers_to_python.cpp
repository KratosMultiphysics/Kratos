//
//   Project Name:        Kratos
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-12-09 20:20:55 $
//   Revision:            $Revision: 1.5 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>

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
#include "includes/define.h"
#include "trilinos_application.h"
#include "trilinos_space.h"
#include "spaces/ublas_space.h"
#include "includes/model_part.h"

//strategies
// #include "solving_strategies/strategies/solving_strategy.h"
// #include "solving_strategies/strategies/residualbased_linear_strategy.h"
// #include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

//schemes
// #include "solving_strategies/schemes/scheme.h"
// #include "custom_strategies/schemes/trilinos_residualbased_incrementalupdate_static_scheme.h"
// #include "custom_strategies/schemes/trilinos_residualbased_lagrangian_monolithic_scheme.h"
// #include "../../incompressible_fluid_application/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme.h"
// #include "custom_strategies/schemes/trilinos_predictorcorrector_velocity_bossak_scheme.h"

//convergence criterias
// #include "solving_strategies/convergencecriterias/convergence_criteria.h"
// #include "solving_strategies/convergencecriterias/displacement_criteria.h"
//
// //Builder And Solver
// // #include "solving_strategies/builder_and_solvers/builder_and_solver.h"
// #include "custom_strategies/builder_and_solvers/trilinos_residualbased_elimination_builder_and_solver.h"
// #include "custom_strategies/convergencecriterias/trilinos_displacement_criteria.h"
// #include "custom_strategies/convergencecriterias/trilinos_up_criteria.h"
// #include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML.h"
// #include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML_vec.h"
// #include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML_mixed.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

//utilities
#include "python/pointer_vector_set_python_interface.h"

//teuchos parameter list
#include "Teuchos_ParameterList.hpp"

#include "external_includes/epetra_default_utility.h"
#include "external_includes/aztec_solver.h"
#include "external_includes/amesos_solver.h"
#include "external_includes/ml_solver.h"

#ifdef AMGCL_INCLUDED
    #include "external_includes/amgcl_solver.h"
#endif

//configuration files
// #include "../../incompressible_fluid_application/custom_strategies/strategies/solver_configuration.h"
// #include "custom_strategies/strategies/trilinos_fractionalstep_configuration.h"
// #include "../../incompressible_fluid_application/custom_strategies/strategies/fractional_step_strategy.h"
// #include "../../incompressible_fluid_application/incompressible_fluid_application.h"



namespace Kratos
{

namespace Python
{

using namespace boost::python;


void  AddLinearSolvers()
{
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
    typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

    class_<TrilinosLinearSolverType, TrilinosLinearSolverType::Pointer > ("TrilinosLinearSolver");

    class_<EpetraDefaultSetter, boost::noncopyable>("EpetraDefaultSetter",init<>())
    .def("SetDefaults", &EpetraDefaultSetter::SetDefaults);

    typedef AztecSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AztecSolverType;
    class_<AztecSolverType, bases<TrilinosLinearSolverType>, boost::noncopyable >
    ("AztecSolver",
     init< Teuchos::ParameterList&, std::string, Teuchos::ParameterList&, double, int, int >())
    .def("SetScalingType", &AztecSolverType::SetScalingType)
    ;

    typedef AmesosSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AmesosSolverType;
    class_<AmesosSolverType, bases<TrilinosLinearSolverType>, boost::noncopyable >
    ("AmesosSolver",
     init<const std::string&, Teuchos::ParameterList& >());

    typedef MultiLevelSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > MLSolverType;
    class_<MLSolverType, bases<TrilinosLinearSolverType>, boost::noncopyable >
    ("MultiLevelSolver",
     init<Teuchos::ParameterList&, Teuchos::ParameterList&, double, int >())
        .def("SetScalingType", &MLSolverType::SetScalingType)
        .def("SetReformPrecAtEachStep", &MLSolverType::SetReformPrecAtEachStep)
        ;

    #ifdef AMGCL_INCLUDED
         enum_<TrilinosAmgclSettings::AMGCLSmoother>("AMGCLSmoother")
            .value("SPAI0", TrilinosAmgclSettings::SPAI0)
            .value("ILU0", TrilinosAmgclSettings::ILU0)
            .value("DAMPED_JACOBI",TrilinosAmgclSettings::DAMPED_JACOBI)
            .value("GAUSS_SEIDEL",TrilinosAmgclSettings::GAUSS_SEIDEL)
            .value("MULTICOLOR_GAUSS_SEIDEL",TrilinosAmgclSettings::MULTICOLOR_GAUSS_SEIDEL)         
            ;
        
        enum_<TrilinosAmgclSettings::AMGCLIterativeSolverType>("AMGCLIterativeSolverType")
            .value("GMRES", TrilinosAmgclSettings::GMRES)
            .value("BICGSTAB", TrilinosAmgclSettings::BICGSTAB)
            .value("CG",TrilinosAmgclSettings::CG)
            .value("BICGSTAB_WITH_GMRES_FALLBACK",TrilinosAmgclSettings::BICGSTAB_WITH_GMRES_FALLBACK)
            .value("BICGSTAB2",TrilinosAmgclSettings::BICGSTAB2)
            ;
            
        enum_<TrilinosAmgclSettings::AMGCLCoarseningType>("AMGCLCoarseningType")
            .value("RUGE_STUBEN", TrilinosAmgclSettings::RUGE_STUBEN)
            .value("AGGREGATION", TrilinosAmgclSettings::AGGREGATION)
            .value("SA",TrilinosAmgclSettings::SA)
            .value("SA_EMIN",TrilinosAmgclSettings::SA_EMIN)
            ;
  
        typedef AmgclMPISolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > AmgclMPISolverType;
        class_<AmgclMPISolverType, bases<TrilinosLinearSolverType>, boost::noncopyable >
            ("AmgclMPISolver",init<double, int,int,bool >())
            .def(init<TrilinosAmgclSettings::AMGCLSmoother,TrilinosAmgclSettings::AMGCLIterativeSolverType,TrilinosAmgclSettings::AMGCLCoarseningType,double,int,int,bool>())
            .def("SetDoubleParameter", &AmgclMPISolverType::SetDoubleParameter)
            .def("SetIntParameter", &AmgclMPISolverType::SetIntParameter)
            ;
            

    #endif
    
    enum_<AztecScalingType>("AztecScalingType")
    .value("NoScaling", NoScaling)
    .value("LeftScaling", LeftScaling)
    .value("SymmetricScaling", SymmetricScaling)
    ;

    enum_<MLSolverType::ScalingType>("MLSolverScalingType")
    .value("NoScaling", MLSolverType::NoScaling)
    .value("LeftScaling", MLSolverType::LeftScaling)
    ;
}


} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
