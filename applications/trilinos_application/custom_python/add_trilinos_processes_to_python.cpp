#if defined(KRATOS_PYTHON)
// External includes
#include "custom_python/add_trilinos_processes_to_python.h"

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
#include "includes/model_part.h"
#include "processes/process.h"
#include "trilinos_application.h"
#include "trilinos_space.h"
#include "spaces/ublas_space.h"

#include "custom_processes/trilinos_levelset_convection_process.h"
#include "custom_processes/trilinos_spalart_allmaras_turbulence_model.h"
#include "custom_processes/trilinos_stokes_initialization_process.h"
#include "custom_processes/trilinos_variational_distance_calculation_process.h"
#include "../FluidDynamicsApplication/custom_processes/spalart_allmaras_turbulence_model.h"
#include "../FluidDynamicsApplication/custom_processes/stokes_initialization_process.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

void AddProcesses(pybind11::module& m)
{
    // Turbulence models
    typedef SpalartAllmarasTurbulenceModel<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> BaseSpAlModelType;

    py::class_<BaseSpAlModelType, BaseSpAlModelType::Pointer, Process>(m, "TrilinosBaseSpAlModel");

    typedef TrilinosSpalartAllmarasTurbulenceModel<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosSpAlModelType;
    py::class_<TrilinosSpAlModelType, TrilinosSpAlModelType::Pointer, BaseSpAlModelType >(m,"TrilinosSpalartAllmarasTurbulenceModel")
        .def(py::init < Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, unsigned int, double, unsigned int, bool, unsigned int>())
        .def("ActivateDES", &SpalartAllmarasTurbulenceModel< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::ActivateDES)
        .def("AdaptForFractionalStep", &SpalartAllmarasTurbulenceModel< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::AdaptForFractionalStep)
        ;

    // Stokes initialization processes
    typedef StokesInitializationProcess<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> BaseStokesInitializationType;

    py::class_<BaseStokesInitializationType, BaseStokesInitializationType::Pointer, Process>(m, "TrilinosBaseStokesInitialization")
        .def("SetConditions",&BaseStokesInitializationType::SetConditions)
        ;

    typedef TrilinosStokesInitializationProcess<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosStokesInitializationType;
    py::class_<TrilinosStokesInitializationType, TrilinosStokesInitializationType::Pointer, BaseStokesInitializationType >
        (m,"TrilinosStokesInitializationProcess")
        .def(py::init<Epetra_MpiComm&, ModelPart&,TrilinosLinearSolverType::Pointer, unsigned int, const Kratos::Variable<int>& >())
        ;

    // Variational distance calculation processes
    typedef VariationalDistanceCalculationProcess<2,TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> BaseDistanceCalculationType2D;
    typedef VariationalDistanceCalculationProcess<3,TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> BaseDistanceCalculationType3D;

    py::class_<BaseDistanceCalculationType2D, BaseDistanceCalculationType2D::Pointer, Process>(m,"BaseDistanceCalculation2D");
    py::class_<BaseDistanceCalculationType3D, BaseDistanceCalculationType3D::Pointer, Process>(m,"BaseDistanceCalculation3D");

    typedef TrilinosVariationalDistanceCalculationProcess<2,TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosDistanceCalculationType2D;
    py::class_<TrilinosDistanceCalculationType2D, TrilinosDistanceCalculationType2D::Pointer, BaseDistanceCalculationType2D >(m,"TrilinosVariationalDistanceCalculationProcess2D")
        .def(py::init<Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, unsigned int>())
        ;

    typedef TrilinosVariationalDistanceCalculationProcess<3,TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosDistanceCalculationType3D;
    py::class_<TrilinosDistanceCalculationType3D, TrilinosDistanceCalculationType3D::Pointer, BaseDistanceCalculationType3D >(m,"TrilinosVariationalDistanceCalculationProcess3D")
        .def(py::init<Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, unsigned int>())
        ;

    // Level set convection processes
    typedef LevelSetConvectionProcess<2, TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> BaseLevelSetConvectionProcess2D;
    typedef LevelSetConvectionProcess<3, TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> BaseLevelSetConvectionProcess3D;

    py::class_<BaseLevelSetConvectionProcess2D, BaseLevelSetConvectionProcess2D::Pointer, Process>(m,"BaseTrilinosLevelSetConvectionProcess2D");
    py::class_<BaseLevelSetConvectionProcess3D, BaseLevelSetConvectionProcess3D::Pointer, Process>(m,"BaseTrilinosLevelSetConvectionProcess3D");

    typedef TrilinosLevelSetConvectionProcess<2, TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosLevelSetConvectionProcess2D;
    py::class_<TrilinosLevelSetConvectionProcess2D, TrilinosLevelSetConvectionProcess2D::Pointer, BaseLevelSetConvectionProcess2D>(m, "TrilinosLevelSetConvectionProcess2D")
        .def(py::init<Epetra_MpiComm&, Variable<double>&, ModelPart&, TrilinosLinearSolverType::Pointer>())
        .def(py::init<Epetra_MpiComm&, Variable<double>&, ModelPart&, TrilinosLinearSolverType::Pointer, const double>())
        .def(py::init<Epetra_MpiComm&, Variable<double>&, ModelPart&, TrilinosLinearSolverType::Pointer, const double, const double>())
        .def(py::init<Epetra_MpiComm&, Variable<double>&, ModelPart&, TrilinosLinearSolverType::Pointer, const double, const double, const unsigned int>())
        ;

    typedef TrilinosLevelSetConvectionProcess<3, TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosLevelSetConvectionProcess3D;
    py::class_<TrilinosLevelSetConvectionProcess3D, TrilinosLevelSetConvectionProcess3D::Pointer, BaseLevelSetConvectionProcess3D>(m, "TrilinosLevelSetConvectionProcess3D")
        .def(py::init<Epetra_MpiComm&, Variable<double>&, ModelPart&, TrilinosLinearSolverType::Pointer>())
        .def(py::init<Epetra_MpiComm&, Variable<double>&, ModelPart&, TrilinosLinearSolverType::Pointer, const double>())
        .def(py::init<Epetra_MpiComm&, Variable<double>&, ModelPart&, TrilinosLinearSolverType::Pointer, const double, const double>())
        .def(py::init<Epetra_MpiComm&, Variable<double>&, ModelPart&, TrilinosLinearSolverType::Pointer, const double, const double, const unsigned int>())
        ;

}

}
}

#endif
