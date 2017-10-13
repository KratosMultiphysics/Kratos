#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>

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

using namespace boost::python;

typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

void AddProcesses()
{
    typedef SpalartAllmarasTurbulenceModel<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> BaseSpAlModelType;

    class_ < BaseSpAlModelType,bases< Process >, boost::noncopyable >
    ( "TrilinosBaseSpAlModel",no_init );

    // Turbulence models
    class_< TrilinosSpalartAllmarasTurbulenceModel< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >, bases<BaseSpAlModelType>, boost::noncopyable >
    ("TrilinosSpalartAllmarasTurbulenceModel", init < Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, unsigned int, double, unsigned int, bool, unsigned int>())
    .def("ActivateDES", &SpalartAllmarasTurbulenceModel< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::ActivateDES)
    .def("AdaptForFractionalStep", &SpalartAllmarasTurbulenceModel< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::AdaptForFractionalStep)
    ;

    typedef StokesInitializationProcess<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> BaseStokesInitializationType;

    class_ < BaseStokesInitializationType,bases< Process >, boost::noncopyable >
            ( "TrilinosBaseStokesInitialization",no_init )
            .def("SetConditions",&BaseStokesInitializationType::SetConditions);

    class_< TrilinosStokesInitializationProcess< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >, bases<BaseStokesInitializationType>, boost::noncopyable >
            ("TrilinosStokesInitializationProcess",init<Epetra_MpiComm&, ModelPart::Pointer,TrilinosLinearSolverType::Pointer, unsigned int, const Kratos::Variable<int>& >())
            ;

    typedef VariationalDistanceCalculationProcess<2,TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> BaseDistanceCalculationType2D;
    typedef VariationalDistanceCalculationProcess<3,TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> BaseDistanceCalculationType3D;

    class_< BaseDistanceCalculationType2D, bases< Process >, boost::noncopyable >("BaseDistanceCalculation2D",no_init);
    class_< BaseDistanceCalculationType3D, bases< Process >, boost::noncopyable >("BaseDistanceCalculation3D",no_init);

    class_< TrilinosVariationalDistanceCalculationProcess<2,TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>,
            bases< BaseDistanceCalculationType2D >, boost::noncopyable >("TrilinosVariationalDistanceCalculationProcess2D",init<Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, unsigned int>() );

    class_< TrilinosVariationalDistanceCalculationProcess<3,TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>,
            bases< BaseDistanceCalculationType3D >, boost::noncopyable >("TrilinosVariationalDistanceCalculationProcess3D",init<Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, unsigned int>() );
}

}
}

#endif
