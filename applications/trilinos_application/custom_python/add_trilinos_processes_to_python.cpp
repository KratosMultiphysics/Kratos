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
#include "../FluidDynamicsApplication/custom_processes/spalart_allmaras_turbulence_model.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos {

    namespace Python {

        using namespace boost::python;

        typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
        typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

        void AddProcesses()
        {
            typedef SpalartAllmarasTurbulenceModel<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> BaseSpAlModelType;

            class_ < BaseSpAlModelType,bases< Process >, boost::noncopyable >
            ( "BaseSpAlModelType",no_init );

            // Turbulence models
            class_< TrilinosSpalartAllmarasTurbulenceModel< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >, bases<BaseSpAlModelType>, boost::noncopyable >
                    ("TrilinosSpalartAllmarasTurbulenceModel", init < Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, unsigned int, double, unsigned int, bool, unsigned int>())
                    .def("ActivateDES", &SpalartAllmarasTurbulenceModel< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::ActivateDES)
                    .def("AdaptForFractionalStep", &SpalartAllmarasTurbulenceModel< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::AdaptForFractionalStep)
                    ;
        }

    }
}

#endif