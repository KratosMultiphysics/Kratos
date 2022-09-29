//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#include "add_trilinos_utilities_to_python.h"

// Trilinos includes
#include "Epetra_FEVector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_MpiComm.h"

// KratosCore dependencies
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"
#include "processes/process.h"
#include "spaces/ublas_space.h"

// TrilinosApplication dependencies
#include "trilinos_space.h"

// FluidDynamics trilinos extensions
#include "custom_processes/trilinos_spalart_allmaras_turbulence_model.h"
#include "custom_processes/trilinos_stokes_initialization_process.h"

// FluidDynamicsApplication dependencies
#include "custom_processes/distance_smoothing_process.h"

namespace Kratos {
namespace Python {

namespace py = pybind11;

typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> UblasLocalSpaceType;
typedef LinearSolver<TrilinosSparseSpaceType, UblasLocalSpaceType> TrilinosLinearSolverType;

// Helper to define Trilinos DistanceSmoothingProcess
template<unsigned int TDim> using TrilinosDistanceSmoothingProcess = DistanceSmoothingProcess<TDim, TrilinosSparseSpaceType, UblasLocalSpaceType, TrilinosLinearSolverType>;

template< class TBinder, unsigned int TDim > void DistanceSmoothingConstructionHelper(TBinder& rBinder)
{
    rBinder.def(py::init([](
        Epetra_MpiComm& rComm, ModelPart& rModelPart, TrilinosLinearSolverType::Pointer pLinearSolver)
        {
            constexpr int RowSizeGuess = (TDim == 2 ? 15 : 40);
            auto p_builder_solver = Kratos::make_shared<TrilinosBlockBuilderAndSolver<TrilinosSparseSpaceType, UblasLocalSpaceType, TrilinosLinearSolverType > >(
                rComm, RowSizeGuess, pLinearSolver);
            return Kratos::make_shared<TrilinosDistanceSmoothingProcess<TDim>>(rModelPart, pLinearSolver, p_builder_solver);
        }));
}

void AddTrilinosProcessesToPython(pybind11::module& m)
{
    // Turbulence models
    using SpalartAllmarasProcess = TrilinosSpalartAllmarasTurbulenceModel<TrilinosSparseSpaceType, UblasLocalSpaceType, TrilinosLinearSolverType>;
    py::class_< SpalartAllmarasProcess, typename SpalartAllmarasProcess::Pointer, Process>(m, "TrilinosSpalartAllmarasTurbulenceModel")
    .def(py::init < Epetra_MpiComm&, ModelPart&, typename TrilinosLinearSolverType::Pointer, unsigned int, double, unsigned int, bool, unsigned int>())
    .def("ActivateDES", &SpalartAllmarasProcess::ActivateDES)
    .def("AdaptForFractionalStep", &SpalartAllmarasProcess::AdaptForFractionalStep)
    ;

    // Stokes initialization processes
    using StokesInitializationProcess = TrilinosStokesInitializationProcess<TrilinosSparseSpaceType, UblasLocalSpaceType, TrilinosLinearSolverType>;
    py::class_< StokesInitializationProcess, typename StokesInitializationProcess::Pointer, Process>(m,"TrilinosStokesInitializationProcess")
    .def(py::init<Epetra_MpiComm&, ModelPart&, typename TrilinosLinearSolverType::Pointer, unsigned int, const Kratos::Variable<int>& >())
    .def("SetConditions",&StokesInitializationProcess::SetConditions)
    ;

    // Distance smoothing processes
    using DistanceSmoothing2DBinderType = py::class_<TrilinosDistanceSmoothingProcess<2>, typename TrilinosDistanceSmoothingProcess<2>::Pointer, Process >;
    using DistanceSmoothing3DBinderType = py::class_<TrilinosDistanceSmoothingProcess<3>, typename TrilinosDistanceSmoothingProcess<3>::Pointer, Process >;

    auto distance_smoothing_2d_binder = DistanceSmoothing2DBinderType(m,"TrilinosDistanceSmoothingProcess2D");
    auto distance_smoothing_3d_binder = DistanceSmoothing3DBinderType(m,"TrilinosDistanceSmoothingProcess3D");

    DistanceSmoothingConstructionHelper<DistanceSmoothing2DBinderType,2>(distance_smoothing_2d_binder);
    DistanceSmoothingConstructionHelper<DistanceSmoothing3DBinderType,3>(distance_smoothing_3d_binder);
}

}
}
