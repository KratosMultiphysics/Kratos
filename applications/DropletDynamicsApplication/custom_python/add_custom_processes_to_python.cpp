//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//
//


// System includes

// External includes
#ifdef KRATOS_USE_AMATRIX
#include "boost/numeric/ublas/matrix.hpp" // for the sparse space dense vector
#endif // KRATOS_USE_AMATRIX

// Project includes
#include "containers/model.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/define_python.h"
#include "includes/model_part.h"
#include "processes/process.h"


#include "custom_processes/ed_reinitialization_process.h"

namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double> > SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;


    py::class_<EDReinitializationProcess, EDReinitializationProcess::Pointer, Process>(m,"EDReinitializationProcess3D")
    .def(py::init< ModelPart&, LinearSolverType::Pointer >())
    .def("Clear", &EDReinitializationProcess::Clear)   
    ;
    // .def(py::init< ModelPart&, Parameters >())
    // .def(py::init< Model&, Parameters >())

    // py::class_<EDReinitializationProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>, EDReinitializationProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"EDReinitializationProcess2D")
    // .def(py::init< ModelPart&, LinearSolverType::Pointer >())
    // .def("Clear", &EDReinitializationProcess<2,SparseSpaceType,LocalSpaceType,LinearSolverType>::Clear)   
    ;
    // py::class_<CalculateElectricForceProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>, CalculateElectricForceProcess<3,SparseSpaceType,LocalSpaceType,LinearSolverType>::Pointer, Process>(m,"CalculateElectricForceProcess3D")
    // .def(py::init< ModelPart&, LinearSolverType::Pointer >())
    // .def(py::init< ModelPart&, Parameters >())
    // .def(py::init< Model&, Parameters >())
    // ;

}

} // namespace Python.

} // Namespace Kratos
