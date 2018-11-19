//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi, Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "spaces/ublas_space.h"
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/FSI_utils.h"
#include "custom_utilities/aitken_utils.h"
#include "custom_utilities/partitioned_fsi_utilities.hpp"
#include "custom_utilities/nodal_update_utilities.h"
namespace Kratos
{

namespace Python
{

void AddCustomUtilitiesToPython(pybind11::module &m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, Matrix, Vector > TSpace;
    typedef NodalUpdateBaseClass< 2 > NodalUpdateBaseClass2DType;
    typedef NodalUpdateBaseClass< 3 > NodalUpdateBaseClass3DType;

    py::class_<FSIUtils>(m,"FSIUtils")
        .def(py::init<>())
        .def("CheckPressureConvergence",&FSIUtils::CheckPressureConvergence)
        .def("StructuralPressurePrediction",&FSIUtils::StructuralPressurePrediction)
        ;

    py::class_<AitkenUtils>(m,"AitkenUtils")
        .def(py::init<>())
        .def("ComputeAitkenFactor",&AitkenUtils::ComputeAitkenFactor)
        .def("ComputeRelaxedDisplacement",&AitkenUtils::ComputeRelaxedDisplacement)
        ;

    py::class_<PartitionedFSIUtilities<TSpace,2>, PartitionedFSIUtilities<TSpace,2>::Pointer>(m,"PartitionedFSIUtilities2D")
        .def(py::init<>())
        .def("CopySkinToElements",&PartitionedFSIUtilities<TSpace,2>::CopySkinToElements)
        .def("GetInterfaceArea",&PartitionedFSIUtilities<TSpace,2>::GetInterfaceArea)
        .def("GetInterfaceResidualSize",&PartitionedFSIUtilities<TSpace,2>::GetInterfaceResidualSize)
        .def("UpdateInterfaceValues",&PartitionedFSIUtilities<TSpace,2>::UpdateInterfaceValues)
        .def("ComputeInterfaceResidualVector",&PartitionedFSIUtilities<TSpace,2>::ComputeInterfaceResidualVector)
        .def("ComputeFluidInterfaceMeshVelocityResidualNorm",&PartitionedFSIUtilities<TSpace,2>::ComputeFluidInterfaceMeshVelocityResidualNorm)
        .def("ComputeAndPrintFluidInterfaceNorms",&PartitionedFSIUtilities<TSpace,2>::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms",&PartitionedFSIUtilities<TSpace,2>::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid",&PartitionedFSIUtilities<TSpace,2>::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure",&PartitionedFSIUtilities<TSpace,2>::CheckCurrentCoordinatesStructure)
        ;

    py::class_<PartitionedFSIUtilities<TSpace,3>, PartitionedFSIUtilities<TSpace,3>::Pointer>(m,"PartitionedFSIUtilities3D")
        .def(py::init<>())
        .def("CopySkinToElements", &PartitionedFSIUtilities<TSpace, 3>::CopySkinToElements)
        .def("GetInterfaceArea", &PartitionedFSIUtilities<TSpace, 3>::GetInterfaceArea)
        .def("GetInterfaceResidualSize", &PartitionedFSIUtilities<TSpace, 3>::GetInterfaceResidualSize)
        .def("UpdateInterfaceValues", &PartitionedFSIUtilities<TSpace, 3>::UpdateInterfaceValues)
        .def("ComputeInterfaceResidualVector", &PartitionedFSIUtilities<TSpace, 3>::ComputeInterfaceResidualVector)
        .def("ComputeFluidInterfaceMeshVelocityResidualNorm", &PartitionedFSIUtilities<TSpace, 3>::ComputeFluidInterfaceMeshVelocityResidualNorm)
        .def("ComputeAndPrintFluidInterfaceNorms", &PartitionedFSIUtilities<TSpace, 3>::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms", &PartitionedFSIUtilities<TSpace, 3>::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid", &PartitionedFSIUtilities<TSpace, 3>::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure", &PartitionedFSIUtilities<TSpace, 3>::CheckCurrentCoordinatesStructure);

    py::class_<NodalUpdateBaseClass<2>>(m,"BaseNodalUpdate2D")
        .def(py::init<>())
        .def("UpdateMeshTimeDerivatives", &NodalUpdateBaseClass<2>::UpdateMeshTimeDerivatives)
        .def("SetMeshTimeDerivativesOnInterface", &NodalUpdateBaseClass<2>::SetMeshTimeDerivativesOnInterface);

    py::class_<NodalUpdateBaseClass<3>>(m,"BaseNodalUpdate3D")
        .def(py::init<>())
        .def("UpdateMeshTimeDerivatives", &NodalUpdateBaseClass<3>::UpdateMeshTimeDerivatives)
        .def("SetMeshTimeDerivativesOnInterface", &NodalUpdateBaseClass<3>::SetMeshTimeDerivativesOnInterface);

    py::class_<NodalUpdateNewmark<2>, NodalUpdateBaseClass2DType>(m,"NodalUpdateNewmark2D")
        .def(py::init<const double>())
        .def("UpdateMeshTimeDerivatives", &NodalUpdateNewmark<2>::UpdateMeshTimeDerivatives)
        .def("SetMeshTimeDerivativesOnInterface", &NodalUpdateNewmark<2>::SetMeshTimeDerivativesOnInterface);

    py::class_<NodalUpdateNewmark<3>, NodalUpdateBaseClass3DType>(m,"NodalUpdateNewmark3D")
        .def(py::init<const double>())
        .def("UpdateMeshTimeDerivatives", &NodalUpdateNewmark<3>::UpdateMeshTimeDerivatives)
        .def("SetMeshTimeDerivativesOnInterface", &NodalUpdateNewmark<3>::SetMeshTimeDerivativesOnInterface);
}

}  // namespace Python.

} // Namespace Kratos
