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

    py::class_<PartitionedFSIUtilities<TSpace,double,2>, PartitionedFSIUtilities<TSpace,double,2>::Pointer>(m,"PartitionedFSIUtilitiesDouble2D")
        .def(py::init<>())
        .def("GetInterfaceArea",&PartitionedFSIUtilities<TSpace,double,2>::GetInterfaceArea)
        .def("GetInterfaceResidualSize",&PartitionedFSIUtilities<TSpace,double,2>::GetInterfaceResidualSize)
        .def("UpdateInterfaceValues",&PartitionedFSIUtilities<TSpace,double,2>::UpdateInterfaceValues)
        .def("ComputeInterfaceResidualNorm",&PartitionedFSIUtilities<TSpace,double,2>::ComputeInterfaceResidualNorm)
        .def("ComputeInterfaceResidualVector",&PartitionedFSIUtilities<TSpace,double,2>::ComputeInterfaceResidualVector)
        .def("ComputeAndPrintFluidInterfaceNorms",&PartitionedFSIUtilities<TSpace,double,2>::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms",&PartitionedFSIUtilities<TSpace,double,2>::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid",&PartitionedFSIUtilities<TSpace,double,2>::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure",&PartitionedFSIUtilities<TSpace,double,2>::CheckCurrentCoordinatesStructure)
        ;

    py::class_<PartitionedFSIUtilities<TSpace,array_1d<double,3>,2>, PartitionedFSIUtilities<TSpace,array_1d<double,3>,2>::Pointer>(m,"PartitionedFSIUtilitiesArray2D")
        .def(py::init<>())
        .def("GetInterfaceArea",&PartitionedFSIUtilities<TSpace,array_1d<double,3>,2>::GetInterfaceArea)
        .def("GetInterfaceResidualSize",&PartitionedFSIUtilities<TSpace,array_1d<double,3>,2>::GetInterfaceResidualSize)
        .def("UpdateInterfaceValues",&PartitionedFSIUtilities<TSpace,array_1d<double,3>,2>::UpdateInterfaceValues)
        .def("ComputeInterfaceResidualNorm",&PartitionedFSIUtilities<TSpace,array_1d<double,3>,2>::ComputeInterfaceResidualNorm)
        .def("ComputeInterfaceResidualVector",&PartitionedFSIUtilities<TSpace,array_1d<double,3>,2>::ComputeInterfaceResidualVector)
        .def("ComputeAndPrintFluidInterfaceNorms",&PartitionedFSIUtilities<TSpace,array_1d<double,3>,2>::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms",&PartitionedFSIUtilities<TSpace,array_1d<double,3>,2>::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid",&PartitionedFSIUtilities<TSpace,array_1d<double,3>,2>::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure",&PartitionedFSIUtilities<TSpace,array_1d<double,3>,2>::CheckCurrentCoordinatesStructure)
        ;

    py::class_<PartitionedFSIUtilities<TSpace,double,3>, PartitionedFSIUtilities<TSpace,double,3>::Pointer>(m,"PartitionedFSIUtilitiesDouble3D")
        .def(py::init<>())
        .def("GetInterfaceArea", &PartitionedFSIUtilities<TSpace,double,3>::GetInterfaceArea)
        .def("GetInterfaceResidualSize", &PartitionedFSIUtilities<TSpace,double,3>::GetInterfaceResidualSize)
        .def("UpdateInterfaceValues", &PartitionedFSIUtilities<TSpace,double,3>::UpdateInterfaceValues)
        .def("ComputeInterfaceResidualNorm", &PartitionedFSIUtilities<TSpace,double,3>::ComputeInterfaceResidualNorm)
        .def("ComputeInterfaceResidualVector", &PartitionedFSIUtilities<TSpace,double,3>::ComputeInterfaceResidualVector)
        .def("ComputeAndPrintFluidInterfaceNorms", &PartitionedFSIUtilities<TSpace,double,3>::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms", &PartitionedFSIUtilities<TSpace,double,3>::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid", &PartitionedFSIUtilities<TSpace,double,3>::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure", &PartitionedFSIUtilities<TSpace,double,3>::CheckCurrentCoordinatesStructure);

    py::class_<PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>, PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::Pointer>(m,"PartitionedFSIUtilitiesArray3D")
        .def(py::init<>())
        .def("GetInterfaceArea", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::GetInterfaceArea)
        .def("GetInterfaceResidualSize", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::GetInterfaceResidualSize)
        .def("UpdateInterfaceValues", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::UpdateInterfaceValues)
        .def("ComputeInterfaceResidualNorm", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::ComputeInterfaceResidualNorm)
        .def("ComputeInterfaceResidualVector", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::ComputeInterfaceResidualVector)
        .def("ComputeAndPrintFluidInterfaceNorms", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure", &PartitionedFSIUtilities<TSpace,array_1d<double,3>,3>::CheckCurrentCoordinatesStructure);

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
