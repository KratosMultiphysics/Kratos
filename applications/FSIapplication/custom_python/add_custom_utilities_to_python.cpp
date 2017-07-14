//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Original author:  Riccardo Rossi, Ruben Zorrilla
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "spaces/ublas_space.h"
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/FSI_utils.h"
#include "custom_utilities/aitken_utils.h"
#include "custom_utilities/partitioned_fsi_utilities.hpp"
#include "custom_utilities/variable_redistribution_utility.h"

namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;
    typedef UblasSpace<double, Matrix, Vector > TSpace;

    class_<FSIUtils>("FSIUtils", init<>())
//		.def("FSIUtils",&FSIUtils::GenerateCouplingElements)
    .def("CheckPressureConvergence",&FSIUtils::CheckPressureConvergence)
    .def("StructuralPressurePrediction",&FSIUtils::StructuralPressurePrediction)
    ;

    class_<AitkenUtils>("AitkenUtils", init<>())
    .def("ComputeAitkenFactor",&AitkenUtils::ComputeAitkenFactor)
    .def("ComputeRelaxedDisplacement",&AitkenUtils::ComputeRelaxedDisplacement)
    ;

    class_<PartitionedFSIUtilities <TSpace,2>, boost::noncopyable >("PartitionedFSIUtilities2D", init< >())
    .def("GetInterfaceArea",&PartitionedFSIUtilities<TSpace,2>::GetInterfaceArea)
    .def("GetInterfaceResidualSize",&PartitionedFSIUtilities<TSpace,2>::GetInterfaceResidualSize)
    .def("SetInterfaceVectorVariable",&PartitionedFSIUtilities<TSpace,2>::SetInterfaceVectorVariable)
    .def("SetAndFixInterfaceVectorVariable",&PartitionedFSIUtilities<TSpace,2>::SetAndFixInterfaceVectorVariable)
    .def("ComputeInterfaceVectorResidual",&PartitionedFSIUtilities<TSpace,2>::ComputeInterfaceVectorResidual)
    .def("ComputeFluidInterfaceMeshVelocityResidualNorm",&PartitionedFSIUtilities<TSpace,2>::ComputeFluidInterfaceMeshVelocityResidualNorm)
    ;

    class_<PartitionedFSIUtilities <TSpace,3>, boost::noncopyable >("PartitionedFSIUtilities3D", init< >())
    .def("GetInterfaceArea",&PartitionedFSIUtilities<TSpace,3>::GetInterfaceArea)
    .def("GetInterfaceResidualSize",&PartitionedFSIUtilities<TSpace,3>::GetInterfaceResidualSize)
    .def("SetInterfaceVectorVariable",&PartitionedFSIUtilities<TSpace,3>::SetInterfaceVectorVariable)
    .def("SetAndFixInterfaceVectorVariable",&PartitionedFSIUtilities<TSpace,3>::SetAndFixInterfaceVectorVariable)
    .def("ComputeInterfaceVectorResidual",&PartitionedFSIUtilities<TSpace,3>::ComputeInterfaceVectorResidual)
    .def("ComputeFluidInterfaceMeshVelocityResidualNorm",&PartitionedFSIUtilities<TSpace,3>::ComputeFluidInterfaceMeshVelocityResidualNorm)
    ;

    typedef void (*DistributePointDoubleType)(ModelPart&, const Variable< double >&, const Variable< double >&, double, unsigned int);
    typedef void (*DistributePointArrayType)(ModelPart&, const Variable< array_1d<double,3> >&, const Variable< array_1d<double,3> >&,double, unsigned int);

    DistributePointDoubleType DistributePointDouble = &VariableRedistributionUtility::DistributePointValues;
    DistributePointArrayType  DistributePointArray  = &VariableRedistributionUtility::DistributePointValues;
    
    typedef void (*ConvertDistributedDoubleType)(ModelPart&, const Variable< double >&, const Variable< double >&);
    typedef void (*ConvertDistributedArrayType)(ModelPart&, const Variable< array_1d<double,3> >&, const Variable< array_1d<double,3> >&);

    ConvertDistributedDoubleType ConvertDistributedDouble = &VariableRedistributionUtility::ConvertDistributedValuesToPoint;
    ConvertDistributedArrayType  ConvertDistributedArray  = &VariableRedistributionUtility::ConvertDistributedValuesToPoint;

    // Note: The StaticMethod thing should be done only once for each set of overloads
    class_< VariableRedistributionUtility, boost::noncopyable >("VariableRedistributionUtilities", no_init)
    .def("DistributePointValues",DistributePointDouble)
    .def("DistributePointValues",DistributePointArray).staticmethod("DistributePointValues")
    .def("ConvertDistributedValuesToPoint",ConvertDistributedDouble)
    .def("ConvertDistributedValuesToPoint",ConvertDistributedArray).staticmethod("ConvertDistributedValuesToPoint")
    ;
}

}  // namespace Python.

} // Namespace Kratos
