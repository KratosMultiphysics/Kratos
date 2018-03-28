// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/ublas_interface.h"

//Utilities
#include "custom_utilities/tree_contact_search.h"
#include "custom_utilities/process_factory_utility.h"
#include "custom_utilities/sparse_matrix_multiplication_utility.h"

namespace Kratos
{
namespace Python
{
void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;

    // Tree contact search
    class_<TreeContactSearch<2, 2>>("TreeContactSearch2D2N", init<ModelPart&>())
    .def(init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&TreeContactSearch<2, 2>::InitializeMortarConditions)
    .def("ClearMortarConditions",&TreeContactSearch<2, 2>::ClearMortarConditions)
    .def("CreatePointListMortar",&TreeContactSearch<2, 2>::CreatePointListMortar)
    .def("UpdatePointListMortar",&TreeContactSearch<2, 2>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&TreeContactSearch<2, 2>::UpdateMortarConditions)
    .def("ResetContactOperators",&TreeContactSearch<2, 2>::ResetContactOperators)
    .def("CheckMortarConditions",&TreeContactSearch<2, 2>::CheckMortarConditions)
    .def("InvertSearch",&TreeContactSearch<2, 2>::InvertSearch)
    ;
    class_<TreeContactSearch<3, 3>>("TreeContactSearch3D3N", init<ModelPart&>())
    .def(init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&TreeContactSearch<3, 3>::InitializeMortarConditions)
    .def("ClearMortarConditions",&TreeContactSearch<3, 3>::ClearMortarConditions)
    .def("CreatePointListMortar",&TreeContactSearch<3, 3>::CreatePointListMortar)
    .def("UpdatePointListMortar",&TreeContactSearch<3, 3>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&TreeContactSearch<3, 3>::UpdateMortarConditions)
    .def("ResetContactOperators",&TreeContactSearch<3, 3>::ResetContactOperators)
    .def("CheckMortarConditions",&TreeContactSearch<3, 3>::CheckMortarConditions)
    .def("InvertSearch",&TreeContactSearch<3, 3>::InvertSearch)
    ;
    class_<TreeContactSearch<3, 4>>("TreeContactSearch3D4N", init<ModelPart&>())
    .def(init<ModelPart&, Parameters>())
    .def("InitializeMortarConditions",&TreeContactSearch<3, 4>::InitializeMortarConditions)
    .def("ClearMortarConditions",&TreeContactSearch<3, 4>::ClearMortarConditions)
    .def("CreatePointListMortar",&TreeContactSearch<3, 4>::CreatePointListMortar)
    .def("UpdatePointListMortar",&TreeContactSearch<3, 4>::UpdatePointListMortar)
    .def("UpdateMortarConditions",&TreeContactSearch<3, 4>::UpdateMortarConditions)
    .def("ResetContactOperators",&TreeContactSearch<3, 4>::ResetContactOperators)
    .def("CheckMortarConditions",&TreeContactSearch<3, 4>::CheckMortarConditions)
    .def("InvertSearch",&TreeContactSearch<3, 4>::InvertSearch)
    ;

    // Adding search related enums
    enum_<SearchTreeType>("SearchTreeType")
    .value("KdtreeInRadius", KdtreeInRadius)
    .value("KdtreeInBox", KdtreeInBox)
    .value("Kdop", Kdop)
    ;

    enum_<CheckGap>("CheckGap")
    .value("NoCheck", NoCheck)
    .value("DirectCheck", DirectCheck)
    .value("MappingCheck", MappingCheck)
    ;

    // Process Factory utility
    class_<ProcessFactoryUtility>("ProcessFactoryUtility", init<>())
    .def(init<boost::python::list&>())
    .def(init<boost::python::object&>())
    .def("AddProcess",&ProcessFactoryUtility::AddProcess)
    .def("AddProcesses",&ProcessFactoryUtility::AddProcesses)
    .def("ExecuteMethod",&ProcessFactoryUtility::ExecuteMethod)
    .def("ExecuteInitialize",&ProcessFactoryUtility::ExecuteInitialize)
    .def("ExecuteBeforeSolutionLoop",&ProcessFactoryUtility::ExecuteBeforeSolutionLoop)
    .def("ExecuteInitializeSolutionStep",&ProcessFactoryUtility::ExecuteInitializeSolutionStep)
    .def("ExecuteFinalizeSolutionStep",&ProcessFactoryUtility::ExecuteFinalizeSolutionStep)
    .def("ExecuteBeforeOutputStep",&ProcessFactoryUtility::ExecuteBeforeOutputStep)
    .def("ExecuteAfterOutputStep",&ProcessFactoryUtility::ExecuteAfterOutputStep)
    .def("ExecuteFinalize",&ProcessFactoryUtility::ExecuteFinalize)
    .def("IsOutputStep",&ProcessFactoryUtility::IsOutputStep)
    .def("PrintOutput",&ProcessFactoryUtility::PrintOutput)
    .def("Clear",&ProcessFactoryUtility::Clear)
    ;

    // Sparse matrix multiplication utility
    class_<SparseMatrixMultiplicationUtility>("SparseMatrixMultiplicationUtility", init<>())
    .def("MatrixMultiplicationSaad",&SparseMatrixMultiplicationUtility::MatrixMultiplicationSaad<CompressedMatrix, CompressedMatrix, CompressedMatrix>)
    .def("MatrixMultiplicationRMerge",&SparseMatrixMultiplicationUtility::MatrixMultiplicationRMerge<CompressedMatrix, CompressedMatrix, CompressedMatrix>)
    .def("MatrixAdd",&SparseMatrixMultiplicationUtility::MatrixAdd<CompressedMatrix, CompressedMatrix>)
    ;
}

}  // namespace Python.

} // Namespace Kratos

