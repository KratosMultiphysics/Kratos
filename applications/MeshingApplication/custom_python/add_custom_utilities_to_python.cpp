// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "custom_python/add_custom_utilities_to_python.h"

/* Utilities */
#include "custom_utilities/projection.h"
#include "custom_utilities/binbased_projection.h"
#include "custom_utilities/local_refine_triangle_mesh.hpp"
#include "custom_utilities/local_refine_prism_mesh.hpp"
#include "custom_utilities/local_refine_tetrahedra_mesh.hpp"

#ifdef  USE_TETGEN_NONFREE_TPL
    #include "custom_utilities/tetgen_volume_mesher.h"
    #include "custom_utilities/tetrahedra_reconnect_utility.h"
#endif 

#include "custom_utilities/cutting_iso_app.h"

#include "utilities/split_tetrahedra.h"

#ifdef PRAGMATIC_ACTIVATED
    #include "external_includes/pragmatic_adapt_3d.h"
#endif

namespace Kratos
{

namespace Python
{
namespace py = pybind11;

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    py::class_<MeshTransfer < 2 > >(m,"MeshTransfer2D")
    .def(py::init< >())
    .def("DirectModelPartInterpolation", &MeshTransfer < 2 > ::DirectInterpolation)
    .def("DirectScalarVarInterpolation", &MeshTransfer < 2 > ::DirectVariableInterpolation<double>)
    .def("DirectVectorialVarInterpolation", &MeshTransfer < 2 > ::DirectVariableInterpolation< array_1d < double, 3 > >)
    ;

    py::class_<MeshTransfer < 3 > >(m,"MeshTransfer3D")
    .def(py::init< >())
    .def("DirectModelPartInterpolation", &MeshTransfer < 3 > ::DirectInterpolation)
    .def("DirectScalarVarInterpolation", &MeshTransfer < 3 > ::DirectVariableInterpolation<double>)
    .def("DirectVectorialVarInterpolation", &MeshTransfer < 3 > ::DirectVariableInterpolation< array_1d < double, 3 > >)
    ;

    py::class_<BinBasedMeshTransfer < 2 > >(m,"BinBasedMeshTransfer2D")
    .def(py::init< >())
    //.def("DirectModelPartInterpolation", &MeshTransfer < 2 > ::DirectInterpolation)
    .def("DirectScalarVarInterpolation", &BinBasedMeshTransfer < 2 > ::DirectVariableInterpolation<double>)
    .def("DirectVectorialVarInterpolation", &BinBasedMeshTransfer < 2 > ::DirectVariableInterpolation< array_1d < double, 3 > >)
    .def("MappingFromMovingMesh_ScalarVar", &BinBasedMeshTransfer < 2 > ::MappingFromMovingMesh<double>)
    .def("MappingFromMovingMesh_VectorialVar", &BinBasedMeshTransfer < 2 > ::MappingFromMovingMesh< array_1d < double, 3 > >)
    .def("MappingFromMovingMesh_VariableMeshes_ScalarVar", &BinBasedMeshTransfer < 2 > ::MappingFromMovingMesh_VariableMeshes<double>)
    .def("MappingFromMovingMesh_VariableMeshes_VectorialVar", &BinBasedMeshTransfer < 2 > ::MappingFromMovingMesh_VariableMeshes< array_1d < double, 3 > >)
    ;

    py::class_<BinBasedMeshTransfer < 3 > >(m,"BinBasedMeshTransfer3D")
    .def(py::init< >())
    //.def("DirectModelPartInterpolation", &MeshTransfer < 3 > ::DirectInterpolation)
    .def("DirectScalarVarInterpolation", &BinBasedMeshTransfer < 3 > ::DirectVariableInterpolation<double>)
    .def("DirectVectorialVarInterpolation", &BinBasedMeshTransfer < 3 > ::DirectVariableInterpolation< array_1d < double, 3 > >)
    .def("MappingFromMovingMesh_ScalarVar", &BinBasedMeshTransfer < 3 > ::MappingFromMovingMesh<double>)
    .def("MappingFromMovingMesh_VectorialVar", &BinBasedMeshTransfer < 3 > ::MappingFromMovingMesh< array_1d < double, 3 > >)
    .def("MappingFromMovingMesh_VariableMeshes_ScalarVar", &BinBasedMeshTransfer < 3 > ::MappingFromMovingMesh_VariableMeshes<double>)
    .def("MappingFromMovingMesh_VariableMeshes_VectorialVar", &BinBasedMeshTransfer < 3 > ::MappingFromMovingMesh_VariableMeshes< array_1d < double, 3 > >)
    ;

    py::class_<LocalRefineTriangleMesh >
    (m,"LocalRefineTriangleMesh")
    .def(py::init<ModelPart&>())
    .def("LocalRefineMesh", &LocalRefineTriangleMesh::LocalRefineMesh)
    ;

    py::class_<LocalRefinePrismMesh >
    (m,"LocalRefinePrismMesh")
    .def(py::init<ModelPart&>())
    .def("LocalRefineMesh", &LocalRefinePrismMesh::LocalRefineMesh)
    ;

    // py::class_<LocalRefineSPrismMesh >
    // (m,"LocalRefineSPrismMesh")
    // .def(py::init<ModelPart&>())
    // .def("LocalRefineMesh", &LocalRefineSPrismMesh::LocalRefineMesh)
    // ;

    py::class_<LocalRefineTetrahedraMesh >
    (m,"LocalRefineTetrahedraMesh")
    .def(py::init<ModelPart&>())
    .def("LocalRefineMesh", &LocalRefineTetrahedraMesh::LocalRefineMesh)
    ;

#ifdef USE_TETGEN_NONFREE_TPL
    py::class_<TetgenVolumeMesher >
    (m,"TetgenVolumeMesher")
    .def(py::init<ModelPart&>())
    .def("AddHole", &TetgenVolumeMesher::AddHole)
    .def("GenerateMesh", &TetgenVolumeMesher::GenerateMesh)
    ;
    
    py::class_<TetrahedraReconnectUtility >(m,"TetrahedraReconnectUtility")
    .def(py::init<ModelPart&>())
    .def("EvaluateQuality", &TetrahedraReconnectUtility::EvaluateQuality)
    .def("TestRemovingElements", &TetrahedraReconnectUtility::TestRemovingElements)
    .def("OptimizeQuality", &TetrahedraReconnectUtility::OptimizeQuality)
    .def("FinalizeOptimization", &TetrahedraReconnectUtility::FinalizeOptimization)
    .def("updateNodesPositions", &TetrahedraReconnectUtility::updateNodesPositions)
    .def("setMaxNumThreads", &TetrahedraReconnectUtility::setMaxNumThreads)
    .def("setBlockSize", &TetrahedraReconnectUtility::setBlockSize)
    .def("isaValidMesh", &TetrahedraReconnectUtility::isaValidMesh)
    ;
#endif
    
#ifdef PRAGMATIC_ACTIVATED
    py::class_<PragmaticAdaptor >(m,"PragmaticAdaptor")
    .def(py::init< >())
    .def("AdaptMesh", &PragmaticAdaptor::AdaptMesh)
    ;
#endif
    
    py::class_<Cutting_Isosurface_Application >(m,"Cutting_Isosurface_Application")
    .def(py::init< >())
    .def("GenerateScalarVarCut", &Cutting_Isosurface_Application::GenerateVariableCut<double>)
    .def("GenerateVectorialComponentVarCut", &Cutting_Isosurface_Application::GenerateVectorialComponentVariableCut<VectorComponentAdaptor< array_1d < double, 3 > > >)
    .def("GenerateVectorialVarCut", &Cutting_Isosurface_Application::GenerateVariableCut< array_1d < double, 3 > >)
    .def("AddModelPartElements", &Cutting_Isosurface_Application::AddModelPartElements)
    .def("AddSkinConditions", &Cutting_Isosurface_Application::AddSkinConditions)
    .def("UpdateCutData", &Cutting_Isosurface_Application::UpdateCutData)
    .def("DeleteCutData", &Cutting_Isosurface_Application::DeleteCutData)
    ;


}
} // namespace Python.

} // Namespace Kratos
