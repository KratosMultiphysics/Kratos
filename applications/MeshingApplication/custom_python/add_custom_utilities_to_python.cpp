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
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "custom_utilities/projection.h"
#include "custom_utilities/binbased_projection.h"

//#include "custom_utilities/GenerateModelPartUtilities.h"

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
    
using namespace pybind11;

//        void GenerateModelTemperaturePart(GenerateModelPartUtilities& GM, ModelPart& origin_model_part, ModelPart& destination_model_part, unsigned int domain_size)
//        {
//            if (domain_size == 2)
//            {
//                GM.GenerateModelPart(origin_model_part, destination_model_part,
//                        KratosComponents<Element>::Get("ConvDiff2D"),
//                        KratosComponents<Condition>::Get("ThermalFace2D"));
//            } else if (domain_size == 3)
//            {
//                GM.GenerateModelPart(origin_model_part, destination_model_part,
//                        KratosComponents<Element>::Get("ConvDiff3D"),
//                        KratosComponents<Condition>::Get("ThermalFace3D"));
//            }
//        }



void AddCustomUtilitiesToPython(pybind11::module& m)
{
    class_<MeshTransfer < 2 > >(m,"MeshTransfer2D")
    .def(init< >())
    .def("DirectModelPartInterpolation", &MeshTransfer < 2 > ::DirectInterpolation)
    .def("DirectScalarVarInterpolation", &MeshTransfer < 2 > ::DirectVariableInterpolation<double>)
    .def("DirectVectorialVarInterpolation", &MeshTransfer < 2 > ::DirectVariableInterpolation< array_1d < double, 3 > >)
    ;

    class_<MeshTransfer < 3 > >(m,"MeshTransfer3D")
    .def(init< >())
    .def("DirectModelPartInterpolation", &MeshTransfer < 3 > ::DirectInterpolation)
    .def("DirectScalarVarInterpolation", &MeshTransfer < 3 > ::DirectVariableInterpolation<double>)
    .def("DirectVectorialVarInterpolation", &MeshTransfer < 3 > ::DirectVariableInterpolation< array_1d < double, 3 > >)
    ;

    class_<BinBasedMeshTransfer < 2 > >(m,"BinBasedMeshTransfer2D")
    .def(init< >())
    //.def("DirectModelPartInterpolation", &MeshTransfer < 2 > ::DirectInterpolation)
    .def("DirectScalarVarInterpolation", &BinBasedMeshTransfer < 2 > ::DirectVariableInterpolation<double>)
    .def("DirectVectorialVarInterpolation", &BinBasedMeshTransfer < 2 > ::DirectVariableInterpolation< array_1d < double, 3 > >)
    .def("MappingFromMovingMesh_ScalarVar", &BinBasedMeshTransfer < 2 > ::MappingFromMovingMesh<double>)
    .def("MappingFromMovingMesh_VectorialVar", &BinBasedMeshTransfer < 2 > ::MappingFromMovingMesh< array_1d < double, 3 > >)
    .def("MappingFromMovingMesh_VariableMeshes_ScalarVar", &BinBasedMeshTransfer < 2 > ::MappingFromMovingMesh_VariableMeshes<double>)
    .def("MappingFromMovingMesh_VariableMeshes_VectorialVar", &BinBasedMeshTransfer < 2 > ::MappingFromMovingMesh_VariableMeshes< array_1d < double, 3 > >)
    ;

    class_<BinBasedMeshTransfer < 3 > >(m,"BinBasedMeshTransfer3D")
    .def(init< >())
    //.def("DirectModelPartInterpolation", &MeshTransfer < 3 > ::DirectInterpolation)
    .def("DirectScalarVarInterpolation", &BinBasedMeshTransfer < 3 > ::DirectVariableInterpolation<double>)
    .def("DirectVectorialVarInterpolation", &BinBasedMeshTransfer < 3 > ::DirectVariableInterpolation< array_1d < double, 3 > >)
    .def("MappingFromMovingMesh_ScalarVar", &BinBasedMeshTransfer < 3 > ::MappingFromMovingMesh<double>)
    .def("MappingFromMovingMesh_VectorialVar", &BinBasedMeshTransfer < 3 > ::MappingFromMovingMesh< array_1d < double, 3 > >)
    .def("MappingFromMovingMesh_VariableMeshes_ScalarVar", &BinBasedMeshTransfer < 3 > ::MappingFromMovingMesh_VariableMeshes<double>)
    .def("MappingFromMovingMesh_VariableMeshes_VectorialVar", &BinBasedMeshTransfer < 3 > ::MappingFromMovingMesh_VariableMeshes< array_1d < double, 3 > >)
    ;



//         class_<GenerateModelPartUtilities > ("GenerateModelPartUtilities")
//         .def(init< >())
//         .def("GenerateModelTemperaturePart", GenerateModelTemperaturePart);


    class_<LocalRefineTriangleMesh >
    (m,"LocalRefineTriangleMesh")
    .def(init<ModelPart&>())
    .def("LocalRefineMesh", &LocalRefineTriangleMesh::LocalRefineMesh)
    ;

    class_<LocalRefinePrismMesh >
    (m,"LocalRefinePrismMesh")
    .def(init<ModelPart&>())
    .def("LocalRefineMesh", &LocalRefinePrismMesh::LocalRefineMesh)
    ;

    // class_<LocalRefineSPrismMesh >
    // (m,"LocalRefineSPrismMesh")
    // .def(init<ModelPart&>())
    // .def("LocalRefineMesh", &LocalRefineSPrismMesh::LocalRefineMesh)
    // ;

    class_<LocalRefineTetrahedraMesh >
    (m,"LocalRefineTetrahedraMesh")
    .def(init<ModelPart&>())
    .def("LocalRefineMesh", &LocalRefineTetrahedraMesh::LocalRefineMesh)
    ;

#ifdef USE_TETGEN_NONFREE_TPL
    class_<TetgenVolumeMesher >
    (m,"TetgenVolumeMesher")
    .def(init<ModelPart&>())
    .def("AddHole", &TetgenVolumeMesher::AddHole)
    .def("GenerateMesh", &TetgenVolumeMesher::GenerateMesh)
    ;
    
    class_<TetrahedraReconnectUtility >(m,"TetrahedraReconnectUtility")
    .def(init<ModelPart&>())
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
    class_<PragmaticAdaptor >(m,"PragmaticAdaptor")
    .def(init< >())
    .def("AdaptMesh", &PragmaticAdaptor::AdaptMesh)
    ;
#endif
    
    class_<Cutting_Isosurface_Application >(m,"Cutting_Isosurface_Application")
    .def(init< >())
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
