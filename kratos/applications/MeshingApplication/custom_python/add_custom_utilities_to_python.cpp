//
//   Project Name:        Kratos
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-15 11:11:35 $
//   Revision:            $Revision: 1.5 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "custom_utilities/projection.h"
#include "custom_utilities/binbased_projection.h"

//#include "custom_utilities/GenerateModelPartUtilities.h"

#include "custom_utilities/local_refine_triangle_mesh.hpp"
#include "custom_utilities/local_refine_prism_mesh.hpp"
#include "custom_utilities/local_refine_sprism_mesh.hpp"
#include "custom_utilities/local_refine_tetrahedra_mesh.hpp"
#include "custom_utilities/tetgen_volume_mesher.h"
#include "custom_utilities/cutting_iso_app.h"
#include "custom_utilities/tetrahedra_reconnect_utility.h"
#include "utilities/split_tetrahedra.h"

#ifdef PRAGMATIC_ACTIVATED
    #include "external_includes/pragmatic_adapt_3d.h"
#endif


namespace Kratos
{

namespace Python
{

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



void AddCustomUtilitiesToPython()
{

    using namespace boost::python;

    class_<MeshTransfer < 2 > >("MeshTransfer2D", init< >())
    .def("DirectModelPartInterpolation", &MeshTransfer < 2 > ::DirectInterpolation)
    .def("DirectScalarVarInterpolation", &MeshTransfer < 2 > ::DirectVariableInterpolation<double>)
    .def("DirectVectorialVarInterpolation", &MeshTransfer < 2 > ::DirectVariableInterpolation< array_1d < double, 3 > >)
    ;

    class_<MeshTransfer < 3 > >("MeshTransfer3D", init< >())
    .def("DirectModelPartInterpolation", &MeshTransfer < 3 > ::DirectInterpolation)
    .def("DirectScalarVarInterpolation", &MeshTransfer < 3 > ::DirectVariableInterpolation<double>)
    .def("DirectVectorialVarInterpolation", &MeshTransfer < 3 > ::DirectVariableInterpolation< array_1d < double, 3 > >)
    ;

    class_<BinBasedMeshTransfer < 2 > >("BinBasedMeshTransfer2D", init< >())
    //.def("DirectModelPartInterpolation", &MeshTransfer < 2 > ::DirectInterpolation)
    .def("DirectScalarVarInterpolation", &BinBasedMeshTransfer < 2 > ::DirectVariableInterpolation<double>)
    .def("DirectVectorialVarInterpolation", &BinBasedMeshTransfer < 2 > ::DirectVariableInterpolation< array_1d < double, 3 > >)
    .def("MappingFromMovingMesh_ScalarVar", &BinBasedMeshTransfer < 2 > ::MappingFromMovingMesh<double>)
    .def("MappingFromMovingMesh_VectorialVar", &BinBasedMeshTransfer < 2 > ::MappingFromMovingMesh< array_1d < double, 3 > >)
    .def("MappingFromMovingMesh_VariableMeshes_ScalarVar", &BinBasedMeshTransfer < 2 > ::MappingFromMovingMesh_VariableMeshes<double>)
    .def("MappingFromMovingMesh_VariableMeshes_VectorialVar", &BinBasedMeshTransfer < 2 > ::MappingFromMovingMesh_VariableMeshes< array_1d < double, 3 > >)
    ;

    class_<BinBasedMeshTransfer < 3 > >("BinBasedMeshTransfer3D", init< >())
    //.def("DirectModelPartInterpolation", &MeshTransfer < 3 > ::DirectInterpolation)
    .def("DirectScalarVarInterpolation", &BinBasedMeshTransfer < 3 > ::DirectVariableInterpolation<double>)
    .def("DirectVectorialVarInterpolation", &BinBasedMeshTransfer < 3 > ::DirectVariableInterpolation< array_1d < double, 3 > >)
    .def("MappingFromMovingMesh_ScalarVar", &BinBasedMeshTransfer < 3 > ::MappingFromMovingMesh<double>)
    .def("MappingFromMovingMesh_VectorialVar", &BinBasedMeshTransfer < 3 > ::MappingFromMovingMesh< array_1d < double, 3 > >)
    .def("MappingFromMovingMesh_VariableMeshes_ScalarVar", &BinBasedMeshTransfer < 3 > ::MappingFromMovingMesh_VariableMeshes<double>)
    .def("MappingFromMovingMesh_VariableMeshes_VectorialVar", &BinBasedMeshTransfer < 3 > ::MappingFromMovingMesh_VariableMeshes< array_1d < double, 3 > >)
    ;



//            class_<GenerateModelPartUtilities, boost::noncopyable > ("GenerateModelPartUtilities", init< >())
//                    .def("GenerateModelTemperaturePart", GenerateModelTemperaturePart);



    class_<LocalRefineTriangleMesh, boost::noncopyable >
    ("LocalRefineTriangleMesh", init<ModelPart&>())
    .def("LocalRefineMesh", &LocalRefineTriangleMesh::LocalRefineMesh)
    ;
    
    class_<LocalRefinePrismMesh, boost::noncopyable >
    ("LocalRefinePrismMesh", init<ModelPart&>())
    .def("LocalRefineMesh", &LocalRefinePrismMesh::LocalRefineMesh)
    ;

    class_<LocalRefineSPrismMesh, boost::noncopyable >
    ("LocalRefineSPrismMesh", init<ModelPart&>())
    .def("LocalRefineMesh", &LocalRefineSPrismMesh::LocalRefineMesh)
    ;

    class_<LocalRefineTetrahedraMesh, boost::noncopyable >
    ("LocalRefineTetrahedraMesh", init<ModelPart&>())
    .def("LocalRefineMesh", &LocalRefineTetrahedraMesh::LocalRefineMesh)
    ;

    class_<TetgenVolumeMesher, boost::noncopyable >
    ("TetgenVolumeMesher", init<ModelPart&>())
    .def("AddHole", &TetgenVolumeMesher::AddHole)
    .def("GenerateMesh", &TetgenVolumeMesher::GenerateMesh)
    ;

    #ifdef PRAGMATIC_ACTIVATED
    class_<PragmaticAdaptor >("PragmaticAdaptor", init< >())
    .def("AdaptMesh", &PragmaticAdaptor::AdaptMesh)
    ;        
    #endif
    
    class_<Cutting_Isosurface_Application >("Cutting_Isosurface_Application", init< >())
    .def("GenerateScalarVarCut", &Cutting_Isosurface_Application::GenerateVariableCut<double>)
    .def("GenerateVectorialComponentVarCut", &Cutting_Isosurface_Application::GenerateVectorialComponentVariableCut<VectorComponentAdaptor< array_1d < double, 3 > > >)
    .def("GenerateVectorialVarCut", &Cutting_Isosurface_Application::GenerateVariableCut< array_1d < double, 3 > >)
    .def("AddModelPartElements", &Cutting_Isosurface_Application::AddModelPartElements)
    .def("AddSkinConditions", &Cutting_Isosurface_Application::AddSkinConditions)
    .def("UpdateCutData", &Cutting_Isosurface_Application::UpdateCutData)
    .def("DeleteCutData", &Cutting_Isosurface_Application::DeleteCutData)
    ;

    class_<TetrahedraReconnectUtility, boost::noncopyable >("TetrahedraReconnectUtility", init<ModelPart&>())
    .def("EvaluateQuality", &TetrahedraReconnectUtility::EvaluateQuality)
    .def("TestRemovingElements", &TetrahedraReconnectUtility::TestRemovingElements)
    .def("OptimizeQuality", &TetrahedraReconnectUtility::OptimizeQuality)
    .def("FinalizeOptimization", &TetrahedraReconnectUtility::FinalizeOptimization)
    .def("updateNodesPositions", &TetrahedraReconnectUtility::updateNodesPositions)
    .def("setMaxNumThreads", &TetrahedraReconnectUtility::setMaxNumThreads)
    .def("setBlockSize", &TetrahedraReconnectUtility::setBlockSize)
    .def("isaValidMesh", &TetrahedraReconnectUtility::isaValidMesh)
    ;
}
} // namespace Python.

} // Namespace Kratos

