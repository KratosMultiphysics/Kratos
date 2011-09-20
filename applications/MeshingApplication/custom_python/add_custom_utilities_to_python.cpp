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
#include "custom_utilities/GenerateModelPartUtilities.h"
#include "custom_utilities/connectivity_preserve_modeler.h"
#include "custom_utilities/local_triangle_refine_mesh.h"
#include "custom_utilities/local_tetrahedra_refine_mesh.h"
#include "custom_utilities/tetgen_volume_mesher.h"
#include "custom_utilities/cutting_app.h"


namespace Kratos
{

    namespace Python
    {

        void GenerateModelTemperaturePart(GenerateModelPartUtilities& GM, ModelPart& origin_model_part, ModelPart& destination_model_part, unsigned int domain_size)
        {
            if (domain_size == 2)
            {
                GM.GenerateModelPart(origin_model_part, destination_model_part,
                        KratosComponents<Element>::Get("ConvDiff2D"),
                        KratosComponents<Condition>::Get("ThermalFace2D"));
            } else if (domain_size == 3)
            {
                GM.GenerateModelPart(origin_model_part, destination_model_part,
                        KratosComponents<Element>::Get("ConvDiff3D"),
                        KratosComponents<Condition>::Get("ThermalFace3D"));
            }
        }

        void GenerateModelPart(ConnectivityPreserveModeler& GM, ModelPart& origin_model_part, ModelPart& destination_model_part, const char* ElementName, const char* ConditionName)
        {
            GM.GenerateModelPart(origin_model_part, destination_model_part,
                    KratosComponents<Element>::Get(ElementName),
                    KratosComponents<Condition>::Get(ConditionName));

        }

        void AddCustomUtilitiesToPython()
        {

            using namespace boost::python;

            class_<MeshTransfer < 2 > >("MeshTransfer2D", init< >())
                    .def("DirectModelPartInterpolation", &MeshTransfer < 2 > ::DirectInterpolation)
                    .def("DirectScalarVarInterpolation", &MeshTransfer < 2 > ::DirectVariableInterpolation<double>)
                    .def("DirectVectorialVarInterpolation", &MeshTransfer < 2 > ::DirectVariableInterpolation< array_1d < double, 3 > >);

            class_<MeshTransfer < 3 > >("MeshTransfer3D", init< >())
                    .def("DirectModelPartInterpolation", &MeshTransfer < 3 > ::DirectInterpolation)
                    .def("DirectScalarVarInterpolation", &MeshTransfer < 3 > ::DirectVariableInterpolation<double>)
                    .def("DirectVectorialVarInterpolation", &MeshTransfer < 3 > ::DirectVariableInterpolation< array_1d < double, 3 > >);

            class_<GenerateModelPartUtilities, boost::noncopyable > ("GenerateModelPartUtilities", init< >())
                    .def("GenerateModelTemperaturePart", GenerateModelTemperaturePart);

            class_<ConnectivityPreserveModeler, boost::noncopyable > ("ConnectivityPreserveModeler", init< >())
                    .def("GenerateModelPart", GenerateModelPart);

            class_<Local_Refine_Triangle_Mesh, boost::noncopyable >
                    ("LocalRefineTriangleMesh", init<ModelPart&>())
                    .def("LocalRefineMesh", &Local_Refine_Triangle_Mesh::Local_Refine_Mesh)
                    ;

            class_<Local_Refine_Tetrahedra_Mesh, boost::noncopyable >
                    ("LocalRefineTetrahedraMesh", init<ModelPart&>())
                    .def("LocalRefineMesh", &Local_Refine_Tetrahedra_Mesh::Local_Refine_Mesh)
                    ;

            class_<TetgenVolumeMesher, boost::noncopyable >
                    ("TetgenVolumeMesher", init<ModelPart&>())
                    .def("AddHole", &TetgenVolumeMesher::AddHole)
                    .def("GenerateMesh", &TetgenVolumeMesher::GenerateMesh)
                    ;

            class_<Cutting_Application >("Cutting_Application", init< >())
                    .def("GenerateCut", &Cutting_Application::GenerateCut)
                    .def("UpdateCutData", &Cutting_Application ::UpdateCutData)
                    .def("AddSkinConditions", &Cutting_Application ::AddSkinConditions)
                    .def("FindSmallestEdge", &Cutting_Application ::FindSmallestEdge)
                    ;

        }


    } // namespace Python.

} // Namespace Kratos

