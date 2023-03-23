//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Ariadna Cortes
//
//

// Project includes
#include "containers/model.h"
#include "includes/node.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_3d_4.h"
#include "utilities/qef_utility.h"
#include "includes/checks.h"
#include "testing/testing.h"
#include "custom_utilities/dual_countouring_mesher.h"

#include "input_output/vtk_output.h"

namespace Kratos {
namespace Testing {
namespace {
    typedef Node<3> NodeType;
    typedef Node<3>::Pointer NodePtrType;
    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType::Pointer GeometryPtrType;
    typedef GeometryType::GeometriesArrayType GeometryArrayType;
    typedef GeometryType::PointsArrayType PointsArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;

    void AddCube(ModelPart& rModelPart) {
        rModelPart.CreateNewNode(1,-0.5,-0.55,-0.8);
        rModelPart.CreateNewNode(2,0.55,-0.55,-0.8);
        rModelPart.CreateNewNode(3,0.55,0.5,-0.8); //s3{-0.5,0.5,-0.5}
        rModelPart.CreateNewNode(4,-0.5,0.5,-0.8);
        rModelPart.CreateNewNode(5,-0.5,-0.55,0.4);
        rModelPart.CreateNewNode(6,0.55,-0.55,0.4);
        rModelPart.CreateNewNode(7,0.55,0.5,0.4);
        rModelPart.CreateNewNode(8,-0.5,0.5,0.4);

        for (std::size_t i = 1; i < 9; i++) {
            rModelPart.pGetNode(i)->FastGetSolutionStepValue(DISTANCE) = 1;
        }
    
        Properties::Pointer p_properties_1(new Properties(0)); 
        rModelPart.CreateNewElement("Element3D3N", 1, {1, 2, 3}, p_properties_1);
        rModelPart.CreateNewElement("Element3D3N", 2, {1, 3, 4}, p_properties_1);
        rModelPart.CreateNewElement("Element3D3N", 3, {1, 2, 5}, p_properties_1);
        rModelPart.CreateNewElement("Element3D3N", 4, {2, 5, 6}, p_properties_1);
        rModelPart.CreateNewElement("Element3D3N", 5, {2, 3, 6}, p_properties_1);
        rModelPart.CreateNewElement("Element3D3N", 6, {3, 6, 7}, p_properties_1);
        rModelPart.CreateNewElement("Element3D3N", 7, {3, 4, 7}, p_properties_1);
        rModelPart.CreateNewElement("Element3D3N", 8, {4, 7, 8}, p_properties_1);
        rModelPart.CreateNewElement("Element3D3N", 9, {4, 1, 8}, p_properties_1);
        rModelPart.CreateNewElement("Element3D3N", 10, {1, 5, 8}, p_properties_1);
        rModelPart.CreateNewElement("Element3D3N", 11, {5, 6, 7}, p_properties_1);
        rModelPart.CreateNewElement("Element3D3N", 12, {5, 7, 8}, p_properties_1);
    }

    void Output(ModelPart& rModelPart2, std::string Name) {
        Parameters this_parameters_skin_part = Parameters(R"({
            "model_part_name"        : "skin_part",
            "output_sub_model_parts" : false,
            "nodal_solution_step_data_variables" : ["DISTANCE"]
        })");

        VtkOutput vtk_skin_part(rModelPart2,this_parameters_skin_part);
        
        vtk_skin_part.PrintOutput("vtk_" + Name);
    }

    void WriteCubeSkinMeshMdpaFile()
    {
        Kratos::shared_ptr<std::iostream> p_input(new std::stringstream(
            R"input(
        Begin ModelPartData
    End ModelPartData
    Begin Properties 0
    End Properties
    Begin Properties 1
    End Properties
    Begin Nodes
            1     0.35     0.4    -0.91    
            2     0.95    -0.3    -0.91    
            3     0.35    -0.3    0.7    
            4     0.45     0.25    0.3    
            5    -0.35     0.2    -0.81    
            6    -0.5     0.32    0.3    
            7    -0.35    -0.3    -0.91    
            8    -0.5    -0.6    0.3    
    End Nodes
    Begin Elements Element3D3N
            1     1    1    2    3    
            2     1    3    4    1    
            3     1    5    1    4    
            4     1    4    6    5    
            5     1    7    5    6    
            6     1    6    8    7    
            7     1    7    8    3    
            8     1    3    2    7    
            9     1    7    2    1    
            10    1    1    5    7    
            11    1    8    6    4    
            12    1    4    3    8    
    End Elements
    Begin SubModelPart workpiece
    Begin SubModelPartNodes
            1
            2
            3
            4
            5
            6
            7
            8
    End SubModelPartNodes
    Begin SubModelPartElements
            1
            2
            3
            4
            5
            6
            7
            8
            9
            10
            11
            12
    End SubModelPartElements
    End SubModelPart
        )input"));

      

        Model current_model;
        ModelPartIO model_part_io(p_input);

        ModelPart& model_part_0 =
            current_model.CreateModelPart("tmp_main_model_part");
        model_part_io.ReadModelPart(model_part_0);

        // Create the output .mdpa file
        std::string output_file_name = "cube_skin_mesh";
        std::fstream output_file;
        output_file.open(output_file_name + ".mdpa", std::fstream::out);
        output_file.close();

        // Fill the output .mdpa file
        ModelPartIO* model_part_io_write = new ModelPartIO(output_file_name, IO::WRITE);
        model_part_io_write->WriteModelPart(model_part_0);
    }
}

    KRATOS_TEST_CASE_IN_SUITE(DualCountouringRemesherCube, KratosMeshingApplicationFastSuite) {

        WriteCubeSkinMeshMdpaFile();

        Parameters mesher_parameters(R"(
        {
            "output_model_part_name" : "voxel_model_part",
            "input_model_part_name" : "skin_model_part",
            "mdpa_file_name" : "cube_skin_mesh",
            "key_plane_generator": {
                "Parameters" : {
                    "voxel_sizes" : [0.14, 0.14, 0.14],
                    "min_point" : [-1, -1, -1],
                    "max_point" : [1, 1, 1]
                }
            },
            "coloring_settings_list": [
            {
                "type" : "cells_in_touch",
                "model_part_name": "skin_model_part.workpiece",
                "color": -1
            },
            {
                "type" : "cells_with_inside_center",
                "model_part_name": "skin_model_part.workpiece",
                "color": -1
            }
            ],
            "entities_generator_list": [
            {
                "type" : "elements_with_cell_color",
                "model_part_name": "voxel_model_part.workpiece",
                "color": -1,
                "properties_id": 1
            } 
            ]
        })");

        Model my_model;
        ModelPart& voxels_part = my_model.CreateModelPart("voxel_model_part");
        voxels_part.AddNodalSolutionStepVariable(DISTANCE);
        ModelPart& skin_model_part = my_model.CreateModelPart("skin_model_part");     
        ModelPart& fited_mesh = my_model.CreateModelPart("fited_mesh"); 
        fited_mesh.AddNodalSolutionStepVariable(DISTANCE);  

        //AddCube(skin_model_part);

        DualCountouringMesher modeler(my_model,mesher_parameters); 
        modeler.SetupGeometryModel();
        modeler.PrepareGeometryModel();
        modeler.SetupModelPart(); 

        Output(skin_model_part,"cube_pre");

        KRATOS_CHECK_EQUAL(skin_model_part.Nodes().size(),8);
        KRATOS_CHECK_EQUAL(skin_model_part.Elements().size(),12); 
    
        modeler.DualCountourAdaptativeRemesh(fited_mesh);

        Output(fited_mesh,"cube_post");

        KRATOS_CHECK_EQUAL(fited_mesh.Elements().size(),269); 
    }
} //Namespace Testing
} //Namespace Kratos