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

    void AddCubeLimit(ModelPart& rModelPart) {
        rModelPart.CreateNewNode(1,-0.95,-0.99,-0.98);
        rModelPart.CreateNewNode(2,0.55,-0.99,-0.98);
        rModelPart.CreateNewNode(3,0.55,0.5,-0.98); 
        rModelPart.CreateNewNode(4,-0.95,0.5,-0.98);
        rModelPart.CreateNewNode(5,-0.95,-0.99,0.4);
        rModelPart.CreateNewNode(6,0.55,-0.99,0.4);
        rModelPart.CreateNewNode(7,0.55,0.5,0.4);
        rModelPart.CreateNewNode(8,-0.95,0.5,0.4);

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

    void AddRandskin_part(ModelPart& rModelPart) {
        rModelPart.CreateNewNode(1,-0.13,-0.55,-0.8);
        rModelPart.CreateNewNode(2,0.55,-0.55,-0.8);
        rModelPart.CreateNewNode(3,0.55,0.5,-0.9); 
        rModelPart.CreateNewNode(4,-0.4,0.5,-0.6);
        rModelPart.CreateNewNode(5,-0.8,-0.65,0.4);
        rModelPart.CreateNewNode(6,0.55,-0.3,0.666);
        rModelPart.CreateNewNode(7,0.55,0.5,0.3);
        rModelPart.CreateNewNode(8,-0.7,0.5,0.6);

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

    void AddTetra(ModelPart& rModelPart) {

        rModelPart.CreateNewNode(1,-0.7,-0.5,-0.5);
        rModelPart.CreateNewNode(2,0.5,-0.25,-0.5);
        rModelPart.CreateNewNode(3,0.5,0.5,-0.5); //s3{-0.5,0.5,-0.5}
        rModelPart.CreateNewNode(4,-0.5,0.7,0.8);

        for (std::size_t i = 1; i < 5; i++) {
            rModelPart.pGetNode(i)->FastGetSolutionStepValue(DISTANCE) = 1;
        }
    
        Properties::Pointer p_properties_1(new Properties(0)); 
        rModelPart.CreateNewElement("Element3D3N", 1, {1, 2, 3}, p_properties_1);
        rModelPart.CreateNewElement("Element3D3N", 2, {1, 3, 4}, p_properties_1);
        rModelPart.CreateNewElement("Element3D3N", 3, {1, 2, 4}, p_properties_1);
        rModelPart.CreateNewElement("Element3D3N", 4, {2, 3, 4}, p_properties_1);
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
}

     KRATOS_TEST_CASE_IN_SUITE(DualCountouringRemesherCube, KratosMeshingApplicationFastSuite) {

        Kratos::shared_ptr<std::iostream> p_input(new std::stringstream(
R"input(
Begin ModelPartData
End ModelPartData
Begin Properties 0
End Properties
Begin Properties 1
End Properties
Begin Nodes
		1     0.03     0.03    0.09    
		2     0.03    -0.03    0.09    
		3     0.03    -0.03    0.03    
		4     0.03     0.03    0.03    
		5    -0.03     0.03    0.09    
		6    -0.03     0.03    0.03    
		7    -0.03    -0.03    0.09    
		8    -0.03    -0.03    0.03    
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

        Model my_model;
        ModelPart& skin_model_part = my_model.CreateModelPart("skin_model_part");
        ModelPart& voxels_part = my_model.CreateModelPart("voxels_part");
        ModelPart& fited_mesh = my_model.CreateModelPart("fited_mesh");
        skin_model_part.AddNodalSolutionStepVariable(DISTANCE);
        voxels_part.AddNodalSolutionStepVariable(DISTANCE);
        fited_mesh.AddNodalSolutionStepVariable(DISTANCE);  

        //ModelPartIO model_part_io(p_input);
        //model_part_io.ReadModelPart(skin_model_part);

        //KRATOS_WATCH(voxels_part.Elements().size());

        AddCube(skin_model_part);
        KRATOS_CHECK_EQUAL(skin_model_part.Nodes().size(),8);
        KRATOS_CHECK_EQUAL(skin_model_part.Elements().size(),12);

        //Output(skin_model_part,"cube_pre");

        DualCountouringMesher mesher; 
        mesher.DualCountourAdaptativeRemesh(skin_model_part, voxels_part, fited_mesh); 
        KRATOS_WATCH('yas');
        
        //Output(skin_model_part,"cube_post");
        
        //KRATOS_CHECK_EQUAL(voxels_part.Elements().size(),1);
    }

    KRATOS_TEST_CASE_IN_SUITE(DualCountouringRemesherCubeLimit, KratosMeshingApplicationFastSuite) {

        Model my_model;
        ModelPart& voxels = my_model.CreateModelPart("Voxels");
        ModelPart& skin_part = my_model.CreateModelPart("skin_part");
        voxels.AddNodalSolutionStepVariable(DISTANCE); 
        skin_part.AddNodalSolutionStepVariable(DISTANCE); 
        
        AddCubeLimit(skin_part);
        Output(skin_part,"cube_limit_pre");

        //DualCountouringMesher mesher; 
        //mesher.DualCountourAdaptativeRemesh(voxels, skin_part); 

        Output(skin_part,"cube_limit_post");
    }

    KRATOS_TEST_CASE_IN_SUITE(DualCountouringRemesherRanskin_part, KratosMeshingApplicationFastSuite) {

        Model my_model;
        ModelPart& voxels = my_model.CreateModelPart("Voxels");
        ModelPart& skin_part = my_model.CreateModelPart("skin_part");
        voxels.AddNodalSolutionStepVariable(DISTANCE); 
        skin_part.AddNodalSolutionStepVariable(DISTANCE); 
        
        AddRandskin_part(skin_part);

        Output(skin_part,"randskin_part_pre");

        array_1d<double,3> cell_size{0.1,0.1,0.1};
        //GeometricalObjectsBins bins(skin_part.ElementsBegin(), skin_part.ElementsEnd(), cell_size);

        DualCountouringMesher mesher; 
        //mesher.DualCountourAdaptativeRemesh(bins, skin_part); 

        Output(skin_part,"randskin_part_post");
    }
} //Namespace Testing
} //Namespace Kratos