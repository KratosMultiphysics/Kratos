//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pablo Agustin Becker
//

// System includes
#include <limits>

// External includes

// kratos includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/element.h"
#include "includes/global_pointer_variables.h"
#include "geometries/geometry.h"
#include "processes/generic_find_elements_neighbours_process.h"
#include "includes/model_part_io.h"
#include "testing/testing.h"

namespace Kratos::Testing {


KRATOS_TEST_CASE_IN_SUITE(HexasGenericFindElementsNeighbourProcessTest,
                          KratosCoreFastSuite)
{
    using namespace Kratos;

    Model model;
    ModelPart& model_part = model.CreateModelPart("test_model_part");
    model_part.AddNodalSolutionStepVariable(DISTANCE);

    // read mdpa
    //the model part has 3 connected hexas and one disconnected one
    Kratos::shared_ptr<std::iostream> p_input(new std::stringstream(
    R"input(
	Begin Properties 1
    End Properties
    Begin Nodes
        1	0	0	0
        2	1	0	0
        3	2	0	0
        4	0	1	0
        5	1	1	0
        6	2	1	0
        7	0	0	1
        8	1	0	1
        9	2	0	1
        10	0	1	1
        11	1	1	1
        12	2	1	1
        13	1	0	2
        14	2	0	2
        15	2	1	2
        16	1	1	2
        21	10	0	0
        22	11	0	0
        23	11	11	0
        24	0	11	0
        25	10	0	1
        26	11	0	1
        27	11	11	1
        28	0	11	1
    End Nodes
    Begin Elements	Element3D8N
        1	1	1	2	5	4	7	8	11	10
        2	1	2	3	6	5	8	9	12	11
        3	1	8	9	12	11	13	14	15	16
        4	1	21	22	23	24	25	26	27	28
    End Elements
    )input"));
    ModelPartIO(p_input).ReadModelPart(model_part);

    //call process to create vector of neighbours
    auto neigh_process = GenericFindElementalNeighboursProcess(model_part);
    neigh_process.Execute();

    //element 1 must have a single neighbour in third face
    const auto& elem1_neighs = model_part.Elements()[1].GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_EXPECT_EQ(elem1_neighs.size(),6); // 6 faces
    for(unsigned int i_face = 0; i_face<elem1_neighs.size(); i_face++){
        if(i_face!=2){
            KRATOS_EXPECT_EQ( elem1_neighs(i_face).get(),nullptr); //no neigh in these faces
        } else{
            KRATOS_EXPECT_NE( elem1_neighs(i_face).get(),nullptr); //real pointer in this face
            KRATOS_EXPECT_EQ( elem1_neighs[i_face].Id(),2); //element #2 is in this face
        }
    }

    //element 2 must have neighbours in two faces
    const auto& elem2_neighs = model_part.Elements()[2].GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_EXPECT_EQ(elem2_neighs.size(),6); // 6 faces
    for(unsigned int i_face = 0; i_face<elem2_neighs.size(); i_face++){
        if(i_face==4){
            KRATOS_EXPECT_NE( elem2_neighs(i_face).get(),nullptr); //real pointer in this face
            KRATOS_EXPECT_EQ( elem2_neighs[i_face].Id(),1); //element #1 is in this face
        } else if(i_face==5) {
            KRATOS_EXPECT_NE( elem2_neighs(i_face).get(),nullptr); //real pointer in this face
            KRATOS_EXPECT_EQ( elem2_neighs[i_face].Id(),3); //element #3 is in this face
        } else {
            KRATOS_EXPECT_EQ( elem2_neighs(i_face).get(),nullptr); //no neigh in these faces
        }
    }

    //element 3 must have a single neighbour in first face
    const auto& elem3_neighs = model_part.Elements()[3].GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_EXPECT_EQ(elem3_neighs.size(),6); // 6 faces
    for(unsigned int i_face = 0; i_face<elem3_neighs.size(); i_face++){
        if(i_face!=0){
            KRATOS_EXPECT_EQ( elem3_neighs(i_face).get(),nullptr); //no neigh in these faces
        } else{
            KRATOS_EXPECT_NE( elem3_neighs(i_face).get(),nullptr); //real pointer in this face
            KRATOS_EXPECT_EQ( elem3_neighs[i_face].Id(),2); //element #2 is in this face
        }
    }

    //element 4 must have no neighbours
    const auto& elem4_neighs = model_part.Elements()[4].GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_EXPECT_EQ(elem4_neighs.size(),6); // 6 faces
    for(unsigned int i_face = 0; i_face<elem4_neighs.size(); i_face++){
        KRATOS_EXPECT_EQ( elem4_neighs(i_face).get(),nullptr); //no neigh in these faces
    }

}

KRATOS_TEST_CASE_IN_SUITE(TetrahedraGenericFindElementsNeighbourProcessTest,
                          KratosCoreFastSuite)
{
    using namespace Kratos;

    Model model;
    ModelPart& model_part = model.CreateModelPart("test_model_part");
    model_part.AddNodalSolutionStepVariable(DISTANCE);

    // read mdpa
    //the model part is a cube formed by 6 tetras
    Kratos::shared_ptr<std::iostream> p_input(new std::stringstream(
    R"input(
	Begin Properties 1
    End Properties
    Begin Nodes
        1	0	0	0
        2	1	0	0
        3	0	1	0
        4	1	1	0
        5	0	0	1
        6	1	0	1
        7	0	1	1
        8	1	1	1
    End Nodes
    Begin Elements	Element3D4N
        1	1	2	6	8	1
        2	1	1	2	4	8
        3	1	5	7	8	1
        4	1	1	3	7	8
        5	1	1	5	6	8
        6	1	3	4	8	1
    End Elements
    )input"));
    ModelPartIO(p_input).ReadModelPart(model_part);

    //call process to create vector of neighbours
    auto neigh_process = GenericFindElementalNeighboursProcess(model_part);
    neigh_process.Execute();

    //checking single element, complete model part test was done for hexas
    const auto& elem1_neighs = model_part.Elements()[1].GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_EXPECT_EQ(elem1_neighs.size(),4); // 4 faces
    for(unsigned int i_face = 0; i_face<elem1_neighs.size(); i_face++){
        if(i_face==0){
            KRATOS_EXPECT_NE( elem1_neighs(i_face).get(),nullptr); //real pointer in this face
            KRATOS_EXPECT_EQ( elem1_neighs[i_face].Id(),5); //element #1 is in this face
        } else if(i_face==1) {
            KRATOS_EXPECT_NE( elem1_neighs(i_face).get(),nullptr); //real pointer in this face
            KRATOS_EXPECT_EQ( elem1_neighs[i_face].Id(),2); //element #3 is in this face
        } else {
            KRATOS_EXPECT_EQ( elem1_neighs(i_face).get(),nullptr); //no neigh in these faces
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrianglesQuadilateralsGenericFindElementsNeighbourProcessTest,
                          KratosCoreFastSuite)
{
    using namespace Kratos;

    Model model;
    ModelPart& model_part = model.CreateModelPart("test_model_part");
    model_part.AddNodalSolutionStepVariable(DISTANCE);

    // read mdpa
    //the model part has 3 connected hexas and one disconnected one
    Kratos::shared_ptr<std::iostream> p_input(new std::stringstream(
    R"input(
	Begin Properties 1
    End Properties
    Begin Nodes
        1	0	0	0
        2	1	0	0
        3	0	1	0
        4   1   1   1
        5   -0.5  -1  1
        6   0	0	2
        7	1	0	2
        8	0	1	2
        9  -1   0   0

    End Nodes
    Begin Elements	Element3D3N
        1	1	1	2   3
        2	1	3	2   4
        3   1   1   5   2
        4   1   6   7   8
    End Elements

    Begin Elements	Element2D4N
        5	1	5 1 3 9
    End Elements


    )input"));
    ModelPartIO(p_input).ReadModelPart(model_part);

    //call process to create vector of neighbours
    auto neigh_process = GenericFindElementalNeighboursProcess(model_part);
    neigh_process.Execute();

    //element 1 must have a two neighs
    const auto& elem1_neighs = model_part.Elements()[1].GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_EXPECT_EQ(elem1_neighs.size(),3); // 3 edges
    for(unsigned int i_face = 0; i_face<elem1_neighs.size(); i_face++){
        if(i_face==0){
            KRATOS_EXPECT_NE( elem1_neighs(i_face).get(),nullptr); //real pointer in this face
            KRATOS_EXPECT_EQ( elem1_neighs[i_face].Id(),2); //element #2 is in this face
        } else if(i_face==2) {
            KRATOS_EXPECT_NE( elem1_neighs(i_face).get(),nullptr); //real pointer in this face
            KRATOS_EXPECT_EQ( elem1_neighs[i_face].Id(),3); //element #3 is in this face
        } else {
            KRATOS_EXPECT_NE( elem1_neighs(i_face).get(),nullptr); //real pointer in this face
            KRATOS_EXPECT_EQ( elem1_neighs[i_face].Id(),5); //element #5 is in this face
        }
    }

    //element 2 must have a single neigh
    const auto& elem2_neighs = model_part.Elements()[2].GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_EXPECT_EQ(elem2_neighs.size(),3); // 3 edges
    for(unsigned int i_face = 0; i_face<elem2_neighs.size(); i_face++){
        if(i_face!=2){
            KRATOS_EXPECT_EQ( elem2_neighs(i_face).get(),nullptr); //no neigh in these faces
        } else{
            KRATOS_EXPECT_NE( elem2_neighs(i_face).get(),nullptr); //real pointer in this face
            KRATOS_EXPECT_EQ( elem2_neighs[i_face].Id(),1); //element #1 is in this face
        }
    }

    //element 3 must have two neighs
    const auto& elem3_neighs = model_part.Elements()[3].GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_EXPECT_EQ(elem3_neighs.size(),3); // 3 edges
    for(unsigned int i_face = 0; i_face<elem3_neighs.size(); i_face++){
        if(i_face==1){
            KRATOS_EXPECT_NE( elem3_neighs(i_face).get(),nullptr); //real pointer in this face
            KRATOS_EXPECT_EQ( elem3_neighs[i_face].Id(),1); //element #5 is in this face
        } else if(i_face==2) {
            KRATOS_EXPECT_NE( elem3_neighs(i_face).get(),nullptr); //real pointer in this face
            KRATOS_EXPECT_EQ( elem3_neighs[i_face].Id(),5); //element #5 is in this face
        } else {
            KRATOS_EXPECT_EQ( elem3_neighs(i_face).get(),nullptr); //no neigh in these faces
        }
    }

    //element 4 must have no neighbours
    const auto& elem4_neighs = model_part.Elements()[4].GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_EXPECT_EQ(elem4_neighs.size(),3); // 3 edges
    for(unsigned int i_face = 0; i_face<elem4_neighs.size(); i_face++){
        KRATOS_EXPECT_EQ( elem4_neighs(i_face).get(),nullptr); //no neigh in these faces
    }

    //element 5 is a quadrilateral that has first and second edges filled
    const auto& elem5_neighs = model_part.Elements()[5].GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_EXPECT_EQ(elem5_neighs.size(),4); // 4 edges
    for(unsigned int i_face = 0; i_face<elem5_neighs.size(); i_face++){
        if(i_face==0){
            KRATOS_EXPECT_NE( elem5_neighs(i_face).get(),nullptr); //real pointer in this face
            KRATOS_EXPECT_EQ( elem5_neighs[i_face].Id(),3); //element #3 is in this face
        } else if(i_face==1) {
            KRATOS_EXPECT_NE( elem5_neighs(i_face).get(),nullptr); //real pointer in this face
            KRATOS_EXPECT_EQ( elem5_neighs[i_face].Id(),1); //element #1 is in this face
        } else {
            KRATOS_EXPECT_EQ( elem5_neighs(i_face).get(),nullptr); //no neigh in these faces
        }
    }

}

} // namespace Kratos::Testing
