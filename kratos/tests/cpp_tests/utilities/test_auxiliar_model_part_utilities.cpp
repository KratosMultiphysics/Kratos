//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"
#include "geometries/triangle_3d_3.h"

// Utilities
#include "utilities/auxiliar_model_part_utilities.h"
#include "utilities/cpp_tests_utilities.h"

namespace Kratos::Testing {

using DataLocation = Globals::DataLocation;

/******************************************************************************************/
/* Helper Functions */
/******************************************************************************************/

void Initialize(ModelPart& this_model_part)
{
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);
}

template<class TContainerType=std::vector<double>>
auto ComputationGetData(const DataLocation DataLoc, ModelPart& this_model_part, int dim)
{
    //Assign random values to Displacement of each node and also save it in test_values for comparision
    int i=0; //To create various values
    double disp_test_value = 1.26;
    TContainerType test_values;

    switch (DataLoc)
        {
            case (DataLocation::NodeHistorical):{
                test_values.resize(this_model_part.NumberOfNodes() * dim);
                for (auto& i_node : this_model_part.Nodes())
                {
                    for(int j=0; j<dim; j++)
		            {
                        i_node.FastGetSolutionStepValue(DISPLACEMENT)[j] = (disp_test_value * (i+j));
                        test_values[i] = (disp_test_value * (i+j));
                        i++;
                    }
                }
                break;
            }
            case (DataLocation::NodeNonHistorical):{
                test_values.resize(this_model_part.NumberOfNodes() * dim);
                for (auto& i_node : this_model_part.Nodes())
                {
                    //array_1d<double, 3> test_vals; // can be done like this
		            for(int j=0; j<dim; j++)
		            {
                        //i_node.SetValue(DISPLACEMENT, test_vals);
                        i_node.GetValue(DISPLACEMENT)[j] = (disp_test_value * (i+j));
                        test_values[i] = (disp_test_value * (i+j));
                        i++;
                    }
                }
                break;
            }
            case (DataLocation::Element):{
                test_values.resize(this_model_part.NumberOfElements() * dim);
                for (auto& i_elem : this_model_part.Elements())
                {
                    for(int j=0; j<dim; j++)
		            {
                        i_elem.GetValue(DISPLACEMENT)[j] = (disp_test_value * (i+j));
                        test_values[i] = (disp_test_value * (i+j));
                        i++;
                    }
                }
                break;
            }
            case (DataLocation::Condition):{
                test_values.resize(this_model_part.NumberOfConditions() * dim);
                for (auto i_cond : this_model_part.Conditions())
                {
                    for(int j=0; j<dim; j++)
		            {
                        i_cond.GetValue(DISPLACEMENT)[j] = (disp_test_value * (i+j));
                        test_values[i] = (disp_test_value * (i+j));
                        i++;
                    }
                }
                break;
            }
            case (DataLocation::ModelPart):{
                test_values.resize(1 * dim);
                for(int j=0; j<dim; j++)
		            {
                        this_model_part[DISPLACEMENT][j] = (disp_test_value * (i+j));
                        test_values[i] = (disp_test_value * (i+j));
                        i++;
                    }
                break;
            }
            case (DataLocation::ProcessInfo):{
                test_values.resize(1 * dim);
                for(int j=0; j<dim; j++)
		            {
                        this_model_part.GetProcessInfo()[DISPLACEMENT][j] = (disp_test_value * (i+j));
                        test_values[i] = (disp_test_value * (i+j));
                        i++;
                    }
                break;
            }
            default:{
                //Throw an error about invalid DataLocation
                KRATOS_ERROR << "unknown Datalocation" << std::endl;
                break;
            }
        }
    return test_values;
}

template<class TContainerType=std::vector<double>>
auto PreComputeSetData(int size, int dim)
{
    TContainerType rData;
    rData.resize(size * dim);
    double disp_test_value = 2.55;
    int counter = 0 ;

    for(int i=0; i<size; i++)
    {
        for(int j=0; j<dim; j++)
        {
            rData[counter++] = disp_test_value * (i+1+j);  //Addition of 1 for the case where size=1 && dim=1, to have non zero random value
        }
    }

    return rData;
}

template<class TContainerType=std::vector<double>>
auto PostComputeSetData(const DataLocation DataLoc, ModelPart& this_model_part, int dim)
{
    TContainerType output_values;
    int counter=0;

    switch (DataLoc)
        {
            case (DataLocation::NodeHistorical):{
                output_values.resize(this_model_part.NumberOfNodes()* dim);
                for (auto& i_node : this_model_part.Nodes())
                {
                    for(int j=0; j<dim; j++)
		            {
                        output_values[counter++] = i_node.FastGetSolutionStepValue(DISPLACEMENT)[j];
                    }
                }
                break;
            }
            case (DataLocation::NodeNonHistorical):{
                output_values.resize(this_model_part.NumberOfNodes()* dim);
                for (auto& i_node : this_model_part.Nodes())
                {
                    //array_1d<double, 3> test_vals; // can be done like this
                    for(int j=0; j<dim; j++)
		            {
                        //i_node.SetValue(DISPLACEMENT, test_vals);
                        output_values[counter++] = i_node.GetValue(DISPLACEMENT)[j];
                    }
                }
                break;
            }
            case (DataLocation::Element):{
                output_values.resize(this_model_part.NumberOfElements()* dim);
                for (auto& i_elem : this_model_part.Elements())
                {
                    for(int j=0; j<dim; j++)
		            {
                        output_values[counter++] = i_elem.GetValue(DISPLACEMENT)[j];
                    }
                }
                break;
            }
            case (DataLocation::Condition):{
                output_values.resize(this_model_part.NumberOfConditions()* dim);
                for (auto& i_cond : this_model_part.Conditions())
                {
                    for(int j=0; j<dim; j++)
		            {
                        output_values[counter++] = i_cond.GetValue(DISPLACEMENT)[j];
                    }
                }
                break;
            }
            case (DataLocation::ModelPart):{
                output_values.resize(1 * dim);
                for(int j=0; j<dim; j++)
		            {
                        output_values[counter++] = this_model_part(DISPLACEMENT)[j];
                    }
                break;
            }
            case (DataLocation::ProcessInfo):{
                output_values.resize(1 * dim);
                for(int j=0; j<dim; j++)
		            {
                        output_values[counter++] = this_model_part.GetProcessInfo()(DISPLACEMENT)[j];
                    }
                break;
            }
            default:{
                //Throw an error about invalid DataLocation
                KRATOS_ERROR << "unknown Datalocation" << std::endl;
                break;
            }
        }
        return output_values;
}

KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_RemoveNodesFromSubModePartsWithoutCorrespondingEntities, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    Initialize(r_model_part);
    auto& r_sub = r_model_part.CreateSubModelPart("SubModelPart");
    r_sub.CreateSubModelPart("SubSubModel");

    r_sub.AddNode(r_model_part.pGetNode(1));
    r_sub.AddNode(r_model_part.pGetNode(2));

    KRATOS_EXPECT_EQ(r_sub.NumberOfNodes(), 2);

    auto utilities = AuxiliarModelPartUtilities(r_model_part);
    utilities.RemoveOrphanNodesFromSubModelParts();

    KRATOS_EXPECT_EQ(r_sub.NumberOfNodes(), 0);

    r_sub.AddNode(r_model_part.pGetNode(1));
    r_sub.AddNode(r_model_part.pGetNode(2));
    r_sub.AddNode(r_model_part.pGetNode(3));
    r_sub.AddNode(r_model_part.pGetNode(4));
    r_sub.AddElement(r_model_part.pGetElement(1));

    KRATOS_EXPECT_EQ(r_sub.NumberOfNodes(), 4);
    KRATOS_EXPECT_EQ(r_sub.NumberOfElements(), 1);
    KRATOS_EXPECT_EQ(r_sub.NumberOfGeometries(), 0);

    utilities.RemoveOrphanNodesFromSubModelParts();

    KRATOS_EXPECT_EQ(r_sub.NumberOfNodes(), 3);
    KRATOS_EXPECT_EQ(r_sub.NumberOfElements(), 1);
    KRATOS_EXPECT_EQ(r_sub.NumberOfGeometries(), 0);

    // Replace element by geometry
    r_sub.AddNode(r_model_part.pGetNode(4));
    r_sub.AddGeometry(r_model_part.pGetElement(1)->pGetGeometry());
    r_sub.RemoveElement(1);

    KRATOS_EXPECT_EQ(r_sub.NumberOfNodes(), 4);
    KRATOS_EXPECT_EQ(r_sub.NumberOfElements(), 0);
    KRATOS_EXPECT_EQ(r_sub.NumberOfGeometries(), 1);

    utilities.RemoveOrphanNodesFromSubModelParts();

    KRATOS_EXPECT_EQ(r_sub.NumberOfNodes(), 3);
    KRATOS_EXPECT_EQ(r_sub.NumberOfElements(), 0);
    KRATOS_EXPECT_EQ(r_sub.NumberOfGeometries(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_CopySubModelPartStructure, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    ModelPart& this_copy_model_part = current_model.CreateModelPart("MainCopied");
    Initialize(this_model_part);
    auto& r_sub = this_model_part.CreateSubModelPart("SubModel");
    r_sub.CreateSubModelPart("SubSubModel");
    AuxiliarModelPartUtilities::CopySubModelPartStructure(this_model_part, this_copy_model_part);

    KRATOS_EXPECT_EQ(this_model_part.HasSubModelPart("Pikachu,pika,pika,pi"), this_copy_model_part.HasSubModelPart("Pikachu,pika,pika,pi"));
    KRATOS_EXPECT_EQ(this_model_part.HasSubModelPart("SubModel"), this_copy_model_part.HasSubModelPart("SubModel"));
    auto& r_sub_copy = this_copy_model_part.GetSubModelPart("SubModel");
    KRATOS_EXPECT_EQ(r_sub.HasSubModelPart("SubSubModel"), r_sub_copy.HasSubModelPart("SubSubModel"));
}

KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_DeepCopyModelPart, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_origin_model_part = current_model.CreateModelPart("Main");
    auto& r_sub = r_origin_model_part.CreateSubModelPart("SubModel");
    r_sub.CreateSubModelPart("SubSubModel");

    // Set in the ProcessInfo a value
    auto& r_process_info = r_origin_model_part.GetProcessInfo();
    r_process_info[STEP] = 1;

    // Adding variables to the model part
    r_origin_model_part.SetBufferSize(2);
    r_origin_model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    Properties::Pointer p_prop = r_origin_model_part.CreateNewProperties(0);
    p_prop->SetValue(DENSITY, 1.0);

    using NodeType = Node;

    // First we create the nodes
    auto p_node_1 = r_origin_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.00);
    auto p_node_2 = r_origin_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.00);
    auto p_node_3 = r_origin_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.01);
    std::vector<NodeType::Pointer> nodes_0 = {p_node_3, p_node_2, p_node_1};

    // Set values
    p_node_1->Set(RIGID, true);
    p_node_1->SetValue(PRESSURE, 15.0);

    auto p_node_4 = r_origin_model_part.CreateNewNode(4, 0.0 , 0.0 , 0.01);
    auto p_node_5 = r_origin_model_part.CreateNewNode(5, 1.0 , 0.0 , 0.01);
    auto p_node_6 = r_origin_model_part.CreateNewNode(6, 0.0 , 1.0 , 0.02);
    std::vector<NodeType::Pointer> nodes_1 = { p_node_4, p_node_5, p_node_6};

    // Set temperature to the nodes
    for (auto& r_node : r_origin_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) = static_cast<double>(r_node.Id());
    }

    // Now we create the "geometries"
    Triangle3D3<NodeType> triangle_0( PointerVector<NodeType>{nodes_0} );
    Triangle3D3<NodeType> triangle_1( PointerVector<NodeType>{nodes_1} );
    auto p_geom_1 = r_origin_model_part.CreateNewGeometry("Triangle3D3", 1, triangle_0);
    auto p_geom_2 = r_origin_model_part.CreateNewGeometry("Triangle3D3", 2, triangle_1);

    // Now we create the "elements"
    auto p_elem_1 = r_origin_model_part.CreateNewElement("Element3D3N", 1, triangle_0, p_prop);
    auto p_elem_2 = r_origin_model_part.CreateNewElement("Element3D3N", 2, triangle_1, p_prop);

    // Set the variables to the elements
    p_elem_1->SetValue(TEMPERATURE, 1.0);

    // Adding elements to the submodelpart
    std::vector<IndexType> element_ids = {2};
    r_sub.AddElements(element_ids);

    // Now we create the "conditions"
    auto p_cond_1 = r_origin_model_part.CreateNewCondition("SurfaceCondition3D3N", 1, triangle_0, p_prop);
    auto p_cond_2 = r_origin_model_part.CreateNewCondition("SurfaceCondition3D3N", 2, triangle_1, p_prop);

    ModelPart& r_copy_model_part = AuxiliarModelPartUtilities(r_origin_model_part).DeepCopyModelPart("MainCopied");

    // Check the structure of the copied model part
    KRATOS_EXPECT_TRUE(r_copy_model_part.HasSubModelPart("SubModel"));
    auto& r_sub_copy = r_copy_model_part.GetSubModelPart("SubModel");
    KRATOS_EXPECT_TRUE(r_sub_copy.HasSubModelPart("SubSubModel"));
    KRATOS_EXPECT_EQ(r_sub_copy.NumberOfNodes(), 0);
    KRATOS_EXPECT_EQ(r_sub_copy.NumberOfGeometries(), 0);
    KRATOS_EXPECT_EQ(r_sub_copy.NumberOfElements(), 1);
    KRATOS_EXPECT_EQ(r_sub_copy.Elements().begin()->Id(), 2);
    KRATOS_EXPECT_EQ(r_sub_copy.NumberOfConditions(), 0);

    // Verify it is the same pointer
    KRATOS_EXPECT_EQ(p_prop.get(), r_origin_model_part.pGetProperties(0).get());
    KRATOS_EXPECT_EQ(p_node_1.get(), r_origin_model_part.pGetNode(1).get());
    KRATOS_EXPECT_EQ(p_node_2.get(), r_origin_model_part.pGetNode(2).get());
    KRATOS_EXPECT_EQ(p_node_3.get(), r_origin_model_part.pGetNode(3).get());
    KRATOS_EXPECT_EQ(p_node_4.get(), r_origin_model_part.pGetNode(4).get());
    KRATOS_EXPECT_EQ(p_node_5.get(), r_origin_model_part.pGetNode(5).get());
    KRATOS_EXPECT_EQ(p_node_6.get(), r_origin_model_part.pGetNode(6).get());

    KRATOS_EXPECT_EQ(p_geom_1.get(), r_origin_model_part.pGetGeometry(1).get());
    KRATOS_EXPECT_EQ(p_geom_2.get(), r_origin_model_part.pGetGeometry(2).get());

    KRATOS_EXPECT_EQ(p_elem_1.get(), r_origin_model_part.pGetElement(1).get());
    KRATOS_EXPECT_EQ(p_elem_2.get(), r_origin_model_part.pGetElement(2).get());

    KRATOS_EXPECT_EQ(p_cond_1.get(), r_origin_model_part.pGetCondition(1).get());
    KRATOS_EXPECT_EQ(p_cond_2.get(), r_origin_model_part.pGetCondition(2).get());

    // Check it is a different pointer
    KRATOS_EXPECT_NE(p_prop.get(), r_copy_model_part.pGetProperties(0).get());
    KRATOS_EXPECT_NE(p_node_1.get(), r_copy_model_part.pGetNode(1).get());
    KRATOS_EXPECT_NE(p_node_2.get(), r_copy_model_part.pGetNode(2).get());
    KRATOS_EXPECT_NE(p_node_3.get(), r_copy_model_part.pGetNode(3).get());
    KRATOS_EXPECT_NE(p_node_4.get(), r_copy_model_part.pGetNode(4).get());
    KRATOS_EXPECT_NE(p_node_5.get(), r_copy_model_part.pGetNode(5).get());
    KRATOS_EXPECT_NE(p_node_6.get(), r_copy_model_part.pGetNode(6).get());

    KRATOS_EXPECT_NE(p_geom_1.get(), r_copy_model_part.pGetGeometry(1).get());
    KRATOS_EXPECT_NE(p_geom_2.get(), r_copy_model_part.pGetGeometry(2).get());

    KRATOS_EXPECT_NE(p_elem_1.get(), r_copy_model_part.pGetElement(1).get());
    KRATOS_EXPECT_NE(p_elem_2.get(), r_copy_model_part.pGetElement(2).get());

    KRATOS_EXPECT_NE(p_cond_1.get(), r_copy_model_part.pGetCondition(1).get());
    KRATOS_EXPECT_NE(p_cond_2.get(), r_copy_model_part.pGetCondition(2).get());

    // Verify values set
    auto& r_new_properties = r_copy_model_part.GetProperties(0);
    KRATOS_EXPECT_TRUE(p_prop->Has(DENSITY));
    KRATOS_EXPECT_TRUE(r_new_properties.Has(DENSITY));
    KRATOS_EXPECT_EQ(r_new_properties.GetValue(DENSITY), p_prop->GetValue(DENSITY));
    auto& r_copy_process_info = r_copy_model_part.GetProcessInfo();
    KRATOS_EXPECT_TRUE(r_process_info.Has(STEP));
    KRATOS_EXPECT_TRUE(r_copy_process_info.Has(STEP));
    KRATOS_EXPECT_EQ(r_copy_process_info.GetValue(STEP),r_process_info.GetValue(STEP));
    KRATOS_EXPECT_FALSE(r_copy_process_info.Has(NL_ITERATION_NUMBER));
    KRATOS_EXPECT_EQ(p_node_1->Is(RIGID), r_copy_model_part.pGetNode(1)->Is(RIGID));
    KRATOS_EXPECT_EQ(p_node_1->Is(ACTIVE), r_copy_model_part.pGetNode(1)->Is(ACTIVE));
    KRATOS_EXPECT_TRUE(r_copy_model_part.pGetNode(1)->Has(PRESSURE));
    KRATOS_EXPECT_EQ(p_node_1->GetValue(PRESSURE), r_copy_model_part.pGetNode(1)->GetValue(PRESSURE));
    KRATOS_EXPECT_FALSE(r_copy_model_part.pGetNode(1)->Has(TEMPERATURE));
    for (auto& r_node : r_copy_model_part.Nodes()) {
        KRATOS_EXPECT_DOUBLE_EQ(r_node.FastGetSolutionStepValue(TEMPERATURE), static_cast<double>(r_node.Id()));
    }

    KRATOS_EXPECT_DOUBLE_EQ(p_elem_1->GetValue(TEMPERATURE), r_copy_model_part.pGetElement(1)->GetValue(TEMPERATURE));
    KRATOS_EXPECT_FALSE(r_copy_model_part.pGetElement(2)->Has(TEMPERATURE));
}

/******************************************************************************************/
/* Testing for GetData and SetData for 6 different DataLocations */
/******************************************************************************************/

//1. Checks the correct work of the Auxiliar model parts utility GetData for scalar data on Node_historicalDatalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetScalarData_Node_historical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);
    std::vector<double> data;

    auto test_values = ComputationGetData(DataLocation::NodeHistorical, this_model_part, 1);
    AuxiliarModelPartUtilities(this_model_part).GetScalarData(DISPLACEMENT_X, DataLocation::NodeHistorical, data);
    KRATOS_EXPECT_VECTOR_NEAR(test_values, data, 1e-15);

}
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetScalarData_Node_historical_ublas_vector, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);
    Vector data;

    auto test_values = ComputationGetData<Vector>(DataLocation::NodeHistorical, this_model_part, 1);
    AuxiliarModelPartUtilities(this_model_part).GetScalarData(DISPLACEMENT_X, DataLocation::NodeHistorical, data);
    KRATOS_EXPECT_VECTOR_NEAR(test_values, data, 1e-15);
}

//2. Checks the correct work of the Auxiliar model parts utility GetData for scalar data on Node_NonHistorical Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetScalarData_Node_Nonhistorical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);
    std::vector<double> data;

    auto test_values = ComputationGetData(DataLocation::NodeNonHistorical, this_model_part, 1);
    AuxiliarModelPartUtilities(this_model_part).GetScalarData(DISPLACEMENT_X, DataLocation::NodeNonHistorical, data);
    KRATOS_EXPECT_VECTOR_NEAR(test_values, data, 1e-15);
}

//3. Checks the correct work of the Auxiliar model parts utility GetData for scalar data on Element Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetScalarData_Element, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);
    std::vector<double> data;

    auto test_values = ComputationGetData(DataLocation::Element, this_model_part, 1);
    AuxiliarModelPartUtilities(this_model_part).GetScalarData(DISPLACEMENT_X, DataLocation::Element, data);
    KRATOS_EXPECT_VECTOR_NEAR(test_values, data, 1e-15);
}

//4. Checks the correct work of the Auxiliar model parts utility GetData for scalar data on Condition Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetScalarData_Condition, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);
    std::vector<double> data;

    auto test_values = ComputationGetData(DataLocation::Condition, this_model_part, 1);
    AuxiliarModelPartUtilities(this_model_part).GetScalarData(DISPLACEMENT_X, DataLocation::Condition, data);
    KRATOS_EXPECT_VECTOR_NEAR(test_values, data, 1e-15);
}

//5. Checks the correct work of the Auxiliar model parts utility GetData for scalar data on ModelPart Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetScalarData_ModelPart, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);
    std::vector<double> data;

    auto test_values = ComputationGetData(DataLocation::ModelPart, this_model_part, 1);
    AuxiliarModelPartUtilities(this_model_part).GetScalarData(DISPLACEMENT_X, DataLocation::ModelPart, data);
    KRATOS_EXPECT_VECTOR_NEAR(test_values, data, 1e-15);
}

//6. Checks the correct work of the Auxiliar model parts utility GetData for scalar data on ProcessInfo Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetScalarData_ProcessInfo, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);
    std::vector<double> data;

    auto test_values = ComputationGetData(DataLocation::ProcessInfo, this_model_part, 1);
    AuxiliarModelPartUtilities(this_model_part).GetScalarData(DISPLACEMENT_X, DataLocation::ProcessInfo, data);
    KRATOS_EXPECT_VECTOR_NEAR(test_values, data, 1e-15);
}

//7. Checks the correct work of the Auxiliar model parts utility GetData for vector data on Node_historicalDatalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetVectorData_Node_historical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);
    std::vector<double> data;

    auto test_values = ComputationGetData(DataLocation::NodeHistorical, this_model_part, 3);
    AuxiliarModelPartUtilities(this_model_part).GetVectorData(DISPLACEMENT, DataLocation::NodeHistorical, data);
    KRATOS_EXPECT_VECTOR_NEAR(test_values, data, 1e-15);
}

KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetVectorData_Node_historical_ublas_vector, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);
    Vector data;

    auto test_values = ComputationGetData<Vector>(DataLocation::NodeHistorical, this_model_part, 3);
    AuxiliarModelPartUtilities(this_model_part).GetVectorData(DISPLACEMENT, DataLocation::NodeHistorical, data);
    KRATOS_EXPECT_VECTOR_NEAR(test_values, data, 1e-15);
}

//8. Checks the correct work of the Auxiliar model parts utility GetData for Vector data on Node_NonHistorical Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetVectorData_Node_Nonhistorical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);
    std::vector<double> data;

    auto test_values = ComputationGetData(DataLocation::NodeNonHistorical, this_model_part, 3);
    AuxiliarModelPartUtilities(this_model_part).GetVectorData(DISPLACEMENT, DataLocation::NodeNonHistorical, data);
    KRATOS_EXPECT_VECTOR_NEAR(test_values, data, 1e-15);
}

//9. Checks the correct work of the Auxiliar model parts utility GetData for scalar data on Element Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetVectorData_Element, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);
    std::vector<double> data;

    auto test_values = ComputationGetData(DataLocation::Element, this_model_part, 3);
    AuxiliarModelPartUtilities(this_model_part).GetVectorData(DISPLACEMENT, DataLocation::Element, data);
    KRATOS_EXPECT_VECTOR_NEAR(test_values, data, 1e-15);
}

//10. Checks the correct work of the Auxiliar model parts utility GetData for Vector data on Condition Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetVectorData_Condition, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);
    std::vector<double> data;

    auto test_values = ComputationGetData(DataLocation::Condition, this_model_part, 3);
    AuxiliarModelPartUtilities(this_model_part).GetVectorData(DISPLACEMENT, DataLocation::Condition, data);
    KRATOS_EXPECT_VECTOR_NEAR(test_values, data, 1e-15);
}

//11. Checks the correct work of the Auxiliar model parts utility GetData for Vector data on ModelPart Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetVectorData_ModelPart, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);
    std::vector<double> data;

    auto test_values = ComputationGetData(DataLocation::ModelPart, this_model_part, 3);
    AuxiliarModelPartUtilities(this_model_part).GetVectorData(DISPLACEMENT, DataLocation::ModelPart, data);
    KRATOS_EXPECT_VECTOR_NEAR(test_values, data, 1e-15);
}

//12. Checks the correct work of the Auxiliar model parts utility GetData for Vector data on ProcessInfo Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetVectorData_ProcessInfo, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);
    std::vector<double> data;

    auto test_values = ComputationGetData(DataLocation::ProcessInfo, this_model_part, 3);
    AuxiliarModelPartUtilities(this_model_part).GetVectorData(DISPLACEMENT, DataLocation::ProcessInfo, data);
    KRATOS_EXPECT_VECTOR_NEAR(test_values, data, 1e-15);
}

//13. Checks the correct work of the Auxiliar model parts utility SetData for Scalar data on NodeHistorical Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetScalarData_Node_historical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);

    auto rData = PreComputeSetData(this_model_part.NumberOfNodes(), 1); //To create an input Data "rData" to feed to SetScalarData()
    AuxiliarModelPartUtilities(this_model_part).SetScalarData(DISPLACEMENT_X, DataLocation::NodeHistorical, rData); //Run the Function SetVariable to Import a "rData" into the Model
    auto output_values = PostComputeSetData(DataLocation::NodeHistorical, this_model_part, 1);

    KRATOS_EXPECT_VECTOR_NEAR(rData, output_values, 1e-15);
}

KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetScalarData_Node_historical_ublas_vector, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);

    Vector rData = PreComputeSetData<Vector>(this_model_part.NumberOfNodes(), 1); //To create an input Data "rData" to feed to SetScalarData()
    AuxiliarModelPartUtilities(this_model_part).SetScalarData(DISPLACEMENT_X, DataLocation::NodeHistorical, rData); //Run the Function SetVariable to Import a "rData" into the Model
    auto output_values = PostComputeSetData(DataLocation::NodeHistorical, this_model_part, 1);

    KRATOS_EXPECT_VECTOR_NEAR(rData, output_values, 1e-15);
}

//14. Checks the correct work of the Auxiliar model parts utility SetData for Scalar data on Node_NonHistorical Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetScalarData_Node_Nonhistorical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);

    auto rData = PreComputeSetData(this_model_part.NumberOfNodes(), 1); //To create an input Data "rData" to feed to SetScalarData()
    AuxiliarModelPartUtilities(this_model_part).SetScalarData(DISPLACEMENT_X, DataLocation::NodeNonHistorical, rData); //Run the Function SetVariable to Import a "rData" into the Model
    auto output_values = PostComputeSetData(DataLocation::NodeNonHistorical, this_model_part, 1);

    KRATOS_EXPECT_VECTOR_NEAR(rData, output_values, 1e-15);
}

//15. Checks the correct work of the Auxiliar model parts utility SetData for Scalar data on Element Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetScalarData_Element, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);

    auto rData = PreComputeSetData(this_model_part.NumberOfElements(), 1); //To create an input Data "rData" to feed to SetScalarData()
    AuxiliarModelPartUtilities(this_model_part).SetScalarData(DISPLACEMENT_X, DataLocation::Element, rData); //Run the Function SetVariable to Import a "rData" into the Model
    auto output_values = PostComputeSetData(DataLocation::Element, this_model_part, 1);

    KRATOS_EXPECT_VECTOR_NEAR(rData, output_values, 1e-15);
}

//16. Checks the correct work of the Auxiliar model parts utility SetData for Scalar data on Condition Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetScalarData_Condition, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);

    auto rData = PreComputeSetData(this_model_part.NumberOfConditions(), 1); //To create an input Data "rData" to feed to SetScalarData()
    AuxiliarModelPartUtilities(this_model_part).SetScalarData(DISPLACEMENT_X, DataLocation::Condition, rData); //Run the Function SetVariable to Import a "rData" into the Model
    auto output_values = PostComputeSetData(DataLocation::Condition, this_model_part, 1);

    KRATOS_EXPECT_VECTOR_NEAR(rData, output_values, 1e-15);
}

//17. Checks the correct work of the Auxiliar model parts utility SetData for Scalar data on ModelPart Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetScalarData_ModelPart, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);

    auto rData = PreComputeSetData(1, 1); //To create an input Data "rData" to feed to SetScalarData()
    AuxiliarModelPartUtilities(this_model_part).SetScalarData(DISPLACEMENT_X, DataLocation::ModelPart, rData); //Run the Function SetVariable to Import a "rData" into the Model
    auto output_values = PostComputeSetData(DataLocation::ModelPart, this_model_part, 1);

    KRATOS_EXPECT_VECTOR_NEAR(rData, output_values, 1e-15);
}

//18. Checks the correct work of the Auxiliar model parts utility SetData for Scalar data on ProcessInfo Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetScalarData_ProcessInfo, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);

    auto rData = PreComputeSetData(1, 1); //To create an input Data "rData" to feed to SetScalarData()
    AuxiliarModelPartUtilities(this_model_part).SetScalarData(DISPLACEMENT_X, DataLocation::ProcessInfo, rData); //Run the Function SetVariable to Import a "rData" into the Model
    auto output_values = PostComputeSetData(DataLocation::ProcessInfo, this_model_part, 1);

    KRATOS_EXPECT_VECTOR_NEAR(rData, output_values, 1e-15);
}

//19. Checks the correct work of the Auxiliar model parts utility SetData for Vector data on NodeHistorical Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetVectorData_Node_historical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);

    auto rData = PreComputeSetData(this_model_part.NumberOfNodes(), 3); //To create an input Data "rData" to feed to SetVectorData()
    AuxiliarModelPartUtilities(this_model_part).SetVectorData(DISPLACEMENT, DataLocation::NodeHistorical, rData); //Run the Function SetVariable to Import a "rData" into the Model
    auto output_values = PostComputeSetData(DataLocation::NodeHistorical, this_model_part, 3);

    KRATOS_EXPECT_VECTOR_NEAR(rData, output_values, 1e-15);
}

KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetVectorData_Node_historical_ublas_vector, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);

    Vector rData = PreComputeSetData<Vector>(this_model_part.NumberOfNodes(), 3); //To create an input Data "rData" to feed to SetVectorData()
    AuxiliarModelPartUtilities(this_model_part).SetVectorData(DISPLACEMENT, DataLocation::NodeHistorical, rData); //Run the Function SetVariable to Import a "rData" into the Model
    auto output_values = PostComputeSetData(DataLocation::NodeHistorical, this_model_part, 3);

    KRATOS_EXPECT_VECTOR_NEAR(rData, output_values, 1e-15);
}

//20. Checks the correct work of the Auxiliar model parts utility SetData for Vector data on Node_NonHistorical Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetVectorData_Node_Nonhistorical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);

    auto rData = PreComputeSetData(this_model_part.NumberOfNodes(), 3); //To create an input Data "rData" to feed to SetVectorData()
    AuxiliarModelPartUtilities(this_model_part).SetVectorData(DISPLACEMENT, DataLocation::NodeNonHistorical, rData); //Run the Function SetVariable to Import a "rData" into the Model
    auto output_values = PostComputeSetData(DataLocation::NodeNonHistorical, this_model_part, 3);

    KRATOS_EXPECT_VECTOR_NEAR(rData, output_values, 1e-15);
}

//21. Checks the correct work of the Auxiliar model parts utility SetData for Vector data on Element Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetVectorData_Element, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);

    auto rData = PreComputeSetData(this_model_part.NumberOfElements(), 3); //To create an input Data "rData" to feed to SetVectorData()
    AuxiliarModelPartUtilities(this_model_part).SetVectorData(DISPLACEMENT, DataLocation::Element, rData); //Run the Function SetVariable to Import a "rData" into the Model
    auto output_values = PostComputeSetData(DataLocation::Element, this_model_part, 3);

    KRATOS_EXPECT_VECTOR_NEAR(rData, output_values, 1e-15);
}

//22. Checks the correct work of the Auxiliar model parts utility SetData for Vector data on Condition Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetVectorData_Condition, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);

    auto rData = PreComputeSetData(this_model_part.NumberOfConditions(), 3); //To create an input Data "rData" to feed to SetVectorData()
    AuxiliarModelPartUtilities(this_model_part).SetVectorData(DISPLACEMENT, DataLocation::Condition, rData); //Run the Function SetVariable to Import a "rData" into the Model
    auto output_values = PostComputeSetData(DataLocation::Condition, this_model_part, 3);

    KRATOS_EXPECT_VECTOR_NEAR(rData, output_values, 1e-15);
}

//23. Checks the correct work of the Auxiliar model parts utility SetData for Vector data on ModelPart Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetVectorData_ModelPart, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);

    auto rData = PreComputeSetData(1, 3); //To create an input Data "rData" to feed to SetVectorData()
    AuxiliarModelPartUtilities(this_model_part).SetVectorData(DISPLACEMENT, DataLocation::ModelPart, rData); //Run the Function SetVariable to Import a "rData" into the Model
    auto output_values = PostComputeSetData(DataLocation::ModelPart, this_model_part, 3);

    KRATOS_EXPECT_VECTOR_NEAR(rData, output_values, 1e-15);
}

//24. Checks the correct work of the Auxiliar model parts utility SetData for Vector data on ProcessInfo Datalocation
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetVectorData_ProcessInfo, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    Initialize(this_model_part);

    auto rData = PreComputeSetData(1, 3); //To create an input Data "rData" to feed to SetVectorData()
    AuxiliarModelPartUtilities(this_model_part).SetVectorData(DISPLACEMENT, DataLocation::ProcessInfo, rData); //Run the Function SetVariable to Import a "rData" into the Model
    auto output_values = PostComputeSetData(DataLocation::ProcessInfo, this_model_part, 3);

    KRATOS_EXPECT_VECTOR_NEAR(rData, output_values, 1e-15);
}

} // namespace Kratos::Testing.
