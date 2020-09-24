//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:   BSD License
//      Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"
#include "utilities/auxiliar_model_part_utilities.h"

// Utilities
#include "utilities/cpp_tests_utilities.h"

namespace Kratos {
namespace Testing {

/**
* 1. Checks the correct work of the Auxiliar model parts utility GetData for scalar data on Node_historicalDatalocation 
*/
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_double_Node_historical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    
    // Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N");

    //Assign random values to x_direction displacement to each node, also save it in test_values to compare later
    int i=0;
    double disp_test_value = 1.26;    
    std::vector<double> test_values;
    test_values.resize(this_model_part.NumberOfNodes());

    for (auto i_node = this_model_part.NodesBegin(); i_node != this_model_part.NodesEnd(); i_node++) 
    {
		i_node->FastGetSolutionStepValue(DISPLACEMENT_X) = static_cast<double>(disp_test_value * i);
        test_values[i] = (disp_test_value * i);
        i++;
    }

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT_X, DataLocation::NodeHistorical);

    KRATOS_CHECK_VECTOR_NEAR(data, test_values, 0.004);

}

/**
* 2. Checks the correct work of the Auxiliar model parts utility GetData for scalar data on Node_NonHistorical Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_double_Node_Nonhistorical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    
    // Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N");

    //Assign random values to x_direction displacement to each node, also save it in test_values to compare later
    int i=0;
    double disp_test_value = 1.26;    
    std::vector<double> test_values;
    test_values.resize(this_model_part.NumberOfNodes());

    for (auto i_node = this_model_part.NodesBegin(); i_node != this_model_part.NodesEnd(); i_node++) 
    {
		i_node->GetValue(DISPLACEMENT_X) = static_cast<double>(disp_test_value * i);
        test_values[i] = (disp_test_value * i);
        i++;
    }

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT_X, DataLocation::NodeNonHistorical);

    KRATOS_CHECK_VECTOR_NEAR(data, test_values, 0.004);

}

/**
* 3. Checks the correct work of the Auxiliar model parts utility GetData for scalar data on Element Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_double_Element, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    
    // Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Assign random values to x_direction displacement to each element, also save it in test_values to compare later
    int i=0;
    double disp_test_value = 2.55;    
    std::vector<double> test_values;
    test_values.resize(this_model_part.NumberOfElements());

    for (auto i_elem = this_model_part.ElementsBegin(); i_elem != this_model_part.ElementsEnd(); i_elem++) 
    {
		i_elem->GetValue(DISPLACEMENT_X) = static_cast<double>(disp_test_value * i);
        test_values[i] = (disp_test_value * i);
        i++;
    }

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT_X, DataLocation::Element);

    KRATOS_CHECK_VECTOR_NEAR(data, test_values, 0.004);

}

/**
* 4. Checks the correct work of the Auxiliar model parts utility GetData for scalar data on Condition Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_double_Condition, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    
    // Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Assign random values to x_direction displacement to each condition, also save it in test_values to compare later
    int i=0;
    double disp_test_value = 2.55;    
    std::vector<double> test_values;
    test_values.resize(this_model_part.NumberOfConditions());

    for (auto i_cond = this_model_part.ConditionsBegin(); i_cond != this_model_part.ConditionsEnd(); i_cond++) 
    {
		i_cond->GetValue(DISPLACEMENT_X) = static_cast<double>(disp_test_value * i);
        test_values[i] = (disp_test_value * i);
        i++;
    }

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT_X, DataLocation::Condition);

    KRATOS_CHECK_VECTOR_NEAR(data, test_values, 0.004);

}

/**
* 5. Checks the correct work of the Auxiliar model parts utility GetData for scalar data on ModelPart Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_double_ModelPart, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    
    // Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Assign random values to x_direction displacement to each modelpart, also save it in test_values to compare later    
    std::vector<double> test_values;
    test_values.resize(1);
    test_values[0] = 3.55;

    this_model_part[DISPLACEMENT_X] = 3.55;

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT_X, DataLocation::ModelPart);

    KRATOS_CHECK_VECTOR_NEAR(data, test_values, 0.004);

}

/**
* 6. Checks the correct work of the Auxiliar model parts utility GetData for scalar data on ProcessInfo Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_double_ProcessInfo, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    
    // Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Assign random values to x_direction displacement to each ProcessInfo, also save it in test_values to compare later    
    std::vector<double> test_values;
    test_values.resize(1);
    test_values[0] = 3.55;

    this_model_part.GetProcessInfo()[DISPLACEMENT_X] = 3.55;

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT_X, DataLocation::ProcessInfo);

    KRATOS_CHECK_VECTOR_NEAR(data, test_values, 0.004);

}

/**
* 7. Checks the correct work of the Auxiliar model parts utility GetData for vector data on Node_historicalDatalocation 
*/
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_Vector_Node_historical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    
    // Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N");

    //Assign random values to x_direction displacement to each node, also save it in test_values to compare later
    int i=0; //To create random values
    double disp_test_value = 1.26;    
    std::vector<double> test_values;
    test_values.resize(this_model_part.NumberOfNodes() * 3); 

    for (auto i_node = this_model_part.NodesBegin(); i_node != this_model_part.NodesEnd(); i_node++) 
    {
        for(int j=0; j<3; j++)
        {
            i_node->FastGetSolutionStepValue(DISPLACEMENT)[j] = static_cast<double>(disp_test_value * i * j);
            test_values[i++] = (disp_test_value * i *j);
        }   
    }

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT, DataLocation::NodeHistorical);

    KRATOS_CHECK_VECTOR_NEAR(data, test_values, 0.004);

}

/**
* 8. Checks the correct work of the Auxiliar model parts utility GetData for Vector data on Node_NonHistorical Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_Vector_Node_Nonhistorical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    
    // Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N");

    //Assign random values to x_direction displacement to each node, also save it in test_values to compare later
    int i=0; //To create random values
    double disp_test_value = 1.26;    
    std::vector<double> test_values;
    test_values.resize(this_model_part.NumberOfNodes() * 3); 

    for (auto i_node = this_model_part.NodesBegin(); i_node != this_model_part.NodesEnd(); i_node++) 
    {
        for(int j=0; j<3; j++)
        {
            i_node->GetValue(DISPLACEMENT)[j] = static_cast<double>(disp_test_value * i * j);
            test_values[i++] = (disp_test_value * i *j);
        }   
    }

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT, DataLocation::NodeNonHistorical);

    KRATOS_CHECK_VECTOR_NEAR(data, test_values, 0.004);

}

/**
* 9. Checks the correct work of the Auxiliar model parts utility GetData for scalar data on Element Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_Vector_Element, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    
    // Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Assign random values to x_direction displacement to each element, also save it in test_values to compare later
    int i=0; //To create random values
    double disp_test_value = 2.55;    
    std::vector<double> test_values;
    test_values.resize(this_model_part.NumberOfElements() * 3); 

    for (auto i_elem = this_model_part.ElementsBegin(); i_elem != this_model_part.ElementsEnd(); i_elem++) 
    {
        for(int j=0; j<3; j++)
        {
            i_elem->GetValue(DISPLACEMENT)[j] = static_cast<double>(disp_test_value * i * j);
            test_values[i++] = (disp_test_value * i *j);
        }   
    }

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT, DataLocation::Element);

    KRATOS_CHECK_VECTOR_NEAR(data, test_values, 0.004);

}

/**
* 10. Checks the correct work of the Auxiliar model parts utility GetData for Vector data on Condition Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_Vector_Condition, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    
    // Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Assign random values to x_direction displacement to each condition, also save it in test_values to compare later
    int i=0;
    double disp_test_value = 2.55;    
    std::vector<double> test_values;
    test_values.resize(this_model_part.NumberOfConditions() * 3);

    for (auto i_cond = this_model_part.ConditionsBegin(); i_cond != this_model_part.ConditionsEnd(); i_cond++) 
    {
        for(int j=0; j<3; j++)
        {
		    i_cond->GetValue(DISPLACEMENT)[j] = static_cast<double>(disp_test_value * i * j);
            test_values[i++] = (disp_test_value * i *j);
        }
    }

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT, DataLocation::Condition);

    KRATOS_CHECK_VECTOR_NEAR(data, test_values, 0.004);

}

/**
* 11. Checks the correct work of the Auxiliar model parts utility GetData for Vector data on ModelPart Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_Vector_ModelPart, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    
    // Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Assign random values to x_direction displacement to each modelpart, also save it in test_values to compare later    
    std::vector<double> test_values;
    test_values.resize(3);

    for(int j=0; j<3; j++)
    {
        this_model_part[DISPLACEMENT][j] = 3.55;
        test_values[j] = 3.55;
	}

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT, DataLocation::ModelPart);

    KRATOS_CHECK_VECTOR_NEAR(data, test_values, 0.004);

}

/**
* 12. Checks the correct work of the Auxiliar model parts utility GetData for Vector data on ProcessInfo Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_GetData_Vector_ProcessInfo, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    
    // Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Assign random values to x_direction displacement to each ProcessInfo, also save it in test_values to compare later    
    std::vector<double> test_values;
    test_values.resize(3);

    for(int j=0; j<3; j++)
    {
        this_model_part.GetProcessInfo()[DISPLACEMENT][j] = 3.55;
        test_values[j] = 3.55;
	}

    auto data = AuxiliarModelPartUtilities(this_model_part).GetVariableData(DISPLACEMENT, DataLocation::ProcessInfo);

    KRATOS_CHECK_VECTOR_NEAR(data, test_values, 0.004);

}

/**
* 13. Checks the correct work of the Auxiliar model parts utility SetData for Scalar data on NodeHistorical Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetData_double_Node_historical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    
    //Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N");

    //Create some values to x_direction displacement to each node
    std::vector<double> rData;
    rData.resize(this_model_part.NumberOfNodes());
    double disp_test_value = 1.26;  

    for(int i=0; i<this_model_part.NumberOfNodes(); i++)
    {
        rData[i] = disp_test_value * i ;
    }
    
    //Run the Function SetVariable to Import a "rData" into the Model
    AuxiliarModelPartUtilities(this_model_part).SetVariableData(DISPLACEMENT_X, DataLocation::NodeHistorical, rData);

    std::vector<double> output_values;
    int i=0;
    for (auto i_node = this_model_part.NodesBegin(); i_node != this_model_part.NodesEnd(); i_node++) 
    {
		output_values[i++] = i_node->FastGetSolutionStepValue(DISPLACEMENT_X);
    }
    
    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 0.004);

}

/**
* 14. Checks the correct work of the Auxiliar model parts utility SetData for Scalar data on Node_NonHistorical Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetData_double_Node_Nonhistorical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    
    //Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N");

    //Create some values to x_direction displacement to each node
    std::vector<double> rData;
    rData.resize(this_model_part.NumberOfNodes());
    double disp_test_value = 1.26;  

    for(int i=0; i<this_model_part.NumberOfNodes(); i++)
    {
        rData[i] = disp_test_value * i ;
    }
    
    //Run the Function SetVariable to Import a "rData" into the Model
    AuxiliarModelPartUtilities(this_model_part).SetVariableData(DISPLACEMENT_X, DataLocation::NodeNonHistorical, rData);

    std::vector<double> output_values;
    int i=0;
    for (auto i_node = this_model_part.NodesBegin(); i_node != this_model_part.NodesEnd(); i_node++) 
    {
		output_values[i++] = i_node->GetValue(DISPLACEMENT_X);
    }
    
    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 0.004);

}

/**
* 15. Checks the correct work of the Auxiliar model parts utility SetData for Scalar data on Element Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetData_double_element, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    
    //Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Create some values to x_direction displacement to each element
    std::vector<double> rData;
    rData.resize(this_model_part.NumberOfElements());
    double disp_test_value = 2.55;  

    for(int i=0; i<this_model_part.NumberOfElements(); i++)
    {
        rData[i] = disp_test_value * i ;
    }
    
    //Run the Function SetVariable to Import a "rData" into the Model
    AuxiliarModelPartUtilities(this_model_part).SetVariableData(DISPLACEMENT_X, DataLocation::Element, rData);

    std::vector<double> output_values;
    int i=0;
    for (auto i_elem = this_model_part.ElementsBegin(); i_elem != this_model_part.ElementsEnd(); i_elem++) 
    {
		output_values[i++] = i_elem->GetValue(DISPLACEMENT_X);
    }
    
    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 0.004);

}

/**
* 16. Checks the correct work of the Auxiliar model parts utility SetData for Scalar data on Condition Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetData_double_condition, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    
    //Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Create some values to x_direction displacement to each condition
    std::vector<double> rData;
    rData.resize(this_model_part.NumberOfConditions());
    double disp_test_value = 2.55;  

    for(int i=0; i<this_model_part.NumberOfConditions(); i++)
    {
        rData[i] = disp_test_value * i ;
    }
    
    //Run the Function SetVariable to Import a "rData" into the Model
    AuxiliarModelPartUtilities(this_model_part).SetVariableData(DISPLACEMENT_X, DataLocation::Condition, rData);

    std::vector<double> output_values;
    int i=0;
    for (auto i_cond = this_model_part.ConditionsBegin(); i_cond != this_model_part.ConditionsEnd(); i_cond++) 
    {
		output_values[i++] = i_cond->GetValue(DISPLACEMENT_X);
    }
    
    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 0.004);

}

/**
* 17. Checks the correct work of the Auxiliar model parts utility SetData for Scalar data on ModelPart Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetData_double_ModelPart, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    
    //Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Create some values to x_direction displacement to each ModelPart
    std::vector<double> rData;
    rData.resize(1);
    rData[0] = 3.55;
    
    //Run the Function SetVariable to Import a "rData" into the Model
    AuxiliarModelPartUtilities(this_model_part).SetVariableData(DISPLACEMENT_X, DataLocation::ModelPart, rData);

    std::vector<double> output_values;
    output_values[0] = this_model_part[DISPLACEMENT_X];
    
    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 0.004);

}

/**
* 18. Checks the correct work of the Auxiliar model parts utility SetData for Scalar data on ProcessInfo Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetData_double_ProcessInfo, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    
    //Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Create some values to x_direction displacement to each process info
    std::vector<double> rData;
    rData.resize(1);
    rData[0] = 3.55;
    
    //Run the Function SetVariable to Import a "rData" into the Model
    AuxiliarModelPartUtilities(this_model_part).SetVariableData(DISPLACEMENT_X, DataLocation::ProcessInfo, rData);

    std::vector<double> output_values;
    output_values[0] = this_model_part.GetProcessInfo()[DISPLACEMENT_X];
    
    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 0.004);

}

/**
* 19. Checks the correct work of the Auxiliar model parts utility SetData for Vector data on NodeHistorical Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetData_Vector_Node_historical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    
    //Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N");

    //Create some random values to 3 dimensional displacement to each node
    std::vector<double> rData;
    rData.resize(this_model_part.NumberOfNodes() * 3);
    double disp_test_value = 1.26;  
    int counter=0;

    for(int i=0; i<this_model_part.NumberOfNodes(); i++)
    {
        for(int j=1; j<=3; j++)//for each node ther are 3 values
        {   
            rData[counter++] = disp_test_value * i * j ;
        }
    }
   
    //Run the Function SetVariable to Import a "rData" into the Model
    AuxiliarModelPartUtilities(this_model_part).SetVariableData(DISPLACEMENT, DataLocation::NodeHistorical, rData);

    std::vector<double> output_values;
    counter=0;

    for (auto i_node = this_model_part.NodesBegin(); i_node != this_model_part.NodesEnd(); i_node++) 
    {
        for(int j=0; j<3 ;j++)
        {
		    output_values[counter++] = i_node->FastGetSolutionStepValue(DISPLACEMENT)[j];
        }
    }
    
    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 0.004);

}

/**
* 20. Checks the correct work of the Auxiliar model parts utility SetData for Vector data on Node_NonHistorical Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetData_Vector_Node_Nonhistorical, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    
    //Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N");

    //Create some random values to 3 Dimensional displacement to each node
    std::vector<double> rData;
    rData.resize(this_model_part.NumberOfNodes() * 3);
    double disp_test_value = 1.26;  
    int counter=0;

    for(int i=0; i<this_model_part.NumberOfNodes(); i++)
    {
        for(int j=1; j<=3; j++)//for each node ther are 3 values
        {   
            rData[counter++] = disp_test_value * i * j ;
        }
    }
    
    //Run the Function SetVariable to Import a "rData" into the Model
    AuxiliarModelPartUtilities(this_model_part).SetVariableData(DISPLACEMENT, DataLocation::NodeNonHistorical, rData);

    std::vector<double> output_values;
    counter=0;

    for (auto i_node = this_model_part.NodesBegin(); i_node != this_model_part.NodesEnd(); i_node++) 
    {
        for(int j=0; j<3; j++)
        {
		    output_values[counter++] = i_node->GetValue(DISPLACEMENT)[j];
        }
    }
    
    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 0.004);

}

/**
* 21. Checks the correct work of the Auxiliar model parts utility SetData for Vector data on Element Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetData_Vector_element, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    
    //Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Create some values to x_direction displacement to each element
    std::vector<double> rData;
    rData.resize(this_model_part.NumberOfElements() * 3);
    double disp_test_value = 2.55;
    int counter=0;  

    for(int i=0; i<this_model_part.NumberOfElements(); i++)
    {
        for(int j=1; j<=3; j++)//for each element ther are 3 values
        {   
            rData[counter++] = disp_test_value * i * j ;
        }
    }
    
    //Run the Function SetVariable to Import a "rData" into the Model
    AuxiliarModelPartUtilities(this_model_part).SetVariableData(DISPLACEMENT, DataLocation::Element, rData);

    std::vector<double> output_values;
    counter=0;

    for (auto i_elem = this_model_part.ElementsBegin(); i_elem != this_model_part.ElementsEnd(); i_elem++) 
    {
        for(int j=0; j<3; j++)
        {
		    output_values[counter++] = i_elem->GetValue(DISPLACEMENT)[j];
        }
    }
    
    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 0.004);

}

/**
* 22. Checks the correct work of the Auxiliar model parts utility SetData for Vector data on Condition Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetData_Vector_condition, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    
    //Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Create some values to x_direction displacement to each condition
    std::vector<double> rData;
    rData.resize(this_model_part.NumberOfConditions() * 3);
    double disp_test_value = 2.55;  
    int counter=0;

    for(int i=0; i<this_model_part.NumberOfConditions(); i++)
    {
        for(int j=1; j<=3; j++)//for each condition ther are 3 values
        {   
            rData[counter++] = disp_test_value * i * j ;
        }
    }
    
    //Run the Function SetVariable to Import a "rData" into the Model
    AuxiliarModelPartUtilities(this_model_part).SetVariableData(DISPLACEMENT, DataLocation::Condition, rData);

    std::vector<double> output_values;
    counter=0;

    for (auto i_cond = this_model_part.ConditionsBegin(); i_cond != this_model_part.ConditionsEnd(); i_cond++) 
    {
        for(int j=0; j<3; j++)
        {
		    output_values[counter++] = i_cond->GetValue(DISPLACEMENT)[j];
        }
    }
    
    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 0.004);

}

/**
* 23. Checks the correct work of the Auxiliar model parts utility SetData for Vector data on ModelPart Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetData_Vector_ModelPart, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    
    //Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Create some values to x_direction displacement to each ModelPart
    std::vector<double> rData;
    rData.resize(3);
    int counter = 0;
    
    for(int j=1; j<=3; j++)//for one ModelPart ther are 3 values
    {   
        rData[counter++] = (3.55) * j ;
    }
    
    //Run the Function SetVariable to Import a "rData" into the Model
    AuxiliarModelPartUtilities(this_model_part).SetVariableData(DISPLACEMENT, DataLocation::ModelPart, rData);

    std::vector<double> output_values;
    
    for(int j=0; j<3; j++)
    {
		output_values[j] = this_model_part[DISPLACEMENT][j];
    }
      
    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 0.004);

}

/**
* 24. Checks the correct work of the Auxiliar model parts utility SetData for Vector data on ProcessInfo Datalocation
*/ 
KRATOS_TEST_CASE_IN_SUITE(AuxiliarModelPartUtilities_SetData_Vector_ProcessInfo, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    
    //Variable addition
    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);

    //Create some values to x_direction displacement to each process info
    std::vector<double> rData;
    rData.resize(3);
    int counter = 0;
    
    for(int j=1; j<=3; j++)//for one Process Info ther are 3 values
    {   
        rData[counter++] = (3.55) * j ;
    }
    
    //Run the Function SetVariable to Import a "rData" into the Model
    AuxiliarModelPartUtilities(this_model_part).SetVariableData(DISPLACEMENT, DataLocation::ProcessInfo, rData);

    std::vector<double> output_values;

    for(int j=0; j<3; j++)
    {
		output_values[j] = this_model_part.GetProcessInfo()[DISPLACEMENT][j];
    }
    
    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 0.004);

}

} // namespace Testing
} // namespace Kratos.
