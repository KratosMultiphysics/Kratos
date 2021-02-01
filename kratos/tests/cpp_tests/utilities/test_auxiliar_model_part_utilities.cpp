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

/******************************************************************************************/
/* Helper Functions */
/******************************************************************************************/

void Initialize(ModelPart& this_model_part)
{
    this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    CppTestsUtilities::Create2DGeometry(this_model_part, "Element2D3N", true, true);
}

std::vector<double> ComputationGetData(const DataLocation DataLoc, ModelPart& this_model_part, int dim)
{
    //Assign random values to Displacement of each node and also save it in test_values for comparision
    int i=0; //To create various values
    double disp_test_value = 1.26;
    std::vector<double> test_values;

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

std::vector<double> PreComputeSetData(int size, int dim)
{
    std::vector<double> rData;
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

std::vector<double> PostComputeSetData(const DataLocation DataLoc, ModelPart& this_model_part, int dim)
{
    std::vector<double> output_values;
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
    KRATOS_CHECK_VECTOR_NEAR(test_values, data, 1e-15);
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
    KRATOS_CHECK_VECTOR_NEAR(test_values, data, 1e-15);
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
    KRATOS_CHECK_VECTOR_NEAR(test_values, data, 1e-15);
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
    KRATOS_CHECK_VECTOR_NEAR(test_values, data, 1e-15);
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
    KRATOS_CHECK_VECTOR_NEAR(test_values, data, 1e-15);
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
    KRATOS_CHECK_VECTOR_NEAR(test_values, data, 1e-15);
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
    KRATOS_CHECK_VECTOR_NEAR(test_values, data, 1e-15);
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
    KRATOS_CHECK_VECTOR_NEAR(test_values, data, 1e-15);
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
    KRATOS_CHECK_VECTOR_NEAR(test_values, data, 1e-15);
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
    KRATOS_CHECK_VECTOR_NEAR(test_values, data, 1e-15);
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
    KRATOS_CHECK_VECTOR_NEAR(test_values, data, 1e-15);
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
    KRATOS_CHECK_VECTOR_NEAR(test_values, data, 1e-15);
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

    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 1e-15);
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

    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 1e-15);
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

    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 1e-15);
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

    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 1e-15);
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

    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 1e-15);
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

    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 1e-15);
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

    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 1e-15);
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

    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 1e-15);
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

    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 1e-15);
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

    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 1e-15);
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

    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 1e-15);
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

    KRATOS_CHECK_VECTOR_NEAR(rData, output_values, 1e-15);
}

} // namespace Testing
} // namespace Kratos.
