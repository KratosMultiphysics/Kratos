// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "containers/model.h"
#include "custom_elements/transient_thermal_element.h"
#include "geo_mechanics_application_variables.h"
#include "testing/testing.h"

using namespace Kratos;

namespace Kratos::Testing {

void GenerateTransientThermalElement2D3N(ModelPart& rModelPart)
{
    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(TEMPERATURE);

    // Set the element properties
    Properties::Pointer pElemProp = rModelPart.CreateNewProperties(0);

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3};
    rModelPart.CreateNewElement("TransientThermalElement2D3N", 1, elemNodes, pElemProp);
}

KRATOS_TEST_CASE_IN_SUITE(EquationIdVectorTransientThermalElement2D3N, KratosGeoMechanicsFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateTransientThermalElement2D3N(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    for (unsigned int i = 0; i < model_part.NumberOfNodes(); i++) {
        p_element->GetGeometry()[i].AddDof(TEMPERATURE);
    }

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
   // p_element->Check(r_current_process_info);

    Element::DofsVectorType ElementalDofList;
    p_element->GetDofList(ElementalDofList, r_current_process_info);

    for (int i = 0; i < ElementalDofList.size(); i++) {
        ElementalDofList[i]->SetEquationId(i);
    }

    Element::EquationIdVectorType EquationIdVector;
    p_element->EquationIdVector(EquationIdVector, r_current_process_info);

    // Check the EquationIdVector values
    for (unsigned int i = 0; i < EquationIdVector.size(); i++) {
        KRATOS_EXPECT_TRUE(EquationIdVector[i] == i);
    }
}

} // namespace Kratos::Testing
