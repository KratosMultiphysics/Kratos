//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/variables.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(SubpropertiesInterface, KratosCoreFastSuite) {
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    Properties::Pointer p1 = r_model_part.CreateNewProperties(1);
    Properties::Pointer p11 = r_model_part.CreateNewProperties(11);
    Properties::Pointer p12 = r_model_part.CreateNewProperties(12);
    Properties::Pointer p111 = r_model_part.CreateNewProperties(111);

    //adding existing subproperties
    p1->AddSubProperties(p11);
    p1->AddSubProperties(p12);
    p11->AddSubProperties(p111);

    //add a new one
    //note that i am creating a new property "1", belonging to 1.12 IT IS DIFFERENT FROM p1
    r_model_part.GetSubProperties("1.12").AddSubProperties(Kratos::make_shared<Properties>(1)); 
    


    p1->SetValue(YOUNG_MODULUS,1.0);
    p11->SetValue(YOUNG_MODULUS,11.0);
    KRATOS_CHECK_EQUAL(p1->GetSubProperties(11).GetValue(YOUNG_MODULUS), 11.0);
    KRATOS_CHECK_EQUAL(r_model_part.GetSubProperties("1.11").GetValue(YOUNG_MODULUS), 11.0);  //SUGGESTION TODO: this should be called GetProperties
    //r_model_part.GetSubProperties("1.11.111").SetValue(YOUNG_MODULUS,111.0); //ERROR TODO: this should work but is failing
    //KRATOS_CHECK_EQUAL(p111->GetValue(YOUNG_MODULUS), 11.0); //ERROR TODO: this should work but is failing

    //1.12.1 is different from 1, even though it has the same index
    //r_model_part.GetSubProperties("1.12.1").SetValue(YOUNG_MODULUS,12345.0);  //ERROR TODO: this should work but is failing
    //KRATOS_CHECK_EQUAL(r_model_part.GetSubProperties("1.12.1").GetValue(YOUNG_MODULUS), 12345.0); //ERROR TODO: this should work but is failing

    r_model_part.GetSubProperties("1.12").GetSubProperties(1).SetValue(YOUNG_MODULUS,12345.0);
    KRATOS_CHECK_EQUAL(r_model_part.GetSubProperties("1.12").GetSubProperties(1).GetValue(YOUNG_MODULUS), 12345.0);

    KRATOS_CHECK_EQUAL(p1->GetValue(YOUNG_MODULUS), 1.0);

}


}
} // namespace Kratos.
