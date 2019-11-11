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
//                   Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/variables.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(SubpropertiesInterface, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    Properties::Pointer p1 = r_model_part.CreateNewProperties(1);
    Properties::Pointer p11 = r_model_part.CreateNewProperties(11);
    Properties::Pointer p12 = r_model_part.CreateNewProperties(12);
    Properties::Pointer p111 = r_model_part.CreateNewProperties(111);

    // Adding existing subproperties
    p1->AddSubProperties(p11);
    p1->AddSubProperties(p12);
    p11->AddSubProperties(p111);

    KRATOS_CHECK(p1->HasSubProperties(11));
    KRATOS_CHECK(p11->HasSubProperties(111));
    KRATOS_CHECK(r_model_part.HasProperties("1.11"));
    KRATOS_CHECK(r_model_part.HasProperties("1.11.111"));

    // Add a new one
    // Note that i am creating a new property "1", belonging to 1.12 IT IS DIFFERENT FROM p1
    r_model_part.GetProperties("1.12").AddSubProperties(Kratos::make_shared<Properties>(1));

    p1->SetValue(YOUNG_MODULUS,1.0);
    p11->SetValue(YOUNG_MODULUS,11.0);
    KRATOS_CHECK_EQUAL(p1->GetSubProperties(11).GetValue(YOUNG_MODULUS), 11.0);
    KRATOS_CHECK_EQUAL(r_model_part.GetProperties("1.11").GetValue(YOUNG_MODULUS), 11.0);
    
    r_model_part.GetProperties("1.11.111").SetValue(YOUNG_MODULUS,111.0);
    KRATOS_CHECK_EQUAL(p111->GetValue(YOUNG_MODULUS), 111.0);

    //1.12.1 is different from 1, even though it has the same index
    r_model_part.GetProperties("1.12.1").SetValue(YOUNG_MODULUS,12345.0);
    KRATOS_CHECK_EQUAL(r_model_part.GetProperties("1.12.1").GetValue(YOUNG_MODULUS), 12345.0);

    r_model_part.GetProperties("1.12").GetSubProperties(1).SetValue(YOUNG_MODULUS,12345.0);
    KRATOS_CHECK_EQUAL(r_model_part.GetProperties("1.12").GetSubProperties(1).GetValue(YOUNG_MODULUS), 12345.0);

    KRATOS_CHECK_EQUAL(p1->GetValue(YOUNG_MODULUS), 1.0);

    std::size_t found = 0;
    for(const auto& prop_it : p1->GetSubProperties())
        if(prop_it.Id() == 11 || prop_it.Id() == 12)
            found++;
        else
            KRATOS_ERROR << "the property with Id " << prop_it.Id() << " should not exist in theh first layer";

    KRATOS_CHECK_EQUAL(found,2);
        
}


}
} // namespace Kratos.
