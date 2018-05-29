//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//  			 Kratos default license: kratos/license.txt
//
//  Main authors:    Marcelo Raschi


// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/read_materials_utility.hpp"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtility, KratosCoreFastSuite)
{
    std::string properties_str = R"(
{
    "properties": [
        {
            "model_part_name": "MATERIAL_HOMOGENEOUS",
            "properties_id": 1,
            "Material": {
                "name": "homog",
                "constitutive_law": {
                    "name": "LinearElastic3DLaw"
                },
                "Variables": {
                    "DENSITY": 7850.0
                },
                "Tables": {
                    "Table1": {
                        "input_variable": "TEMPERATURE",
                        "output_variable": "YOUNG_MODULUS",
                        "data": [
                            [0, 100],
                            [1, 200],
                            [2, 300]
                        ]
                    }
                }
            }
        },
        {
            "model_part_name": "MATERIAL_MULTISCALE",
            "properties_id": 2,
            "Material": {
                "name": "homog",
                "constitutive_law": {
                    "name": "LinearElastic3DLaw"
                },
                "Variables": {},
                "Tables": {}
            }
        }
    ]
}
)";

    ModelPart model_part("dummy");
    ReadMaterialsUtility(model_part, properties_str);

    KRATOS_CHECK_EQUAL(model_part.NumberOfProperties(), 2);
    Properties& r_property = model_part.GetProperties(1, 0);
	KRATOS_CHECK(r_property.HasVariables());
    KRATOS_CHECK(r_property.HasTables());
    KRATOS_CHECK(r_property.Has(CONSTITUTIVE_LAW));
    KRATOS_CHECK(r_property.HasTable(TEMPERATURE, YOUNG_MODULUS));
    KRATOS_CHECK_EQUAL(r_property(DENSITY), 7850.0);
    KRATOS_CHECK_EQUAL(r_property(CONSTITUTIVE_LAW)->Info(), "ConstitutiveLaw");
}

} // namespace Testing.
} // namespace Kratos.
