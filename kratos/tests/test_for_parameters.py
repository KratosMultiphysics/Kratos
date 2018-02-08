import KratosMultiphysics

test_settings = KratosMultiphysics.Parameters("""
    {
            "Parameters": {
                    "materials_filename": "materials.json"
            }
    }
    """)

print(test_settings["Parameters"].IsArray())

test_settings = KratosMultiphysics.Parameters("""
    {
            "Parameters": {
                    "some_value": [1,2,3,4,5]
            }
    }
    """)

for value in test_settings["Parameters"]["some_value"]:
    print(value)

for i in range(test_settings["Parameters"]["some_value"].size()):
    print(test_settings["Parameters"]["some_value"][i])

for value in test_settings["Parameters"]["some_value"]:
    print(value.GetDouble())