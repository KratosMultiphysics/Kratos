//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

// System includes


// External includes


// Project includes
#include "custom_utilities/read_materials.hpp"
//#include "structural_mechanics_application_variables.h"

namespace Kratos
{

    ReadMaterialProcess::ReadMaterialProcess(ModelPart &rModelPart,
                                             Parameters parameters)
            : mrModelPart(rModelPart)
    {
        Parameters default_parameters(R"(
            {
                "materials_filename" : "please specify the file to be opened"
            }  )"
        );

        parameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        KRATOS_INFO("::[Reading materials process DEBUG]::") << "Started" << std::endl;

        // read json string in materials file, create Parameters
        std::string materials_filename = parameters["materials_filename"].GetString();
        std::ifstream infile(materials_filename);
        std::stringstream buffer;
        buffer << infile.rdbuf();
        Parameters materials(buffer.str());

        for (auto i = 0; i < materials["properties"].size(); ++i)
        {
            Parameters material = materials["properties"].GetArrayItem(i);
            AssignPropertyBlock(material);
        }

        KRATOS_INFO("::[Reading materials process DEBUG]::") << "Finished" << std::endl;
    }


    void ReadMaterialProcess::Execute()
    {
    }

    void ReadMaterialProcess::AssignPropertyBlock(Parameters data)
    {
        /* Set constitutive law and material properties and assign to elements and conditions.
        Arguments:
        data -- a dictionary or json object defining properties for a model part.

        Example:
        data = {
            "model_part_name" : "Plate",
            "properties_id" : 1,
            "Material" : {
                "constitutive_law" : {
                    "name" : "LinearElasticPlaneStress2DLaw"
                },
                "Variables" : {
                    "YOUNG_MODULUS" : 200e9,
                    "POISSON_RATIO" : 0.3,
                    "RESIDUAL_VECTOR" : [1.5,0.3,-2.58],
                    "LOCAL_INERTIA_TENSOR" : [[0.27,0.0],[0.0,0.27]]
                },
                "Tables" : {
                    "Table1" : {
                        "input_variable" : "TEMPERATURE",
                        "output_variable" : "YOUNG_MODULUS",
                        "data" : [
                            [0.0,  100.0],
                            [20.0, 90.0],
                            [30.0, 85.0],
                            [35.0, 80.0]
                        ]
                    }
                }
            }
        }
        */

        //KRATOS_WATCH(data);

        // Get the properties for the specified model part.
        IndexType property_id = data["properties_id"].GetInt();
        IndexType mesh_id = 0;
        Properties prop = mrModelPart.GetProperties(property_id, mesh_id);

        //TODO(marcelo): Complete the "keys()" part
        //if (len(data["Material"]["Variables"].keys()) > 0 and prop.HasVariables())
        if (prop.HasVariables())
            KRATOS_INFO("::[Reading materials process DEBUG]::") << "Property " << property_id << " already has variables." << std::endl;
        //if (len(data["Material"]["Tables"].keys()) > 0 && prop.HasTables())
        if (prop.HasTables())
            KRATOS_INFO("::[Reading materials process DEBUG]::") << "Property " << property_id << " already has tables." << std::endl;

        // Assign the properties to the model part's elements and conditions.
        /*
        for (auto elem = mrModelPart.ElementsBegin(); elem != mrModelPart.ElementsEnd(); elem++)
            //elem-> Propertie = prop
            //KRATOS_WATCH(elem->pGetProperties());
            elem->SetProperties(prop);

        for cond in model_part.Conditions:
            cond.Properties = prop

        Parameters mat = data["Material"];

        //Set the CONSTITUTIVE_LAW for the current properties.
        constitutive_law = self._GetConstitutiveLaw( mat["constitutive_law"]["name"].GetString() )

        prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, constitutive_law.Clone())

        # Add / override the values of material parameters in the properties
        for key, value in mat["Variables"].items():
            var = self._GetVariable(key)
            if value.IsDouble():
                prop.SetValue( var, value.GetDouble() )
            elif value.IsInt():
                prop.SetValue( var, value.GetInt() )
            elif value.IsBool():
                prop.SetValue( var, value.GetBool() )
            elif value.IsString():
                prop.SetValue( var, value.GetString() )
            elif value.IsMatrix():
                prop.SetValue( var, value.GetMatrix() )
            elif value.IsVector():
                prop.SetValue( var, value.GetVector() )
            else:
                raise ValueError("Type of value is not available")

        # Add / override tables in the properties
        for key, table in mat["Tables"].items():
            table_name = key

            input_var = self._GetVariable(table["input_variable"].GetString())
            output_var = self._GetVariable(table["output_variable"].GetString())

            new_table = KratosMultiphysics.PiecewiseLinearTable()

            for i in range(table["data"].size()):
                new_table.AddRow(table["data"][i][0].GetDouble(), table["data"][i][1].GetDouble())

            prop.SetTable(input_var,output_var,new_table)
        */
    }

    /*
    ConstitutiveLaw ReadMaterialProces::GetConstitutiveLaw(my_string)
    {
        //Return the Constitutive Law object named by the string argument.
        //
        //Example:
        //recommended usage:
        //constitutive_law = self._GetConstitutiveLaw('LinearElastic3DLaw')
        //deprecated:
        //constitutive_law = self._GetConstitutiveLaw('KratosMultiphysics.StructuralMechanicsApplication.LinearElastic3DLaw')
        //constitutive_law = self._GetConstitutiveLaw('StructuralMechanicsApplication.LinearElastic3DLaw')

        //return KratosGlobals<ConstitutiveLaw>("constitutive_law_name");
    }
    */


}  // namespace Kratos.
