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
#include "custom_processes/read_materials_process.hpp"
#include "structural_mechanics_application_variables.h"
#include "input_output/logger.h"
#include "includes/kratos_parameters.h"
#include "includes/checks.h"

namespace Kratos
{

    ReadMaterialProcess::ReadMaterialProcess(ModelPart &rModelPart,
                                             Parameters parameters)
            : mrModelPart(rModelPart),
              mParametersFilename(parameters)
    {
        Parameters default_parameters(R"(
            {
                "materials_filename" : "please specify the file to be opened"
            }  )"
        );

        mParametersFilename.RecursivelyValidateAndAssignDefaults(default_parameters);
    }


    void ReadMaterialProcess::Execute()
    {
        LoggerMessage logger("::[Reading materials process DEBUG]::");

        logger << "Started";

        Parameters materials = mParametersFilename.GetValue("materials_filename");
        //for(auto prop = materials.begin(); prop != materials.end(); prop++)
        //    AssignPropertyBlock(prop->GetValue() );
        logger << materials;

        logger << "Finished";
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

        /*
        // Get the properties for the specified model part.
        IndexType property_id = data["properties_id"].GetInt()
        IndexType mesh_id = 0;
        Properties prop = mrModelPart.GetProperties(property_id, mesh_id);

        if len(data["Material"]["Variables"].keys()) > 0 and prop.HasVariables():
            KratosMultiphysics.Logger.PrintInfo("::[Reading materials process]:: ", "Property", str(property_id), "already has variables." )
        if len(data["Material"]["Tables"].keys()) > 0 and prop.HasTables():
            KratosMultiphysics.Logger.PrintInfo("::[Reading materials process]:: ", "Property", str(property_id), "already has tables." )

        # Assign the properties to the model part's elements and conditions.
        for elem in model_part.Elements:
            elem.Properties = prop

        for cond in model_part.Conditions:
            cond.Properties = prop

        mat = data["Material"]

        # Set the CONSTITUTIVE_LAW for the current properties.
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



}  // namespace Kratos.
