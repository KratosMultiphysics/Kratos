//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes


// External includes


// Project includes
#include "custom_processes/print_after_formfinding_process.h"
#include "structural_mechanics_application_variables.h"
#include "includes/model_part_io.h"


namespace Kratos
{

    PrintAfterFormfindingProcess::PrintAfterFormfindingProcess(ModelPart &rModelPart,
                                                                 Parameters OutputParameters) 
                                                                 : mModelPart(rModelPart), 
                                                                   mOutputParameters(OutputParameters)
    {
        Parameters default_parameters(R"(
            {
                "python_module"   : "print_after_formfinding_process",
                "kratos_module"   : "KratosMultiphysics.StructuralMechanicsApplication",
                "help"                  : "This process writes the mesh resulting from the formfinding in a .mdpa-file",
                "process_name"          : "PrintAfterFormfindingProcess",
                "Parameters"            : {
                    "model_part_name"   : "Structure"
                }
            }  )"
        );

        mOutputParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    }




    void PrintAfterFormfindingProcess::Execute()
    {
        //// erase nodal data
        //mModelPart.GetNodalSolutionStepVariablesList().clear();
       // 
        //// erase elemental data
        //for( auto& ele: mModelPart.Elements()){
        //    const Variable<Matrix> variable = KratosComponents<Variable<Matrix>>::Get("MEMBRANE_PRESTRESS");
        //    if(ele.Has(variable)){
        //        const Matrix membrane_prestress(ele.GetValue(MEMBRANE_PRESTRESS));
        //        ele.Data().Clear();
        //        ele.SetValue(MEMBRANE_PRESTRESS, membrane_prestress);
        //    }
        //    else
        //        ele.Data().Clear();

        //}

        //// erase conditional data
        //for( auto& cond: mModelPart.Conditions()){
        //        cond.Data().Clear();

        //}

        //ModelPartIO model_part_io("formfinding_out", IO::WRITE);
        //model_part_io.WriteModelPart(mModelPart);
        //const TVariableType variable = KratosComponents<Matrix>::Get(MEMBRANE_PRESTRESS);
        //(*mpStream) << "Begin "<<rObjectName<<"alData "<<variable.Name()<<std::endl;
        //for(auto& object : rThisObjectContainer){
        //    if(object.Has(variable)){
        //        (*mpStream)<<object.Id()<<"\t"<<object.GetValue(variable)<<std::endl;
        //    }
        //}
        //(*mpStream)<<"End "<<rObjectName<<"alData\n"<<std::endl;
        
    }
}  // namespace Kratos.
