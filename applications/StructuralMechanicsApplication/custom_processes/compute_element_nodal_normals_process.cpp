// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//

// System includes

// External includes

// Project includes
#include "custom_processes/compute_element_nodal_normals_process.h"      // include the header of the process
#include "containers/model.h"                                            // include Model and ModelPart definitions to use them in the Create() method. needed for rModel.GetModelPart()
#include "utilities/normal_calculation_utils.h"                          // include utility to compute normals
#include "includes/variables.h"                                          // include core variable "NORMAL" 

namespace Kratos
{

Process::Pointer ComputeElementNodalNormalsProcess::Create(
    Model& rModel,              // retrieve the ModelPart on which the process will operate
    Parameters ThisParameters
    )
{
    ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());       // Fills missing optional parameters with defaults. If user omits a key, it gets added from GetDefaultParameters(). In Create() it ensures model_part_name exists and is valid before GetModelPart()
    return Kratos::make_shared<ComputeElementNodalNormalsProcess>(          // Elements/conditions follow intrusive-pointer patterns, whereas processes follow shared_ptr.
        rModel.GetModelPart(ThisParameters["model_part_name"].GetString()),
        ThisParameters
    );
}

ComputeElementNodalNormalsProcess::ComputeElementNodalNormalsProcess(
    ModelPart& rModelPart,
    Parameters ThisParameters
    )
    : mpModelPart(&rModelPart)
{
    ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());       // In constructor: ensures runtime flags exist and are valid before reading booleans. Catches typos or wrong types early ( recompute_each_step: "yes" instead of bool)
    mComputeOnInitialize = ThisParameters["compute_on_initialize"].GetBool();
    mRecomputeEachStep = ThisParameters["recompute_each_step"].GetBool();
}

void ComputeElementNodalNormalsProcess::ExecuteInitialize()
{
    if (mComputeOnInitialize) {
        ComputeNormals();     
    }
}

void ComputeElementNodalNormalsProcess::ExecuteInitializeSolutionStep()
{
    if (mRecomputeEachStep) {
        ComputeNormals();           // call ComputeNormals() only if recompute_each_step is true
    }
}

const Parameters ComputeElementNodalNormalsProcess::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "model_part_name"       : "",
        "compute_on_initialize" : true,     
        "recompute_each_step"   : false
    })");       // by default, compute nodal normals on nodes from the reference configuration only 
}

void ComputeElementNodalNormalsProcess::ComputeNormals()
{
    KRATOS_TRY;

    KRATOS_ERROR_IF_NOT(mpModelPart)
        << "ComputeElementNodalNormalsProcess was called without a valid ModelPart." << std::endl;

    auto& r_model_part = *mpModelPart;

    KRATOS_ERROR_IF(r_model_part.NumberOfNodes() == 0)
        << "ModelPart has no nodes: " << r_model_part.Name() << std::endl;

    NormalCalculationUtils normal_utils;
    normal_utils.CalculateUnitNormals<ModelPart::ElementsContainerType, false>(
        r_model_part,           // compute and populate nodal NORMAL for the selected model part nodes
        true,
        NORMAL
    );

    KRATOS_CATCH("");
}

} // namespace Kratos
