//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/apply_constant_vectorvalue_process.h"

namespace Kratos
{

KRATOS_CREATE_LOCAL_FLAG(ApplyConstantVectorValueProcess,X_COMPONENT_FIXED, 0);
KRATOS_CREATE_LOCAL_FLAG(ApplyConstantVectorValueProcess,Y_COMPONENT_FIXED, 1);
KRATOS_CREATE_LOCAL_FLAG(ApplyConstantVectorValueProcess,Z_COMPONENT_FIXED, 2);

/***********************************************************************************/
/***********************************************************************************/

ApplyConstantVectorValueProcess::ApplyConstantVectorValueProcess(ModelPart& rModelPart,
                                Parameters parameters
                                ) : Process(Flags()),
                                    mrModelPart(rModelPart)
{
    KRATOS_TRY

    // Get defaults parameters
    Parameters default_parameters = GetDefaultParameters();

    // Some values need to be mandatorily prescribed since no meaningful default Value exist. For this reason try accessing to them
    // So that an error is thrown if they don't exist
    KRATOS_ERROR_IF(parameters["direction"].IsArray() && parameters["direction"].size() != 3) << "direction vector is not a vector or it does not have size 3. Direction vector currently passed" << parameters.PrettyPrintJsonString() << std::endl;
    KRATOS_ERROR_IF_NOT(parameters["modulus"].IsNumber()) << "modulus shall be a number. Parameter list in which is included is :" <<  parameters.PrettyPrintJsonString() << std::endl;
    KRATOS_ERROR_IF_NOT(parameters["variable_name"].IsString()) << "variable_name shall be a String. Parameter list in which is included is :" << parameters.PrettyPrintJsonString() << std::endl;
    KRATOS_ERROR_IF_NOT(parameters["model_part_name"].IsString()) << "model_part_name shall be a String. Parameter list in which is included is :" << parameters.PrettyPrintJsonString() << std::endl;

    // Validate against defaults -- this also ensures no type mismatch
    parameters.ValidateAndAssignDefaults(default_parameters);

    // Read from the parameters and assign to the values
    mMeshId = parameters["mesh_id"].GetInt();

    this->Set(X_COMPONENT_FIXED, parameters["is_fixed_x"].GetBool());
    this->Set(Y_COMPONENT_FIXED, parameters["is_fixed_y"].GetBool());
    this->Set(Z_COMPONENT_FIXED, parameters["is_fixed_z"].GetBool());

    // Get the modulus and variable name
    mVariableName = parameters["variable_name"].GetString();
    mModulus = parameters["modulus"].GetDouble();

    mDirection.resize(3,false);
    mDirection[0] = parameters["direction"][0].GetDouble();
    mDirection[1] = parameters["direction"][1].GetDouble();
    mDirection[2] = parameters["direction"][2].GetDouble();

    const double dim_norm = norm_2(mDirection);
    KRATOS_ERROR_IF(dim_norm < 1e-20) << " Norm of direction given is approximately zero. Please give a direction vector with a non zero norm : current Value of direction vector = " << mDirection << std::endl;

    // Normalize the direction
    mDirection /= dim_norm;

    const bool has_variable = KratosComponents<Variable<array_1d<double,3>>>::Has(mVariableName);
    KRATOS_ERROR_IF_NOT(has_variable) << "Not defined the variable " << mVariableName << std::endl;
    const Variable<array_1d<double,3> >& r_variable = KratosComponents< Variable<array_1d<double,3>>>::Get(mVariableName);

    KRATOS_ERROR_IF(mMeshId >= rModelPart.NumberOfMeshes()) << "Mesh does not exist in rModelPart: mesh id is --> " << mMeshId << std::endl;

    KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(r_variable)) << "Trying to fix a variable that is not in the rModelPart - variable: "  << mVariableName << std::endl;

    KRATOS_ERROR_IF(mDirection.size() != 3) << "Direction vector is expected to have size 3. Direction vector currently passed " << mDirection << std::endl;

    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(mVariableName+ "_X")) << "Not defined the variable " << mVariableName+ "_X" << std::endl;
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(mVariableName+ "_Y")) << "Not defined the variable " << mVariableName+ "_Y" << std::endl;
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(mVariableName+ "_Z")) << "Not defined the variable " << mVariableName+ "_Z" << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

ApplyConstantVectorValueProcess::ApplyConstantVectorValueProcess(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const double Modulus,
    const Vector& Direction,
    std::size_t MeshId,
    const Flags Options
    ) : Process(Options) ,
        mrModelPart(rModelPart),
        mModulus(Modulus),
        mDirection(Direction),
        mMeshId(MeshId)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF_NOT(this->IsDefined(X_COMPONENT_FIXED) ) << "Please specify if component x is to be fixed or not  (flag X_COMPONENT_FIXED)" << std::endl;
    KRATOS_ERROR_IF_NOT(this->IsDefined(Y_COMPONENT_FIXED) ) << "Please specify if component y is to be fixed or not  (flag Y_COMPONENT_FIXED)" << std::endl;
    KRATOS_ERROR_IF_NOT(this->IsDefined(Z_COMPONENT_FIXED) ) << "Please specify if the variable is to be fixed or not (flag Z_COMPONENT_FIXED)" << std::endl;

    mVariableName = rVariable.Name();

    KRATOS_ERROR_IF(mMeshId >= rModelPart.NumberOfMeshes()) << "Mesh does not exist in rModelPart: mesh id is --> " << mMeshId << std::endl;

    KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(rVariable)) << "Trying to fix a variable that is not in the rModelPart - variable: " << mVariableName << std::endl;

    KRATOS_ERROR_IF(mDirection.size() != 3) << "Direction vector is expected to have size 3. Direction vector currently passed " << mDirection << std::endl;

    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(mVariableName+ "_X")) << "Not defined the variable " << mVariableName+ "_X" << std::endl;
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(mVariableName+ "_Y")) << "Not defined the variable " << mVariableName+ "_Y" << std::endl;
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(mVariableName+ "_Z")) << "Not defined the variable " << mVariableName+ "_Z" << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void ApplyConstantVectorValueProcess::ExecuteInitialize()
{
    // Compute the Value to be applied
    const array_1d<double,3> values = mModulus * mDirection;

    // Retrieve the variables
    const Variable<double>& r_variable_x = KratosComponents<Variable<double>>::Get(mVariableName + std::string("_X"));
    const Variable<double>& r_variable_y = KratosComponents<Variable<double>>::Get(mVariableName + std::string("_Y"));
    const Variable<double>& r_variable_z = KratosComponents<Variable<double>>::Get(mVariableName + std::string("_Z"));

    // Apply the values
    InternalApplyValue(r_variable_x, this->Is(X_COMPONENT_FIXED),  values[0]);
    InternalApplyValue(r_variable_y, this->Is(Y_COMPONENT_FIXED),  values[1]);
    InternalApplyValue(r_variable_z, this->Is(Z_COMPONENT_FIXED),  values[2]);
}

/***********************************************************************************/
/***********************************************************************************/

void ApplyConstantVectorValueProcess::InternalApplyValue(
    const Variable<double>& rVariable,
    const bool ToBeFixed,
    const double Value
    )
{
    // Get the number of nodes in the model part
    const std::size_t number_of_nodes = mrModelPart.GetMesh(mMeshId).Nodes().size();

    // Check if there are nodes in the model part
    if(number_of_nodes != 0) {
        // Get an iterator to the beginning of the nodes
        auto it_begin = mrModelPart.GetMesh(mMeshId).NodesBegin();

        // Check if the dofs are there (on the first node). Throw an error if trying to fix a variable that was not allocated
        KRATOS_ERROR_IF(ToBeFixed && (it_begin->HasDofFor(rVariable) == false)) << "Trying to fix a dofs which was not allocated. Variable is --> " << rVariable.Name() << std::endl;

        // Iterate over all nodes in the model part
        block_for_each(mrModelPart.GetMesh(mMeshId).Nodes(), [&](Node& rNode){
            // Fix the variable if needed
            if(ToBeFixed) {
                rNode.Fix(rVariable);
            }
            // Assign the constant value to the variable for the node
            rNode.FastGetSolutionStepValue(rVariable) = Value;
        });
    }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters ApplyConstantVectorValueProcess::GetDefaultParameters() const
{
    Parameters default_settings(R"({
        "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
        "mesh_id"        : 0,
        "variable_name"  : "PLEASE_PRESCRIBE_VARIABLE_NAME",
        "is_fixed_x"     : false,
        "is_fixed_y"     : false,
        "is_fixed_z"     : false,
        "modulus"        : 1.0,
        "direction"      : [1.0, 0.0, 0.0]
    })");
    return default_settings;
}

}  // namespace Kratos.


