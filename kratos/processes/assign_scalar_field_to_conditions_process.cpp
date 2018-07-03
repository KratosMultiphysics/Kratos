//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "processes/assign_scalar_field_to_conditions_process.h"

namespace Kratos
{
AssignScalarFieldToConditionsProcess::AssignScalarFieldToConditionsProcess(
    ModelPart& rModelPart,
    Parameters rParameters
    ) : Process() ,
        mrModelPart(rModelPart)
{
    KRATOS_TRY

    Parameters default_parameters( R"(
    {
        "model_part_name":"MODEL_PART_NAME",
        "mesh_id": 0,
        "variable_name": "VARIABLE_NAME",
        "interval"        : [0.0, 1e30],
        "value"           : "please give an expression in terms of the variable x, y, z, t",
        "local_axes" : {}
    }  )" );

    // Validate against defaults -- this ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mMeshId       = rParameters["mesh_id"].GetInt();
    mVariableName = rParameters["variable_name"].GetString();

    mpFunction = PythonGenericFunctionUtility::Pointer( new PythonGenericFunctionUtility(rParameters["value"].GetString(),  rParameters["local_axes"]));

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void AssignScalarFieldToConditionsProcess::Execute()
{
    KRATOS_TRY;

    ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();

    const double rCurrentTime = rCurrentProcessInfo[TIME];

    if( KratosComponents< Variable<double> >::Has( mVariableName ) ) { //case of scalar variable
        InternalAssignValueScalar<>(KratosComponents< Variable<double> >::Get(mVariableName), rCurrentTime);
    } else if( KratosComponents< array_1d_component_type >::Has( mVariableName ) ) { //case of component variable
        InternalAssignValueScalar<>(KratosComponents< array_1d_component_type >::Get(mVariableName), rCurrentTime);
    } else if( KratosComponents< Variable<Vector> >::Has( mVariableName ) ) { //case of vector variable
        InternalAssignValueVector<>(KratosComponents< Variable<Vector> >::Get(mVariableName), rCurrentTime);
    } else {
        KRATOS_ERROR << "Not able to set the variable. Attempting to set variable:" << mVariableName << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/


void AssignScalarFieldToConditionsProcess::CallFunction(
    const Condition::Pointer& pCondition,
    const double Time,
    Vector& rValue
    )
{
    GeometryType& rConditionGeometry = pCondition->GetGeometry();
    const SizeType size = rConditionGeometry.size();

    if(rValue.size() !=  size) {
        rValue.resize(size,false);
    }

    for (IndexType i=0; i<size; ++i) {
        rValue[i] = mpFunction->CallFunction(rConditionGeometry[i].X(),rConditionGeometry[i].Y(),rConditionGeometry[i].Z(),Time  );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AssignScalarFieldToConditionsProcess::CallFunctionComponents(
    const Condition::Pointer& pCondition,
    const double Time,
    double& rValue
    )
{
    GeometryType& rConditionGeometry = pCondition->GetGeometry();
    const array_1d<double,3>& center = rConditionGeometry.Center();

    rValue = mpFunction->CallFunction(center[0],center[1],center[2], Time  );
}

/***********************************************************************************/
/***********************************************************************************/

void AssignScalarFieldToConditionsProcess::CallFunctionLocalSystem(
    const Condition::Pointer& pCondition,
    const double Time,
    Vector& rValue
    )
{

    GeometryType& rConditionGeometry = pCondition->GetGeometry();
    const SizeType size = rConditionGeometry.size();

    if (rValue.size() !=  size) {
        rValue.resize(size,false);
    }

    for (IndexType i=0; i<size; ++i) {
        rValue[i] = mpFunction->RotateAndCallFunction(rConditionGeometry[i].X(),rConditionGeometry[i].Y(),rConditionGeometry[i].Z(), Time  );
    }
}

/***********************************************************************************/
/***********************************************************************************/

void AssignScalarFieldToConditionsProcess::CallFunctionLocalSystemComponents(
    const Condition::Pointer& pCondition,
    const double Time,
    double& rValue
    )
{
    GeometryType& rConditionGeometry = pCondition->GetGeometry();

    const array_1d<double,3>& center = rConditionGeometry.Center();

    rValue = mpFunction->RotateAndCallFunction(center[0],center[1],center[2], Time  );
}

/***********************************************************************************/
/***********************************************************************************/

void AssignScalarFieldToConditionsProcess::AssignTimeDependentValue(
    const Condition::Pointer& pCondition,
    const double Time,
    Vector& rValue,
    const double Value
    )
{
    GeometryType& rConditionGeometry = pCondition->GetGeometry();
    const SizeType size = rConditionGeometry.size();

    if(rValue.size() !=  size) {
        rValue.resize(size,false);
    }

    for(IndexType i=0; i<size; ++i) {
        rValue[i] = Value;
    }
}
}
