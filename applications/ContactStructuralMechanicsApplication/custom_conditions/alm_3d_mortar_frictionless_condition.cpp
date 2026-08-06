// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// External includes

// Project includes

// #include "includes/global_variables.h"
#include "custom_conditions/alm_3d_mortar_frictionless_condition.h"

/* Utilities */
#include "utilities/geometrical_projection_utilities.h"
// #include "utilities/math_utils.h"
#include "custom_utilities/mortar_explicit_contribution_utilities.h"

namespace Kratos
{

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

Condition::Pointer ALM3dMortarFrictionlessCondition::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ALM3dMortarFrictionlessCondition>(NewId, this->GetParentGeometry().Create(rThisNodes), pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer ALM3dMortarFrictionlessCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ALM3dMortarFrictionlessCondition>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer ALM3dMortarFrictionlessCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom) const
{
    return Kratos::make_intrusive<ALM3dMortarFrictionlessCondition>(NewId, pGeom, pProperties, pMasterGeom);
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

ALM3dMortarFrictionlessCondition::~ALM3dMortarFrictionlessCondition() = default;

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::Initialize(rCurrentProcessInfo);

    // We reset the ISOLATED flag
    this->Set(ISOLATED, false);

    KRATOS_CATCH("Initialize");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::InitializeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH("InitializeSolutionStep");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::InitializeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH("InitializeNonLinearIteration");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH("FinalizeSolutionStep");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::FinalizeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH("FinalizeNonLinearIteration");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    // Calculate condition system
    CalculateConditionSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("CalculateLocalSystem");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Calculate condition system
    VectorType aux_right_hand_side_vector = Vector();
    CalculateConditionSystem(rLeftHandSideMatrix, aux_right_hand_side_vector, rCurrentProcessInfo, true, false);
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Calculate condition system
    MatrixType aux_left_hand_side_matrix = Matrix();
    CalculateConditionSystem(aux_left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo, false, true);
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateConditionSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLHS,
    const bool ComputeRHS
    )
{
    KRATOS_TRY

    // TODOOOOO

    KRATOS_CATCH("CalculateConditionSystem");
}

/***********************************************************************************/
/***********************************************************************************/

bool ALM3dMortarFrictionlessCondition::CheckIsolatedElement(
    const double DeltaTime,
    const bool HalfJump
    )
{
    if (this->Is(ISOLATED))
        return true;

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& CurrentProcessInfo
    ) const
{
    // KRATOS_ERROR << "You are calling to the base class method EquationIdVector, check your condition definition" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    // KRATOS_ERROR << "You are calling to the base class method GetDofList, check your condition definition" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    KRATOS_CATCH("CalculateOnIntegrationPoints");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>> &rVariable,
    std::vector<array_1d<double, 3>> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("CalculateOnIntegrationPoints");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    KRATOS_CATCH("CalculateOnIntegrationPoints");
}

/***********************************************************************************/
/***********************************************************************************/

int ALM3dMortarFrictionlessCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int check_result = BaseType::Check(rCurrentProcessInfo);


    return check_result;

    KRATOS_CATCH("Check")
}

/***********************************************************************************/
/***********************************************************************************/


} // Namespace Kratos
