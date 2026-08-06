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

// #include <algorithm>
// #ifdef KRATOS_DEBUG
// #include <iomanip>
// #endif

// External includes

// Project includes
#include "includes/global_variables.h"
// #include "custom_conditions/mortar_contact_condition.h"

/* Utilities */
#include "utilities/geometrical_projection_utilities.h"
#include "utilities/math_utils.h"
#include "custom_utilities/alm_3d_mortar_frictionless_condition.h"

namespace Kratos
{

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

Condition::Pointer ALM3dMortarFrictionlessCondition::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties ) const
{
    KRATOS_ERROR << "You are calling to the base class method Create, check your condition declaration" << std::endl;

    return Kratos::make_intrusive< ALM3dMortarFrictionlessCondition >( NewId, this->GetParentGeometry().Create( rThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer ALM3dMortarFrictionlessCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_ERROR << "You are calling to the base class method Create, check your condition declaration" << std::endl;

    return Kratos::make_intrusive< ALM3dMortarFrictionlessCondition >( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer ALM3dMortarFrictionlessCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom) const
{
    KRATOS_ERROR << "You are calling to the base class method Create, check your condition declaration" << std::endl;

    return Kratos::make_intrusive< ALM3dMortarFrictionlessCondition >( NewId, pGeom, pProperties, pMasterGeom );
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

ALM3dMortarFrictionlessCondition::~ALM3dMortarFrictionlessCondition( )
= default;

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::Initialize(rCurrentProcessInfo);

    // We reset the ISOLATED flag
    this->Set(ISOLATED, false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::InitializeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::InitializeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    BaseType::FinalizeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Resizing as needed
    // ResizeLHS(rLeftHandSideMatrix);
    // ResizeRHS(rRightHandSideVector);

    // Calculate condition system
    CalculateConditionSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo );

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Resizing as needed
    // ResizeLHS(rLeftHandSideMatrix);

    // Creating an auxiliary vector
    VectorType aux_right_hand_side_vector = Vector();

    // Calculate condition system
    CalculateConditionSystem(rLeftHandSideMatrix, aux_right_hand_side_vector, rCurrentProcessInfo, true, false );
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Creating an auxiliary matrix
    MatrixType aux_left_hand_side_matrix = Matrix();

    // Resizing as needed
    ResizeRHS(rRightHandSideVector);

    // Calculate condition system
    CalculateConditionSystem(aux_left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo, false );
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // rMassMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // rDampingMatrix.resize(0, 0, false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::AddExplicitContribution(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // const IndexType integration_order = GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
    // MortarExplicitContributionUtilities<TDim, TNumNodes, FrictionalCase::FRICTIONLESS_PENALTY, TNormalVariation, TNumNodesMaster>::AddExplicitContributionOfMortarCondition(this, rCurrentProcessInfo, integration_order, IsAxisymmetric(), false);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<double >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // const IndexType integration_order = this->GetProperties().Has(INTEGRATION_ORDER_CONTACT) ? this->GetProperties().GetValue(INTEGRATION_ORDER_CONTACT) : 2;
    // MortarExplicitContributionUtilities<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ComputeNodalArea(this, rCurrentProcessInfo, rDestinationVariable, integration_order, false);
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double, 3> >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "You are calling to the base class method AddExplicitContribution, check your condition definition" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::CalculateConditionSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLHS,
    const bool ComputeRHS
    )
{
    KRATOS_TRY;

    // TODOOOOO

    KRATOS_CATCH( "" );
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

// void ALM3dMortarFrictionlessCondition::CalculateLocalLHS(
//     Matrix& rLocalLHS,
//     const MortarConditionMatrices& rMortarConditionMatrices,
//     const DerivativeDataType& rDerivativeData,
//     const IndexType rActiveInactive,
//     const ProcessInfo& rCurrentProcessInfo
//     )
// {
//     KRATOS_ERROR << "You are calling to the base class method CalculateLocalLHS, check your condition definition" << std::endl;
// }

// /***********************************************************************************/
// /***********************************************************************************/

// void ALM3dMortarFrictionlessCondition::CalculateLocalRHS(
//     Vector& rLocalRHS,
//     const MortarConditionMatrices& rMortarConditionMatrices,
//     const DerivativeDataType& rDerivativeData,
//     const IndexType rActiveInactive,
//     const ProcessInfo& rCurrentProcessInfo
//     )
// {
//     KRATOS_ERROR << "You are calling to the base class method CalculateLocalRHS, check your condition definition" << std::endl;
// }

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

void ALM3dMortarFrictionlessCondition<TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::GetDofList(
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
    KRATOS_TRY;

    // const GeometryType::IntegrationPointsArrayType &integration_points = GetParentGeometry().IntegrationPoints();

    // if ( rOutput.size() != integration_points.size() )
    //     rOutput.resize( integration_points.size() );

    // for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
    //     rOutput[point_number] = 0.0;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // const GeometryType::IntegrationPointsArrayType &integration_points = GetParentGeometry().IntegrationPoints();

    // if ( rOutput.size() != integration_points.size() )
    //     rOutput.resize( integration_points.size() );

    // for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
    //     rOutput[point_number] = ZeroVector(3);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // const GeometryType::IntegrationPointsArrayType &integration_points = GetParentGeometry().IntegrationPoints();

    // if ( rOutput.size() != integration_points.size() )
    //     rOutput.resize( integration_points.size() );

    // for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number)
    //     rOutput[point_number] = ZeroVector(3);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

int ALM3dMortarFrictionlessCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // // Base class checks for positive Jacobian and Id > 0
    // int ierr = Condition::Check(rCurrentProcessInfo);
    // if(ierr != 0) return ierr;

    // const GeometryType& r_current_slave = this->GetParentGeometry();
    // KRATOS_ERROR_IF(r_current_slave.NumberOfGeometryParts() == 0) << "YOU HAVE NOT INITIALIZED THE PAIR GEOMETRY IN THE ALM3dMortarFrictionlessCondition" << std::endl;

    // // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    // for ( IndexType i = 0; i < TNumNodes; ++i ) {
    //     const auto& r_node = r_current_slave[i];
    //     KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,r_node)
    //     KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(WEIGHTED_GAP,r_node)
    //     KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NORMAL,r_node)

    //     KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node)
    //     KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node)
    //     KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node)
    // }

    // return ierr;

    return 0;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

// bool ALM3dMortarFrictionlessCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::IsAxisymmetric() const
// {
//     return false;
// }

// /***********************************************************************************/
// /***********************************************************************************/

// double ALM3dMortarFrictionlessCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::GetAxisymmetricCoefficient(const GeneralVariables& rVariables) const
// {
//     return 1.0;
// }

/***********************************************************************************/
/***********************************************************************************/

// template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
// void ALM3dMortarFrictionlessCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ResizeLHS(MatrixType& rLeftHandSideMatrix)
// {
//     // Resizing as needed the LHS
//     if ( rLeftHandSideMatrix.size1() != MatrixSize || rLeftHandSideMatrix.size2() != MatrixSize )
//             rLeftHandSideMatrix.resize( MatrixSize, MatrixSize, false );
// }

/***********************************************************************************/
/***********************************************************************************/

// template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
// void ALM3dMortarFrictionlessCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ResizeRHS(VectorType& rRightHandSideVector)
// {
//     // Resizing as needed the RHS
//     if ( rRightHandSideVector.size() != MatrixSize )
//         rRightHandSideVector.resize( MatrixSize, false );
// }

/***********************************************************************************/
/***********************************************************************************/

// template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
// void ALM3dMortarFrictionlessCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ZeroLHS(MatrixType& rLeftHandSideMatrix)
// {
//     rLeftHandSideMatrix = ZeroMatrix( MatrixSize, MatrixSize );
// }

// /***********************************************************************************/
// /***********************************************************************************/

// template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
// void ALM3dMortarFrictionlessCondition< TDim, TNumNodes, TFrictional, TNormalVariation, TNumNodesMaster>::ZeroRHS(VectorType& rRightHandSideVector)
// {
//     rRightHandSideVector = ZeroVector( MatrixSize );
// }

/***********************************************************************************/
/***********************************************************************************/

// Frictionless cases
// template class ALM3dMortarFrictionlessCondition<2, 2, FrictionalCase::FRICTIONLESS, false, 2>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONLESS, false, 3>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONLESS, false, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONLESS, false, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONLESS, false, 3>;
// template class ALM3dMortarFrictionlessCondition<2, 2, FrictionalCase::FRICTIONLESS, true, 2>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONLESS, true, 3>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONLESS, true, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONLESS, true, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONLESS, true, 3>;

// // Frictionless components cases
// template class ALM3dMortarFrictionlessCondition<2, 2, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 2>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 3>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 3>;
// template class ALM3dMortarFrictionlessCondition<2, 2, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 2>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 3>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 3>;

// // Frictional cases
// template class ALM3dMortarFrictionlessCondition<2, 2, FrictionalCase::FRICTIONAL, false, 2>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONAL, false, 3>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONAL, false, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONAL, false, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONAL, false, 3>;
// template class ALM3dMortarFrictionlessCondition<2, 2, FrictionalCase::FRICTIONAL, true, 2>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONAL, true, 3>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONAL, true, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONAL, true, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONAL, true, 3>;

// // Frictionless penalty cases
// template class ALM3dMortarFrictionlessCondition<2, 2, FrictionalCase::FRICTIONLESS_PENALTY, false, 2>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONLESS_PENALTY, false, 3>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONLESS_PENALTY, false, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONLESS_PENALTY, false, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONLESS_PENALTY, false, 3>;
// template class ALM3dMortarFrictionlessCondition<2, 2, FrictionalCase::FRICTIONLESS_PENALTY, true, 2>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONLESS_PENALTY, true, 3>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONLESS_PENALTY, true, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONLESS_PENALTY, true, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONLESS_PENALTY, true, 3>;

// // Frictional penalty cases
// template class ALM3dMortarFrictionlessCondition<2, 2, FrictionalCase::FRICTIONAL_PENALTY, false, 2>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONAL_PENALTY, false, 3>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONAL_PENALTY, false, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONAL_PENALTY, false, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONAL_PENALTY, false, 3>;
// template class ALM3dMortarFrictionlessCondition<2, 2, FrictionalCase::FRICTIONAL_PENALTY, true, 2>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONAL_PENALTY, true, 3>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONAL_PENALTY, true, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 3, FrictionalCase::FRICTIONAL_PENALTY, true, 4>;
// template class ALM3dMortarFrictionlessCondition<3, 4, FrictionalCase::FRICTIONAL_PENALTY, true, 3>;

} // Namespace Kratos
