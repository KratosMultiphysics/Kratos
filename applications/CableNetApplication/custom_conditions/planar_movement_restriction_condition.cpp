//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/planar_movement_restriction_condition.h"
#include "includes/variables.h"
#include "includes/checks.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{
PlanarMovementRestrictionCondition3D1N::PlanarMovementRestrictionCondition3D1N(PlanarMovementRestrictionCondition3D1N const& rOther)
    : BaseType(rOther)
{
}

/***********************************************************************************/
/***********************************************************************************/

PlanarMovementRestrictionCondition3D1N& PlanarMovementRestrictionCondition3D1N::operator=(PlanarMovementRestrictionCondition3D1N const& rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Condition::operator=(rOther);

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PlanarMovementRestrictionCondition3D1N::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<PlanarMovementRestrictionCondition3D1N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PlanarMovementRestrictionCondition3D1N::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<PlanarMovementRestrictionCondition3D1N>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PlanarMovementRestrictionCondition3D1N::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<PlanarMovementRestrictionCondition3D1N>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void PlanarMovementRestrictionCondition3D1N::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY
    if (rResult.size() != msLocalSize) {
        rResult.resize(msLocalSize);
    }

    rResult[0] = GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[1] = GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[2] = GetGeometry()[0].GetDof(DISPLACEMENT_Z).EquationId();
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
void PlanarMovementRestrictionCondition3D1N::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY
    if (rElementalDofList.size() != msLocalSize) {
        rElementalDofList.resize(msLocalSize);
    }

    rElementalDofList[0] = GetGeometry()[0].pGetDof(DISPLACEMENT_X);
    rElementalDofList[1] = GetGeometry()[0].pGetDof(DISPLACEMENT_Y);
    rElementalDofList[2] = GetGeometry()[0].pGetDof(DISPLACEMENT_Z);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void PlanarMovementRestrictionCondition3D1N::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    if (rValues.size() != msLocalSize) {
        rValues.resize(msLocalSize);
    }

    const auto& disp = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, Step);
    rValues[0] = disp[0];
    rValues[1] = disp[1];
    rValues[2] = disp[2];
}

/***********************************************************************************/
/***********************************************************************************/

void PlanarMovementRestrictionCondition3D1N::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    if (rValues.size() != msLocalSize) {
        rValues.resize(msLocalSize, false);
    }

    const auto& vel = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, Step);
    rValues[0] = vel[0];
    rValues[1] = vel[1];
    rValues[2] = vel[2];
}

/***********************************************************************************/
/***********************************************************************************/

void PlanarMovementRestrictionCondition3D1N::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    if (rValues.size() != msLocalSize) {
        rValues.resize(msLocalSize, false);
    }

    const auto& acc = GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION, Step);

    rValues[0] = acc[0];
    rValues[1] = acc[1];
    rValues[2] = acc[2];
}

/***********************************************************************************/
/***********************************************************************************/

void PlanarMovementRestrictionCondition3D1N::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    Vector current_nodal_displacements = ZeroVector(msLocalSize);
    GetValuesVector(current_nodal_displacements, 0);

    rRightHandSideVector = ZeroVector(msLocalSize);

    // only activated when penetrating
    if (CalculateNormalDistance(current_nodal_displacements) < 0.0) {
        const double alpha = GetProperties()[YOUNG_MODULUS]; // simplified "spring stiffness"
        Vector plane_normal = GetProperties()[NORMAL_TO_WALL]; 
        plane_normal /= norm_2(plane_normal);
        const Vector plane_position = GetProperties()[WALL_POSITION]; 

        const double scaling_factor = alpha*(plane_normal[0]*(-plane_position[0] + this->GetGeometry()[0].X0() 
            + current_nodal_displacements[0]) + plane_normal[1]*(-plane_position[1] + this->GetGeometry()[0].Y0() 
            + current_nodal_displacements[1]) + plane_normal[2]*(-plane_position[2] + this->GetGeometry()[0].Z0() 
            + current_nodal_displacements[2]) );

        noalias(rRightHandSideVector) -= scaling_factor*plane_normal;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void PlanarMovementRestrictionCondition3D1N::CalculateLeftHandSide(
            MatrixType& rLeftHandSideMatrix,
            const ProcessInfo& rCurrentProcessInfo)
{

    Vector current_nodal_displacements = ZeroVector(msLocalSize);
    GetValuesVector(current_nodal_displacements, 0);

    rLeftHandSideMatrix = ZeroMatrix(msLocalSize,msLocalSize);

    // only activated when penetrating
    if (CalculateNormalDistance(current_nodal_displacements) < 0.0) {

        const double alpha = GetProperties()[YOUNG_MODULUS]; // simplified "spring stiffness"
        Vector plane_normal = GetProperties()[NORMAL_TO_WALL]; 
        plane_normal /= norm_2(plane_normal);

        rLeftHandSideMatrix(0, 0) = std::pow(plane_normal[0], 2);
        rLeftHandSideMatrix(0, 1) = plane_normal[0]*plane_normal[1];
        rLeftHandSideMatrix(0, 2) = plane_normal[0]*plane_normal[2];

        rLeftHandSideMatrix(1, 0) = plane_normal[0]*plane_normal[1];
        rLeftHandSideMatrix(1, 1) = std::pow(plane_normal[1], 2);
        rLeftHandSideMatrix(1, 2) = plane_normal[1]*plane_normal[2];

        rLeftHandSideMatrix(2, 0) = plane_normal[0]*plane_normal[2];
        rLeftHandSideMatrix(2, 1) = plane_normal[1]*plane_normal[2];
        rLeftHandSideMatrix(2, 2) = std::pow(plane_normal[2], 2);

        rLeftHandSideMatrix *= alpha;

    }
}

/***********************************************************************************/
/***********************************************************************************/
void PlanarMovementRestrictionCondition3D1N::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void PlanarMovementRestrictionCondition3D1N::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if(rMassMatrix.size1() != 0) {
        rMassMatrix.resize(0, 0, false);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void PlanarMovementRestrictionCondition3D1N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if(rDampingMatrix.size1() != 0) {
        rDampingMatrix.resize(0, 0, false);
    }
}


double PlanarMovementRestrictionCondition3D1N::CalculateNormalDistance(const Vector& rCurrentDisplacements) const
{
    Vector plane_normal = GetProperties()[NORMAL_TO_WALL]; 
    plane_normal /= norm_2(plane_normal);
    const Vector plane_position = GetProperties()[WALL_POSITION]; 

    Vector initial_pos = ZeroVector(msLocalSize);
    initial_pos[0] = this->GetGeometry()[0].X0();
    initial_pos[1] = this->GetGeometry()[0].Y0();
    initial_pos[2] = this->GetGeometry()[0].Z0();

    return inner_prod(initial_pos + rCurrentDisplacements - plane_position, plane_normal);
}

/***********************************************************************************/
/***********************************************************************************/

int PlanarMovementRestrictionCondition3D1N::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    // Base check
    Condition::Check(rCurrentProcessInfo);

    // Check that the condition's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry().Points()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,r_node)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node)
    }


    KRATOS_ERROR_IF((GetGeometry().WorkingSpaceDimension() != 3) || (GetGeometry().size() != 1))
        << "The planar movement condition works only in 3D and with 1 noded condition" << std::endl;
    
    KRATOS_ERROR_IF(!GetProperties().Has(NORMAL_TO_WALL))
        << "No NORMAL_TO_WALL given" << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(WALL_POSITION))
        << "No WALL_POSITION given" << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(YOUNG_MODULUS))
        << "No YOUNG_MODULUS given" << std::endl;

    KRATOS_ERROR_IF(CalculateNormalDistance(ZeroVector(3)) < 0)
        << "NORMAL_TO_WALL points in wrong direction" << std::endl;

    return 0;
}



/***********************************************************************************/
/***********************************************************************************/

void PlanarMovementRestrictionCondition3D1N::AddExplicitContribution(
    const VectorType& rRHS,
    const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double,3> >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension       = GetGeometry().WorkingSpaceDimension();

    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL ) {
        for(SizeType i=0; i< number_of_nodes; ++i) {
            SizeType index = dimension * i;

            array_1d<double, 3 >& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            for(SizeType j=0; j<dimension; ++j) {
                AtomicAdd(r_force_residual[j], rRHS[index + j]);
            }
        }
    }

    KRATOS_CATCH( "" )
}


/***********************************************************************************/
/***********************************************************************************/

void PlanarMovementRestrictionCondition3D1N::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
}

/***********************************************************************************/
/***********************************************************************************/

void PlanarMovementRestrictionCondition3D1N::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
}

} // Namespace Kratos
