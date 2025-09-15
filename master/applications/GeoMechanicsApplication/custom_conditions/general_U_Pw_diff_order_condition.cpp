// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_2d_4.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/triangle_3d_3.h"

// Project includes
#include "custom_conditions/general_U_Pw_diff_order_condition.hpp"
#include "custom_utilities/dof_utilities.h"

namespace Kratos
{

Condition::Pointer GeneralUPwDiffOrderCondition::Create(IndexType               NewId,
                                                        NodesArrayType const&   ThisNodes,
                                                        PropertiesType::Pointer pProperties) const
{
    return Create(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Condition::Pointer GeneralUPwDiffOrderCondition::Create(IndexType               NewId,
                                                        GeometryType::Pointer   pGeom,
                                                        PropertiesType::Pointer pProperties) const
{
    return make_intrusive<GeneralUPwDiffOrderCondition>(NewId, pGeom, pProperties);
}

void GeneralUPwDiffOrderCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& r_geometry = GetGeometry();

    switch (r_geometry.PointsNumber()) {
    case 3: // 2D L3P2
        mpPressureGeometry = make_shared<Line2D2<Node>>(r_geometry(0), r_geometry(1));
        break;
    case 4: // 2D L4P3
        mpPressureGeometry = make_shared<Line2D3<Node>>(r_geometry(0), r_geometry(1), r_geometry(2));
        break;
    case 5: // 2D L5P4
        mpPressureGeometry =
            make_shared<Line2D4<Node>>(r_geometry(0), r_geometry(1), r_geometry(2), r_geometry(3));
        break;
    case 6: // 3D T6P3
        mpPressureGeometry = make_shared<Triangle3D3<Node>>(r_geometry(0), r_geometry(1), r_geometry(2));
        break;
    case 8: // 3D Q8P4
        mpPressureGeometry = make_shared<Quadrilateral3D4<Node>>(r_geometry(0), r_geometry(1),
                                                                 r_geometry(2), r_geometry(3));
        break;
    case 9: // 3D Q9P4
        mpPressureGeometry = make_shared<Quadrilateral3D4<Node>>(r_geometry(0), r_geometry(1),
                                                                 r_geometry(2), r_geometry(3));
        break;
    default:
        KRATOS_ERROR << "Unexpected geometry type for different order interpolation element" << std::endl;
    }

    KRATOS_CATCH("")
}

void GeneralUPwDiffOrderCondition::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo&) const
{
    rConditionDofList = GetDofs();
}

void GeneralUPwDiffOrderCondition::CalculateLocalSystem(Matrix&            rLeftHandSideMatrix,
                                                        Vector&            rRightHandSideVector,
                                                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto&    r_geom = GetGeometry();
    const SizeType condition_size =
        r_geom.PointsNumber() * r_geom.WorkingSpaceDimension() + mpPressureGeometry->PointsNumber();

    // Resetting the LHS
    rLeftHandSideMatrix = ZeroMatrix(condition_size, condition_size);

    // Resetting the RHS
    rRightHandSideVector = ZeroVector(condition_size);

    // calculation flags
    bool CalculateLHSMatrixFlag      = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                 CalculateLHSMatrixFlag, CalculateResidualVectorFlag);

    KRATOS_CATCH("")
}

void GeneralUPwDiffOrderCondition::CalculateRightHandSide(Vector&            rRightHandSideVector,
                                                          const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geom = GetGeometry();
    const auto  condition_size =
        r_geom.PointsNumber() * r_geom.WorkingSpaceDimension() + mpPressureGeometry->PointsNumber();

    // Resetting the RHS
    rRightHandSideVector = ZeroVector(condition_size);

    // calculation flags
    bool CalculateLHSMatrixFlag      = false;
    bool CalculateResidualVectorFlag = true;
    auto temp                        = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateLHSMatrixFlag, CalculateResidualVectorFlag);
}

void GeneralUPwDiffOrderCondition::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

void GeneralUPwDiffOrderCondition::CalculateAll(const Matrix&,
                                                Vector&            rRightHandSideVector,
                                                const ProcessInfo& rCurrentProcessInfo,
                                                bool,
                                                bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Definition of variables
    ConditionVariables Variables;
    this->InitializeConditionVariables(Variables, rCurrentProcessInfo);

    // Loop over integration points
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
        GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    for (unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++) {
        // compute element kinematics (Np)
        this->CalculateKinematics(Variables, PointNumber);

        // Compute Condition Vector
        this->CalculateConditionVector(Variables, PointNumber);

        // Calculating weighting coefficient for integration
        Variables.IntegrationCoefficient =
            this->CalculateIntegrationCoefficient(PointNumber, Variables.JContainer, IntegrationPoints);

        // Contributions to the right hand side
        if (CalculateResidualVectorFlag) this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }

    KRATOS_CATCH("")
}

void GeneralUPwDiffOrderCondition::InitializeConditionVariables(ConditionVariables& rVariables,
                                                                const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geom      = GetGeometry();
    const auto  num_u_nodes = r_geom.PointsNumber();
    const auto  num_p_nodes = mpPressureGeometry->PointsNumber();
    const auto number_of_integration_points = r_geom.IntegrationPointsNumber(this->GetIntegrationMethod());

    rVariables.NuContainer = r_geom.ShapeFunctionsValues(this->GetIntegrationMethod());

    rVariables.NpContainer = mpPressureGeometry->ShapeFunctionsValues(this->GetIntegrationMethod());

    rVariables.Nu.resize(num_u_nodes, false);
    rVariables.Np.resize(num_p_nodes, false);

    rVariables.JContainer.resize(number_of_integration_points, false);
    for (auto& j : rVariables.JContainer)
        j.resize(r_geom.WorkingSpaceDimension(), r_geom.LocalSpaceDimension(), false);
    r_geom.Jacobian(rVariables.JContainer, this->GetIntegrationMethod());
}

void GeneralUPwDiffOrderCondition::CalculateKinematics(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    // Setting the shape function vector
    noalias(rVariables.Nu) = row(rVariables.NuContainer, PointNumber);
    noalias(rVariables.Np) = row(rVariables.NpContainer, PointNumber);

    KRATOS_CATCH("")
}

void GeneralUPwDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateConditionVector method for a particular "
                    "condition ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH("")
}

double GeneralUPwDiffOrderCondition::CalculateIntegrationCoefficient(
    IndexType                                       PointNumber,
    const GeometryType::JacobiansType&              JContainer,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateIntegrationCoefficient method for a particular "
                    "condition ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH("")
}

void GeneralUPwDiffOrderCondition::CalculateAndAddRHS(Vector& rRightHandSideVector, ConditionVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateAndAddConditionForce(rRightHandSideVector, rVariables);

    KRATOS_CATCH("")
}

void GeneralUPwDiffOrderCondition::CalculateAndAddConditionForce(Vector& rRightHandSideVector,
                                                                 ConditionVariables& rVariables)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateAndAddConditionForce method for a particular "
                    "condition ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH("")
}

Condition::DofsVectorType GeneralUPwDiffOrderCondition::GetDofs() const
{
    const auto& r_geometry = this->GetGeometry();
    return Geo::DofUtilities::ExtractUPwDofsFromNodes(r_geometry, *mpPressureGeometry,
                                                      r_geometry.WorkingSpaceDimension());
}

std::string GeneralUPwDiffOrderCondition::Info() const { return "GeneralUPwDiffOrderCondition"; }

} // Namespace Kratos.
