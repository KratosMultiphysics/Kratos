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
#include "geometries/line_2d_4.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"

// Project includes
#include "custom_conditions/general_U_Pw_diff_order_condition.hpp"
#include "custom_utilities/dof_utilities.h"

namespace Kratos
{

Condition::Pointer GeneralUPwDiffOrderCondition::
    Create(IndexType NewId,
           NodesArrayType const& ThisNodes,
           PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new GeneralUPwDiffOrderCondition(NewId,
                                                               GetGeometry().Create(ThisNodes),
                                                               pProperties));
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const SizeType NumUNodes = rGeom.PointsNumber();

    switch(NumUNodes) {
        case 3: //2D L3P2
            mpPressureGeometry = make_shared<Line2D2<Node>>(rGeom(0), rGeom(1));
            break;
        case 4: //2D L4P3
            mpPressureGeometry = make_shared<Line2D3<Node>>(rGeom(0), rGeom(1), rGeom(2));
            break;
        case 5: //2D L5P4
            mpPressureGeometry = make_shared<Line2D4<Node>>(rGeom(0), rGeom(1), rGeom(2), rGeom(3));
            break;
        case 6: //3D T6P3
            mpPressureGeometry = make_shared<Triangle3D3<Node>>(rGeom(0), rGeom(1), rGeom(2));
            break;
        case 8: //3D Q8P4
            mpPressureGeometry = make_shared<Quadrilateral3D4<Node>>(rGeom(0), rGeom(1), rGeom(2), rGeom(3));
            break;
        case 9: //3D Q9P4
            mpPressureGeometry = make_shared<Quadrilateral3D4<Node>>(rGeom(0), rGeom(1), rGeom(2), rGeom(3));
            break;
        default:
            KRATOS_ERROR << "Unexpected geometry type for different order interpolation element" << std::endl;
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo&) const
{
    rConditionDofList = GetDofs();
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::
    CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                         VectorType& rRightHandSideVector,
                         const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType ConditionSize = NumUNodes * Dim + NumPNodes;

    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != ConditionSize )
        rLeftHandSideMatrix.resize( ConditionSize, ConditionSize, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( ConditionSize, ConditionSize );

    //Resetting the RHS
    if ( rRightHandSideVector.size() != ConditionSize )
        rRightHandSideVector.resize( ConditionSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( ConditionSize );

    //calculation flags
    bool CalculateLHSMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll(rLeftHandSideMatrix,
                 rRightHandSideVector,
                 rCurrentProcessInfo,
                 CalculateLHSMatrixFlag,
                 CalculateResidualVectorFlag);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::
    CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix,
                           const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_ERROR << "GeneralUPwDiffOrderCondition::CalculateLeftHandSide is not implemented" << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::
    CalculateRightHandSide( VectorType& rRightHandSideVector,
                            const ProcessInfo& rCurrentProcessInfo )
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType ConditionSize = NumUNodes * Dim + NumPNodes;

    //Resetting the RHS
    if ( rRightHandSideVector.size() != ConditionSize )
        rRightHandSideVector.resize( ConditionSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( ConditionSize );

    //calculation flags
    bool CalculateLHSMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateLHSMatrixFlag, CalculateResidualVectorFlag);
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::
    CalculateAll(MatrixType& rLeftHandSideMatrix,
                 VectorType& rRightHandSideVector,
                 const ProcessInfo& rCurrentProcessInfo,
                 bool CalculateLHSMatrixFlag,
                 bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    //Definition of variables
    ConditionVariables Variables;
    this->InitializeConditionVariables(Variables,rCurrentProcessInfo);

    //Loop over integration points
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    for ( unsigned int PointNumber = 0; PointNumber < IntegrationPoints.size(); PointNumber++ )
    {
        //compute element kinematics (Np)
        this->CalculateKinematics(Variables,PointNumber);

        //Compute Condition Vector
        this->CalculateConditionVector(Variables,PointNumber);

        //Calculating weighting coefficient for integration
        Variables.IntegrationCoefficient = 
            this->CalculateIntegrationCoefficient(PointNumber,
                                                  Variables.JContainer,
                                                  IntegrationPoints);

        //Contributions to the left hand side
        if ( CalculateLHSMatrixFlag )
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        //Contributions to the right hand side
        if ( CalculateResidualVectorFlag )
            this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
void GeneralUPwDiffOrderCondition::
    InitializeConditionVariables(ConditionVariables& rVariables,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType NumGPoints = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());
    const SizeType WorkingDim = rGeom.WorkingSpaceDimension();
    const SizeType LocalDim = rGeom.LocalSpaceDimension();

    (rVariables.NuContainer).resize(NumGPoints,NumUNodes,false);
    rVariables.NuContainer = rGeom.ShapeFunctionsValues(this->GetIntegrationMethod());

    (rVariables.NpContainer).resize(NumGPoints,NumPNodes,false);
    rVariables.NpContainer = mpPressureGeometry->ShapeFunctionsValues(this->GetIntegrationMethod());

    (rVariables.Nu).resize(NumUNodes,false);
    (rVariables.Np).resize(NumPNodes,false);

    (rVariables.JContainer).resize(NumGPoints,false);
    for(SizeType i = 0; i<NumGPoints; ++i)
        ((rVariables.JContainer)[i]).resize(WorkingDim,LocalDim,false);
    rGeom.Jacobian(rVariables.JContainer, this->GetIntegrationMethod());
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::CalculateKinematics(ConditionVariables& rVariables,unsigned int PointNumber)
{
    KRATOS_TRY

    //Setting the shape function vector
    noalias(rVariables.Nu) = row( rVariables.NuContainer, PointNumber);
    noalias(rVariables.Np) = row( rVariables.NpContainer, PointNumber);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateConditionVector method for a particular condition ... illegal operation!!" << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
double GeneralUPwDiffOrderCondition::
    CalculateIntegrationCoefficient(const IndexType PointNumber,
                                    const GeometryType::JacobiansType& JContainer,
                                    const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const

{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateIntegrationCoefficient method for a particular condition ... illegal operation!!" << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables)
{

}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::CalculateAndAddRHS(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateAndAddConditionForce(rRightHandSideVector, rVariables);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateAndAddConditionForce method for a particular condition ... illegal operation!!" << std::endl;

    KRATOS_CATCH( "" )
}

Condition::DofsVectorType GeneralUPwDiffOrderCondition::GetDofs() const
{
    Condition::DofsVectorType result;

    for (const auto& r_node : GetGeometry()) {
        result.push_back(r_node.pGetDof(DISPLACEMENT_X));
        result.push_back(r_node.pGetDof(DISPLACEMENT_Y));
        if (GetGeometry().WorkingSpaceDimension() == 3) {
            result.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }
    }

    const auto water_pressure_dofs =
        Geo::DofUtilities::ExtractDofsFromNodes(*mpPressureGeometry, WATER_PRESSURE);
    result.insert(result.end(), water_pressure_dofs.begin(), water_pressure_dofs.end());

    return result;
}

} // Namespace Kratos.
