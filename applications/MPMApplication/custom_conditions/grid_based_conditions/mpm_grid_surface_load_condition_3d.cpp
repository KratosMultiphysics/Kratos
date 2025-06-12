//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/grid_based_conditions/mpm_grid_surface_load_condition_3d.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

MPMGridSurfaceLoadCondition3D::MPMGridSurfaceLoadCondition3D()
{
}

//***********************************************************************************
//***********************************************************************************

MPMGridSurfaceLoadCondition3D::MPMGridSurfaceLoadCondition3D(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    )
    : MPMGridBaseLoadCondition(NewId, pGeometry)
{
}

//***********************************************************************************
//***********************************************************************************

MPMGridSurfaceLoadCondition3D::MPMGridSurfaceLoadCondition3D(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    )
    : MPMGridBaseLoadCondition(NewId, pGeometry, pProperties)
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer MPMGridSurfaceLoadCondition3D::Create(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<MPMGridSurfaceLoadCondition3D>(NewId, pGeometry, pProperties);
}

//***********************************************************************************
//***********************************************************************************

Condition::Pointer MPMGridSurfaceLoadCondition3D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<MPMGridSurfaceLoadCondition3D>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

MPMGridSurfaceLoadCondition3D::~MPMGridSurfaceLoadCondition3D()
{
}

//***********************************************************************************
//***********************************************************************************

void MPMGridSurfaceLoadCondition3D::CalculateAndSubKp(
    Matrix& rK,
    const array_1d<double, 3 >& rge,
    const array_1d<double, 3 >& rgn,
    const Matrix& rDN_De,
    const RowMatrix& rN,
    const double Pressure,
    const double Weight)
{
    KRATOS_TRY

    Matrix Kij(3, 3);
    BoundedMatrix<double, 3, 3 > Cross_ge;
    BoundedMatrix<double, 3, 3 > Cross_gn;
    double coeff;
    const unsigned int number_of_nodes = GetGeometry().size();

    MakeCrossMatrix(Cross_ge, rge);
    MakeCrossMatrix(Cross_gn, rgn);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const unsigned int RowIndex = i * 3;
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
            const unsigned int ColIndex = j * 3;

            coeff = Pressure * rN[i] * rDN_De(j, 1) * Weight;
            noalias(Kij) = coeff * Cross_ge;

            coeff = Pressure * rN[i] * rDN_De(j, 0) * Weight;

            noalias(Kij) -= coeff * Cross_gn;

            //TAKE CARE: the load correction matrix should be SUBTRACTED not added
            MathUtils<double>::SubtractMatrix(rK, Kij, RowIndex, ColIndex);
        }
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void MPMGridSurfaceLoadCondition3D::MakeCrossMatrix(
    BoundedMatrix<double, 3, 3 > & M,
    const array_1d<double, 3 > & U)
{
    M(0, 0) = 0.0;
    M(0, 1) = -U[2];
    M(0, 2) = U[1];
    M(1, 0) = U[2];
    M(1, 1) = 0.0;
    M(1, 2) = -U[0];
    M(2, 0) = -U[1];
    M(2, 1) = U[0];
    M(2, 2) = 0.0;
}

//***********************************************************************************
//***********************************************************************************

void MPMGridSurfaceLoadCondition3D::CalculateAndAddPressureForce(
    VectorType& rResidualVector,
    const RowMatrix& N,
    const array_1d<double, 3 >& Normal,
    const double Pressure,
    const double Weight,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const unsigned int number_of_nodes = GetGeometry().size();

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const int index = 3 * i;
        const double coeff = Pressure * N[i] * Weight;
        rResidualVector[index    ] -= coeff * Normal[0];
        rResidualVector[index + 1] -= coeff * Normal[1];
        rResidualVector[index + 2] -= coeff * Normal[2];
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void MPMGridSurfaceLoadCondition3D::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int matrix_size = number_of_nodes * 3;

    //Resizing as needed the LHS
    if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
    {
        if (rLeftHandSideMatrix.size1() != matrix_size)
        {
            rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
        }

        noalias(rLeftHandSideMatrix) = ZeroMatrix(matrix_size, matrix_size); //resetting LHS
    }

    // Resizing as needed the RHS
    if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
    {
        if (rRightHandSideVector.size() != matrix_size)
        {
            rRightHandSideVector.resize(matrix_size, false);
        }

        rRightHandSideVector = ZeroVector(matrix_size); //resetting RHS
    }

    // Reading integration points and local gradients
    IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geometry);
    const auto& integration_points = r_geometry.IntegrationPoints(integration_method);
    const auto& r_DN_De = r_geometry.ShapeFunctionsLocalGradients(integration_method);
    const Matrix& r_N = r_geometry.ShapeFunctionsValues(integration_method);

    // Calculating actual jacobian
    GeometryType::JacobiansType J;
    J = r_geometry.Jacobian(J,integration_method);

    // Vector with a loading applied to the elemnt
    array_1d<double, 3 > surface_load = ZeroVector(3);
    if( this->Has( SURFACE_LOAD ) )
    {
        noalias(surface_load) = this->GetValue( SURFACE_LOAD );
    }

    // Pressure applied to the element itself
    double pressure_on_condition = 0.0;
    if( this->Has( PRESSURE ) )
    {
        pressure_on_condition += this->GetValue( PRESSURE );
    }
    if( this->Has( NEGATIVE_FACE_PRESSURE ) )
    {
        pressure_on_condition += this->GetValue( NEGATIVE_FACE_PRESSURE );
    }
    if( this->Has( POSITIVE_FACE_PRESSURE ) )
    {
        pressure_on_condition -= this->GetValue( POSITIVE_FACE_PRESSURE );
    }

    // Auxiliary terms
    Vector pressure_on_nodes(number_of_nodes, pressure_on_condition); //note that here we initialize from the value applied to the condition

    for (unsigned int i = 0; i < pressure_on_nodes.size(); i++)
    {
        if( r_geometry[i].SolutionStepsDataHas( NEGATIVE_FACE_PRESSURE) )
        {
            pressure_on_nodes[i] += r_geometry[i].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE );
        }
        if( r_geometry[i].SolutionStepsDataHas( POSITIVE_FACE_PRESSURE) )
        {
            pressure_on_nodes[i] -= r_geometry[i].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE );
        }
    }

    array_1d<double, 3 > ge, gn;

    for (unsigned int point_number = 0; point_number < integration_points.size(); point_number++)
    {
        const double det_j = MathUtils<double>::GeneralizedDet(J[point_number]);
        const double integration_weight = GetIntegrationWeight(integration_points, point_number, det_j);
        const auto& N = row(r_N, point_number);

        ge[0] = J[point_number](0, 0);
        gn[0] = J[point_number](0, 1);
        ge[1] = J[point_number](1, 0);
        gn[1] = J[point_number](1, 1);
        ge[2] = J[point_number](2, 0);
        gn[2] = J[point_number](2, 1);

        array_1d<double, 3 > normal;
        MathUtils<double>::UnitCrossProduct(normal, gn, ge);

        // Calculating the pressure on the gauss point
        double pressure = 0.0;
        for (unsigned int ii = 0; ii < number_of_nodes; ii++)
        {
            pressure += N[ii] * pressure_on_nodes[ii];
        }

        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true)
        {
            if (std::abs(pressure) > std::numeric_limits<double>::epsilon())
            {
                CalculateAndSubKp(rLeftHandSideMatrix, ge, gn, r_DN_De[point_number], N, pressure, integration_weight);
            }
        }

        // RIGHT HAND SIDE VECTOR
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            if (std::abs(pressure) > std::numeric_limits<double>::epsilon())
            {
                CalculateAndAddPressureForce(rRightHandSideVector, N, normal, pressure, integration_weight, rCurrentProcessInfo);
            }
        }

        //generic load on gauss point
        array_1d<double, 3> gauss_load = surface_load;
        for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
        {
            if( r_geometry[ii].SolutionStepsDataHas( SURFACE_LOAD ) )
            {
                noalias(gauss_load) += N[ii]*r_geometry[ii].FastGetSolutionStepValue( SURFACE_LOAD );
            }
        }

        for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
        {
            const unsigned int base = ii * 3;
            for(unsigned int k = 0; k < 3; ++k)
            {
                rRightHandSideVector[base+k] += integration_weight * N[ii] * gauss_load[k];
            }
        }
    }


    KRATOS_CATCH("")
}

} // Namespace Kratos.
