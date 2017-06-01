// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/surface_load_condition_3d.h"
#include "utilities/math_utils.h"
// #include "custom_utilities/sd_math_utils.h"

namespace Kratos
{

SurfaceLoadCondition3D::SurfaceLoadCondition3D()
{
}

// Constructor

SurfaceLoadCondition3D::SurfaceLoadCondition3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : BaseLoadCondition(NewId, pGeometry)
{
}

// Constructor

SurfaceLoadCondition3D::SurfaceLoadCondition3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : BaseLoadCondition(NewId, pGeometry, pProperties)
{
}

//***********************************************************************************
//***********************************************************************************

Condition::Pointer SurfaceLoadCondition3D::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
{
    return boost::make_shared<SurfaceLoadCondition3D>(NewId, pGeom, pProperties);
}

Condition::Pointer SurfaceLoadCondition3D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return boost::make_shared<SurfaceLoadCondition3D>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

//***********************************************************************************
//***********************************************************************************
// Destructor

SurfaceLoadCondition3D::~SurfaceLoadCondition3D()
{
}

//***********************************************************************************
//***********************************************************************************


//***********************************************************************************
//***********************************************************************************

void SurfaceLoadCondition3D::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoadCondition3D::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

}


//***********************************************************************************
//***********************************************************************************
// --------- //
//  PRIVATE  //
// --------- //


//***********************************************************************************
//***********************************************************************************

void SurfaceLoadCondition3D::CalculateAndSubKp(
    Matrix& K,
    array_1d<double, 3 > & ge,
    array_1d<double, 3 > & gn,
    const Matrix& DN_De,
    const Vector& N,
    double pressure,
    double weight)
{
    KRATOS_TRY

    boost::numeric::ublas::bounded_matrix<double, 3, 3 > Kij;
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > Cross_ge;
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > Cross_gn;
    double coeff;
    unsigned int number_of_nodes = GetGeometry().size();

    MakeCrossMatrix(Cross_ge, ge);
    MakeCrossMatrix(Cross_gn, gn);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        unsigned int RowIndex = i * 3;
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
            unsigned int ColIndex = j * 3;

            coeff = pressure * N[i] * DN_De(j, 1) * weight;
            noalias(Kij) = coeff * Cross_ge;

            coeff = pressure * N[i] * DN_De(j, 0) * weight;

            noalias(Kij) -= coeff * Cross_gn;

            //TAKE CARE: the load correction matrix should be SUBTRACTED not added
            SubtractMatrix(K, Kij, RowIndex, ColIndex);
        }
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoadCondition3D::MakeCrossMatrix(
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > & M,
    array_1d<double, 3 > & U)
{
    M(0, 0) = 0.00;
    M(0, 1) = -U[2];
    M(0, 2) = U[1];
    M(1, 0) = U[2];
    M(1, 1) = 0.00;
    M(1, 2) = -U[0];
    M(2, 0) = -U[1];
    M(2, 1) = U[0];
    M(2, 2) = 0.00;
}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoadCondition3D::CrossProduct(
    array_1d<double, 3 > & cross,
    array_1d<double, 3 > & a,
    array_1d<double, 3 > & b)
{
    cross[0] = a[1] * b[2] - a[2] * b[1];
    cross[1] = a[2] * b[0] - a[0] * b[2];
    cross[2] = a[0] * b[1] - a[1] * b[0];
}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoadCondition3D::ExpandReducedMatrix(
    Matrix& Destination,
    Matrix& ReducedMatrix)
{
    KRATOS_TRY

    unsigned int size = ReducedMatrix.size2();

    for (unsigned int i = 0; i < size; i++)
    {
        int rowindex = i * 3;
        for (unsigned int j = 0; j < size; j++)
        {
            unsigned int colindex = j * 3;
            for (unsigned int ii = 0; ii < 3; ii++)
                Destination(rowindex + ii, colindex + ii) += ReducedMatrix(i, j);
        }
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoadCondition3D::SubtractMatrix(
    MatrixType& Destination,
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > & InputMatrix,
    int InitialRow,
    int InitialCol)
{
    KRATOS_TRY

    for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++)
            Destination(InitialRow + i, InitialCol + j) -= InputMatrix(i, j);

    KRATOS_CATCH("")
}




//***********************************************************************************
//***********************************************************************************

void SurfaceLoadCondition3D::CalculateAndAdd_PressureForce(
    VectorType& residualvector,
    const Vector& N,
    const array_1d<double, 3 > & v3,
    double pressure,
    double weight,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().size();

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        int index = 3 * i;
        double coeff = pressure * N[i] * weight;
        residualvector[index] += coeff * v3[0];
        residualvector[index + 1] += coeff * v3[1];
        residualvector[index + 2] += coeff * v3[2];
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoadCondition3D::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;
    //const unsigned int dim= 3;
    //resizing as needed the LHS
    if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
    {
        if (rLeftHandSideMatrix.size1() != MatSize)
            rLeftHandSideMatrix.resize(MatSize, MatSize, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize); //resetting LHS
    }

    //resizing as needed the RHS
    if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
    {
        if (rRightHandSideVector.size() != MatSize)
            rRightHandSideVector.resize(MatSize, false);
        rRightHandSideVector = ZeroVector(MatSize); //resetting RHS
    }

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = GetGeometry().ShapeFunctionsLocalGradients();
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

    //calculating actual jacobian
    GeometryType::JacobiansType J;
    J = GetGeometry().Jacobian(J);

    //vector with a loading applied to the elemnt
    array_1d<double, 3 > SurfaceLoad = ZeroVector(3);
    if( this->Has( SURFACE_LOAD ) )
        noalias(SurfaceLoad) = this->GetValue( SURFACE_LOAD );

    //pressure applied to the element itself
    double pressure_on_condition = 0.0;
    if( this->Has( NEGATIVE_FACE_PRESSURE ) )
        pressure_on_condition += this->GetValue( NEGATIVE_FACE_PRESSURE );
    if( this->Has( POSITIVE_FACE_PRESSURE ) )
        pressure_on_condition -= this->GetValue( POSITIVE_FACE_PRESSURE );



    //auxiliary terms
    array_1d<double, 3 > BodyForce;
    Vector PressureOnNodes(number_of_nodes, pressure_on_condition); //note that here we initialize from the value applied to the condition

    for (unsigned int i = 0; i < PressureOnNodes.size(); i++)
    {
        if( GetGeometry()[i].SolutionStepsDataHas( NEGATIVE_FACE_PRESSURE) )
            PressureOnNodes[i] += GetGeometry()[i].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE );
        if( GetGeometry()[i].SolutionStepsDataHas( POSITIVE_FACE_PRESSURE) )
            PressureOnNodes[i] -= GetGeometry()[i].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE );
    }

    array_1d<double, 3 > ge;
    array_1d<double, 3 > gn;
    array_1d<double, 3 > v3;
    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
    {
        double IntegrationWeight = GetGeometry().IntegrationPoints()[PointNumber].Weight();
        auto& N = row(Ncontainer, PointNumber);

        ge[0] = J[PointNumber](0, 0);
        gn[0] = J[PointNumber](0, 1);
        ge[1] = J[PointNumber](1, 0);
        gn[1] = J[PointNumber](1, 1);
        ge[2] = J[PointNumber](2, 0);
        gn[2] = J[PointNumber](2, 1);

        CrossProduct(v3, ge, gn);

        // calculating the pressure on the gauss point
        double pressure = 0.00;
        for (unsigned int ii = 0; ii < number_of_nodes; ii++)
            pressure += N[ii] * PressureOnNodes[ii];


        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true)
        {
            if (pressure != 0.00)
                CalculateAndSubKp(rLeftHandSideMatrix, ge, gn, DN_DeContainer[PointNumber], N, pressure, IntegrationWeight);
        }

        // RIGHT HAND SIDE VECTOR
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            if (pressure != 0.00)
                CalculateAndAdd_PressureForce(rRightHandSideVector, N, v3, pressure, IntegrationWeight, rCurrentProcessInfo);
        }

        //generic load on gauss point
        array_1d<double,3> gauss_load = SurfaceLoad;
        for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
        {
            if( GetGeometry()[ii].SolutionStepsDataHas( SURFACE_LOAD ) )
            {
                noalias(gauss_load) += N[ii]*GetGeometry()[ii].FastGetSolutionStepValue( SURFACE_LOAD );
            }
        }

        for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
        {
            unsigned int base = ii*3;
            for(unsigned int k=0; k<3; ++k)
                rRightHandSideVector[base+k] += N[ii]*gauss_load[k];
        }
    }


    KRATOS_CATCH("")
}


//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
