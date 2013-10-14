//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/surface_load_3D_condition.hpp"
#include "custom_utilities/sd_math_utils.hpp"

#include "solid_mechanics_application.h"

namespace Kratos
{

//***********************************************************************************
//***********************************************************************************
SurfaceLoad3DCondition::SurfaceLoad3DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//***********************************************************************************
//***********************************************************************************
SurfaceLoad3DCondition::SurfaceLoad3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
SurfaceLoad3DCondition::SurfaceLoad3DCondition( SurfaceLoad3DCondition const& rOther )
    : Condition(rOther)
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer SurfaceLoad3DCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new SurfaceLoad3DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//***********************************************************************************
//***********************************************************************************
SurfaceLoad3DCondition::~SurfaceLoad3DCondition()
{
}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int conditiondimension = number_of_nodes * 3;

    if (rResult.size() != conditiondimension)
        rResult.resize(conditiondimension);

    unsigned int index = 0;
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        index = i * 3;
        rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
    }
    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    ElementalDofList.resize(0);

    for (unsigned int i = 0; i < GetGeometry().size(); i++)
    {
        ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
    }
}

//***********************************************************************************
//***********************************************************************************


//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateConditionalSystem(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateConditionalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::MassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rMassMatrix.resize(0, 0, false);

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::DampMatrix(
    MatrixType& rDampMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rDampMatrix.resize(0, 0, false);

    KRATOS_CATCH("")
}


//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::GetValuesVector(
    Vector& values,
    int Step)
{
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;
    if (values.size() != MatSize)
        values.resize(MatSize, false);

    unsigned int index = 0;
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        index = i * 3;
        values[index] = disp[0];
        values[index + 1] = disp[1];
        values[index + 2] = disp[2];
    }
}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::GetFirstDerivativesVector(
    Vector& values,
    int Step)
{
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;
    if (values.size() != MatSize)
        values.resize(MatSize, false);
    unsigned int index = 0;
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        index = i * 3;
        values[index] = vel[0];
        values[index + 1] = vel[1];
        values[index + 2] = vel[2];
    }

}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::GetSecondDerivativesVector(
    Vector& values,
    int Step)
{
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;
    if (values.size() != MatSize)
        values.resize(MatSize, false);
    unsigned int index = 0;
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3 > & acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        index = i * 3;
        values[index] = acc[0];
        values[index + 1] = acc[1];
        values[index + 2] = acc[2];
    }
}

//***********************************************************************************
//***********************************************************************************
// --------- //
//  PRIVATE  //
// --------- //


//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::CalculateAndSubKp(
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
    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        RowIndex = i * 3;
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
            ColIndex = j * 3;

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

void SurfaceLoad3DCondition::MakeCrossMatrix(
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

void SurfaceLoad3DCondition::CrossProduct(
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

void SurfaceLoad3DCondition::ExpandReducedMatrix(
    Matrix& Destination,
    Matrix& ReducedMatrix)
{
    KRATOS_TRY

    unsigned int size = ReducedMatrix.size2();
    unsigned int rowindex = 0;
    unsigned int colindex = 0;

    for (unsigned int i = 0; i < size; i++)
    {
        rowindex = i * 3;
        for (unsigned int j = 0; j < size; j++)
        {
            colindex = j * 3;
            for (unsigned int ii = 0; ii < 3; ii++)
                Destination(rowindex + ii, colindex + ii) += ReducedMatrix(i, j);
        }
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::SubtractMatrix(
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

void SurfaceLoad3DCondition::CalculateAndAddFacePressure(
    VectorType& rF,
    const Vector& rN,
    const array_1d<double, 3 > & rNormal,
    double rPressure,
    double rIntegrationWeight)
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension  = 3;

    unsigned int index = 0;
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        array_1d<double, 3 > & ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(FORCE_EXTERNAL);
        index = dimension * i;
        double  DiscretePressure = rPressure * rN[i] * rIntegrationWeight;
        rF[index]     += DiscretePressure * rNormal[0];
        rF[index + 1] += DiscretePressure * rNormal[1];
        rF[index + 2] += DiscretePressure * rNormal[2];

	GetGeometry()[i].SetLock();
        ExternalForce[0] += DiscretePressure * rNormal[0];
        ExternalForce[1] += DiscretePressure * rNormal[1];
        ExternalForce[2] += DiscretePressure * rNormal[2];
	GetGeometry()[i].UnSetLock();
    }

    KRATOS_CATCH("")
}


//***********************************************************************
//***********************************************************************

void SurfaceLoad3DCondition::CalculateAndAddSurfaceLoad(Vector& rF,
        const Vector& rN,
        Vector& rForce,
        double  rIntegrationWeight )
{

    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension  = 3;
    unsigned int index = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        index = dimension * i;

        array_1d<double, 3 > & ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(FORCE_EXTERNAL);

	GetGeometry()[i].SetLock();
        for ( unsigned int idim = 0; idim < number_of_nodes; idim++ )
        {
            rF[index+idim] += rN[i] * rForce[idim] * rIntegrationWeight;

            ExternalForce[idim] += rN[i] * rForce[idim] * rIntegrationWeight;
        }
	GetGeometry()[i].UnSetLock();
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::CalculateConditionalSystem(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = 3;
    unsigned int MatSize = number_of_nodes * dimension;
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

    //PRESSURE CONDITION:
    Vector PressureOnNodes(number_of_nodes);
    PressureOnNodes = ZeroVector(number_of_nodes);

    double PressureCondition = 0;
    //PressureCondition = GetValue( PRESSURE );

    for (unsigned int i = 0; i < PressureOnNodes.size(); i++)
    {
        PressureOnNodes[i] = PressureCondition;
        PressureOnNodes[i]+= GetGeometry()[i].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE ) - GetGeometry()[i].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE );
    }


    //FORCE CONDITION:
    std::vector<array_1d<double,3> > ForceArray (number_of_nodes);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        ForceArray[i]=GetGeometry()[i].FastGetSolutionStepValue( FACE_LOAD );
    }

    //Plane definition:
    array_1d<double, 3 > ge; //tangent 1
    array_1d<double, 3 > gn; //tangent 2

    array_1d<double, 3 > NormalVector;

    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
    {
        double IntegrationWeight = GetGeometry().IntegrationPoints()[PointNumber].Weight();

        ge[0] = J[PointNumber](0, 0);
        gn[0] = J[PointNumber](0, 1);
        ge[1] = J[PointNumber](1, 0);
        gn[1] = J[PointNumber](1, 1);
        ge[2] = J[PointNumber](2, 0);
        gn[2] = J[PointNumber](2, 1);

        CrossProduct(NormalVector, ge, gn);
       // NormalVector /= norm_2(NormalVector);
//         KRATOS_WATCH(NormalVector);

        // calculating the pressure and force on the gauss point
        double gauss_pressure = 0.00;
        Vector ForceLoad = ZeroVector(dimension);

        for (unsigned int ii = 0; ii < number_of_nodes; ii++)
        {
            gauss_pressure += Ncontainer(PointNumber, ii) * PressureOnNodes[ii];
            for ( unsigned int j = 0; j < dimension; j++ )
            {
                ForceLoad[j] += Ncontainer( PointNumber, ii ) * ForceArray[ii][j];
            }
        }

        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true)
        {
            if (gauss_pressure != 0.00)
                CalculateAndSubKp(rLeftHandSideMatrix, ge, gn, DN_DeContainer[PointNumber], row(Ncontainer, PointNumber), gauss_pressure, IntegrationWeight);
        }

        // RIGHT HAND SIDE VECTOR
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            if (gauss_pressure != 0.00)
                CalculateAndAddFacePressure(rRightHandSideVector, row(Ncontainer, PointNumber),
                                            NormalVector, gauss_pressure, IntegrationWeight);

            if (norm_2(ForceLoad) != 0.00)
                CalculateAndAddSurfaceLoad(rRightHandSideVector, row(Ncontainer, PointNumber),
                                           ForceLoad, IntegrationWeight);
        }
    }


    KRATOS_CATCH("")
}

int SurfaceLoad3DCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
