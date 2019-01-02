// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Long Chen, Anna Rehr
//


// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "custom_elements/prestress_membrane_element.hpp"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

// Constructor
PrestressMembraneElement::PrestressMembraneElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{

}

// Constructor
PrestressMembraneElement::PrestressMembraneElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element( NewId, pGeometry, pProperties )
{

}

//***********************************************************************************
//***********************************************************************************

Element::Pointer PrestressMembraneElement::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties) const

{
    return Kratos::make_shared< PrestressMembraneElement >(NewId, GetGeometry().Create(rThisNodes), pProperties);
}

//***********************************************************************************
//***********************************************************************************

Element::Pointer PrestressMembraneElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const

{
    return Kratos::make_shared< PrestressMembraneElement >(NewId, pGeom, pProperties);
}

//***********************************************************************************
//***********************************************************************************
// Destructor
PrestressMembraneElement::~PrestressMembraneElement()
{
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)

{
  KRATOS_TRY;

  unsigned int num_nodes, local_size;
  unsigned int local_index = 0;

  num_nodes = GetGeometry().size();
  local_size = num_nodes * 3;

  const unsigned int d_pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

  if (rResult.size() != local_size)
      rResult.resize(local_size, false);

  for (unsigned int i_node = 0; i_node < num_nodes; ++i_node)
  {
      rResult[local_index++] = this->GetGeometry()[i_node].GetDof(DISPLACEMENT_X, d_pos).EquationId();
      rResult[local_index++] = this->GetGeometry()[i_node].GetDof(DISPLACEMENT_Y, d_pos + 1).EquationId();
      rResult[local_index++] = this->GetGeometry()[i_node].GetDof(DISPLACEMENT_Z, d_pos + 2).EquationId();
  }

  KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)

{
    unsigned int num_nodes, local_size;
    num_nodes = GetGeometry().size();
    local_size = num_nodes * 3;

    if (rElementalDofList.size() != local_size)
        rElementalDofList.resize(local_size);

    unsigned int local_index = 0;

    for (unsigned int i_node = 0; i_node < num_nodes; ++i_node)
    {
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DISPLACEMENT_X);
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DISPLACEMENT_Z);
    }
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::Initialize()

{
    KRATOS_TRY

    // reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    // define container sizes
    mDetJ0.resize(integration_points.size());
    mGab0.resize(integration_points.size());
    mGVector.resize(integration_points.size(), ZeroMatrix(2, 2));

    // compute base vectors in reference configuration, metrics
    ComputeBaseVectors(integration_points);

    // initialize prestress
    ComputePrestress(integration_points.size());

    // initialize formfinding
    if(this->Has(IS_FORMFINDING))
        InitializeFormfinding(integration_points.size());

    // Initialize Material
    InitializeMaterial(integration_points.size());

    KRATOS_CATCH( "" )
}



//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)

{
    //calculation flags
    const bool calculate_stiffness_matrix_flag = false;
    const bool calculate_residual_vector_flag = true;
    MatrixType temp = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, calculate_stiffness_matrix_flag, calculate_residual_vector_flag);
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)

{
    //calculation flags
    const bool calculate_stiffness_matrix_flag = true;
    const bool calculate_residual_vector_flag = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, calculate_stiffness_matrix_flag, calculate_residual_vector_flag);
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& Output,
    const ProcessInfo& rCurrentProcessInfo)

{
    KRATOS_ERROR << "CalculateOnIntegrationPoints not implemented yet" << std::endl;
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo)

{
    KRATOS_TRY

    // LUMPED MASS MATRIX
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int mat_size = number_of_nodes * 3;

    if (rMassMatrix.size1() != mat_size)
    {
        rMassMatrix.resize(mat_size, mat_size);
    }

    noalias(rMassMatrix) = ZeroMatrix(mat_size, mat_size);

    double total_mass = mTotalDomainInitialSize * GetProperties()[THICKNESS] * GetProperties()[DENSITY];

    Vector lump_fact;

    lump_fact = GetGeometry().LumpingFactors(lump_fact);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        double temp = lump_fact[i] * total_mass;

        for (unsigned int j = 0; j < 3; j++)
        {
            unsigned int index = i * 3 + j;
            rMassMatrix(index, index) = temp;
        }
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo)

{
    KRATOS_TRY

    // LUMPED DAMPING MATRIX

    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int mat_size = number_of_nodes * 3;

    if (rDampingMatrix.size1() != mat_size)
        rDampingMatrix.resize(mat_size, mat_size);

    rDampingMatrix = ZeroMatrix(mat_size, mat_size);

    double total_mass = mTotalDomainInitialSize * GetProperties()[THICKNESS] * GetProperties()[DENSITY];

    Vector lump_fact;

    lump_fact = GetGeometry().LumpingFactors(lump_fact);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        double temp = lump_fact[i] * total_mass;

        for (unsigned int j = 0; j < 3; j++)
        {
            unsigned int index = i * 3 + j;
            rDampingMatrix(index, index) = temp;
        }
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::FinalizeSolutionStep(
    ProcessInfo& rCurrentProcessInfo)
{
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
        mConstitutiveLawVector[i]->FinalizeSolutionStep(GetProperties(),
            GetGeometry(),
            row(GetGeometry().ShapeFunctionsValues(), i),
            rCurrentProcessInfo);
    }
}

//##### From here, old code

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::GetValuesVector(
    Vector& rValues,
    int Step)

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int mat_size = number_of_nodes * 3;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3>& disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        unsigned int index = i * 3;
        rValues[index] = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];
    }
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::GetFirstDerivativesVector(
    Vector& rValues,
    int Step)

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int mat_size = number_of_nodes * 3;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        unsigned int index = i * 3;
        rValues[index] = vel[0];
        rValues[index + 1] = vel[1];
        rValues[index + 2] = vel[2];
    }

}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::GetSecondDerivativesVector(
    Vector& rValues,
    int Step)

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int mat_size = number_of_nodes * 3;

    if (rValues.size() != mat_size)
        rValues.resize(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3>& acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        unsigned int index = i * 3;
        rValues[index] = acc[0];
        rValues[index + 1] = acc[1];
        rValues[index + 2] = acc[2];
    }
}

//***********************************************************************************
//***********************************************************************************
// --------- //
//  PRIVATE  //
// --------- //
//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateAndAddKm(
    Matrix& rK,
    Matrix& rB,
    Matrix& rD,
    const double& rWeight)

{
    KRATOS_TRY

    unsigned int dim = rB.size2();
    Matrix temp(3, dim);
    noalias(temp) = prod(rD, rB);
    temp *= rWeight;
    Matrix Km(dim, dim);
    noalias(Km) = prod(trans(rB), temp);
    noalias(rK) += Km;

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateAndAddNonlinearKm(
    Matrix& rK,
    Matrix& rB11,
    Matrix& rB22,
    Matrix& rB12,
    Vector& rSD,
    const double& rWeight)

{
    KRATOS_TRY;

    const unsigned int number_of_nodes = GetGeometry().size();

    for (unsigned int n = 0; n < number_of_nodes; n++)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            for (unsigned int m = 0; m <= n; m++)
            {
                unsigned int check = 3;
                if (m == n)
                    check = i + 1;
                for (unsigned int j = 0; j < check; j++)
                {
                    rK(3 * n + i, 3 * m + j) += (rSD[0] * rB11(3 * n + i, 3 * m + j) + rSD[1] * rB22(3 * n + i, 3 * m + j) + rSD[2] * rB12(3 * n + i, 3 * m + j))*rWeight;
                    rK(3 * m + i, 3 * n + j) = rK(3 * n + i, 3 * m + j);
                }
            }
        }
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateQ(
    BoundedMatrix<double, 3, 3>& Q,
    const unsigned int rPointNumber)

{
    KRATOS_TRY
    Q(0, 0) = std::pow(mGVector[rPointNumber](0, 0), 2);
    Q(0, 1) = std::pow(mGVector[rPointNumber](0, 1), 2);
    Q(0, 2) = 2.00*mGVector[rPointNumber](0, 0)*mGVector[rPointNumber](0, 1);

    Q(1, 0) = std::pow(mGVector[rPointNumber](1, 0), 2);
    Q(1, 1) = std::pow(mGVector[rPointNumber](1, 1), 2);
    Q(1, 2) = 2.00*mGVector[rPointNumber](1, 0) * mGVector[rPointNumber](1, 1);

    Q(2, 0) = 2.00 * mGVector[rPointNumber](0, 0) * mGVector[rPointNumber](1, 0);
    Q(2, 1) = 2.00 * mGVector[rPointNumber](0, 1)*mGVector[rPointNumber](1, 1);
    Q(2, 2) = 2.00 * (mGVector[rPointNumber](0, 0) * mGVector[rPointNumber](1, 1) + mGVector[rPointNumber](0, 1)*mGVector[rPointNumber](1, 0));

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateB(
    Matrix& rB,
    const BoundedMatrix<double, 3, 3>& rQ,
    const Matrix& DN_De,
    const array_1d<double, 3>& g1,
    const array_1d<double, 3>& g2)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    Matrix b(3, number_of_nodes * 3);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        unsigned int index = 3 * i;

        //first line
        b(0, index) = DN_De(i, 0) * g1[0];
        b(0, index + 1) = DN_De(i, 0) * g1[1];
        b(0, index + 2) = DN_De(i, 0) * g1[2];

        //second line
        b(1, index) = DN_De(i, 1) * g2[0];
        b(1, index + 1) = DN_De(i, 1) * g2[1];
        b(1, index + 2) = DN_De(i, 1) * g2[2];

        //third line
        b(2, index) = 0.5*(DN_De(i, 1) * g1[0] + DN_De(i, 0) * g2[0]);
        b(2, index + 1) = 0.5*(DN_De(i, 1) * g1[1] + DN_De(i, 0) * g2[1]);
        b(2, index + 2) = 0.5*(DN_De(i, 1) * g1[2] + DN_De(i, 0) * g2[2]);
    }

    rB = prod(rQ, b);

    KRATOS_CATCH("")
}


//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateStrain(
    Vector& rStrainVector,
    array_1d<double, 3>& rgab,
    array_1d<double, 3>& rGab)

{
    KRATOS_TRY

    rStrainVector[0] = 0.5 * (rgab[0] - rGab[0]);
    rStrainVector[1] = 0.5 * (rgab[1] - rGab[1]);
    rStrainVector[2] = 0.5 * (rgab[2] - rGab[2]);

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************
void PrestressMembraneElement::CalculateAndAdd_BodyForce(
    const Vector& rN,
    const ProcessInfo& rCurrentProcessInfo,
    array_1d<double, 3>& BodyForce,
    VectorType& rRightHandSideVector,
    const double& rWeight)

{
    KRATOS_TRY

        unsigned int number_of_nodes = this->GetGeometry().size();
    const double density = this->GetProperties()[DENSITY];


    noalias(BodyForce) = ZeroVector(3);
    for (unsigned int i = 0; i<number_of_nodes; i++)
        BodyForce += rN[i] * this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    BodyForce *= density;

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        int index = 3 * i;

        for (unsigned int j = 0; j < 3; j++)
            rRightHandSideVector[index + j] += rWeight * rN[i] * BodyForce[j];
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************
void PrestressMembraneElement::CalculateAndAdd_PressureForce(
    VectorType& rResidualVector,
    const Vector& rN,
    const array_1d<double, 3>& rv3,
    const double& rPressure,
    const double& rWeight,
    const ProcessInfo& rCurrentProcessInfo)

{
    KRATOS_TRY

        unsigned int number_of_nodes = GetGeometry().size();

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        int index = 3 * i;
        double coeff = rPressure * rN[i] * rWeight;
        rResidualVector[index] += coeff * rv3[0];
        rResidualVector[index + 1] += coeff * rv3[1];
        rResidualVector[index + 2] += coeff * rv3[2];
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool rCalculateStiffnessMatrixFlag,
    const bool rCalculateResidualVectorFlag)

{
    KRATOS_TRY

    // Initializing all needed variables
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int mat_size = number_of_nodes * 3;

    Vector strain_vector(3);
    Vector stress_vector(3);

    // set up Constitutive Law
    ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags &constitutive_law_options=Values.GetOptions();
    constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    // for formfinding: Constitutive Tensor is 0 Tensor
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    if(this->Has(IS_FORMFINDING)){
        if(this->GetValue(IS_FORMFINDING))
            constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    }

    Values.SetStrainVector(strain_vector);       // this is the input parameter
    Values.SetStressVector(stress_vector);       // this is an output parameter

    // resizing as needed the LHS
    if (rCalculateStiffnessMatrixFlag == true)    // calculation of the matrix is required
    {
        if (rLeftHandSideMatrix.size1() != mat_size)
        {
            rLeftHandSideMatrix.resize(mat_size, mat_size);
        }

        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
    }

    //resizing as needed the RHS
    if (rCalculateResidualVectorFlag == true) //calculation of the matrix is required
    {
        if (rRightHandSideVector.size() != mat_size)
        {
            rRightHandSideVector.resize(mat_size);
        }

        rRightHandSideVector = ZeroVector(mat_size); //resetting RHS
    }

    // reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    const GeometryType::ShapeFunctionsGradientsType& DN_De_container = GetGeometry().ShapeFunctionsLocalGradients();

    // calculating actual jacobian
    GeometryType::JacobiansType J;

    J = GetGeometry().Jacobian(J);

    for (unsigned int point_number = 0; point_number < integration_points.size(); point_number++)
    {
        // reading integration weight, shape function value and its gradients at this integration point
        const double integration_weight = integration_points[point_number].Weight();

        Matrix DN_De = DN_De_container[point_number];

        // covariant metric in deformed system
        array_1d<double, 3> gab;

        // Transformation Matrix Q
        BoundedMatrix<double, 3, 3> Q;
        noalias(Q) = ZeroMatrix(3, 3);
        // basis vectors in deformed system
        array_1d<double, 3> g1, g2, g3;

        CalculateQ(Q, point_number);
        CalculateMetricDeformed(point_number, DN_De, gab, g1, g2);
        CalculateStrain(strain_vector, gab, mGab0[point_number]);

        Vector cartesian_strain_vector = prod(Q, strain_vector); //in refence configuration

        // Constitutive Matrices D
        Matrix D(3, 3, 0);
        Values.SetConstitutiveMatrix(D); // this is an output parameter

        //Deformation Gradient
        Values.SetDeformationGradientF(CalculateDeformationGradient(point_number));

        mConstitutiveLawVector[point_number]->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

        // Deformations for Non-linear force vector
        Vector strain_deformation;

        strain_deformation = prod(trans(D), cartesian_strain_vector);

        // Adding the pre-stress values as forces over length
        for (int i = 0; i < 3; ++i)
            strain_deformation[i] += GetValue(MEMBRANE_PRESTRESS)(i,point_number);

        // calculate B matrices
        Matrix B(3, mat_size);
        noalias(B) = ZeroMatrix(3, mat_size);
        CalculateB(B, Q, DN_De, g1, g2);

        // integration on the REFERENCE CONFIGURATION
        double DetJ0 = mDetJ0[point_number];
        double int_reference_weight = integration_weight * DetJ0 * GetProperties()[THICKNESS];

        // Nonlinear Deformation
        Matrix strain_local_cart_11(number_of_nodes * 3, number_of_nodes * 3);
        noalias(strain_local_cart_11) = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);
        Matrix strain_local_cart_22(number_of_nodes * 3, number_of_nodes * 3);
        noalias(strain_local_cart_22) = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);
        Matrix strain_local_cart_12(number_of_nodes * 3, number_of_nodes * 3);
        noalias(strain_local_cart_12) = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);

        CalculateSecondVariationStrain(DN_De, strain_local_cart_11, strain_local_cart_22, strain_local_cart_12, Q);

        // LEFT HAND SIDE MATRIX
        if (rCalculateStiffnessMatrixFlag == true)
        {
            // adding membrane contribution to the stiffness matrix
            CalculateAndAddKm(rLeftHandSideMatrix, B, D, int_reference_weight);

            // Adding non-linear-contribution to stiffness matrix
            CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
                strain_local_cart_11, strain_local_cart_22, strain_local_cart_12,
                strain_deformation,
                int_reference_weight);
        }

        // RIGHT HAND SIDE VECTOR
        if (rCalculateResidualVectorFlag == true)         // calculation of the matrix is required
        {
            // operation performed: rRighthandSideVector -= Weight* IntForce
            noalias(rRightHandSideVector) -= int_reference_weight* prod(trans(B), strain_deformation);
        }
    } // end loop over integration points

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************
void PrestressMembraneElement::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}
//***********************************************************************************
//***********************************************************************************
void PrestressMembraneElement::CalculateMetricDeformed(const unsigned int rPointNumber, Matrix DN_De,
    array_1d<double, 3>& rgab,
    array_1d<double, 3>& rg1,
    array_1d<double, 3>& rg2)
{
    GeometryType::JacobiansType J_act;
    J_act = GetGeometry().Jacobian(J_act);

    array_1d<double, 3> g3;

    rg1[0] = J_act[rPointNumber](0, 0);
    rg2[0] = J_act[rPointNumber](0, 1);
    rg1[1] = J_act[rPointNumber](1, 0);
    rg2[1] = J_act[rPointNumber](1, 1);
    rg1[2] = J_act[rPointNumber](2, 0);
    rg2[2] = J_act[rPointNumber](2, 1);

    //basis vector g3
    MathUtils<double>::CrossProduct(g3, rg1, rg2);
    //differential area dA
    double dA = norm_2(g3);
    //normal vector _n
    array_1d<double, 3> n = g3 / dA;

    //GetCovariantMetric
    rgab[0] = inner_prod(rg1,rg1);
    rgab[1] = inner_prod(rg2,rg2);
    rgab[2] = inner_prod(rg1,rg2);

}

//***********************************************************************************
//***********************************************************************************
void PrestressMembraneElement::CalculateSecondVariationStrain(Matrix DN_De,
    Matrix & rStrainLocalCart11,
    Matrix & rStrainLocalCart22,
    Matrix & rStrainLocalCart12,
    BoundedMatrix<double, 3, 3>& rQ)
{
    const unsigned int number_of_nodes = GetGeometry().size();

    Vector dd_strain_curvilinear(3);

    for (unsigned int n = 0; n < number_of_nodes; ++n)
    {
        for (unsigned int i = 0; i < 3; ++i)
        {
            for (unsigned int m = 0; m <= n; m++)
            {

                unsigned int limit = i + 1;
                if (m < n)
                    limit = 3;
                for (unsigned int j = 0; j < limit; j++)
                {
                    noalias(dd_strain_curvilinear) = ZeroVector(3);
                    if (j == i)
                    {
                        dd_strain_curvilinear[0] = DN_De(n, 0)*DN_De(m, 0);
                        dd_strain_curvilinear[1] = DN_De(n, 1)*DN_De(m, 1);
                        dd_strain_curvilinear[2] = 0.5*(DN_De(n, 0)*DN_De(m, 1) + DN_De(n, 1)*DN_De(m, 0));

                        rStrainLocalCart11(3 * n + i, 3 * m + j) = rQ(0, 0)*dd_strain_curvilinear[0] + rQ(0, 1)*dd_strain_curvilinear[1] + rQ(0, 2)*dd_strain_curvilinear[2];
                        rStrainLocalCart22(3 * n + i, 3 * m + j) = rQ(1, 0)*dd_strain_curvilinear[0] + rQ(1, 1)*dd_strain_curvilinear[1] + rQ(1, 2)*dd_strain_curvilinear[2];
                        rStrainLocalCart12(3 * n + i, 3 * m + j) = rQ(2, 0)*dd_strain_curvilinear[0] + rQ(2, 1)*dd_strain_curvilinear[1] + rQ(2, 2)*dd_strain_curvilinear[2];

                    }
                }
            }
        }
    }
}

//***********************************************************************************
//***********************************************************************************
void PrestressMembraneElement::ProjectPrestress(
    const unsigned int rPointNumber){

    // Compute local cartesian basis E1, E2, E3
    array_1d<double,3> E1 = column( GetValue(BASE_REF_1),rPointNumber )/norm_2(column( GetValue(BASE_REF_1),rPointNumber ));
    array_1d<double,3> E2, E3;
    MathUtils<double>::CrossProduct(E3,E1,column( GetValue(BASE_REF_2),rPointNumber ));
    E3 /= norm_2(E3);
    MathUtils<double>::CrossProduct(E2,E3,E1);
    E2 /= norm_2(E2);

    // definition in-plane prestress vectors
    array_1d<double, 3> T1, T2, T3;

    // Alternative 1: planar projection according to dissertation Roland Wuechner
    const std::string this_proyection_type = GetProperties().Has(PROJECTION_TYPE_COMBO) ? GetProperties()[PROJECTION_TYPE_COMBO] : "planar";
    if(this_proyection_type == "planar"){
        // definition prestress axes
        array_1d<double,3> global_prestress_axis1, global_prestress_axis2, global_prestress_axis3;
        for(unsigned int i=0; i<3;i++){
            global_prestress_axis1[i] = GetValue(PRESTRESS_AXIS_1)(i,rPointNumber);
            global_prestress_axis2[i] = GetValue(PRESTRESS_AXIS_2)(i,rPointNumber);
        }

        KRATOS_DEBUG_ERROR_IF(norm_2(global_prestress_axis1) < std::numeric_limits<double>::epsilon() && norm_2(global_prestress_axis1) > -std::numeric_limits<double>::epsilon()) << "division by zero!" << std::endl;
        KRATOS_DEBUG_ERROR_IF(norm_2(global_prestress_axis2) < std::numeric_limits<double>::epsilon() && norm_2(global_prestress_axis2) > -std::numeric_limits<double>::epsilon()) << "division by zero!" << std::endl;

        // normalization global prestress axes
        global_prestress_axis1 /= norm_2(global_prestress_axis1);
        global_prestress_axis2 /= norm_2(global_prestress_axis2);

        // Compute global_prestress_axis3
        MathUtils<double>::CrossProduct(global_prestress_axis3, global_prestress_axis1, global_prestress_axis2);
        global_prestress_axis3 /= norm_2(global_prestress_axis3);

        // Compute T1, T2, T3
        array_1d<double, 3> normal_projection_plane;
        MathUtils<double>::CrossProduct(normal_projection_plane, global_prestress_axis3, global_prestress_axis1);
        MathUtils<double>::CrossProduct(T1, normal_projection_plane, E3);
        T1 /= norm_2(T1);
        MathUtils<double>::CrossProduct(T2,E3,T1);
        T2 /= norm_2(T2);
        T3 = E3;
    } else if(this_proyection_type == "radial"){ // Alternative 2: radial projection (T1 = radial direction, T3 = out-of-plane direction, T2=T3xT1)
        // definition prestress axes
        array_1d<double,3> global_prestress_axis1; // = direction of rotational axis
        for(unsigned int i=0; i<3;i++)
            global_prestress_axis1[i] = GetValue(PRESTRESS_AXIS_1)(i,rPointNumber);

        KRATOS_DEBUG_ERROR_IF(norm_2(global_prestress_axis1) < std::numeric_limits<double>::epsilon() && norm_2(global_prestress_axis1) > -std::numeric_limits<double>::epsilon()) << "division by zero!" << std::endl;

        // compute T1, T2, T3
        MathUtils<double>::CrossProduct(T1, global_prestress_axis1, E3);
        T1 /= norm_2(T1);
        MathUtils<double>::CrossProduct(T2, E3, T1);
        T2 /= norm_2(T2);
        T3 = E3;
    }

    // Transform prestresses in the local cartesian cosy in reference configuration
    BoundedMatrix<double,3,3> origin, target, tensor;
    noalias(origin) = ZeroMatrix(3,3);
    noalias(target) = ZeroMatrix(3,3);
    noalias(tensor) = ZeroMatrix(3,3);

    for (int i=0;i<3;i++){
        origin(i,0) = T1(i);
        origin(i,1) = T2(i);
        origin(i,2) = T3(i);
        target(i,0) = E1(i);
        target(i,1) = E2(i);
        target(i,2) = E3(i);
    }

    tensor.clear();
    tensor(0,0) = this->GetValue(MEMBRANE_PRESTRESS)(0,rPointNumber);
    tensor(1,1) = this->GetValue(MEMBRANE_PRESTRESS)(1,rPointNumber);
    tensor(1,0) = this->GetValue(MEMBRANE_PRESTRESS)(2,rPointNumber);
    tensor(0,1) = this->GetValue(MEMBRANE_PRESTRESS)(2,rPointNumber);

    // Transformation Stress Tensor
    StructuralMechanicsMathUtilities::TensorTransformation<3>(origin,origin,target,target,tensor);

    // store prestress values in the elemental data
    Matrix& prestress_matrix = GetValue(MEMBRANE_PRESTRESS);
    prestress_matrix(0,rPointNumber) = tensor(0,0);
    prestress_matrix(1,rPointNumber) = tensor(1,1);
    prestress_matrix(2,rPointNumber) = tensor(1,0);
}
//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo){

}
//***********************************************************************************
//***********************************************************************************
void PrestressMembraneElement::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo){
    // for formfinding: update basevectors (and prestress in case of anisotropy)
    if(this->Has(IS_FORMFINDING)){
        if(this->GetValue(IS_FORMFINDING)){
            const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
            if(mAnisotropicPrestress == true){
                // update prestress in the anisotropic case
                for (unsigned int point_number = 0; point_number < integration_points.size(); ++point_number){
                    if(mAnisotropicPrestress == true)
                        UpdatePrestress(point_number);
                }
            }
            //update base vectors in reference configuration, metrics
            ComputeBaseVectors(integration_points);
        }
    }
}


//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::InitializeMaterial(const unsigned int NumberIntegrationPoints){
    if (mConstitutiveLawVector.size() == 0)
    {
        mConstitutiveLawVector.resize(NumberIntegrationPoints);

        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
        {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();

            mConstitutiveLawVector[i]->InitializeMaterial(GetProperties(), GetGeometry(), row(GetGeometry().ShapeFunctionsValues(), i));
        }
    }
}

//***********************************************************************************
//***********************************************************************************
void PrestressMembraneElement::UpdatePrestress(const unsigned int rPointNumber){
    // --1--computation relevant CoSys
    array_1d<double, 3> g1, g2, g3, gab; //base vectors/metric actual config
    array_1d<double, 3> G3;  // base vector reference config
    array_1d<double, 3> E1_tot,E2_tot, E3_tot; //local cartesian base vectors initial reference config
    array_1d<double, 3> E1, E2, E3; // local cartesian base vectors initial config
    array_1d<double, 3> base_ref_contra_tot_1, base_ref_contra_tot_2; //contravariant base vectors initial reference config

    ComputeRelevantCoSys(rPointNumber, g1, g2, g3, gab, G3, E1_tot, E2_tot, E3_tot, E1, E2, E3, base_ref_contra_tot_1, base_ref_contra_tot_2);

    // --2-- computation total deformation gradient
    BoundedMatrix<double,3,3> deformation_gradient_total;
    for(unsigned int i=0; i<3; i++){
        for(unsigned int j=0; j<3; j++){
            deformation_gradient_total(i,j) = base_ref_contra_tot_1(j)*g1(i) + base_ref_contra_tot_2(j)*g2(i) + mG3Initial[rPointNumber](j)*g3(i);
        }
    }

    //--3--Compute the eigenvalues of the total deformation gradient
    BoundedMatrix<double,3,3> origin, target, tensor;
    noalias(origin) = ZeroMatrix(3,3);
    noalias(target) = ZeroMatrix(3,3);
    noalias(tensor) = ZeroMatrix(3,3);
    double lambda_1, lambda_2;

    ComputeEigenvaluesDeformationGradient(rPointNumber,
                    origin, target, tensor,
                    base_ref_contra_tot_1,base_ref_contra_tot_2,
                    E1_tot, E2_tot, E3_tot,
                    gab,
                    lambda_1, lambda_2);

    //--4--Compute the eigenvectors in the reference and actual configuration
    BoundedMatrix<double,3,3> N_act; // eigenvectors in actual configuration
    noalias(N_act) = ZeroMatrix(3,3);
    ComputeEigenvectorsDeformationGradient(rPointNumber,
                                tensor, origin,
                                deformation_gradient_total,
                                E1_tot, E2_tot,
                                lambda_1, lambda_2,
                                N_act);

    //--5--Compute the modified prestress
    ModifyPrestress(rPointNumber,
                    origin, target,tensor,
                    E1, E2, E3, G3,
                    g1, g2, g3, N_act,
                    lambda_1, lambda_2);
}
//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::ComputeRelevantCoSys(const unsigned int rPointNumber,
             array_1d<double, 3>& rg1,array_1d<double, 3>& rg2,array_1d<double, 3>& rg3, array_1d<double, 3>& rgab,
             array_1d<double, 3>& rG3,
             array_1d<double, 3>& rE1Tot, array_1d<double, 3>& rE2Tot,array_1d<double, 3>& rE3Tot,
             array_1d<double, 3>& rE1,array_1d<double, 3>& rE2,array_1d<double, 3>& rE3,
             array_1d<double, 3>& rBaseRefContraTot1,array_1d<double, 3>& rBaseRefContraTot2){

    const GeometryType::ShapeFunctionsGradientsType& DN_De_container = GetGeometry().ShapeFunctionsLocalGradients();
    Matrix DN_De = DN_De_container[rPointNumber];
    CalculateMetricDeformed(rPointNumber, DN_De, rgab, rg1, rg2);

    // compute the out-of-plane direction and normalize it
    MathUtils<double>::CrossProduct(rg3,rg1,rg2);
    rg3 /= norm_2(rg3);

    // computation of the out-of-plane direction in the reference configuration
    MathUtils<double>::CrossProduct(rG3,column( GetValue(BASE_REF_1),rPointNumber ),column( GetValue(BASE_REF_2),rPointNumber ));
    rG3 /= norm_2(rG3);

    //Compute cartesian basis in total reference configuration
    // rE1Tot = mG1Initial/|mG1Initial|, rE3Tot = mG3Initial
    rE1Tot = mG1Initial[rPointNumber]/norm_2(mG1Initial[rPointNumber]);
    rE3Tot = mG3Initial[rPointNumber];
    MathUtils<double>::CrossProduct(rE2Tot,rE3Tot,rE1Tot);
    rE2Tot /= norm_2(rE2Tot);

    //Compute cartesian basis in (updated) reference configuration
    // rE1 = BASE_REF_1/|BASE_REF_1|, rE3 = rG3
    rE1 = column( GetValue(BASE_REF_1),rPointNumber )/norm_2(column( GetValue(BASE_REF_1),rPointNumber));
    rE3 = rG3;
    MathUtils<double>::CrossProduct(rE2,rE3,rE1);
    rE2 /= norm_2(rE2);

    // Compute contravariant basis in total Reference configuration
    // -->Compute the metric in total reference configuration
    array_1d<double, 3> metric_reference_tot;
    metric_reference_tot(0)=inner_prod(mG1Initial[rPointNumber],mG1Initial[rPointNumber]);
    metric_reference_tot(1)=inner_prod(mG2Initial[rPointNumber],mG2Initial[rPointNumber]);
    metric_reference_tot(2)=inner_prod(mG2Initial[rPointNumber],mG1Initial[rPointNumber]);

    // -->Invert metric in (total) reference configuration
    double det_metric_tot = metric_reference_tot(0)*metric_reference_tot(1)-metric_reference_tot(2)*metric_reference_tot(2);
    array_1d<double, 3> inv_metric_tot;
    KRATOS_DEBUG_ERROR_IF(det_metric_tot < std::numeric_limits<double>::epsilon() && det_metric_tot > -std::numeric_limits<double>::epsilon()) << "division by zero!" << std::endl;
    inv_metric_tot(0) = 1.0/det_metric_tot * metric_reference_tot(1);
    inv_metric_tot(1) = 1.0/det_metric_tot * metric_reference_tot(0);
    inv_metric_tot(2) = -1.0/det_metric_tot * metric_reference_tot(2);

    // -->Compute contravariant basis (in total reference configuration)
    rBaseRefContraTot1 = inv_metric_tot(0)*mG1Initial[rPointNumber] + inv_metric_tot(2)*mG2Initial[rPointNumber];
    rBaseRefContraTot2 = inv_metric_tot(2)*mG1Initial[rPointNumber] + inv_metric_tot(1)*mG2Initial[rPointNumber];

}
//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::ComputeEigenvaluesDeformationGradient(const unsigned int rPointNumber,
                    BoundedMatrix<double,3,3>& rOrigin, BoundedMatrix<double,3,3>& rTarget, BoundedMatrix<double,3,3>& rTensor,
                    const array_1d<double, 3>& rBaseRefContraTot1, const array_1d<double, 3>& rBaseRefContraTot2,
                    const array_1d<double, 3>& rE1Tot, const array_1d<double, 3>& rE2Tot, const array_1d<double, 3>& rE3Tot,
                    const array_1d<double, 3>& rgab,
                    double& rLambda1, double& rLambda2){
    for (int i=0;i<3;i++){
        rOrigin(i,0) = rBaseRefContraTot1(i);
        rOrigin(i,1) = rBaseRefContraTot2(i);
        rOrigin(i,2) = mG3Initial[rPointNumber](i);
        rTarget(i,0) = rE1Tot(i);
        rTarget(i,1) = rE2Tot(i);
        rTarget(i,2) = rE3Tot(i);
    }
    rTensor.clear();
    rTensor(0,0) = rgab(0);
    rTensor(1,1) = rgab(1);
    rTensor(1,0) = rgab(2);
    rTensor(0,1) = rgab(2);

    // Transformation C rTensor
    StructuralMechanicsMathUtilities::TensorTransformation<3>(rOrigin,rOrigin,rTarget,rTarget,rTensor);

    // Compute the Eigenvalues
    // (C-lambda_i*I)N=0
    // rTensor:Eigenvalue in out-of-plane direction lambda_3 = 1 (no strain in this direction)

    //solution quadratic formula
    const double root = std::sqrt(std::abs((rTensor(0,0)+rTensor(1,1))*(rTensor(0,0)+rTensor(1,1))-4.0*(rTensor(0,0)*rTensor(1,1)-rTensor(0,1)*rTensor(1,0))));
    rLambda1 = (rTensor(0,0) + rTensor(1,1) + root)/2.0;
    rLambda2 = (rTensor(0,0) + rTensor(1,1) - root)/2.0;

    // Eigenvectors of C: lambda^2 --> root needed for "real" eigenvalues
    rLambda1 = std::sqrt(rLambda1);
    rLambda2 = std::sqrt(rLambda2);
}
//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::ComputeEigenvectorsDeformationGradient(const unsigned int rPointNumber,
                                BoundedMatrix<double,3,3>& rTensor, BoundedMatrix<double,3,3>& rOrigin,
                                const BoundedMatrix<double,3,3>& rDeformationGradientTotal,
                                const array_1d<double, 3>& rE1Tot, const array_1d<double, 3>& rE2Tot,
                                const double Lambda1, const double Lambda2,
                                BoundedMatrix<double,3,3>& rNAct){
    //Check if rTensor is the Identity Matrix -> no deformation
    //(rTensor-lambda^2*I)=ZeroMatrix
    bool principal_strain_state = false;
    rTensor(0,0) -= Lambda1*Lambda1;
    rTensor(1,1) -= Lambda2*Lambda2;
    if (std::abs(rTensor(0,0))<1.0e-6 && std::abs(rTensor(1,0))<1.0e-6 && std::abs(rTensor(0,1))<1.0e-6 && std::abs(rTensor(1,1))<1.0e-6)
        principal_strain_state = true;

    // compute C (=rTensor)
    noalias(rTensor) = prod(trans(rDeformationGradientTotal),rDeformationGradientTotal);

    // Eigenvectors in reference configuration: N_ref
    array_1d<double, 3> N_ref_1, N_ref_2, N_ref_3;
    N_ref_3 = mG3Initial[rPointNumber];

    if (principal_strain_state) {
        N_ref_1 = rE1Tot;
        N_ref_2 = rE2Tot;
    }
    else {
        //Compute the first principal direction
        //Eigenvector lies in the plan spanned up by G1 and G2
        //-> (C-lambda^2*I)N =(C-lambda^2*I)(alpha1*G1+alpha2*G2)=B1*alpha1+B2*alpha2 =0

        rTensor(0,0) -= Lambda1*Lambda1;
        rTensor(1,1) -= Lambda1*Lambda1;
        rTensor(2,2) -= Lambda1*Lambda1;

        column(rOrigin,0) = mG1Initial[rPointNumber];
        column(rOrigin,1) = mG2Initial[rPointNumber];

        BoundedMatrix<double,3,3> B;
        noalias(B) = prec_prod(rTensor, rOrigin);

        // compute alpha1 and alpha2
        double alpha_1,alpha_2;
        unsigned int imax=0;  unsigned int jmax=0;
        for (unsigned int i=0;i<3;i++){
            for (unsigned int j=0;j<2;j++){
                if(std::abs(B(i,j))>std::abs(B(imax,jmax))){
                    imax=i;
                    jmax=j;
                }
            }
        }
        if (jmax==0){
          alpha_2 = 1.0;
          alpha_1 = -B(imax,1)/B(imax,0);
        }
        else{
          alpha_1 = 1.0;
          alpha_2 = -B(imax,0)/B(imax,1);
        }
        N_ref_1 = alpha_1*mG1Initial[rPointNumber] + alpha_2*mG2Initial[rPointNumber];
        MathUtils<double>::CrossProduct(N_ref_2,N_ref_3,N_ref_1);
    }
    //Eigenvectors in the actual configuration from n=F*N
    for (unsigned int i=0;i<3;i++){
            for (unsigned int j=0;j<3;j++){
                rNAct(i,0) += rDeformationGradientTotal(i,j)*N_ref_1(j);
                rNAct(i,1) += rDeformationGradientTotal(i,j)*N_ref_2(j);
                rNAct(i,2) += rDeformationGradientTotal(i,j)*N_ref_3(j);
            }
        }

    // normalize eigenvectors
    double length;
    for (unsigned int j=0;j<3;j++){
        length = norm_2(column(rNAct,j));
        for(unsigned int i=0;i<3;i++)
            rNAct(i,j) /= length;
        }
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::ModifyPrestress(const unsigned int rPointNumber,
                    BoundedMatrix<double,3,3>& rOrigin, BoundedMatrix<double,3,3>& rTarget,BoundedMatrix<double,3,3>& rTensor,
                    const array_1d<double, 3>& rE1, const array_1d<double, 3>& rE2, const array_1d<double, 3>& rE3, const array_1d<double, 3>& rG3,
                    const array_1d<double, 3>& rg1, const array_1d<double, 3>& rg2, const array_1d<double, 3>& rg3, const BoundedMatrix<double,3,3>& rNAct,
                    const double Lambda1, const double Lambda2){
    // Transform prestresses from the local cosy in reference config to the curvilinear cosy in reference config, covariant
    for (int i=0;i<3;i++){
        rOrigin(i,0) = rE1(i);
        rOrigin(i,1) = rE2(i);
        rOrigin(i,2) = rE3(i);
        rTarget(i,0) = GetValue(BASE_REF_1)(i,rPointNumber);
        rTarget(i,1) = GetValue(BASE_REF_2)(i,rPointNumber);
        rTarget(i,2) = rG3(i);
    }
    rTensor.clear();
    rTensor(0,0) = GetValue(MEMBRANE_PRESTRESS)(0,rPointNumber);
    rTensor(1,1) = GetValue(MEMBRANE_PRESTRESS)(1,rPointNumber);
    rTensor(1,0) = GetValue(MEMBRANE_PRESTRESS)(2,rPointNumber);
    rTensor(0,1) = GetValue(MEMBRANE_PRESTRESS)(2,rPointNumber);

    StructuralMechanicsMathUtilities::TensorTransformation<3>(rOrigin,rOrigin,rTarget,rTarget,rTensor);

    // Computation actual stress in the covariant basis
    double detF;
    array_1d<double, 3> G1G2, g1g2;
    MathUtils<double>::CrossProduct(G1G2,column( GetValue(BASE_REF_1),rPointNumber ),column( GetValue(BASE_REF_2),rPointNumber ));
    MathUtils<double>::CrossProduct(g1g2,rg1,rg2);
    detF = norm_2(g1g2)/norm_2(G1G2);
    rTensor(0,0) /= detF;
    rTensor(1,1) /= detF;
    rTensor(1,0) /= detF;
    rTensor(0,1) /= detF;

    // Transform the prestress in the principal directions
    for (int i=0;i<3;i++){
        rOrigin(i,0) = rg1(i);
        rOrigin(i,1) = rg2(i);
        rOrigin(i,2) = rg3(i);
        rTarget(i,0) = rNAct(i,0);
        rTarget(i,1) = rNAct(i,1);
        rTarget(i,2) = rNAct(i,2);
    }

    StructuralMechanicsMathUtilities::TensorTransformation<3>(rOrigin,rOrigin,rTarget,rTarget,rTensor);

    // compute lambda_mod
    double lambda_mod_1, lambda_mod_2;
    double lambda_max = this->GetValue(LAMBDA_MAX);
    if(Lambda1 > lambda_max)
        lambda_mod_1 = lambda_max;
    else if(Lambda1 < 1.0/lambda_max)
        lambda_mod_1 = 1.0/lambda_max;
    else
        lambda_mod_1 = Lambda1;

    if(Lambda2 > lambda_max)
        lambda_mod_2 = lambda_max;
    else if(Lambda2 < 1.0/lambda_max)
        lambda_mod_2 = 1.0/lambda_max;
    else
        lambda_mod_2 = Lambda2;

    // compute modified prestress
    rTensor(0,0) *= Lambda1*lambda_mod_2 / (Lambda2* lambda_mod_1);
    rTensor(1,1) *= Lambda2*lambda_mod_1 / (Lambda1* lambda_mod_2);

    //transform the prestress into the actual local cartesian basis
    array_1d<double, 3> e1, e2, e3;
    e1 = rg1/norm_2(rg1);
    MathUtils<double>::CrossProduct(e3,e1,rg2);
    MathUtils<double>::CrossProduct(e2,e3,e1);
    e2 /= norm_2(e2);
    e3 /= norm_2(e3);
    for (int i=0;i<3;i++){
        rOrigin(i,0) = rNAct(i,0);
        rOrigin(i,1) = rNAct(i,1);
        rOrigin(i,2) = rNAct(i,2);
        rTarget(i,0) = e1(i);
        rTarget(i,1) = e2(i);
        rTarget(i,2) = e3(i);
    }

    StructuralMechanicsMathUtilities::TensorTransformation<3>(rOrigin,rOrigin,rTarget,rTarget,rTensor);

    Matrix& prestress_modified = GetValue(MEMBRANE_PRESTRESS);
    prestress_modified(0,rPointNumber) = rTensor(0,0);
    prestress_modified(1,rPointNumber) = rTensor(1,1);
    prestress_modified(2,rPointNumber) = rTensor(1,0);

}
void PrestressMembraneElement::ComputePrestress(const unsigned int rIntegrationPointSize){
    // initialize prestress matrix and prestress directions (with dummy zero values)
    unsigned int strain_size = this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    Matrix prestress_matrix(strain_size,rIntegrationPointSize,0), prestress_direction1(3,rIntegrationPointSize,0), prestress_direction2(3,rIntegrationPointSize,0);
    if(this->Has(MEMBRANE_PRESTRESS) == false)
        this->SetValue(MEMBRANE_PRESTRESS,prestress_matrix);
    this->SetValue(PRESTRESS_AXIS_1, prestress_direction1);
    this->SetValue(PRESTRESS_AXIS_2, prestress_direction2);

    const bool to_compute = GetProperties().Has(PROJECTION_TYPE_COMBO) ? GetProperties()[PROJECTION_TYPE_COMBO] != "file": true;

    if(to_compute) {
        if(GetProperties().Has(PRESTRESS_VECTOR)){
            Matrix& prestress_variable = this->GetValue(MEMBRANE_PRESTRESS);
            Matrix& prestress_axis_1 = this->GetValue(PRESTRESS_AXIS_1);
            Matrix& prestress_axis_2 = this->GetValue(PRESTRESS_AXIS_2);

            if(std::abs(GetProperties()[PRESTRESS_VECTOR](0)- GetProperties()[PRESTRESS_VECTOR](1)) <std::numeric_limits<double>::epsilon() && GetProperties()[PRESTRESS_VECTOR](2) == 0)
                mAnisotropicPrestress = false;

            else
                mAnisotropicPrestress = true;

            // read prestress from the material properties
            for(unsigned int point_number = 0; point_number < rIntegrationPointSize;point_number ++){
                // anisotropic case: prestress projection in the membrane
                if(mAnisotropicPrestress) {
                    // Initialize Prestress
                    for (unsigned int i_strain=0; i_strain<strain_size; i_strain++){
                        prestress_variable(i_strain,point_number) = GetProperties()[PRESTRESS_VECTOR](i_strain);
                        if(GetProperties().Has(PRESTRESS_AXIS_1_GLOBAL))
                            prestress_axis_1(i_strain,point_number) = GetProperties()[PRESTRESS_AXIS_1_GLOBAL](i_strain);
                        if(GetProperties().Has(PRESTRESS_AXIS_2_GLOBAL))
                            prestress_axis_2(i_strain,point_number) = GetProperties()[PRESTRESS_AXIS_2_GLOBAL](i_strain);
                    }
                    // in case that no prestress directions are prescribed: hardcode the prestress direction
                    if (GetProperties().Has(PRESTRESS_AXIS_1_GLOBAL) == false){
                        prestress_axis_1(0,point_number) = 1;
                        prestress_axis_1(1,point_number) = 0;
                        prestress_axis_1(2,point_number) = 0;
                    }
                    if (GetProperties().Has(PRESTRESS_AXIS_2_GLOBAL) == false){
                        prestress_axis_2(0,point_number) = 0;
                        prestress_axis_2(1,point_number) = 1;
                        prestress_axis_2(2,point_number) = 0;
                    }

                    // Project prestress in membrane plane
                    ProjectPrestress(point_number);
                }

                // in case of isotropic prestress: set prestress in the first step (no transformation necessary)
                else {
                    for (unsigned int i_strain=0; i_strain<strain_size; i_strain++){
                        prestress_variable(i_strain,point_number) = GetProperties()[PRESTRESS_VECTOR](i_strain);
                    }
                }
            }
        }
    }
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::ComputeBaseVectors(const GeometryType::IntegrationPointsArrayType& rIntegrationPoints){
    // compute Jacobian
    GeometryType::JacobiansType J0;
    J0 = GetGeometry().Jacobian(J0);

    mTotalDomainInitialSize = 0.00;
    Matrix dummy;
    SetValue(BASE_REF_1,dummy);
    SetValue(BASE_REF_2, dummy);

    Matrix& base_1 = GetValue(BASE_REF_1);
    Matrix& base_2 = GetValue(BASE_REF_2);
    base_1.resize(3,rIntegrationPoints.size());
    base_2.resize(3,rIntegrationPoints.size());

    // Calculating geometry tensors in reference configuration on Integration points
    for (unsigned int point_number = 0; point_number < rIntegrationPoints.size(); point_number++)
    {
        // getting information for integration
        double integration_weight = rIntegrationPoints[point_number].Weight();

        // base vectors in reference configuration
        array_1d<double, 3> G1, G2, G3;

        G1[0] = J0[point_number](0, 0);
        G2[0] = J0[point_number](0, 1);
        G1[1] = J0[point_number](1, 0);
        G2[1] = J0[point_number](1, 1);
        G1[2] = J0[point_number](2, 0);
        G2[2] = J0[point_number](2, 1);

        // Store base vectors in reference configuration
        column(base_1,point_number) = G1;
        column(base_2,point_number) = G2;
        // base vector G3
        MathUtils<double>::CrossProduct(G3, G1, G2);
        // differential area dA
        const double dA = norm_2(G3);

        KRATOS_DEBUG_ERROR_IF(dA < 1E-15) << "PrestressMembraneElement with id" << this->GetId()
                                          << "has ZERO differential area!" << std::endl;

        // normal vector n
        array_1d<double, 3> n = G3 / dA;

        // Get CovariantMetric
        mGab0[point_number][0] = std::pow(G1[0], 2) + std::pow(G1[1], 2) + std::pow(G1[2], 2);
        mGab0[point_number][1] = std::pow(G2[0], 2) + std::pow(G2[1], 2) + std::pow(G2[2], 2);
        mGab0[point_number][2] = G1[0] * G2[0] + G1[1] * G2[1] + G1[2] * G2[2];

        array_1d<double, 3> Gab_contravariant_1, Gab_contravariant_2;
        ComputeContravariantBaseVectors(Gab_contravariant_1,Gab_contravariant_2, point_number);

        // build local cartesian coordinate system
        double lg1 = norm_2(G1);
        array_1d<double, 3> E1 = G1 / lg1;
        double lg_contravariant_2 = norm_2(Gab_contravariant_2);
        array_1d<double, 3> E2 = Gab_contravariant_2 / lg_contravariant_2;

        BoundedMatrix<double, 2, 2> G;
        G(0, 0) = inner_prod(E1, Gab_contravariant_1);
        G(0, 1) = inner_prod(E1, Gab_contravariant_2);
        G(1, 0) = inner_prod(E2, Gab_contravariant_1);
        G(1, 1) = inner_prod(E2, Gab_contravariant_2);

        noalias(mGVector[point_number]) = ZeroMatrix(2, 2);
        // saving the G matrix for this point number
        noalias(mGVector[point_number]) = G;

        // Calculate the reduced mass matrix
        mDetJ0[point_number] = norm_2(G3);

        // Calculating the total area
        mTotalDomainInitialSize += mDetJ0[point_number] * integration_weight;
    }
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::ComputeContravariantBaseVectors(
                        array_1d<double, 3>& rG1Contra,
                        array_1d<double, 3>& rG2Contra,
                        const unsigned int rPointNumber){
    // determinant metric
    double det_metric = mGab0[rPointNumber][0]*mGab0[rPointNumber][1]- mGab0[rPointNumber][2]*mGab0[rPointNumber][2];

    // contravariant metric
    array_1d<double, 3> metric_contra;
    metric_contra[0]=  1.0/det_metric * mGab0[rPointNumber][1];
    metric_contra[1]=  1.0/det_metric * mGab0[rPointNumber][0];
    metric_contra[2]= -1.0/det_metric * mGab0[rPointNumber][2];

    // contravariant base vectors
    rG1Contra = metric_contra[0]*column( GetValue(BASE_REF_1),rPointNumber ) + metric_contra[2]*column( GetValue(BASE_REF_2),rPointNumber );
    rG2Contra = metric_contra[2]*column( GetValue(BASE_REF_1),rPointNumber ) + metric_contra[1]*column( GetValue(BASE_REF_2),rPointNumber );
}
//***********************************************************************************
//***********************************************************************************

const Matrix PrestressMembraneElement::CalculateDeformationGradient(const unsigned int rPointNumber){
    // Compute contravariant base vectors in reference configuration
    array_1d<double, 3> G1_contra, G2_contra;
    ComputeContravariantBaseVectors(G1_contra, G2_contra, rPointNumber);

    // Compute G3
    array_1d<double, 3> G3;
    MathUtils<double>::CrossProduct(G3, column( GetValue(BASE_REF_1),rPointNumber ), column( GetValue(BASE_REF_2),rPointNumber ));
    G3 = G3/ norm_2(G3);

    // Compute g1, g2, g3
    const GeometryType::ShapeFunctionsGradientsType& DN_De_container = GetGeometry().ShapeFunctionsLocalGradients();
    Matrix DN_De = DN_De_container[rPointNumber];
    array_1d<double, 3> g1, g2, g3, gab;
    CalculateMetricDeformed(rPointNumber, DN_De, gab, g1, g2);

    MathUtils<double>::CrossProduct(g3,g1,g2);
    g3 /= norm_2(g3);

    BoundedMatrix<double,3,3> deformation_gradient;
    for(unsigned int i=0; i<3; i++){
        for(unsigned int j=0; j<3; j++){
            deformation_gradient(i,j) = G1_contra(j)*g1(i) + G2_contra(j)*g2(i) + G3(j)*g3(i);
        }
    }
    return deformation_gradient;
}
//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::InitializeFormfinding(const unsigned int rIntegrationPointSize){

    if(mAnisotropicPrestress == true){
        // store base vectors of the initial configuration
        mG1Initial.resize(rIntegrationPointSize);
        mG2Initial.resize(rIntegrationPointSize);
        mG3Initial.resize(rIntegrationPointSize);

        for(unsigned int point_number = 0; point_number < rIntegrationPointSize;point_number ++){
            mG1Initial[point_number] = column( GetValue(BASE_REF_1),point_number );
            mG2Initial[point_number] = column( GetValue(BASE_REF_2),point_number );

            // out-of-plane direction
            MathUtils<double>::CrossProduct(mG3Initial[point_number],mG1Initial[point_number],mG2Initial[point_number]);
            mG3Initial[point_number] /= norm_2(mG3Initial[point_number]);
        }

        // set lambda_max
        if(GetProperties().Has(LAMBDA_MAX))
            this->SetValue(LAMBDA_MAX,GetProperties()[LAMBDA_MAX]);
        else
            this->SetValue(LAMBDA_MAX,1.2);
    }
}
//***********************************************************************************
//***********************************************************************************
int PrestressMembraneElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const unsigned int number_of_nodes = this->GetGeometry().size();
    // const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    // Verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(DENSITY)
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(THICKNESS)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( unsigned int i = 0; i < number_of_nodes; i++ ) {
        Node<3> &r_node = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,r_node)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node)
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(this->GetProperties().Has( CONSTITUTIVE_LAW ))
        << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const unsigned int strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    KRATOS_ERROR_IF( strain_size != 3) << "Wrong constitutive law used. This is a membrane element! "
        << "Expected strain size is 3 (el id = " << this->Id() << ")" << std::endl;

    //check constitutive law
    if(GetValue(IS_FORMFINDING)== false){
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
        {
            mConstitutiveLawVector[i]->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);

            ConstitutiveLaw::Features LawFeatures;
            mConstitutiveLawVector[i]->GetLawFeatures(LawFeatures);

            KRATOS_ERROR_IF(LawFeatures.mOptions.IsNot(ConstitutiveLaw::PLANE_STRESS_LAW))
                << "Constitutive law is compatible only with a plane stress 2D law for "
                << "membrane element with Id " << this->Id() << std::endl;

            KRATOS_ERROR_IF(LawFeatures.mOptions.IsNot(ConstitutiveLaw::INFINITESIMAL_STRAINS))
                << "Constitutive law is compatible only with a law using infinitessimal "
                << "strains for membrane element with Id " << this->Id() << std::endl;

            KRATOS_ERROR_IF(LawFeatures.mStrainSize != 3) << "Constitutive law expects a strain "
                << "size different from 3 for membrane element with Id "<< this->Id() <<std::endl;
        }
    }
    return 0;

    KRATOS_CATCH("");
}

void PrestressMembraneElement::save(Serializer& rSerializer) const
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
      rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
      rSerializer.save("DetJ0", mDetJ0);
      rSerializer.save("TotalDomainInitialSize", mTotalDomainInitialSize);
      rSerializer.save("G_ab", mGab0);
      rSerializer.save("G_Vector", mGVector);
      rSerializer.save("AnisotropicPrestress", mAnisotropicPrestress);
      rSerializer.save("G1Initial", mG1Initial);
      rSerializer.save("G2Initial", mG2Initial);
      rSerializer.save("G3Initial", mG3Initial);
    }

    void PrestressMembraneElement::load(Serializer& rSerializer)
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
      rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
      rSerializer.load("DetJ0", mDetJ0);
      rSerializer.load("TotalDomainInitialSize", mTotalDomainInitialSize);
      rSerializer.load("G_ab", mGab0);
      rSerializer.load("G_Vector", mGVector);
      rSerializer.load("AnisotropicPrestress", mAnisotropicPrestress);
      rSerializer.load("G1Initial", mG1Initial);
      rSerializer.load("G2Initial", mG2Initial);
      rSerializer.load("G3Initial", mG3Initial);
    }



//***********************************************************************************
//***********************************************************************************
} // Namespace Kratos.
