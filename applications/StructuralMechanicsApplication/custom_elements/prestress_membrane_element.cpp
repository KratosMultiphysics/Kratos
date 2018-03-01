// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Long Chen
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
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const

{
    return Kratos::make_shared< PrestressMembraneElement >(NewId, GetGeometry().Create(ThisNodes), pProperties);
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

  unsigned int NumNodes, LocalSize;
  unsigned int LocalIndex = 0;

  NumNodes = GetGeometry().size();
  LocalSize = NumNodes * 3;
  //const unsigned int NumNodes = GetGeometry().size();
  //const unsigned int LocalSize = NumNodes * 3;

  unsigned int dpos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

  if (rResult.size() != LocalSize)
      rResult.resize(LocalSize, false);

  for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
  {
      rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(DISPLACEMENT_X, dpos).EquationId();
      rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(DISPLACEMENT_Y, dpos + 1).EquationId();
      rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(DISPLACEMENT_Z, dpos + 2).EquationId();
  }

  KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo)

{
    unsigned int NumNodes, LocalSize;
    NumNodes = GetGeometry().size();
    LocalSize = NumNodes * 3;

    if (ElementalDofList.size() != LocalSize)
        ElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        ElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(DISPLACEMENT_X);
        ElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(DISPLACEMENT_Y);
        ElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(DISPLACEMENT_Z);
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
    mG1.resize(integration_points.size());
    mG2.resize(integration_points.size());
    mG1Initial.resize(integration_points.size());
    mG2Initial.resize(integration_points.size());
    mG3Initial.resize(integration_points.size());
    mStrainsVector.resize(integration_points.size());
    mStressesVector.resize(integration_points.size());
    mCauchyStressesVector.resize(integration_points.size()); //VM
    mGVector.resize(integration_points.size(), ZeroMatrix(2, 2));

    // compute base vectors in reference configuration, metrics
    ComputeBaseVectors(integration_points);

    // determine if the prestress state is isotropic or anisotropic
    if(GetProperties().Has(MEMBRANE_PRESTRESS)){
        if(GetProperties()[MEMBRANE_PRESTRESS](0,0)== GetProperties()[MEMBRANE_PRESTRESS](0,1) && GetProperties()[MEMBRANE_PRESTRESS](0,2) == 0)
            mAnisotropicPrestress = false;

        else
            mAnisotropicPrestress = true;
    }
    // initialize prestress matrix and prestress directions (with dummy zero values)
    unsigned int strain_size = this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    Matrix prestress_matrix(strain_size,integration_points.size(),0), prestress_direction1(3,integration_points.size(),0), prestress_direction2(3,integration_points.size(),0);
    this->SetValue(MEMBRANE_PRESTRESS,prestress_matrix);
    this->SetValue(PRESTRESS_AXIS_1_GLOBAL, prestress_direction1);
    this->SetValue(PRESTRESS_AXIS_2_GLOBAL, prestress_direction2);

    // set lambda_max
    if(GetProperties().Has(LAMBDA_MAX))
        this->SetValue(LAMBDA_MAX,GetProperties()[LAMBDA_MAX]);
    else
        this->SetValue(LAMBDA_MAX,1.2);
    mStep = 0;

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
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateLocalSystem(
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

    //rMassMatrix.resize(0,0);
    // LUMPED MASS MATRIX
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int mat_size = number_of_nodes * 3;

    if (rMassMatrix.size1() != mat_size)
    {
        rMassMatrix.resize(mat_size, mat_size);
    }

    rMassMatrix = ZeroMatrix(mat_size, mat_size);

    double TotalMass = mTotalDomainInitialSize * GetProperties()[THICKNESS] * GetProperties()[DENSITY];

    Vector LumpFact;

    LumpFact = GetGeometry().LumpingFactors(LumpFact);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        double temp = LumpFact[i] * TotalMass;

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

    double TotalMass = mTotalDomainInitialSize * GetProperties()[THICKNESS] * GetProperties()[DENSITY];

    Vector LumpFact;

    LumpFact = GetGeometry().LumpingFactors(LumpFact);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        double temp = LumpFact[i] * TotalMass;

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
    Vector& values,
    int Step)

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int mat_size = number_of_nodes * 3;

    if (values.size() != mat_size)
        values.resize(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3>& disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        unsigned int index = i * 3;
        values[index] = disp[0];
        values[index + 1] = disp[1];
        values[index + 2] = disp[2];
    }
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::GetFirstDerivativesVector(
    Vector& values,
    int Step)

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int mat_size = number_of_nodes * 3;

    if (values.size() != mat_size)
        values.resize(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        unsigned int index = i * 3;
        values[index] = vel[0];
        values[index + 1] = vel[1];
        values[index + 2] = vel[2];
    }

}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::GetSecondDerivativesVector(
    Vector& values,
    int Step)

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int mat_size = number_of_nodes * 3;

    if (values.size() != mat_size)
        values.resize(mat_size);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        const array_1d<double, 3>& acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        unsigned int index = i * 3;
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

void PrestressMembraneElement::CalculateAndAddKm(
    Matrix& K,
    Matrix& B,
    Matrix& D,
    double weight)

{
    KRATOS_TRY

    unsigned int dim = B.size2();
    Matrix temp(3, dim);
    noalias(temp) = prod(D, B);
    temp *= weight;
    Matrix Km(dim, dim);
    noalias(Km) = prod(trans(B), temp);
    noalias(K) += Km;

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateAndAddNonlinearKm(
    Matrix& K,
    Matrix& B11,
    Matrix& B22,
    Matrix& B12,
    Vector& SD,
    double weight)

{
    KRATOS_TRY;

    unsigned int number_of_nodes = GetGeometry().size();
    Matrix TestK = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);
    //TestK = K;
    //std::cout << " Get to CalculateAndAddNonlinearKm!" << std::endl;

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
                    /*
                    K(3 * n + i, 3 * m + j) += (SD[0] * B11(3 * n + i, 3 * m + j) + SD[1] * B22(3 * n + i, 3 * m + j) + SD[2] * B12(3 * n + i, 3 * m + j))*weight;
                    K(3 * m + j, 3 * n + i) += (SD[0] * B11(3 * n + i, 3 * m + j) + SD[1] * B22(3 * n + i, 3 * m + j) + SD[2] * B12(3 * n + i, 3 * m + j))*weight;


                    TestK(3*n + i, 3*m + j) += (SD[0] * B11(3 * n + i, 3 * m + j) + SD[1] * B22(3 * n + i, 3 * m + j) + SD[2] * B12(3 * n + i, 3 * m + j))*weight;
                    TestK(3 * m + i, 3 * n + j) = TestK(3 * n + i, 3 * m + j);
                    */
                    K(3 * n + i, 3 * m + j) += (SD[0] * B11(3 * n + i, 3 * m + j) + SD[1] * B22(3 * n + i, 3 * m + j) + SD[2] * B12(3 * n + i, 3 * m + j))*weight;
                    K(3 * m + i, 3 * n + j) = K(3 * n + i, 3 * m + j);

                    //if((3*m + i) == (3 * n + j))
                    //    K(3*m + i, 3*n + j) = K(3*n + i, 3*m + j);
                    //else
                    //    K(3 * m + j, 3 * n + i) += (SD[0] * B11(3 * n + i, 3 * m + j) + SD[1] * B22(3 * n + i, 3 * m + j) + SD[2] * B12(3 * n + i, 3 * m + j))*weight;
                }
            }
        }
    }
    //if (this->GetId() == 42) {
    //    std::cout << "############" << std::endl;
    //    KRATOS_WATCH(TestK);
    //    std::cout << "============" << std::endl;
    //    KRATOS_WATCH(K);
    //    std::cout << "############" << std::endl;

    //}

    //K = TestK;
    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************
void PrestressMembraneElement::ClearNodalForces()
{
    KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        if (GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_FORCE) && GetGeometry()[i].SolutionStepsDataHas(INTERNAL_FORCE))
        {

            array_1d<double, 3 > & ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
            array_1d<double, 3 > & InternalForce = GetGeometry()[i].FastGetSolutionStepValue(INTERNAL_FORCE);

            GetGeometry()[i].SetLock();
            ExternalForce.clear();
            InternalForce.clear();
            GetGeometry()[i].UnSetLock();

        }

    }

    KRATOS_CATCH("")
}




//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateQ(
    boost::numeric::ublas::bounded_matrix<double, 3, 3>& Q,
    Matrix& mG)

{
    KRATOS_TRY
    Q(0, 0) = std::pow(mG(0, 0), 2);
    Q(0, 1) = std::pow(mG(0, 1), 2);
    Q(0, 2) = 2.00*mG(0, 0)*mG(0, 1);

    Q(1, 0) = std::pow(mG(1, 0), 2);
    Q(1, 1) = std::pow(mG(1, 1), 2);
    Q(1, 2) = 2.00*mG(1, 0) * mG(1, 1);

    Q(2, 0) = 2.00 * mG(0, 0) * mG(1, 0);
    Q(2, 1) = 2.00 * mG(0, 1)*mG(1, 1);
    Q(2, 2) = 2.00 * (mG(0, 0) * mG(1, 1) + mG(0, 1)*mG(1, 0));

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateB(
    Matrix& B,
    boost::numeric::ublas::bounded_matrix<double, 3, 3>& Q,
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

    B = prod(Q, b);

    KRATOS_CATCH("")
}


//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateStrain(
    Vector& StrainVector,
    array_1d<double, 3>& gab,
    array_1d<double, 3>& gab0)

{
    KRATOS_TRY

    StrainVector[0] = 0.5 * (gab[0] - gab0[0]);
    StrainVector[1] = 0.5 * (gab[1] - gab0[1]);
    StrainVector[2] = 0.5 * (gab[2] - gab0[2]);

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateAndAdd_BodyForce(
    const Vector& N,
    const ProcessInfo& rCurrentProcessInfo,
    array_1d<double, 3>& BodyForce,
    VectorType& rRightHandSideVector,
    double weight)

{
    KRATOS_TRY

        unsigned int number_of_nodes = this->GetGeometry().size();
    const double density = this->GetProperties()[DENSITY];


    noalias(BodyForce) = ZeroVector(3);
    for (unsigned int i = 0; i<number_of_nodes; i++)
        BodyForce += N[i] * this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    BodyForce *= density;

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        int index = 3 * i;

        for (unsigned int j = 0; j < 3; j++)
            rRightHandSideVector[index + j] += weight * N[i] * BodyForce[j];
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::CalculateAndAdd_PressureForce(
    VectorType& residualvector,
    const Vector& N,
    const array_1d<double, 3>& v3,
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

void PrestressMembraneElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag)

{
    KRATOS_TRY

    // Initializing all needed variables
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int mat_size = number_of_nodes * 3;

    // Matrix B(3, mat_size);
    // change to: array_1d<double, 3> StrainVector;
    Vector StrainVector(3);
    Vector StressVector(3);

    // set up Constitutive Law
    ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    Values.SetStrainVector(StrainVector);       // this is the input parameter
    Values.SetStressVector(StressVector);       // this is an output parameter

    // resizing as needed the LHS
    if (CalculateStiffnessMatrixFlag == true)    // calculation of the matrix is required
    {
        if (rLeftHandSideMatrix.size1() != mat_size)
        {
            rLeftHandSideMatrix.resize(mat_size, mat_size);
        }

        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
    }

    //resizing as needed the RHS
    if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
    {
        if (rRightHandSideVector.size() != mat_size)
        {
            rRightHandSideVector.resize(mat_size);
        }

        rRightHandSideVector = ZeroVector(mat_size); //resetting RHS
    }

    // Initializing the Nodal coordinates
    // change to: boost::numeric::ublas::bounded_matrix
    Matrix xyz_reference;   // Nodal coordinates in the reference configuration
    Matrix xyz_actual;      // Nodal coordinates in the actual configuration

    // reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = GetGeometry().ShapeFunctionsLocalGradients();

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

    // calculating actual jacobian
    GeometryType::JacobiansType J;

    J = GetGeometry().Jacobian(J);

    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
    {
        // reading integration weight, shape function value and its gradients at this integration point
        const double IntegrationWeight = integration_points[PointNumber].Weight();

        Vector ShapeFunctionN = row(Ncontainer, PointNumber);
        Matrix DN_De = DN_DeContainer[PointNumber];

        // covariant metric in deformed system
        array_1d<double, 3> gab;

        // Transformation Matrix Q
        boost::numeric::ublas::bounded_matrix<double, 3, 3> Q = ZeroMatrix(3, 3);
        // basis vectors in deformed system
        array_1d<double, 3> g1;
        array_1d<double, 3> g2;
        array_1d<double, 3> g3;

        CalculateQ(Q, mGVector[PointNumber]);
        CalculateMetricDeformed(PointNumber, DN_De, gab, g1, g2);
        CalculateStrain(StrainVector, gab, mGab0[PointNumber]);

        Vector CartesianStrainVector = prod(Q, StrainVector); //in refence configuration

        // Constitutive Matrices D
        Matrix D(3, 3);

        Values.SetConstitutiveMatrix(D); // this is an output parameter

        //Deformation Gradient
        Values.SetDeformationGradientF(CalculateDeformationGradient(PointNumber));

        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);     // Why is the curviliear strains are used here?

        // Deformations for Non-linear force vector
        Vector StrainDeformation;

        StrainDeformation = prod(trans(D), CartesianStrainVector);

        // Adding the pre-stress values as forces over length
        for (int i = 0; i < 3; ++i)
            StrainDeformation[i] += GetValue(MEMBRANE_PRESTRESS)(i,PointNumber);


        // calculate B matrices
        // B matrices:
        Matrix B = ZeroMatrix(3, mat_size);
        //CalculateB(B, Q, DN_De, mG1[PointNumber], mG2[PointNumber]);
        CalculateB(B, Q, DN_De, g1, g2);

        // integration on the REFERENCE CONFIGURATION
        double DetJ0 = mDetJ0[PointNumber];
        double IntToReferenceWeight = IntegrationWeight * DetJ0 * GetProperties()[THICKNESS];

        // Nonlinear Deformation
        Matrix Strain_locCartesian_11 = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);
        Matrix Strain_locCartesian_22 = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);
        Matrix Strain_locCartesian_12 = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);

        CalculateSecondVariationStrain(DN_De, Strain_locCartesian_11, Strain_locCartesian_22, Strain_locCartesian_12, Q, g1, g2);

        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true)
        {
            // adding membrane contribution to the stiffness matrix
            CalculateAndAddKm(rLeftHandSideMatrix, B, D, IntToReferenceWeight);

            // adding non-linear-contribution to stiffness matrix
            CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
                Strain_locCartesian_11, Strain_locCartesian_22, Strain_locCartesian_12,
                StrainDeformation,
                IntToReferenceWeight);
        }

        // RIGHT HAND SIDE VECTOR
        if (CalculateResidualVectorFlag == true)         // calculation of the matrix is required
        {
            // operation performed: rRighthandSideVector -= Weight* IntForce
            noalias(rRightHandSideVector) -= IntToReferenceWeight* prod(trans(B), StrainDeformation);
        }
    } // end loop over integration points

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************
void PrestressMembraneElement::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR)
    {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }

    if (rVariable == PK2_STRESS_TENSOR)
    {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }
    // VM
    if (rVariable == CAUCHY_STRESS_TENSOR)
    {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }
    // VM
}
//***********************************************************************************
//***********************************************************************************
void PrestressMembraneElement::CalculateMetricDeformed(const unsigned int& PointNumber, Matrix DN_De,
    array_1d<double, 3>& gab,
    array_1d<double, 3>& g1,
    array_1d<double, 3>& g2)
{
    //Matrix J0;
    //Jacobian(DN_De, J0);

    GeometryType::JacobiansType J_act;
    J_act = GetGeometry().Jacobian(J_act);

    //KRATOS_WATCH(J_act);

    //auxiliary terms
    //array_1d<double, 3> g1;
    //array_1d<double, 3> g2;
    array_1d<double, 3> g3;

    //double IntegrationWeight = GetGeometry().IntegrationPoints()[0].Weight();

    g1[0] = J_act[PointNumber](0, 0);
    g2[0] = J_act[PointNumber](0, 1);
    g1[1] = J_act[PointNumber](1, 0);
    g2[1] = J_act[PointNumber](1, 1);
    g1[2] = J_act[PointNumber](2, 0);
    g2[2] = J_act[PointNumber](2, 1);

    //basis vector g3
    MathUtils<double>::CrossProduct(g3, g1, g2);
    //differential area dA
    double dA = norm_2(g3);
    //normal vector _n
    array_1d<double, 3> n = g3 / dA;

    //GetCovariantMetric
    gab[0] = std::pow(g1[0], 2) + std::pow(g1[1], 2) + std::pow(g1[2], 2);
    gab[1] = std::pow(g2[0], 2) + std::pow(g2[1], 2) + std::pow(g2[2], 2);
    gab[2] = g1[0] * g2[0] + g1[1] * g2[1] + g1[2] * g2[2];
}

//***********************************************************************************
//***********************************************************************************
void PrestressMembraneElement::CalculateSecondVariationStrain(Matrix DN_De,
    Matrix & Strain_locCartesian11,
    Matrix & Strain_locCartesian22,
    Matrix & Strain_locCartesian12,
    boost::numeric::ublas::bounded_matrix<double, 3, 3>& Q,
    array_1d<double, 3>& g1,
    array_1d<double, 3>& g2)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    //Matrix dg3_n = ZeroMatrix(3, 3);
    //Matrix dg3_m = ZeroMatrix(3, 3);

    Vector ddStrain_curvilinear = ZeroVector(3);

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
                    ddStrain_curvilinear = ZeroVector(3);
                    if (j == i)
                    {
                        ddStrain_curvilinear[0] = DN_De(n, 0)*DN_De(m, 0);
                        ddStrain_curvilinear[1] = DN_De(n, 1)*DN_De(m, 1);
                        ddStrain_curvilinear[2] = 0.5*(DN_De(n, 0)*DN_De(m, 1) + DN_De(n, 1)*DN_De(m, 0));

                        Strain_locCartesian11(3 * n + i, 3 * m + j) = Q(0, 0)*ddStrain_curvilinear[0] + Q(0, 1)*ddStrain_curvilinear[1] + Q(0, 2)*ddStrain_curvilinear[2];
                        Strain_locCartesian22(3 * n + i, 3 * m + j) = Q(1, 0)*ddStrain_curvilinear[0] + Q(1, 1)*ddStrain_curvilinear[1] + Q(1, 2)*ddStrain_curvilinear[2];
                        Strain_locCartesian12(3 * n + i, 3 * m + j) = Q(2, 0)*ddStrain_curvilinear[0] + Q(2, 1)*ddStrain_curvilinear[1] + Q(2, 2)*ddStrain_curvilinear[2];

                    }
                }
            }
        }
    }
}
//************************************************************************************
//************************************************************************************
void PrestressMembraneElement::CalculateMembraneElasticityTensor(
    Matrix& D
    )
{
    double NU = GetProperties()[POISSON_RATIO];
    double E = GetProperties()[YOUNG_MODULUS];
    double coeff = E / (1 - NU*NU);
    D(0, 0) = coeff;
    D(0, 1) = NU*coeff;
    D(0, 2) = 0.0;

    D(1, 0) = NU*coeff;
    D(1, 1) = coeff;
    D(1, 2) = 0.0;

    D(2, 0) = 0.0;
    D(2, 1) = 0.0;
    D(2, 2) = 0.5*(1 - NU)*coeff;
}

//***********************************************************************************
//***********************************************************************************
void PrestressMembraneElement::TransformPrestress(
    const unsigned int PointNumber){

    // definition prestress axes
    array_1d<double,3> global_prestress_axis1, global_prestress_axis2, global_prestress_axis3;
    for(unsigned int i=0; i<3;i++){
        global_prestress_axis1[i] = GetValue(PRESTRESS_AXIS_1_GLOBAL)(i,PointNumber);
        global_prestress_axis2[i] = GetValue(PRESTRESS_AXIS_2_GLOBAL)(i,PointNumber);
    }

    // normalization global prestress axes
    global_prestress_axis1 /= norm_2(global_prestress_axis1);
    global_prestress_axis2 /= norm_2(global_prestress_axis2);

    // Compute global_prestress_axis3
    MathUtils<double>::CrossProduct(global_prestress_axis3, global_prestress_axis1, global_prestress_axis2);
    global_prestress_axis3 /= norm_2(global_prestress_axis3);

    // Compute local cartesian basis E1, E2, E3
    array_1d<double,3> E1 = mG1[PointNumber]/norm_2(mG1[PointNumber]);
    array_1d<double,3> E2, E3;
    MathUtils<double>::CrossProduct(E3,E1,mG2[PointNumber]);
    E3 /= norm_2(E3);
    MathUtils<double>::CrossProduct(E2,E3,E1);
    E2 /= norm_2(E2);


    // Compute T1, T2, T3
    array_1d<double, 3> normal_projection_plane;
    MathUtils<double>::CrossProduct(normal_projection_plane, global_prestress_axis3, global_prestress_axis1);
    array_1d<double, 3> T1, T2, T3;
    MathUtils<double>::CrossProduct(T1, normal_projection_plane, E3);
    T1 /= norm_2(T1);
    MathUtils<double>::CrossProduct(T2,E3,T1);
    T2 /= norm_2(T2);
    T3 = E3;



    // Transform prestresses in the local cartesian cosy in reference configuration
    bounded_matrix<double,3,3> origin = ZeroMatrix(3,3);
    bounded_matrix<double,3,3> target = ZeroMatrix(3,3);
    bounded_matrix<double,3,3> tensor = ZeroMatrix(3,3);

    for (int i=0;i<3;i++){
        origin(i,0) = T1(i);
        origin(i,1) = T2(i);
        origin(i,2) = T3(i);
        target(i,0) = E1(i);
        target(i,1) = E2(i);
        target(i,2) = E3(i);
    }

    tensor.clear();
    tensor(0,0) = this->GetValue(MEMBRANE_PRESTRESS)(0,PointNumber);
    tensor(1,1) = this->GetValue(MEMBRANE_PRESTRESS)(1,PointNumber);
    tensor(1,0) = this->GetValue(MEMBRANE_PRESTRESS)(2,PointNumber);
    tensor(0,1) = this->GetValue(MEMBRANE_PRESTRESS)(2,PointNumber);

    // Transformation Stress Tensor
    StructuralMechanicsMathUtilities::TensorTransformation<3>(origin,origin,target,target,tensor);

    // store prestress values in the elemental data
    Matrix& prestress_matrix = GetValue(MEMBRANE_PRESTRESS);
    prestress_matrix(0,PointNumber) = tensor(0,0);
    prestress_matrix(1,PointNumber) = tensor(1,1);
    prestress_matrix(2,PointNumber) = tensor(1,0);
}
//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo){
// in each solution step: update reference configuration, update prestress (in the anisotropic case)

    // reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    // increase step
    mStep++;

    // Prestress update for anisotropic prestress (in step 1: nothing to update)
    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
    {
        if(mAnisotropicPrestress == true && mStep != 1)
            UpdatePrestress(PointNumber);
    }

    //update base vectors in reference configuration, metrics
    ComputeBaseVectors(integration_points);

    // Calculating geometry tensors in reference configuration on Integration points
    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
    {
        // if step =1: in case of anisotropic remeshing: Store initial reference configuration for later use in the prestress update
        if(mStep ==1 && mAnisotropicPrestress == true){

            mG1Initial[PointNumber] = mG1[PointNumber];
            mG2Initial[PointNumber] = mG2[PointNumber];

            // out-of-plane direction
            MathUtils<double>::CrossProduct(mG3Initial[PointNumber],mG1Initial[PointNumber],mG2Initial[PointNumber]);
            mG3Initial[PointNumber] /= norm_2(mG3Initial[PointNumber]);

        }

        // Compute the prestress (in the local cartesian CoSy in reference configuration)
        PrestressComputation(PointNumber);

    }

    // Initialize Material
    //InitializeMaterial(integration_points.size());
}
//***********************************************************************************
//***********************************************************************************
void PrestressMembraneElement::InitializeNonLinearIteration(){

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
void PrestressMembraneElement::UpdatePrestress(const unsigned int PointNumber){
    // --1--computation relevant CoSys
    array_1d<double, 3> g1, g2, g3, gab; //base vectors/metric actual config
    array_1d<double, 3> G3;  // base vector reference config
    array_1d<double, 3> E1_tot,E2_tot, E3_tot; //local cartesian base vectors initial reference config
    array_1d<double, 3> E1, E2, E3; // local cartesian base vectors initial config
    array_1d<double, 3> base_ref_contra_tot_1, base_ref_contra_tot_2; //contravariant base vectors initial reference config

    ComputeRelevantCoSys(PointNumber, g1, g2, g3, gab, G3, E1_tot, E2_tot, E3_tot, E1, E2, E3, base_ref_contra_tot_1, base_ref_contra_tot_2);

    // --2-- computation total deformation gradient
    bounded_matrix<double,3,3> deformation_gradient_total;
    for(unsigned int i=0; i<3; i++){
        for(unsigned int j=0; j<3; j++){
            deformation_gradient_total(i,j) = base_ref_contra_tot_1(j)*g1(i) + base_ref_contra_tot_2(j)*g2(i) + mG3Initial[PointNumber](j)*g3(i);
        }
    }

    //--3--Compute the eigenvalues of the total deformation gradient
    bounded_matrix<double,3,3> origin = ZeroMatrix(3,3);
    bounded_matrix<double,3,3> target = ZeroMatrix(3,3);
    bounded_matrix<double,3,3> tensor = ZeroMatrix(3,3);
    double lambda_1, lambda_2;

    ComputeEigenvaluesDeformationGradient(PointNumber,
                    origin, target, tensor,
                    base_ref_contra_tot_1,base_ref_contra_tot_2,
                    E1_tot, E2_tot, E3_tot,
                    gab,
                    lambda_1, lambda_2);

    //--4--Compute the eigenvectors in the reference and actual configuration
    bounded_matrix<double,3,3> N_act = ZeroMatrix(3,3); // eigenvectors in actual configuration
    ComputeEigenvectorsDeformationGradient(PointNumber,
                                tensor, origin,
                                deformation_gradient_total,
                                E1_tot, E2_tot,
                                lambda_1, lambda_2,
                                N_act);

    //--5--Compute the modified prestress
    ModifyPrestress(PointNumber,
                    origin, target,tensor,
                    E1, E2, E3, G3,
                    g1, g2, g3, N_act,
                    lambda_1, lambda_2);
}
//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::ComputeRelevantCoSys(const unsigned int PointNumber,
             array_1d<double, 3>& rg1,array_1d<double, 3>& rg2,array_1d<double, 3>& rg3, array_1d<double, 3>& rgab,
             array_1d<double, 3>& rG3,
             array_1d<double, 3>& rE1Tot, array_1d<double, 3>& rE2Tot,array_1d<double, 3>& rE3Tot,
             array_1d<double, 3>& rE1,array_1d<double, 3>& rE2,array_1d<double, 3>& rE3,
             array_1d<double, 3>& rBaseRefContraTot1,array_1d<double, 3>& rBaseRefContraTot2){

    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = GetGeometry().ShapeFunctionsLocalGradients();
    Matrix DN_De = DN_DeContainer[PointNumber];
    CalculateMetricDeformed(PointNumber, DN_De, rgab, rg1, rg2);

    // compute the out-of-plane direction and normalize it
    MathUtils<double>::CrossProduct(rg3,rg1,rg2);
    rg3 /= norm_2(rg3);

    // computation of the out-of-plane direction in the reference configuration
    MathUtils<double>::CrossProduct(rG3,mG1[PointNumber],mG2[PointNumber]);
    rG3 /= norm_2(rG3);

    //Compute cartesian basis in total reference configuration
    // rE1Tot = mG1/|mG1|, rE3Tot = mG3Initial
    rE1Tot = mG1Initial[PointNumber]/norm_2(mG1Initial[PointNumber]);
    rE3Tot = mG3Initial[PointNumber];
    MathUtils<double>::CrossProduct(rE2Tot,rE3Tot,rE1Tot);
    rE2Tot /= norm_2(rE2Tot);

    //Compute cartesian basis in (updated) reference configuration
    // rE1 = mG1/|mG1|, rE3 = mG3Initial
    rE1 = mG1[PointNumber]/norm_2(mG1[PointNumber]);
    rE3 = rG3;
    MathUtils<double>::CrossProduct(rE2,rE3,rE1);
    rE2 /= norm_2(rE2);

    // Compute contravariant basis in total Reference configuration
    // -->Compute the metric in total reference configuration
    array_1d<double, 3> metric_reference_tot;
    metric_reference_tot(0)=inner_prod(mG1Initial[PointNumber],mG1Initial[PointNumber]);
    metric_reference_tot(1)=inner_prod(mG2Initial[PointNumber],mG2Initial[PointNumber]);
    metric_reference_tot(2)=inner_prod(mG2Initial[PointNumber],mG1Initial[PointNumber]);

    // -->Invert metric in (total) reference configuration
    double det_metric_tot = metric_reference_tot(0)*metric_reference_tot(1)-metric_reference_tot(2)*metric_reference_tot(2);
    array_1d<double, 3> inv_metric_tot;
    #ifdef KRATOS_DEBUG
        KRATOS_ERROR_IF(det_metric_tot < std::numeric_limits<double>::epsilon() && det_metric_tot > -std::numeric_limits<double>::epsilon()) << "division by zero!" << std::endl;
    #endif
    inv_metric_tot(0) = 1.0/det_metric_tot * metric_reference_tot(1);
    inv_metric_tot(1) = 1.0/det_metric_tot * metric_reference_tot(0);
    inv_metric_tot(2) = -1.0/det_metric_tot * metric_reference_tot(2);

    // -->Compute contravariant basis (in total reference configuration)
    rBaseRefContraTot1 = inv_metric_tot(0)*mG1Initial[PointNumber] + inv_metric_tot(2)*mG2Initial[PointNumber];
    rBaseRefContraTot2 = inv_metric_tot(2)*mG1Initial[PointNumber] + inv_metric_tot(1)*mG2Initial[PointNumber];

}
//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::ComputeEigenvaluesDeformationGradient(const unsigned int PointNumber,
                    bounded_matrix<double,3,3>& rOrigin, bounded_matrix<double,3,3>& rTarget, bounded_matrix<double,3,3>& rTensor,
                    const array_1d<double, 3>& rBaseRefContraTot1, const array_1d<double, 3>& rBaseRefContraTot2,
                    const array_1d<double, 3>& rE1Tot, const array_1d<double, 3>& rE2Tot, const array_1d<double, 3>& rE3Tot,
                    const array_1d<double, 3>& rgab,
                    double& rLambda1, double& rLambda2){
    for (int i=0;i<3;i++){
        rOrigin(i,0) = rBaseRefContraTot1(i);
        rOrigin(i,1) = rBaseRefContraTot2(i);
        rOrigin(i,2) = mG3Initial[PointNumber](i);
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

void PrestressMembraneElement::ComputeEigenvectorsDeformationGradient(const unsigned int PointNumber,
                                bounded_matrix<double,3,3>& rTensor, bounded_matrix<double,3,3>& rOrigin,
                                const bounded_matrix<double,3,3>& rDeformationGradientTotal,
                                const array_1d<double, 3>& rE1Tot, const array_1d<double, 3>& rE2Tot,
                                const double Lambda1, const double Lambda2,
                                bounded_matrix<double,3,3>& rNAct){
    //Check if rTensor is the Identity Matrix -> no deformation
    //(rTensor-lambda^2*I)=ZeroMatrix
    bool principal_strain_state = false;
    rTensor(0,0) -= Lambda1*Lambda1;
    rTensor(1,1) -= Lambda2*Lambda2;
    if (std::abs(rTensor(0,0))<1.0e-6 && std::abs(rTensor(1,0))<1.0e-6 && std::abs(rTensor(0,1))<1.0e-6 && std::abs(rTensor(1,1))<1.0e-6)
        principal_strain_state = true;

    // compute C (=rTensor)
    rTensor = prec_prod(trans(rDeformationGradientTotal),rDeformationGradientTotal);

    // Eigenvectors in reference configuration: N_ref
    array_1d<double, 3> N_ref_1, N_ref_2, N_ref_3;
    N_ref_3 = mG3Initial[PointNumber];

    if (principal_strain_state){
        N_ref_1 = rE1Tot;
        N_ref_2 = rE2Tot;
    }

    else{
        //Compute the first principal direction
        //Eigenvector lies in the plan spanned up by G1 and G2
        //-> (C-lambda^2*I)N =(C-lambda^2*I)(alpha1*G1+alpha2*G2)=B1*alpha1+B2*alpha2 =0

        rTensor(0,0) -= Lambda1*Lambda1;
        rTensor(1,1) -= Lambda1*Lambda1;
        rTensor(2,2) -= Lambda1*Lambda1;

        column(rOrigin,0) = mG1Initial[PointNumber];
        column(rOrigin,1) = mG2Initial[PointNumber];

        bounded_matrix<double,3,3> B;
        B = prec_prod(rTensor, rOrigin);

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
        N_ref_1 = alpha_1*mG1Initial[PointNumber] + alpha_2*mG2Initial[PointNumber];
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

void PrestressMembraneElement::ModifyPrestress(const unsigned int PointNumber,
                    bounded_matrix<double,3,3>& rOrigin, bounded_matrix<double,3,3>& rTarget,bounded_matrix<double,3,3>& rTensor,
                    const array_1d<double, 3>& rE1, const array_1d<double, 3>& rE2, const array_1d<double, 3>& rE3, const array_1d<double, 3>& rG3,
                    const array_1d<double, 3>& rg1, const array_1d<double, 3>& rg2, const array_1d<double, 3>& rg3, const bounded_matrix<double,3,3>& rNAct,
                    const double Lambda1, const double Lambda2){
    // Transform prestresses from the local cosy in reference config to the curvilinear cosy in reference config, covariant
    for (int i=0;i<3;i++){
        rOrigin(i,0) = rE1(i);
        rOrigin(i,1) = rE2(i);
        rOrigin(i,2) = rE3(i);
        rTarget(i,0) = mG1[PointNumber](i);
        rTarget(i,1) = mG2[PointNumber](i);
        rTarget(i,2) = rG3(i);
    }
    rTensor.clear();
    rTensor(0,0) = GetValue(MEMBRANE_PRESTRESS)(0,PointNumber);
    rTensor(1,1) = GetValue(MEMBRANE_PRESTRESS)(1,PointNumber);
    rTensor(1,0) = GetValue(MEMBRANE_PRESTRESS)(2,PointNumber);
    rTensor(0,1) = GetValue(MEMBRANE_PRESTRESS)(2,PointNumber);

    StructuralMechanicsMathUtilities::TensorTransformation<3>(rOrigin,rOrigin,rTarget,rTarget,rTensor);

    // Computation actual stress in the covariant basis
    double detF;
    array_1d<double, 3> G1G2, g1g2;
    MathUtils<double>::CrossProduct(G1G2,mG1[PointNumber],mG2[PointNumber]);
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
    prestress_modified(0,PointNumber) = rTensor(0,0);
    prestress_modified(1,PointNumber) = rTensor(1,1);
    prestress_modified(2,PointNumber) = rTensor(1,0);

}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::PrestressComputation(const unsigned int PointNumber){
    // initialize prestress matrix and prestress directions
    Matrix& prestress_variable = this->GetValue(MEMBRANE_PRESTRESS);
    Matrix& prestress_axis_1 = this->GetValue(PRESTRESS_AXIS_1_GLOBAL);
    Matrix& prestress_axis_2 = this->GetValue(PRESTRESS_AXIS_2_GLOBAL);
    unsigned int strain_size = this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // anisotropic case: prestress projection in the first step
    if(mAnisotropicPrestress && mStep == 1){
            // Initialize Prestress
            for (unsigned int i_strain=0; i_strain<strain_size; i_strain++){
                prestress_variable(i_strain,PointNumber) = GetProperties()[MEMBRANE_PRESTRESS](0,i_strain);
                if(GetProperties().Has(PRESTRESS_AXIS_2_GLOBAL) && GetProperties().Has(PRESTRESS_AXIS_1_GLOBAL) == true){
                    prestress_axis_1(i_strain,PointNumber) = GetProperties()[PRESTRESS_AXIS_1_GLOBAL](0,i_strain);
                    prestress_axis_2(i_strain,PointNumber) = GetProperties()[PRESTRESS_AXIS_2_GLOBAL](0,i_strain);
                }
            }
            // in case that no prestress directions are prescribed: hardcode the prestress direction
            if ((GetProperties().Has(PRESTRESS_AXIS_2_GLOBAL) == false) || (GetProperties().Has(PRESTRESS_AXIS_1_GLOBAL) == false)){
                prestress_axis_1(0,PointNumber) = 1;
                prestress_axis_1(1,PointNumber) = 0;
                prestress_axis_1(2,PointNumber) = 0;

                prestress_axis_2(0,PointNumber) = 0;
                prestress_axis_2(1,PointNumber) = 1;
                prestress_axis_2(2,PointNumber) = 0;
            }

            // if(this->Id()==1)
            //     std::cout<<"prestress before transformation: "<<GetValue(MEMBRANE_PRESTRESS)<<std::endl;
            TransformPrestress(PointNumber);
            // if(this->Id()==1)
            //     std::cout<<"prestress after transformation: "<<GetValue(MEMBRANE_PRESTRESS)<<std::endl;
        }

        // in case of isotropic prestress: set prestress in the first step (no transformation necessary)
        if((mAnisotropicPrestress == false)&& mStep ==1){
            for (unsigned int i_strain=0; i_strain<strain_size; i_strain++){
                prestress_variable(i_strain,PointNumber) = GetProperties()[MEMBRANE_PRESTRESS](0,i_strain);
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

    // Calculating geometry tensors in reference configuration on Integration points
    for (unsigned int PointNumber = 0; PointNumber < rIntegrationPoints.size(); PointNumber++)
    {
        // getting information for integration
        double IntegrationWeight = rIntegrationPoints[PointNumber].Weight();

        // base vectors in reference configuration
        array_1d<double, 3> G1;
        array_1d<double, 3> G2;
        array_1d<double, 3> G3;

        G1[0] = J0[PointNumber](0, 0);
        G2[0] = J0[PointNumber](0, 1);
        G1[1] = J0[PointNumber](1, 0);
        G2[1] = J0[PointNumber](1, 1);
        G1[2] = J0[PointNumber](2, 0);
        G2[2] = J0[PointNumber](2, 1);


        // Store base vectors in reference configuration
        mG1[PointNumber] = G1;
        mG2[PointNumber] = G2;


        //KRATOS_WATCH(J0);

        // base vector G3
        MathUtils<double>::CrossProduct(G3, G1, G2);
        // differential area dA
        const double dA = norm_2(G3);

        KRATOS_DEBUG_ERROR_IF(dA < 1E-15) << "PrestressMembraneElement with id" << this->GetId()
                                          << "has ZERO differential area!" << std::endl;

        // normal vector n
        array_1d<double, 3> n = G3 / dA;

        // Get CovariantMetric
        mGab0[PointNumber][0] = std::pow(G1[0], 2) + std::pow(G1[1], 2) + std::pow(G1[2], 2);
        mGab0[PointNumber][1] = std::pow(G2[0], 2) + std::pow(G2[1], 2) + std::pow(G2[2], 2);
        mGab0[PointNumber][2] = G1[0] * G2[0] + G1[1] * G2[1] + G1[2] * G2[2];

        array_1d<double, 3> Gab_contravariant_1, Gab_contravariant_2;
        ComputeContravariantBaseVectors(Gab_contravariant_1,Gab_contravariant_2, PointNumber);

        // build local cartesian coordinate system
        double lg1 = norm_2(G1);
        array_1d<double, 3> E1 = G1 / lg1;
        double lg_contravariant_2 = norm_2(Gab_contravariant_2);
        array_1d<double, 3> E2 = Gab_contravariant_2 / lg_contravariant_2;

        boost::numeric::ublas::bounded_matrix<double, 2, 2> mG;
        //boost::numeric::ublas::bounded_matrix<double, 2, 3> mG;
        mG(0, 0) = inner_prod(E1, Gab_contravariant_1);
        mG(0, 1) = inner_prod(E1, Gab_contravariant_2);
        mG(1, 0) = inner_prod(E2, Gab_contravariant_1);
        mG(1, 1) = inner_prod(E2, Gab_contravariant_2);

        mGVector[PointNumber] = ZeroMatrix(2, 2);
        // saving the G matrix for this point number
        noalias(mGVector[PointNumber]) = mG;

        // Calculate the reduced mass matrix
        mDetJ0[PointNumber] = norm_2(G3);

        // Calculating the total area
        mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
    }
}

//***********************************************************************************
//***********************************************************************************

void PrestressMembraneElement::ComputeContravariantBaseVectors(
                        array_1d<double, 3>& rG1Contra, 
                        array_1d<double, 3>& rG2Contra,
                        const unsigned int& rPointNumber){
    // determinant metric
    double det_metric = mGab0[rPointNumber][0]*mGab0[rPointNumber][1]- mGab0[rPointNumber][2]*mGab0[rPointNumber][2];
    
    // contravariant metric
    array_1d<double, 3> metric_contra;
    metric_contra[0]= 1.0/det_metric * mGab0[rPointNumber][1];
    metric_contra[1]= 1.0/det_metric * mGab0[rPointNumber][0];
    metric_contra[2]= -1.0/det_metric * mGab0[rPointNumber][2];

    // contravariant base vectors
    rG1Contra = metric_contra[0]*mG1[rPointNumber] + metric_contra[2]*mG2[rPointNumber];
    rG2Contra = metric_contra[2]*mG1[rPointNumber] + metric_contra[1]*mG2[rPointNumber];
}
//***********************************************************************************
//***********************************************************************************

const Matrix PrestressMembraneElement::CalculateDeformationGradient(const unsigned int& rPointNumber)
{
    // Compute contravariant base vectors in reference configuration
    array_1d<double, 3> G1_contra, G2_contra;
    ComputeContravariantBaseVectors(G1_contra, G2_contra, rPointNumber);

    // Compute G3
    array_1d<double, 3> G3;
    MathUtils<double>::CrossProduct(G3, mG1[rPointNumber], mG2[rPointNumber]);
    G3 = G3/ norm_2(G3);

    // Compute g1, g2, g3
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = GetGeometry().ShapeFunctionsLocalGradients();
    Matrix DN_De = DN_DeContainer[rPointNumber];
    array_1d<double, 3> g1, g2, g3, gab;
    CalculateMetricDeformed(rPointNumber, DN_De, gab, g1, g2);

    MathUtils<double>::CrossProduct(g3,g1,g2);
    g3 /= norm_2(g3);
    
    bounded_matrix<double,3,3> deformation_gradient;
    for(unsigned int i=0; i<3; i++){
        for(unsigned int j=0; j<3; j++){
            deformation_gradient(i,j) = G1_contra(j)*g1(i) + G2_contra(j)*g2(i) + G3(j)*g3(i);
        }
    }
    return deformation_gradient;
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

    /* @Anna please use these checks
    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(this->GetProperties().Has( CONSTITUTIVE_LAW )) << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const unsigned int strain_size = this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    if ( dimension == 2 ) {
        KRATOS_ERROR_IF( strain_size < 3 || strain_size > 4) << "Wrong constitutive law used. This is a 2D element! expected strain size is 3 or 4 (el id = ) " << this->Id() << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT(strain_size == 6) << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) "<<  this->Id() << std::endl;
    }

    // Check constitutive law
    if ( mConstitutiveLawVector.size() > 0 ) {
        return mConstitutiveLawVector[0]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );
    }
    */

    //verify that the constitutive law exists
    /*if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
    {
        KRATOS_THROW_ERROR(std::logic_error, "constitutive law not provided for property ", this->GetProperties().Id())
    }


    //verify that the constitutive law has the correct dimension
    if (this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize() != 3)
        KRATOS_THROW_ERROR(std::logic_error, "wrong constitutive law used. This is a 3D element with expected strain size is 3 (el id = ) ", this->Id())

    //check constitutive law
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
        mConstitutiveLawVector[i]->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);

        ConstitutiveLaw::Features LawFeatures;
        mConstitutiveLawVector[i]->GetLawFeatures(LawFeatures);

        if (LawFeatures.mOptions.IsNot(ConstitutiveLaw::PLANE_STRESS_LAW))
            KRATOS_THROW_ERROR(std::logic_error, "Constitutive law is compatible only with a plane stress 2D law for membrane element with Id", this->Id())

            if (LawFeatures.mOptions.IsNot(ConstitutiveLaw::INFINITESIMAL_STRAINS))
                KRATOS_THROW_ERROR(std::logic_error, "Constitutive law is compatible only with a law using infinitessimal strains for membrane element with Id", this->Id())

                if (LawFeatures.mStrainSize != 3) KRATOS_THROW_ERROR(std::logic_error, "Constitutive law expects a strain size different from 3 for membrane element with Id", this->Id())
    }
    */
    return 0;

    KRATOS_CATCH("");
}

void PrestressMembraneElement::save(Serializer& rSerializer) const
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
      rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
      rSerializer.save("DetJ0", mDetJ0);
      rSerializer.save("TotalDomainInitialSize", mTotalDomainInitialSize);
      rSerializer.save("StrainsVector", mStrainsVector);
      rSerializer.save("StressesVector", mStressesVector);
      rSerializer.save("CauchyStressesVector", mCauchyStressesVector);
      rSerializer.save("G1", mG1);
      rSerializer.save("G2", mG2);
      rSerializer.save("G_ab", mGab0);
      rSerializer.save("G_Vector", mGVector);
      rSerializer.save("Step", mStep);
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
      rSerializer.load("StrainsVector", mStrainsVector);
      rSerializer.load("StressesVector", mStressesVector);
      rSerializer.load("CauchyStressesVector", mCauchyStressesVector);
      rSerializer.load("G1", mG1);
      rSerializer.load("G2", mG2);
      rSerializer.load("G_ab", mGab0);
      rSerializer.load("G_Vector", mGVector);
      rSerializer.load("Step", mStep);
      rSerializer.load("AnisotropicPrestress", mAnisotropicPrestress);
      rSerializer.load("G1Initial", mG1Initial);
      rSerializer.load("G2Initial", mG2Initial);
      rSerializer.load("G3Initial", mG3Initial);
    }



//***********************************************************************************
//***********************************************************************************
} // Namespace Kratos.
