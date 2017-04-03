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
#include "custom_elements/membrane1_element.hpp"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

// Constructor
Membrane1Element::Membrane1Element( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
}

// Constructor
Membrane1Element::Membrane1Element(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element( NewId, pGeometry, pProperties )
{
}

//***********************************************************************************
//***********************************************************************************

Element::Pointer Membrane1Element::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const

{
    return boost::make_shared< Membrane1Element >(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

//***********************************************************************************
//***********************************************************************************
// Destructor
Membrane1Element::~Membrane1Element()
{
}

//***********************************************************************************
//***********************************************************************************

void Membrane1Element::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)

{
  KRATOS_TRY;

  unsigned int number_of_nodes = GetGeometry().size();
  unsigned int dim = number_of_nodes * 3;

  if (rResult.size() != dim)
      rResult.resize(dim);

  for (unsigned int i = 0; i < number_of_nodes; ++i)
  {
      int index = i * 3;
      rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
      rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
      rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
  }

  Vector equationIdVector;
  equationIdVector.resize(rResult.size());
  for (unsigned int i = 0; i < rResult.size(); ++i)
      equationIdVector[i] = rResult[i];
  //KRATOS_WATCH(equationIdVector);

  KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void Membrane1Element::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo)

{

    // TODO refered to vms.cpp in FluidDynamicsApplication

    ElementalDofList.resize(0);

    for (unsigned int i = 0; i < GetGeometry().size(); ++i)
    {
        ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
    }
}

//***********************************************************************************
//***********************************************************************************

void Membrane1Element::Initialize()

{

    //std::cout << "Membrane1Element::Initialize(): start" << this->GetId() << std::endl;
    KRATOS_TRY

    // reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = GetGeometry().ShapeFunctionsLocalGradients();

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

    // Initialize Variables
    mdensity = GetProperties()[DENSITY];
    mThickness0 = GetProperties()[THICKNESS];
    mThickness = 0.00;

    //// temporary set the pre-stress
    mPreStress[0] = GetProperties()[PRESTRESS_11];
    mPreStress[1] = GetProperties()[PRESTRESS_22];
    mPreStress[2] = GetProperties()[PRESTRESS_12];

    if (mPreStress[0] != mPreStress[1] || mPreStress[2] != 0.0)
        KRATOS_THROW_ERROR(std::invalid_argument, "Only Isotropic Pre-stress state is considered in the membrane1 implementation! Further implementation not done yet!", "");

    
    // resizing jacobian inverses containers
    mDetJ0.resize(integration_points.size());
    mGab0.resize(integration_points.size());
    mG1.resize(integration_points.size());
    mG2.resize(integration_points.size());

    GeometryType::JacobiansType J0;
    J0 = GetGeometry().Jacobian(J0);

    //KRATOS_WATCH(J0);

    mStrainsVector.resize(integration_points.size());
    mStressesVector.resize(integration_points.size());
    mCauchyStressesVector.resize(integration_points.size()); //VM

    mG_Vector.resize(integration_points.size(), ZeroMatrix(2, 2));

    
    // Calculating geometry tensors in reference configuration on Integration points
    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
    {
        // getting information for integration
        double IntegrationWeight = integration_points[PointNumber].Weight();
        //std::cout << "Membrane1Element::Initialize(): 2 " << this->GetId() << std::endl;
        // base vectors
        array_1d<double, 3> g1;
        array_1d<double, 3> g2;
        array_1d<double, 3> g3;

        g1[0] = J0[PointNumber](0, 0);
        g2[0] = J0[PointNumber](0, 1);
        g1[1] = J0[PointNumber](1, 0);
        g2[1] = J0[PointNumber](1, 1);
        g1[2] = J0[PointNumber](2, 0);
        g2[2] = J0[PointNumber](2, 1);


        // Store base vectors in reference configuration
        mG1[PointNumber] = g1;
        mG2[PointNumber] = g2;

        //KRATOS_WATCH(J0);

        // base vector g3
        CrossProduct(g3, g1, g2);
        // differential area dA
        double dA = norm_2(g3);
        // normal vector n
        array_1d<double, 3> n = g3 / dA;

        // Get CovariantMetric
        mGab0[PointNumber][0] = pow(g1[0], 2) + pow(g1[1], 2) + pow(g1[2], 2);
        mGab0[PointNumber][1] = pow(g2[0], 2) + pow(g2[1], 2) + pow(g2[2], 2);
        mGab0[PointNumber][2] = g1[0] * g2[0] + g1[1] * g2[1] + g1[2] * g2[2];

        double inverse_determinant_gab = 1.0 / (mGab0[PointNumber][0] * mGab0[PointNumber][1] - mGab0[PointNumber][2] * mGab0[PointNumber][2]);

        array_1d<double, 3> gab_contravariant;
        gab_contravariant[0] = inverse_determinant_gab*mGab0[PointNumber][1];
        gab_contravariant[1] = -inverse_determinant_gab*mGab0[PointNumber][2];
        gab_contravariant[2] = inverse_determinant_gab*mGab0[PointNumber][0];

        array_1d<double, 3> gab_contravariant_1 = g1*gab_contravariant[0] + g2*gab_contravariant[1];
        array_1d<double, 3> gab_contravariant_2 = g1*gab_contravariant[1] + g2*gab_contravariant[2];
        
        // build local cartesian coordinate system
        double lg1 = norm_2(g1);
        array_1d<double, 3> e1 = g1 / lg1;
        double lg_contravariant_2 = norm_2(gab_contravariant_2);
        array_1d<double, 3> e2 = gab_contravariant_2 / lg_contravariant_2;

        boost::numeric::ublas::bounded_matrix<double, 2, 2> mG;
        //boost::numeric::ublas::bounded_matrix<double, 2, 3> mG;
        mG(0, 0) = inner_prod(e1, gab_contravariant_1);
        mG(0, 1) = inner_prod(e1, gab_contravariant_2);
        mG(1, 0) = inner_prod(e2, gab_contravariant_1);
        mG(1, 1) = inner_prod(e2, gab_contravariant_2);

        mG_Vector[PointNumber] = ZeroMatrix(2, 2);
        // saving the G matrix for this point number
        noalias(mG_Vector[PointNumber]) = mG;

        // Calculate the reduced mass matrix
        mDetJ0[PointNumber] = norm_2(g3);

        // Calculating the total area
        mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;

    }

    // Initialize Material

    if (mConstitutiveLawVector.size() == 0)
    {
        mConstitutiveLawVector.resize(integration_points.size());

        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
        {
            /*    ConstitutiveLaw<Node<3> >::Pointer material = ConstitutiveLaw<Node<3> >::Pointer( new Isotropic2D() );
            mConstitutiveLawVector[i] = material;*/


            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();

            mConstitutiveLawVector[i]->InitializeMaterial(GetProperties(), GetGeometry(), row(GetGeometry().ShapeFunctionsValues(), i));
        }
    }

    KRATOS_CATCH( "" )
}



//***********************************************************************************
//***********************************************************************************

void Membrane1Element::CalculateRightHandSide(
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

void Membrane1Element::CalculateLocalSystem(
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

void Membrane1Element::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& Output,
    const ProcessInfo& rCurrentProcessInfo)

{
    KRATOS_ERROR << "CalculateOnIntegrationPoints not implemented yet" << std::endl;
    //std::cout << "CalculateOnIntegrationPoints not implemented yet" << std::endl;
}

//***********************************************************************************
//***********************************************************************************

void Membrane1Element::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo)

{
    KRATOS_TRY

    //std::cout << "Membrane1Element::CalculateMassMatrix : this is the same as MembraneElement" << std::endl;
    //rMassMatrix.resize(0,0);
    // LUMPED MASS MATRIX
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;

    if (rMassMatrix.size1() != MatSize)
        rMassMatrix.resize(MatSize, MatSize);

    rMassMatrix = ZeroMatrix(MatSize, MatSize);

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

void Membrane1Element::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo)

{
    KRATOS_TRY

    // LUMPED DAMPING MATRIX
    std::cout << "Membrane1Element::CalculateDampingMatrix : this is the same as MembraneElement" << std::endl;
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;

    if (rDampingMatrix.size1() != MatSize)
        rDampingMatrix.resize(MatSize, MatSize);

    rDampingMatrix = ZeroMatrix(MatSize, MatSize);

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

void Membrane1Element::FinalizeSolutionStep(
    ProcessInfo& rCurrentProcessInfo)
{
    for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
    {
        //
        //            ConstitutiveLaw::Parameters Values (GetGeometry(),GetProperties(),rCurrentProcessInfo);
        //            Values.GetOptions().Set (ConstitutiveLaw::COMPUTE_STRAIN, false);
        //            Values.GetOptions().Set (ConstitutiveLaw::COMPUTE_STRESS);
        //            Values.GetOptions().Set (ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        //            Matrix dummy = ZeroMatrix ( 0, 0 );
        //            Vector StrainVector = mStrainsVector[i];
        //            Values.SetStrainVector (StrainVector); //this has to be the input parameter
        //            Values.SetStressVector (StressVector);
        //            Values.SetConstitutiveMatrix (dummy);
        //            Values.SetShapeFunctionsValues ( row ( GetGeometry().ShapeFunctionsValues(), PointNumber ) );
        //
        //            mConstitutiveLawVector[PointNumber]->FinalizeMaterialResponse (Values,ConstitutiveLaw::StressMeasure_PK2 );

        mConstitutiveLawVector[i]->FinalizeSolutionStep(GetProperties(),
            GetGeometry(),
            row(GetGeometry().ShapeFunctionsValues(), i),
            rCurrentProcessInfo);
    }
}

//##### From here, old code

//***********************************************************************************
//***********************************************************************************

void Membrane1Element::GetValuesVector(
    Vector& values,
    int Step)

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int MatSize = number_of_nodes * 3;

    if (values.size() != MatSize)
        values.resize(MatSize);

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

void Membrane1Element::GetFirstDerivativesVector(
    Vector& values,
    int Step)

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int MatSize = number_of_nodes * 3;

    if (values.size() != MatSize)
        values.resize(MatSize);

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

void Membrane1Element::GetSecondDerivativesVector(
    Vector& values,
    int Step)

{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int MatSize = number_of_nodes * 3;

    if (values.size() != MatSize)
        values.resize(MatSize);

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


void Membrane1Element::CrossProduct(
    array_1d<double, 3>& cross,
    array_1d<double, 3>& a,
    array_1d<double, 3>& b)

{
    cross[0] = a[1] * b[2] - a[2] * b[1];
    cross[1] = a[2] * b[0] - a[0] * b[2];
    cross[2] = a[0] * b[1] - a[1] * b[0];
}

//***********************************************************************************
//***********************************************************************************

void Membrane1Element::CalculateAndAddKm(
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

void Membrane1Element::CalculateAndAddNonlinearKm(
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
                int check = 3;
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
void Membrane1Element::ClearNodalForces()
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

void Membrane1Element::CalculateQ(
    boost::numeric::ublas::bounded_matrix<double, 3, 3>& Q,
    Matrix& mG)

{
    KRATOS_TRY
    Q(0, 0) = pow(mG(0, 0), 2);
    Q(0, 1) = pow(mG(0, 1), 2);
    Q(0, 2) = 2.00*mG(0, 0)*mG(0, 1);

    Q(1, 0) = pow(mG(1, 0), 2);
    Q(1, 1) = pow(mG(1, 1), 2);
    Q(1, 2) = 2.00*mG(1, 0) * mG(1, 1);

    Q(2, 0) = 2.00 * mG(0, 0) * mG(1, 0);
    Q(2, 1) = 2.00 * mG(0, 1)*mG(1, 1);
    Q(2, 2) = 2.00 * (mG(0, 0) * mG(1, 1) + mG(0, 1)*mG(1, 0));

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void Membrane1Element::CalculateB(
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

void Membrane1Element::CalculateStrain(
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

void Membrane1Element::CalculateAndAdd_BodyForce(
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

void Membrane1Element::CalculateAndAdd_PressureForce(
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

void Membrane1Element::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,
    bool CalculateResidualVectorFlag)

{
    KRATOS_TRY

    //std::cout << "Membrane1Element::CalculateAll" << std::endl;

    
    // Initializing all needed variables
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int MatSize = number_of_nodes * 3;
    
    // set up Constitutive Law
    ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    //Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRAIN, false);        // why set to false?
    //Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
    Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR); 
    

    // Matrix B(3, MatSize);
    // change to: array_1d<double, 3> StrainVector;
    Vector StrainVector(3);
    Vector StressVector(3);

    // resizing as needed the LHS
    if (CalculateStiffnessMatrixFlag == true)    // calculation of the matrix is required
    {
        if (rLeftHandSideMatrix.size1() != MatSize)
            rLeftHandSideMatrix.resize(MatSize, MatSize);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize); //resetting LHS
    }

    //resizing as needed the RHS
    if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
    {
        if (rRightHandSideVector.size() != MatSize)
            rRightHandSideVector.resize(MatSize);

        rRightHandSideVector = ZeroVector(MatSize); //resetting RHS
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
   
    //KRATOS_WATCH(integration_points.size());
    bool is_initialized = false;
    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
    {
        // reading integration weight, shape function value and its gradients at this integration point
        double IntegrationWeight = GetGeometry().IntegrationPoints()[PointNumber].Weight();

        double coord0 = GetGeometry().IntegrationPoints()[PointNumber].Coordinate(0);
        double coord1 = GetGeometry().IntegrationPoints()[PointNumber].Coordinate(1);
        double coord2 = GetGeometry().IntegrationPoints()[PointNumber].Coordinate(2);

        /*
        KRATOS_WATCH(integration_points.size());
        KRATOS_WATCH(PointNumber);
        KRATOS_WATCH(IntegrationWeight);
        KRATOS_WATCH(coord0);
        KRATOS_WATCH(coord1);
        KRATOS_WATCH(coord2);
        */

        Vector ShapeFunctionN = row(Ncontainer, PointNumber);
        Matrix DN_De = DN_DeContainer[PointNumber];

        /*
        KRATOS_WATCH(ShapeFunctionN);
        KRATOS_WATCH(DN_De);
        */

        // covariant metric in deformed system
        array_1d<double, 3> gab;
        
        bool formfingding_flag = true;

        // Transformation Matrix Q
        boost::numeric::ublas::bounded_matrix<double, 3, 3> Q = ZeroMatrix(3, 3);
        // basis vectors in deformed system
        array_1d<double, 3> g1;
        array_1d<double, 3> g2;
        array_1d<double, 3> g3;


        //if (!is_initialized)
        //{
        //    this->Initialize();
        //    is_initialized = true;
        //}

        CalculateQ(Q, mG_Vector[PointNumber]);
        CalculateMetricDeformed(PointNumber, DN_De, gab, g1, g2);
        CalculateStrain(StrainVector, gab, mGab0[PointNumber]);

        //if (this->GetId() == 1)
        //{
        //    KRATOS_WATCH(Q);
        //    KRATOS_WATCH(g1);
        //    KRATOS_WATCH(g2);
        //    KRATOS_WATCH(StrainVector);
        //}

        Vector CartesianStrainVector = prod(Q, StrainVector);

        // Constitutive Matrices D
        Matrix D = ZeroMatrix(3, 3);
        //boost::numeric::ublas::bounded_matrix<double, 3, 3> D;
        //KRATOS_WATCH(StressVector);

        
        Values.SetStrainVector(StrainVector);       // this is the input parameter
        Values.SetStressVector(StressVector);       // this is an output parameter
        Values.SetConstitutiveMatrix(D);            // this is an output parameter

        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);     // Why is the curviliear strains are used here?
        /**/
        //KRATOS_WATCH(StressVector);

        // Deformations for Non-linear force vector
        Vector StrainDeformation;

        //if (this->GetId() == 42)
        //{
        //    std::cout << "before my function!" << std::endl;
        //    KRATOS_WATCH(D);
        //}

        StrainDeformation = prod(trans(D), CartesianStrainVector);

        //KRATOS_WATCH(StrainDeformation);
        array_1d<double,3> pre_stress_tensor;   // Vector with the Cauchy Pre-Stress components in local cartesian frame
        
        // Getting the prestress values        
        pre_stress_tensor(0) = mPreStress[0];
        pre_stress_tensor(1) = mPreStress[1];
        pre_stress_tensor(2) = mPreStress[2];


        // !!!! the thickness is considered at the end by the IntToReferenceWeight !!!
        // pre-integration with thickness
        //pre_stress_tensor *= mThickness0;

        //KRATOS_WATCH(pre_stress_tensor);
        //KRATOS_WATCH(StrainDeformation);


        // taking out the pre-stress
        // Adding the pre-stress values as forces over length
        StrainDeformation += pre_stress_tensor;


        //std::cout << "After adding the pre-stress tensor" << std::endl;
        //KRATOS_WATCH(StrainDeformation);


        // calculate B matrices
        // B matrices:
        Matrix B = ZeroMatrix(3, MatSize);
        //CalculateB(B, Q, DN_De, mG1[PointNumber], mG2[PointNumber]);
        CalculateB(B, Q, DN_De, g1, g2);

        // integration on the REFERENCE CONFIGURATION
        double DetJ0 = mDetJ0[PointNumber];
        double IntToReferenceWeight = IntegrationWeight * DetJ0 * mThickness0;

        // Nonlinear Deformation
        Matrix Strain_locCartesian_11 = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);
        Matrix Strain_locCartesian_22 = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);
        Matrix Strain_locCartesian_12 = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);

        CalculateSecondVariationStrain(DN_De, Strain_locCartesian_11, Strain_locCartesian_22, Strain_locCartesian_12, Q, g1, g2);

        //if (this->GetId() == 1)
        //{

        //    KRATOS_WATCH(Strain_locCartesian_11);
        //    KRATOS_WATCH(Strain_locCartesian_22);
        //    KRATOS_WATCH(Strain_locCartesian_12);
        //}


        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true)
        {
            // adding membrane contribution to the stiffness matrix
            //KRATOS_WATCH(B);
            //KRATOS_WATCH(D);
            //KRATOS_WATCH(IntToReferenceWeight);
            CalculateAndAddKm(rLeftHandSideMatrix, B, D, IntToReferenceWeight);

            //if (this->GetId() == 42)
            //{
            //    std::cout << "Klinear" << std::endl;
            //    KRATOS_WATCH(D);
            //    KRATOS_WATCH(rLeftHandSideMatrix);
            //}
            // adding non-linear-contribution to stiffness matrix
            CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
                Strain_locCartesian_11, Strain_locCartesian_22, Strain_locCartesian_12,
                StrainDeformation,
                IntToReferenceWeight);

            //if (this->GetId() == 42)
            //{
            //    std::cout << "KNonlinear" << std::endl;
            //    KRATOS_WATCH(rLeftHandSideMatrix);
            //}
        }

        // RIGHT HAND SIDE VECTOR
        if (CalculateResidualVectorFlag == true)         // calculation of the matrix is required
        {
            // operation performed: rRighthandSideVector -= Weight* IntForce
            noalias(rRightHandSideVector) -= IntToReferenceWeight* prod(trans(B), StrainDeformation);
        }
    } // end loop over integration points
    //if(CalculateStiffnessMatrixFlag==true)
    //    KRATOS_WATCH(rLeftHandSideMatrix);

    //if(CalculateResidualVectorFlag==true)
    //    KRATOS_WATCH(rRightHandSideVector);

    // DEBUG INFO
    //if (this->GetId() == 1 /*|| this->GetId() == 2*/)
    //{
    //    //KRATOS_WATCH(rLeftHandSideMatrix);
    //    KRATOS_WATCH(rRightHandSideVector);
    //}

    //if (this->GetId() == 42)
    //{
    //    std::cout << "K" << std::endl;
    //    KRATOS_WATCH(rLeftHandSideMatrix);
    //}
    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************
//
//void Membrane1Element::Calculate_GlobalStressVector(
//    array_1d<double, 6>& GlobalVector,
//    Vector& LocalStressVector,
//    array_1d<double, 3>& v1,
//    array_1d<double, 3>& v2)
//
//{
//    KRATOS_TRY
//
//        array_1d<double, 6> temp;
//
//    //adding the component S11
//    noalias(temp) = VoigtTensorComponents(v1, v1);
//    temp *= LocalStressVector[0];
//    noalias(GlobalVector) += temp;
//
//    //adding the component S22
//    noalias(temp) = VoigtTensorComponents(v2, v2);
//    temp *= LocalStressVector[1];
//    noalias(GlobalVector) += temp;
//
//    //adding the component S12 (& S21)
//    noalias(temp) = VoigtTensorComponents(v1, v2);
//    noalias(temp) += VoigtTensorComponents(v2, v1);
//    temp *= LocalStressVector[2];
//    noalias(GlobalVector) += temp;
//
//    KRATOS_CATCH("")
//}
//
////***********************************************************************************
////***********************************************************************************
//
////auxiliary function needed in the calculation of output stresses
//inline array_1d<double, 6> Membrane1Element::VoigtTensorComponents(
//    array_1d<double, 3>& a,
//    array_1d<double, 3>& b)
//
//{
//    array_1d<double, 6> v;
//
//    v[0] = a[0] * b[0];
//    v[1] = a[1] * b[1];
//    v[2] = a[2] * b[2];
//    v[3] = a[0] * b[1];
//    v[4] = a[1] * b[2];
//    v[5] = a[0] * b[2];
//
//    return v;
//}

//************************************************************************************
//************************************************************************************
void Membrane1Element::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
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
void Membrane1Element::CalculateMetricDeformed(unsigned int& PointNumber, Matrix DN_De,
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
    CrossProduct(g3, g1, g2);
    //differential area dA
    double dA = norm_2(g3);
    //normal vector _n
    array_1d<double, 3> n = g3 / dA;

    //GetCovariantMetric
    gab[0] = pow(g1[0], 2) + pow(g1[1], 2) + pow(g1[2], 2);
    gab[1] = pow(g2[0], 2) + pow(g2[1], 2) + pow(g2[2], 2);
    gab[2] = g1[0] * g2[0] + g1[1] * g2[1] + g1[2] * g2[2];
}

//***********************************************************************************
//***********************************************************************************
void Membrane1Element::CalculateSecondVariationStrain(Matrix DN_De,
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

                int limit = i + 1;
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
void Membrane1Element::CalculateMembraneElasticityTensor(
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

int  Membrane1Element::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

        //verify that the variables are correctly initialized

        if (VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "");

            if (DISPLACEMENT.Key() == 0)
                KRATOS_THROW_ERROR(std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "")

                if (ACCELERATION.Key() == 0)
                    KRATOS_THROW_ERROR(std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "")

                    if (DENSITY.Key() == 0)
                        KRATOS_THROW_ERROR(std::invalid_argument, "DENSITY has Key zero! (check if the application is correctly registered", "")

                        if (VOLUME_ACCELERATION.Key() == 0)
                            KRATOS_THROW_ERROR(std::invalid_argument, "BODY_FORCE has Key zero! (check if the application is correctly registered", "")

                            if (THICKNESS.Key() == 0)
                                KRATOS_THROW_ERROR(std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "")

                                //verify that the dofs exist
                                for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
                                {
                                    if (this->GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false)
                                        KRATOS_THROW_ERROR(std::invalid_argument, "missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id())

                                        if (this->GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false || this->GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false)
                                            KRATOS_THROW_ERROR(std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id())
                                }

    //verify that the constitutive law exists
    if (this->GetProperties().Has(CONSTITUTIVE_LAW) == false)
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

    return 0;

    KRATOS_CATCH("");
}



//***********************************************************************************
//***********************************************************************************
} // Namespace Kratos.
