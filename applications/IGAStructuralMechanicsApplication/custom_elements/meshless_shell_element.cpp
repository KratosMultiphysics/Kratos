//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Riccardo Rossi
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/meshless_base_element.h"
#include "custom_elements/meshless_membrane_element.h"

#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"

#include "utilities/math_utils.h"

#include "geometries/geometry.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
MeshlessShellElement::MeshlessShellElement(
	IndexType NewId,
	GeometryType::Pointer pGeometry)
    : MeshlessBaseElement(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//************************************************************************************
//************************************************************************************
MeshlessShellElement::MeshlessShellElement(
	IndexType NewId,
	GeometryType::Pointer pGeometry,  
	PropertiesType::Pointer pProperties)
    : MeshlessBaseElement(NewId, pGeometry, pProperties)
{
}

Element::Pointer MeshlessShellElement::Create(
	IndexType NewId, 
	NodesArrayType const& ThisNodes,  
	PropertiesType::Pointer pProperties) const
{
	return boost::make_shared< MeshlessShellElement >(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

MeshlessShellElement::~MeshlessShellElement()
{
}

//***********************************************************************************
//***********************************************************************************

void MeshlessShellElement::EquationIdVector(
	EquationIdVectorType& rResult,
	ProcessInfo& rCurrentProcessInfo)

{
	KRATOS_TRY
	unsigned int number_of_nodes = GetGeometry().size();
	unsigned int dim = number_of_nodes * 3;

	if (rResult.size() != dim)
		rResult.resize(dim);

	for (unsigned int i = 0; i < number_of_nodes; i++)
	{
		int index = i * 3;
		rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
		rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
		rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
	}

	KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void MeshlessShellElement::GetDofList(
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

void MeshlessShellElement::Initialize()

{
	KRATOS_TRY
	// Get values of shape functions and derivatives. Derivatives are first and sevon column: first dreivatives. Third, fourth and fith are second derivatives
	double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
	Vector ShapeFunctionsN = this->GetValue(SHAPE_FUNCTION_VALUES);
	Matrix DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
	Matrix DDN_DDe = this->GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

	// Initialize Variables
	mdensity = GetProperties()[DENSITY];
	mThickness0 = GetProperties()[THICKNESS];
	//WIRD DAS GEBRAUCHT??

	mThickness = 0.00;

	// Initialize Material
	mConstitutiveLawVector = GetProperties()[CONSTITUTIVE_LAW]->Clone();
	ProcessInfo emptyProcessInfo = ProcessInfo();
	//mConstitutiveLawVector->SetValue(INTEGRATION_WEIGHT, integration_weight, emptyProcessInfo);
	//mConstitutiveLawVector->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES), emptyProcessInfo);
	mConstitutiveLawVector->InitializeMaterial(GetProperties(), GetGeometry(), ShapeFunctionsN);

	//Initialize Jacobian
	Matrix J0;
	//Calculate Jacobian
	Jacobian(DN_De, J0);

	//Basis Vectors
	array_1d<double, 3> g1, g2, g3;

	g1[0] = J0(0, 0);
	g2[0] = J0(0, 1);
	g1[1] = J0(1, 0);
	g2[1] = J0(1, 1);
	g1[2] = J0(2, 0);
	g2[2] = J0(2, 1);

	//Basis Vector g3
	CrossProduct(g3, g1, g2);
	//Differential Area dA
	double dA = norm_2(g3);
	//Normal Vector n
	array_1d<double, 3> n = g3 / dA;

	//Covariant Metric
	mGab0[0] = pow(g1[0], 2) + pow(g1[1], 2) + pow(g1[2], 2);
	mGab0[1] = pow(g2[0], 2) + pow(g2[1], 2) + pow(g2[2], 2);
	mGab0[2] = g1[0] * g2[0] + g1[1] * g2[1] + g1[2] * g2[2];

	double inverse_determinant_gab = 1.0 / (mGab0[0] * mGab0[1] - mGab0[2] * mGab0[2]);

	//Contravariant Metric
	array_1d<double, 3> gab_contravariant;
	gab_contravariant[0] = inverse_determinant_gab*mGab0[1];
	gab_contravariant[1] = -inverse_determinant_gab*mGab0[2];
	gab_contravariant[2] = inverse_determinant_gab*mGab0[0];

	array_1d<double, 3> g_contravariant_1 = g1*gab_contravariant[0] + g2*gab_contravariant[1];
	array_1d<double, 3> g_contravariant_2 = g1*gab_contravariant[1] + g2*gab_contravariant[2];

	//Hessian Matrix
	Matrix H = ZeroMatrix(3, 3);
	Hessian(H, DDN_DDe);

	mCurvature0[0] = H(0, 0)*n[0] + H(1, 0)*n[1] + H(2, 0)*n[2];
	mCurvature0[1] = H(0, 1)*n[0] + H(1, 1)*n[1] + H(2, 1)*n[2];
	mCurvature0[2] = H(0, 2)*n[0] + H(1, 2)*n[1] + H(2, 2)*n[2];

	double length_g1 = norm_2(g1);
	array_1d<double, 3> e1 = g1 / length_g1;
	double lg_contravariant_2 = norm_2(g_contravariant_2);
	array_1d<double, 3> e2 = g_contravariant_2 / lg_contravariant_2;

	Matrix mG = ZeroMatrix(2, 2);
	mG(0, 0) = inner_prod(e1, g_contravariant_1);
	mG(0, 1) = inner_prod(e1, g_contravariant_2);
	mG(1, 0) = inner_prod(e2, g_contravariant_1);
	mG(1, 1) = inner_prod(e2, g_contravariant_2);

	boost::numeric::ublas::bounded_matrix<double, 3, 3> Q = ZeroMatrix(3,3);
	CalculateQ(Q, mG);
	mQ = Q;

	//Calculate the reduced mass matrix
	mDetJ0 = norm_2(g3); //norm_2(V3);


	//calculating the total area
	//mTotalDomainInitialSize = mDetJ0 * integration_weight;

	KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************
void MeshlessShellElement::CalculateRightHandSide(
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

void MeshlessShellElement::CalculateLocalSystem(
	MatrixType& rLeftHandSideMatrix,
	VectorType& rRightHandSideVector,
	ProcessInfo& rCurrentProcessInfo)

{
	//calculation flags
	bool CalculateStiffnessMatrixFlag = true;
	bool CalculateResidualVectorFlag = true;

	CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************
void MeshlessShellElement::CalculateOnIntegrationPoints(
  const Variable<double>& rVariable,
  std::vector<double>& rOutput,
  const ProcessInfo& rCurrentProcessInfo
  )
{
  if (rOutput.size() != 1)
  {
    rOutput.resize(1);
  }

  if (rVariable == VON_MISES_STRESS)
  {
    //reading in of integration weight, shape function values and shape function derivatives
    double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
    Vector ShapeFunctionsN = this->GetValue(SHAPE_FUNCTION_VALUES);
    Matrix DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
    Matrix DDN_DDe = this->GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = 3; // get value from dimension of derivatives
    const unsigned int strain_size = 3;// mConstitutiveLawVector[0]->GetStrainSize();

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions = Values.GetOptions();
	///////////////////////////////////////////////////make back
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    //Values.SetStrainVector(this_constitutive_variables.StrainVector);

// StrainStress
    // covariant metric in deformed system
    // curvature covariant metric in deformed system
    array_1d<double, 3> curvature;

    // Transformation Matrix Q - 
    boost::numeric::ublas::bounded_matrix<double, 3, 3>  Q = mQ;

    Vector StrainVector = ZeroVector(3);
    Vector CurvatureVector = ZeroVector(3);

    // basis vectors in deformed system
    //array_1d<double, 3> g1, g2;
    MetricVariables DeformedMetric(3);
    CalculateMetricDeformed(DN_De, DDN_DDe, DeformedMetric);
    // covariant metric in deformed system
	array_1d<double, 3> gab;
	gab[0] = DeformedMetric.gab[0];
	gab[1] = DeformedMetric.gab[1];
	gab[2] = DeformedMetric.gab[2];


    CalculateStrain(StrainVector, gab, mGab0);

    this_constitutive_variables.StrainVector = prod(Q, StrainVector);
	//////VERMUTLICH FALSCH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CalculateCurvature(CurvatureVector, curvature, mCurvature0);
    this_constitutive_variables.CurvatureVector = prod(Q, CurvatureVector);

    //Constitive Matrices DMembrane and DCurvature
    //Matrix DMembrane = ZeroMatrix(3, 3);
    //Matrix DCurvature = ZeroMatrix(3, 3);

    Values.SetStrainVector(StrainVector); //this is the input parameter
    Vector StressVector;
    Values.SetStressVector(StressVector); //this is an ouput parameter
    Values.SetConstitutiveMatrix(this_constitutive_variables.DMembrane); //this is an ouput parameter
                                             //KRATOS_WATCH(DMembrane)

    mConstitutiveLawVector->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);
    //KRATOS_WATCH(DMembrane)
    double thickness = this->GetProperties().GetValue(THICKNESS);
    //KRATOS_WATCH(this_constitutive_variables.DMembrane)
    this_constitutive_variables.DCurvature = this_constitutive_variables.DMembrane*(pow(thickness, 2) / 12);

    //Local Cartesian Foreces and Moments
    //Vector ForceVector_in_Q_coordinates = prod(trans(DMembrane), StrainVector_in_Q_coordinates);
    //Vector MomentVector_in_Q_coordinates = prod(trans(DCurvature), CurvatureVector_in_Q_coordinates);
    //Vector ForceVector_in_Q_coordinates = ZeroVector(3);
    //Vector MomentVector_in_Q_coordinates = ZeroVector(3);


//end Strain Stress

    double detF = DeformedMetric.dA / mDetJ0;

//elemental stiffness

    double youngs_modulus = this->GetProperties().GetValue(YOUNG_MODULUS);
    double poisson_ratio = this->GetProperties().GetValue(POISSON_RATIO);
    Matrix dm = ZeroMatrix(3, 3);
    Matrix db = ZeroMatrix(3, 3);
    dm(0, 0) = 1.0;
    dm(0, 1) = poisson_ratio;
    dm(1, 0) = poisson_ratio;
    dm(1, 1) = 1.0;
    dm(2, 2) = (1.0 - poisson_ratio) / 2.0;
    db = youngs_modulus*pow(thickness, 3) / (12.0*(1.0 - poisson_ratio*poisson_ratio))*dm;
    dm = youngs_modulus*thickness / (1.0 - poisson_ratio*poisson_ratio)*dm;

    Vector n_pk2_ca = prod(dm, this_constitutive_variables.StrainVector);

    Matrix T_E_G = ZeroMatrix(3, 3);
    T_E_G(0, 0) = mQ(0, 0);
    T_E_G(0, 1) = mQ(1, 0);
    T_E_G(0, 2) = mQ(2, 0);
    T_E_G(1, 0) = mQ(0, 1);
    T_E_G(1, 1) = mQ(1, 1);
    T_E_G(1, 2) = mQ(2, 1);
    T_E_G(2, 0) = mQ(0, 2)*0.5;
    T_E_G(2, 1) = mQ(1, 2)*0.5;
    T_E_G(2, 2) = mQ(2, 2)*0.5;

    Vector n_pk2 = prod(T_E_G, n_pk2_ca);
    //array_1d<double, 3> n_pk2 = prod(T_E_G, n_pk2_ca);
    // Cauchy normal force in g1,g2
    array_1d<double, 3> n_cau = 1.0 / detF*n_pk2; ///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                  //c_vector<cfloat,3> n_cau = n_pk2;
    Vector n = ZeroVector(7);
    // Cauchy normal force in normalized g1,g2
    n[0] = sqrt(DeformedMetric.gab[0] / DeformedMetric.gab_con[0])*n_cau[0];
    n[1] = sqrt(DeformedMetric.gab[1] / DeformedMetric.gab_con[1])*n_cau[1];
    n[2] = sqrt(DeformedMetric.gab[0] / DeformedMetric.gab_con[1])*n_cau[2];
    // Cauchy normal force in local cartesian e1,e2
    array_1d<double, 3> n_e = prod(DeformedMetric.T, n_cau);
    n[0] = n_e[0];
    n[1] = n_e[1];
    n[2] = n_e[2];
    // Principal normal forces
    n[3] = 0.5*(n_e[0] + n_e[1] + sqrt(pow(n_e[0] - n_e[1], 2) + 4.0*pow(n_e[2], 2)));
    n[4] = 0.5*(n_e[0] + n_e[1] - sqrt(pow(n_e[0] - n_e[1], 2) + 4.0*pow(n_e[2], 2)));

    // -------------------  moments -------------------------
    // PK2 moment in local cartesian E1,E2
    array_1d<double, 3> m_pk2_ca = prod(db, this_constitutive_variables.CurvatureVector);
    // PK2 moment in G1,G2
    array_1d<double, 3> m_pk2 = prod(T_E_G, m_pk2_ca);
    // Cauchy moment in g1,g2
    array_1d<double, 3> m_cau = 1.0 / detF*m_pk2; ///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                  //c_vector<cfloat,3> m_cau = m_pk2;
                                                  // Cauchy moment in normalized g1,g2
    Vector m = ZeroVector(7);
    m[0] = sqrt(DeformedMetric.gab[0] / DeformedMetric.gab_con[0])*m_cau[0];
    m[1] = sqrt(DeformedMetric.gab[1] / DeformedMetric.gab_con[1])*m_cau[1];
    m[2] = sqrt(DeformedMetric.gab[0] / DeformedMetric.gab_con[1])*m_cau[2];
    // Cauchy moment in local cartesian e1,e2
    array_1d<double, 3> m_e = prod(DeformedMetric.T, m_cau);
    m[0] = m_e[0];
    m[1] = m_e[1];
    m[2] = m_e[2];
    // principal moment
    m[3] = 0.5*(m_e[0] + m_e[1] + sqrt(pow(m_e[0] - m_e[1], 2) + 4.0*pow(m_e[2], 2)));
    m[4] = 0.5*(m_e[0] + m_e[1] - sqrt(pow(m_e[0] - m_e[1], 2) + 4.0*pow(m_e[2], 2)));

    double W = pow(thickness, 2) / 6.0;
    double sigma_1_top = m[0] / W + n[0] / thickness;
    double sigma_2_top = m[1] / W + n[1] / thickness;
    double sigma_3_top = m[2] / W + n[2] / thickness;
    double vMises = pow(pow(sigma_1_top, 2) + pow(sigma_2_top, 2) - sigma_1_top*sigma_2_top + 3 * pow(sigma_3_top, 2), 0.5);

    rOutput[0] = vMises;
    //// Reading integration points
    //const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    //// Displacements vector
    //Vector displacements;
    //GetValuesVector(displacements);

    ////for (unsigned int point_number = 0; point_number < integration_points.size(); point_number++)
    ////{
    //  // Compute element kinematics B, F, DN_DX ...
    //  CalculateKinematicVariables(this_kinematic_variables, point_number, integration_points);

    //  // Compute material reponse
    //  CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure(), displacements);

    //  const Matrix stress_tensor = MathUtils<double>::StressVectorToTensor(this_constitutive_variables.StressVector);

    //  double sigma_equivalent = 0.0;

    //  if (dimension == 2)
    //  {
    //    sigma_equivalent = std::pow((stress_tensor(0, 0) - stress_tensor(1, 1)), 2.0) +
    //      3 * (stress_tensor(0, 1) * stress_tensor(1, 0));
    //  }
    //  else
    //  {
    //    sigma_equivalent = 0.5*(std::pow((stress_tensor(0, 0) - stress_tensor(1, 1)), 2.0) +
    //      std::pow((stress_tensor(1, 1) - stress_tensor(2, 2)), 2.0) +
    //      std::pow((stress_tensor(2, 2) - stress_tensor(0, 0)), 2.0) +
    //      6 * (stress_tensor(0, 1) * stress_tensor(1, 0) +
    //        stress_tensor(1, 2) * stress_tensor(2, 1) +
    //        stress_tensor(2, 0) * stress_tensor(0, 2)));
    //  }

      //if (sigma_equivalent < 0.0)
      //{
      //  rOutput[0] = 0.0;
      //}
      //else
      //{
        //rOutput[0] = 5;// std::sqrt(sigma_equivalent);
      //}
  }
  else if (rVariable == DAMAGE_T)
  {
	  mConstitutiveLawVector->GetValue(DAMAGE_T, rOutput[0]);
  }
  else if (rVariable == DAMAGE_C)
  {
	  mConstitutiveLawVector->GetValue(DAMAGE_C, rOutput[0]);
  }
  else
  {
    rOutput[0] = 0.0;// mConstitutiveLawVector[0]->GetValue(rVariable, rOutput[0]);
  }
}

//***********************************************************************************
//***********************************************************************************

void MeshlessShellElement::FinalizeSolutionStep(
	ProcessInfo& rCurrentProcessInfo)
{
		//
		//            ConstitutiveLaw::Parameters Values (GetGeometry(),GetProperties(),rCurrentProcessInfo);
		//            Values.GetOptions().Set (ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
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

	mConstitutiveLawVector->FinalizeSolutionStep(GetProperties(),
		GetGeometry(),
		this->GetValue(SHAPE_FUNCTION_VALUES),
		rCurrentProcessInfo);

}

//***********************************************************************************
//***********************************************************************************
// --------- //
//  PRIVATE  //
// --------- //
//***********************************************************************************
//***********************************************************************************
void MeshlessShellElement::CalculateAndAddKm(
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
	//KRATOS_WATCH(Km)
	noalias(K) += Km;

	KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void MeshlessShellElement::CalculateAndAddNonlinearKm(
	Matrix& K,
	Matrix& B11,
	Matrix& B22,
	Matrix& B12,
	Vector& SD,
	double weight)

{
	KRATOS_TRY
		
	unsigned int number_of_nodes = GetGeometry().size();

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
					unsigned int nbase = 3 * n + i;
					unsigned int mbase = 3 * m + j;

					K(nbase, mbase) += (SD[0] * B11(nbase, mbase) 
						+ SD[1] * B22(nbase, mbase) 
						+ SD[2] * B12(nbase, mbase))*weight;
					K(mbase, nbase) += (SD[0] * B11(nbase, mbase) 
						+ SD[1] * B22(nbase, mbase) 
						+ SD[2] * B12(nbase, mbase))*weight;
				}
			}
		}
	}
	KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void MeshlessShellElement::CalculateQ(
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
/**doxygen
*/
void MeshlessShellElement::CalculateBMembrane(
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
		b(0, index)     = DN_De(i, 0) * g1[0];
		b(0, index + 1) = DN_De(i, 0) * g1[1];
		b(0, index + 2) = DN_De(i, 0) * g1[2];

		//second line
		b(1, index)     = DN_De(i, 1) * g2[0];
		b(1, index + 1) = DN_De(i, 1) * g2[1];
		b(1, index + 2) = DN_De(i, 1) * g2[2];

		//third line
		b(2, index)     = 0.5*(DN_De(i, 1) * g1[0] + DN_De(i, 0) * g2[0]);
		b(2, index + 1) = 0.5*(DN_De(i, 1) * g1[1] + DN_De(i, 0) * g2[1]);
		b(2, index + 2) = 0.5*(DN_De(i, 1) * g1[2] + DN_De(i, 0) * g2[2]);
	}

	B = prod(Q, b);

	KRATOS_CATCH("")
}


void MeshlessShellElement::CalculateBCurvature(
  Matrix& B,
  boost::numeric::ublas::bounded_matrix<double, 3, 3>& Q,
  const Matrix& DN_De,
  const Matrix& DDN_DDe,
	const array_1d<double, 3>& g1,
	const array_1d<double, 3>& g2)
{
	KRATOS_TRY
	const unsigned int number_of_nodes = GetGeometry().size();
	Matrix dg3 = ZeroMatrix(3, 3);
	Matrix dn = ZeroMatrix(3, 3);
	Matrix b = ZeroMatrix(3, number_of_nodes * 3);

	//calculation of Hessian H
	Matrix H = ZeroMatrix(3, 3);
	Hessian(H, DDN_DDe);

	KRATOS_WATCH(H)
	KRATOS_WATCH(DDN_DDe)
	KRATOS_WATCH(g1)
	KRATOS_WATCH(g2)
	
	//basis vector g3
	array_1d<double, 3> g3;
	array_1d<double, 3> n;
	CrossProduct(g3, g1, g2);
	KRATOS_WATCH(g3)

	//differential area dA
	double dA = norm_2(g3);

	//normal vector n
	n = g3 / dA;

	double invdA = 1 / dA;
	double inddA3 = 1 / pow(dA, 3);

	for (unsigned int i = 0; i < number_of_nodes; i++)
	{
		unsigned int index = 3 * i;
		//first line
		dg3(0, 0) = 0;
		dg3(0, 1) = -DN_De(i, 0) * g2[2] + DN_De(i, 1)*g1[2];
		dg3(0, 2) = DN_De(i, 0) * g2[1] - DN_De(i, 1)*g1[1];

		//second line
		dg3(1, 0) = DN_De(i, 0) * g2[2] - DN_De(i, 1)*g1[2];
		dg3(1, 1) = 0;
		dg3(1, 2) = -DN_De(i, 0)*g2[0] + DN_De(i, 1)*g1[0];

		//third line
		dg3(2, 0) = - DN_De(i, 0) * g2[1] + DN_De(i, 1) * g1[1];
		dg3(2, 1) = DN_De(i, 0) * g2[0] - DN_De(i, 1) * g1[0];
		dg3(2, 2) = 0;

		//KRATOS_WATCH(dg3)

		for (unsigned int j = 0; j < 3; j++)
		{
			double g3dg3lg3 = (g3[0] * dg3(j, 0) + g3[1] * dg3(j, 1) + g3[2] * dg3(j, 2))*inddA3;

			dn(j, 0) = dg3(j, 0)*invdA - g3[0] * g3dg3lg3;
			dn(j, 1) = dg3(j, 1)*invdA - g3[1] * g3dg3lg3;
			dn(j, 2) = dg3(j, 2)*invdA - g3[2] * g3dg3lg3;
		}

		// curvature vector [K11,K22,K12] referred to curvilinear coordinate system
		b(0, index)		= 0 - (DDN_DDe(i, 0) * n[0] + H(0, 0)*dn(0, 0) + H(1, 0)*dn(0, 1) + H(2, 0)*dn(0, 2));
		b(0, index + 1) = 0 - (DDN_DDe(i, 0) * n[1] + H(0, 0)*dn(1, 0) + H(1, 0)*dn(1, 1) + H(2, 0)*dn(1, 2));
		b(0, index + 2) = 0 - (DDN_DDe(i, 0) * n[2] + H(0, 0)*dn(2, 0) + H(1, 0)*dn(2, 1) + H(2, 0)*dn(2, 2));

		//second line
		b(1, index)		= 0 - (DDN_DDe(i, 1) * n[0] + H(0, 1)*dn(0, 0) + H(1, 1)*dn(0, 1) + H(2, 1)*dn(0, 2));
		b(1, index + 1) = 0 - (DDN_DDe(i, 1) * n[1] + H(0, 1)*dn(1, 0) + H(1, 1)*dn(1, 1) + H(2, 1)*dn(1, 2));
		b(1, index + 2) = 0 - (DDN_DDe(i, 1) * n[2] + H(0, 1)*dn(2, 0) + H(1, 1)*dn(2, 1) + H(2, 1)*dn(2, 2));

		//third line
		b(2, index)		= 0 - (DDN_DDe(i, 2) * n[0] + H(0, 2)*dn(0, 0) + H(1, 2)*dn(0, 1) + H(2, 2)*dn(0, 2));
		b(2, index + 1) = 0 - (DDN_DDe(i, 2) * n[1] + H(0, 2)*dn(1, 0) + H(1, 2)*dn(1, 1) + H(2, 2)*dn(1, 2));
		b(2, index + 2) = 0 - (DDN_DDe(i, 2) * n[2] + H(0, 2)*dn(2, 0) + H(1, 2)*dn(2, 1) + H(2, 2)*dn(2, 2));
	}

	B = prod(Q, b);
	KRATOS_CATCH("")
}
//***********************************************************************************
//***********************************************************************************
void MeshlessShellElement::CalculateStrain(
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

void MeshlessShellElement::CalculateCurvature(
	Vector& CurvatureVector,
	array_1d<double, 3>& bv,
	array_1d<double, 3>& bv_ref)

{
	KRATOS_TRY

	CurvatureVector[0] = (bv[0] - bv_ref[0]);
	CurvatureVector[1] = (bv[1] - bv_ref[1]);
	CurvatureVector[2] = (bv[2] - bv_ref[2]);

	KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************
void MeshlessShellElement::CalculateMetricDeformed(const Matrix& DN_De,
  const Matrix& DDN_DDe,
	array_1d<double, 3>& gab,
	array_1d<double, 3>& curvature_coefficient,
	array_1d<double, 3>& g1,
	array_1d<double, 3>& g2)
{
	Matrix J;
	Jacobian(DN_De, J);

  //KRATOS_WATCH(J)

	//auxiliary terms
	array_1d<double, 3> g3;

	//double IntegrationWeight = GetGeometry().IntegrationPoints()[0].Weight();

	g1[0] = J(0, 0);
	g2[0] = J(0, 1);
	g1[1] = J(1, 0);
	g2[1] = J(1, 1);
	g1[2] = J(2, 0);
	g2[2] = J(2, 1);

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

	Matrix H = ZeroMatrix(3, 3);
	Hessian(H, DDN_DDe);

	curvature_coefficient[0] = H(0, 0)*n[0] + H(1, 0)*n[1] + H(2, 0)*n[2];
	curvature_coefficient[1] = H(0, 1)*n[0] + H(1, 1)*n[1] + H(2, 1)*n[2];
	curvature_coefficient[2] = H(0, 2)*n[0] + H(1, 2)*n[1] + H(2, 2)*n[2];

  ////contravariant metric gab_con and base vectors g_con
  //Vector gab_con = ZeroVector(3);
  //double invdetGab = 1.0 / (gab[0] * gab[1] - gab[2] * gab[2]);
  //gab_con[0] = invdetGab*gab[1];
  //gab_con[2] = -invdetGab*gab[2];
  //gab_con[1] = invdetGab*gab[0];

  //array_1d<double, 3> g_con_1 = g1*gab_con[0] + g2*gab_con[2];
  //array_1d<double, 3> g_con_2 = g1*gab_con[2] + g2*gab_con[1];

  ////local cartesian coordinates
  //double lg1 = norm_2(g1);
  //array_1d<double, 3> e1 = g1 / lg1;
  //double lg_con2 = norm_2(g_con_2);
  //array_1d<double, 3> e2 = g_con_2 / lg_con2;

  //Matrix T_G_E = ZeroMatrix(3, 3);
  ////Transformation matrix T from contravariant to local cartesian basis
  //double eG11 = inner_prod(e1, g1);
  //double eG12 = inner_prod(e1, g2);
  //double eG21 = inner_prod(e2, g1);
  //double eG22 = inner_prod(e2, g2);
  //T_G_E(0, 0) = eG11*eG11;
  //T_G_E(0, 1) = eG12*eG12;
  //T_G_E(0, 2) = 2.0*eG11*eG12;
  //T_G_E(1, 0) = eG21*eG21;
  //T_G_E(1, 1) = eG22*eG22;
  //T_G_E(1, 2) = 2.0*eG21*eG22;
  //T_G_E(2, 0) = eG11*eG21;
  //T_G_E(2, 1) = eG12*eG22;
  //T_G_E(2, 2) = eG11*eG22 + eG12*eG21;
}
//***********************************************************************************
//***********************************************************************************
void MeshlessShellElement::CalculateMetricDeformed(
  const Matrix& DN_De,
  const Matrix& DDN_DDe, 
  MetricVariables& metric)
{
  Matrix J;
  Jacobian(DN_De, J);

  
  metric.g1[0] = J(0, 0);
  metric.g2[0] = J(0, 1);
  metric.g1[1] = J(1, 0);
  metric.g2[1] = J(1, 1);
  metric.g1[2] = J(2, 0);
  metric.g2[2] = J(2, 1);

  //basis vector g3
  CrossProduct2(metric.g3, metric.g1, metric.g2);
  //differential area dA
  metric.dA = norm_2(metric.g3);
  //normal vector _n
  array_1d<double, 3> n = metric.g3 / metric.dA;

  
  //GetCovariantMetric
  metric.gab[0] = pow(metric.g1[0], 2) + pow(metric.g1[1], 2) + pow(metric.g1[2], 2);
  metric.gab[1] = pow(metric.g2[0], 2) + pow(metric.g2[1], 2) + pow(metric.g2[2], 2);
  metric.gab[2] = metric.g1[0] * metric.g2[0] + metric.g1[1] * metric.g2[1] + metric.g1[2] * metric.g2[2];

  
  Hessian(metric.H, DDN_DDe);

  
  metric.curvature[0] = metric.H(0, 0)*n[0] + metric.H(1, 0)*n[1] + metric.H(2, 0)*n[2];
  metric.curvature[1] = metric.H(0, 1)*n[0] + metric.H(1, 1)*n[1] + metric.H(2, 1)*n[2];
  metric.curvature[2] = metric.H(0, 2)*n[0] + metric.H(1, 2)*n[1] + metric.H(2, 2)*n[2];

  
  //contravariant metric gab_con and base vectors g_con
  //Vector gab_con = ZeroVector(3);
  double invdetGab = 1.0 / (metric.gab[0] * metric.gab[1] - metric.gab[2] * metric.gab[2]);
  metric.gab_con[0] = invdetGab*metric.gab[1];
  metric.gab_con[2] = -invdetGab*metric.gab[2];
  metric.gab_con[1] = invdetGab*metric.gab[0];

  
  array_1d<double, 3> g_con_1 = metric.g1*metric.gab_con[0] + metric.g2*metric.gab_con[2];
  array_1d<double, 3> g_con_2 = metric.g1*metric.gab_con[2] + metric.g2*metric.gab_con[1];

  
  //local cartesian coordinates
  double lg1 = norm_2(metric.g1);
  array_1d<double, 3> e1 = metric.g1 / lg1;
  double lg_con2 = norm_2(g_con_2);
  array_1d<double, 3> e2 = g_con_2 / lg_con2;

  //Matrix T_G_E = ZeroMatrix(3, 3);
  //Transformation matrix T from contravariant to local cartesian basis
  double eG11 = inner_prod(e1, metric.g1);
  double eG12 = inner_prod(e1, metric.g2);
  double eG21 = inner_prod(e2, metric.g1);
  double eG22 = inner_prod(e2, metric.g2);

  metric.T = ZeroMatrix(3, 3);
  metric.T(0, 0) = eG11*eG11;
  metric.T(0, 1) = eG12*eG12;
  metric.T(0, 2) = 2.0*eG11*eG12;
  metric.T(1, 0) = eG21*eG21;
  metric.T(1, 1) = eG22*eG22;
  metric.T(1, 2) = 2.0*eG21*eG22;
  metric.T(2, 0) = eG11*eG21;
  metric.T(2, 1) = eG12*eG22;
  metric.T(2, 2) = eG11*eG22 + eG12*eG21;
}
//***********************************************************************************
//***********************************************************************************
void MeshlessShellElement::CalculateSecondVariationStrainCurvature(
  const Matrix& DN_De,
  const Matrix& DDN_DDe,
	Matrix& Strain_in_Q_coordinates11,
	Matrix& Strain_in_Q_coordinates22,
	Matrix& Strain_in_Q_coordinates12,
	Matrix& Curvature_in_Q_coordinates11,
	Matrix& Curvature_in_Q_coordinates22,
	Matrix& Curvature_in_Q_coordinates12,
	boost::numeric::ublas::bounded_matrix<double, 3, 3>& Q,
	array_1d<double, 3>& g1,
	array_1d<double, 3>& g2)
{
	const unsigned int number_of_nodes = GetGeometry().size();
	Matrix dg3_n = ZeroMatrix(3, 3);
	Matrix dg3_m = ZeroMatrix(3, 3);

	Vector ddStrain_curvilinear = ZeroVector(3);
	Vector ddCurvature_curvilinear = ZeroVector(3);

	Matrix dn_n = ZeroMatrix(3, 3);
	Matrix dn_m = ZeroMatrix(3, 3);

	//basis vector g3
	array_1d<double, 3> g3;
	array_1d<double, 3> ddn;


	CrossProduct(g3, g1, g2);

	Matrix H = ZeroMatrix(3, 3);

	this->Hessian(H, DDN_DDe);

	//differential area dA
	double dA = norm_2(g3);

	//normal vector n
	array_1d<double, 3> n = g3 / dA;

	double invdA  = 1 / dA;
	double invdA3 = 1 / pow(dA, 3);
	double invdA5 = 1 / pow(dA, 5);

	for (unsigned int n = 0; n < number_of_nodes; n++)
	{
		//first line --- dg(1,0)*g2(2)-dg(2,0)*g2(1) + g1(1)*dg(2,1)-g1(2)*dg(1,1);
		dg3_n(0, 0) = 0;
		dg3_n(0, 1) = -DN_De(n, 0) * g2[2] + DN_De(n, 1)*g1[2];
		dg3_n(0, 2) = DN_De(n, 0) * g2[1] - DN_De(n, 1)*g1[1];

		//second line --- dg(2,0)*g2(0)-dg(0,0)*g2(2) + g1(2)*dg(0,1)-g1(0)*dg(2,1);
		dg3_n(1, 0) = DN_De(n, 0) * g2[2] - DN_De(n, 1)*g1[2];
		dg3_n(1, 1) = 0;
		dg3_n(1, 2) = -DN_De(n, 0)*g2[0] + DN_De(n, 1)*g1[0];

		//third line --- dg(0,0)*g2(1)-dg(1,0)*g2(0) + g1(0)*dg(1,1)-g1(1)*dg(0,1);
		dg3_n(2, 0) = -DN_De(n, 0) * g2[1] + DN_De(n, 1) * g1[1];
		dg3_n(2, 1) = DN_De(n, 0) * g2[0] - DN_De(n, 1) * g1[0];
		dg3_n(2, 2) = 0;

		//std::cout << "dg3_n: " << dg3_n << std::endl;

		for (unsigned int i = 0; i < 3; i++)
		{
			double g3dg3n = (g3[0] * dg3_n(i, 0) + g3[1] * dg3_n(i, 1) + g3[2] * dg3_n(i, 2));
			double g3dg3lg3n = g3dg3n*invdA3;

			dn_n(i, 0) = dg3_n(i, 0)*invdA - g3[0] * g3dg3lg3n;
			dn_n(i, 1) = dg3_n(i, 1)*invdA - g3[1] * g3dg3lg3n;
			dn_n(i, 2) = dg3_n(i, 2)*invdA - g3[2] * g3dg3lg3n;

			//std::cout << dn_n << std::endl;

			for (unsigned int m = 0; m <= n; m++)
			{
				//first line --- dg(1,0)*g2(2)-dg(2,0)*g2(1) + g1(1)*dg(2,1)-g1(2)*dg(1,1);
				dg3_m(0, 0) = 0;
				dg3_m(0, 1) = -DN_De(m, 0) * g2[2] + DN_De(m, 1)*g1[2];
				dg3_m(0, 2) = DN_De(m, 0) * g2[1] - DN_De(m, 1)*g1[1];

				//second line --- dg(2,0)*g2(0)-dg(0,0)*g2(2) + g1(2)*dg(0,1)-g1(0)*dg(2,1);
				dg3_m(1, 0) = DN_De(m, 0) * g2[2] - DN_De(m, 1)*g1[2];
				dg3_m(1, 1) = 0;
				dg3_m(1, 2) = -DN_De(m, 0)*g2[0] + DN_De(m, 1)*g1[0];

				//third line --- dg(0,0)*g2(1)-dg(1,0)*g2(0) + g1(0)*dg(1,1)-g1(1)*dg(0,1);
				dg3_m(2, 0) = -DN_De(m, 0) * g2[1] + DN_De(m, 1) * g1[1];
				dg3_m(2, 1) = DN_De(m, 0) * g2[0] - DN_De(m, 1) * g1[0];
				dg3_m(2, 2) = 0;


				//std::cout << "dg3_m: " << dg3_m << std::endl;

				int limit = i+1;
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
						 
						Strain_in_Q_coordinates11(3*n + i, 3*m + j) = Q(0, 0)*ddStrain_curvilinear[0] + Q(0, 1)*ddStrain_curvilinear[1] + Q(0, 2)*ddStrain_curvilinear[2];
						Strain_in_Q_coordinates22(3*n + i, 3*m + j) = Q(1, 0)*ddStrain_curvilinear[0] + Q(1, 1)*ddStrain_curvilinear[1] + Q(1, 2)*ddStrain_curvilinear[2];
						Strain_in_Q_coordinates12(3*n + i, 3*m + j) = Q(2, 0)*ddStrain_curvilinear[0] + Q(2, 1)*ddStrain_curvilinear[1] + Q(2, 2)*ddStrain_curvilinear[2];

					}
					// curvature
					array_1d<double, 3> ddg3;
					ddg3[0] = ddg3[1] = ddg3[2] = 0;
					double direction = 4 - i - j;
					double ddirection = i - j;
					if (ddirection == -1)  ddg3(direction - 1) = DN_De(n, 0)*DN_De(m, 1) - DN_De(n, 1)*DN_De(m, 0);
					else if (ddirection == 2) ddg3(direction - 1) = DN_De(n, 0)*DN_De(m, 1) - DN_De(n, 1)*DN_De(m, 0);
					else if (ddirection == 1) ddg3(direction - 1) = -DN_De(n, 0)*DN_De(m, 1) + DN_De(n, 1)*DN_De(m, 0);
					else if (ddirection == -2) ddg3(direction - 1) = -DN_De(n, 0)*DN_De(m, 1) + DN_De(n, 1)*DN_De(m, 0);

					double g3dg3m = (g3[0] * dg3_m(j, 0) + g3[1] * dg3_m(j, 1) + g3[2] * dg3_m(j, 2));
					double g3dg3lg3m = g3dg3m*invdA3;

					dn_m(j, 0) = dg3_m(j, 0)*invdA - g3[0] * g3dg3lg3m;
					dn_m(j, 1) = dg3_m(j, 1)*invdA - g3[1] * g3dg3lg3m;
					dn_m(j, 2) = dg3_m(j, 2)*invdA - g3[2] * g3dg3lg3m;


					double c = -(ddg3[0] * g3[0] + ddg3[1] * g3[1] + ddg3[2] * g3[2]
						+ dg3_n(i, 0)*dg3_m(j,0) + dg3_n(i,1)*dg3_m(j,1) + dg3_n(i,2)*dg3_m(j,2)
						)*invdA3;


					double d = 3.0*g3dg3n * g3dg3m * invdA5;

					ddn[0] = ddg3[0] * invdA3 - g3dg3lg3m * dg3_n(i,0) - g3dg3lg3n * dg3_m(j,0) + (c + d)*g3[0];
					ddn[1] = ddg3[1] * invdA3 - g3dg3lg3m * dg3_n(i,1) - g3dg3lg3n * dg3_m(j,1) + (c + d)*g3[1];
					ddn[2] = ddg3[2] * invdA3 - g3dg3lg3m * dg3_n(i,2) - g3dg3lg3n * dg3_m(j,2) + (c + d)*g3[2];

					ddCurvature_curvilinear[0] = DDN_DDe(n, 0)*dn_m(j,i) + DDN_DDe(m, 0)*dn_n(i,j)
						+ H(0, 0)*ddn[0] + H(1, 0)*ddn[1] + H(2, 0)*ddn[2];
					ddCurvature_curvilinear[1] = DDN_DDe(n, 1)*dn_m(j,i) + DDN_DDe(m, 1)*dn_n(i,j)
						+ H(0, 1)*ddn[0] + H(1, 1)*ddn[1] + H(2, 1)*ddn[2];
					ddCurvature_curvilinear[2] = DDN_DDe(n, 2)*dn_m(j,i) + DDN_DDe(m, 2)*dn_n(i,j)
						+ H(0, 2)*ddn[0] + H(1, 2)*ddn[1] + H(2, 2)*ddn[2];

					Curvature_in_Q_coordinates11(3*n + i, 3*m + j) = Q(0, 0)*ddCurvature_curvilinear[0] + Q(0, 1)*ddCurvature_curvilinear[1] + Q(0, 2)*ddCurvature_curvilinear[2];
					Curvature_in_Q_coordinates22(3*n + i, 3*m + j) = Q(1, 0)*ddCurvature_curvilinear[0] + Q(1, 1)*ddCurvature_curvilinear[1] + Q(1, 2)*ddCurvature_curvilinear[2];
					Curvature_in_Q_coordinates12(3*n + i, 3*m + j) = Q(2, 0)*ddCurvature_curvilinear[0] + Q(2, 1)*ddCurvature_curvilinear[1] + Q(2, 2)*ddCurvature_curvilinear[2];
				}
			}
		}
	}
}
void MeshlessShellElement::CalculateMassMatrix(
	MatrixType& rMassMatrix,
	ProcessInfo& rCurrentProcessInfo
)
{
	KRATOS_TRY;

	//double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
	//Vector ShapeFunctionsN = this->GetValue(SHAPE_FUNCTION_VALUES);
	//Matrix DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

	//double density = this->GetProperties().GetValue(DENSITY);

	//KRATOS_WATCH(DN_De)

	//Vector g1, g2, g3;

	//Matrix J;
	//Jacobian(DN_De, J);
	//KRATOS_WATCH(J)

	//g1[0] = J(0, 0);
	//g2[0] = J(0, 1);
	//g1[1] = J(1, 0);
	//g2[1] = J(1, 1);
	//g1[2] = J(2, 0);
	//g2[2] = J(2, 1);

	//CrossProduct2(g3, g1, g2);
	////differential area dA
	//double dA = norm_2(g3);
	//KRATOS_WATCH(dA)

	//unsigned int dimension = 3;
	//unsigned int number_of_nodes = ShapeFunctionsN.size();
	//unsigned int mat_size = dimension * number_of_nodes;
	//
	//if (rMassMatrix.size1() != mat_size)
	//{
	//	rMassMatrix.resize(mat_size, mat_size, false);
	//}
	//rMassMatrix = ZeroMatrix(mat_size, mat_size);
	//KRATOS_WATCH(rMassMatrix)
	//for (int r = 0; r<number_of_nodes; r++)
	//{
	//	for (int s = 0; s<number_of_nodes; s++)
	//	{
	//		rMassMatrix(3 * s, 3 * r) = ShapeFunctionsN(s)*ShapeFunctionsN(r) * density * dA * integration_weight;
	//		rMassMatrix(3 * s + 1, 3 * r + 1) = rMassMatrix(3 * s, 3 * r);
	//		rMassMatrix(3 * s + 2, 3 * r + 2) = rMassMatrix(3 * s, 3 * r);
	//	}
	//}
	//KRATOS_WATCH(rMassMatrix)

	KRATOS_CATCH("")
}
//***********************************************************************************
//***********************************************************************************
void MeshlessShellElement::CalculateAll(
	MatrixType& rLeftHandSideMatrix,
	VectorType& rRightHandSideVector,
	const ProcessInfo& rCurrentProcessInfo,
	bool CalculateStiffnessMatrixFlag,
	bool CalculateResidualVectorFlag)

{
	KRATOS_TRY
	std::cout << "Calculate All Shell Element!" << std::endl;
	// definition of problem size
	const unsigned int number_of_nodes = GetGeometry().size();
	unsigned int MatSize = number_of_nodes * 3;

	//set up Constitutive Law
	ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

	Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
	Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
	Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

	//resizing as needed the LHS
	if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
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

	//reading in of integration weight, shape function values and shape function derivatives
	double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
	Vector ShapeFunctionsN = this->GetValue(SHAPE_FUNCTION_VALUES);
	Matrix DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
	Matrix DDN_DDe = this->GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

	// covariant metric in deformed system
	array_1d<double, 3> gab;
	// curvature covariant metric in deformed system
	array_1d<double, 3> curvature;

	// Transformation Matrix Q - 
	boost::numeric::ublas::bounded_matrix<double, 3, 3>  Q = mQ;

  Vector StrainVector = ZeroVector(3);
  Vector CurvatureVector = ZeroVector(3);

  // basis vectors in deformed system
  array_1d<double, 3> g1, g2;
  CalculateMetricDeformed(DN_De, DDN_DDe, gab, curvature, g1, g2);
  CalculateStrain(StrainVector, gab, mGab0);

  Vector StrainVector_in_Q_coordinates = prod(Q, StrainVector);
  CalculateCurvature(CurvatureVector, curvature, mCurvature0);
  Vector CurvatureVector_in_Q_coordinates = prod(Q, CurvatureVector);

  //Constitive Matrices DMembrane and DCurvature
  Matrix DMembrane = ZeroMatrix(3, 3);
  Matrix DCurvature = ZeroMatrix(3, 3);

	Values.SetStrainVector(StrainVector_in_Q_coordinates); //this is the input parameter
	Vector StressVector;
	Values.SetStressVector(StressVector);    //this is an ouput parameter
	Values.SetConstitutiveMatrix(DMembrane); //this is an ouput parameter
	mConstitutiveLawVector->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);
	double thickness = this->GetProperties().GetValue(THICKNESS);
	DCurvature = DMembrane*(pow(thickness, 2) / 12);

	//Local Cartesian Foreces and Moments
	//Vector ForceVector_in_Q_coordinates = prod(trans(DMembrane), StrainVector_in_Q_coordinates);
	//Vector MomentVector_in_Q_coordinates = prod(trans(DCurvature), CurvatureVector_in_Q_coordinates);
	Vector ForceVector_in_Q_coordinates = ZeroVector(3);
	Vector MomentVector_in_Q_coordinates = ZeroVector(3);

	ForceVector_in_Q_coordinates = prod(trans(DMembrane), StrainVector_in_Q_coordinates);
	MomentVector_in_Q_coordinates = prod(trans(DCurvature), CurvatureVector_in_Q_coordinates);

	//double damage_t = 0.0;
	//mConstitutiveLawVector->GetValue(DAMAGE_T, damage_t);
	////KRATOS_WATCH(damage_t)
	//double damage_c = 0.0;
	//mConstitutiveLawVector->GetValue(DAMAGE_C, damage_c);
	//KRATOS_WATCH(damage_c)

	KRATOS_WATCH(Q)
	// calculate B MATRICES
	//B matrices:
	Matrix BMembrane = ZeroMatrix(3, MatSize);
	Matrix BCurvature = ZeroMatrix(3, MatSize);
	CalculateBMembrane(BMembrane, Q, DN_De, g1, g2);
    KRATOS_WATCH(BMembrane)
	CalculateBCurvature(BCurvature, Q, DN_De, DDN_DDe, g1, g2);
    KRATOS_WATCH(BCurvature)

	// integration on the REFERENCE CONFIGURATION
	double DetJ0 = mDetJ0;
	double IntToReferenceWeight = integration_weight * DetJ0 * mThickness0;

	// Nonlinear Deformation
	Matrix Strain_in_Q_coordinates11    = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);
	Matrix Strain_in_Q_coordinates22    = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);
	Matrix Strain_in_Q_coordinates12    = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);
	Matrix Curvature_in_Q_coordinates11 = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);
	Matrix Curvature_in_Q_coordinates22 = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);
	Matrix Curvature_in_Q_coordinates12 = ZeroMatrix(number_of_nodes * 3, number_of_nodes * 3);
	CalculateSecondVariationStrainCurvature(DN_De, DDN_DDe,
		Strain_in_Q_coordinates11, Strain_in_Q_coordinates22, Strain_in_Q_coordinates12,
		Curvature_in_Q_coordinates11, Curvature_in_Q_coordinates22, Curvature_in_Q_coordinates12, Q, g1, g2);

	// LEFT HAND SIDE MATRIX
	if (CalculateStiffnessMatrixFlag == true)
	{
		//adding membrane contributions to the stiffness matrix
		CalculateAndAddKm(rLeftHandSideMatrix, BMembrane, DMembrane, IntToReferenceWeight);
		//adding curvature contributions to the stiffness matrix
		CalculateAndAddKm(rLeftHandSideMatrix, BCurvature, DCurvature, IntToReferenceWeight);

		// adding  non-linear-contribution to Stiffness-Matrix
		CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
			Strain_in_Q_coordinates11, Strain_in_Q_coordinates22, Strain_in_Q_coordinates12,
			ForceVector_in_Q_coordinates,
			IntToReferenceWeight);

		CalculateAndAddNonlinearKm(rLeftHandSideMatrix,
			Curvature_in_Q_coordinates11, Curvature_in_Q_coordinates22, Curvature_in_Q_coordinates12,
			MomentVector_in_Q_coordinates,
			IntToReferenceWeight);
	}

	//if(this->Id() == 30) //TODO: remove this! it is just for debugging purposes
	//{
	//	KRATOS_WATCH(StrainVector)
	//	KRATOS_WATCH(StressVector)
	//}

	// RIGHT HAND SIDE VECTOR
	if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
	{
		// operation performed: rRightHandSideVector -= Weight*IntForce
		noalias(rRightHandSideVector) -= IntToReferenceWeight * prod(trans(BMembrane), ForceVector_in_Q_coordinates);
		noalias(rRightHandSideVector) += IntToReferenceWeight * prod(trans(BCurvature), MomentVector_in_Q_coordinates);
	}

  //if (this->Id() == 1) //TODO: remove this! it is just for debugging purposes
  //{

    //KRATOS_WATCH(ShapeFunctionsN)
    //KRATOS_WATCH(DN_De)
    //KRATOS_WATCH(DDN_DDe)

    //KRATOS_WATCH(IntToReferenceWeight)

    //KRATOS_WATCH(BMembrane)
    //KRATOS_WATCH(BCurvature)

    //KRATOS_WATCH(rLeftHandSideMatrix)
    //KRATOS_WATCH(rRightHandSideVector)
  //}
  //		Vector displacements;
	//	this->GetValuesVector(displacements, 0);
	//	KRATOS_WATCH(displacements);

	//	KRATOS_WATCH("coordintates")
	//		for (unsigned int i = 0; i < GetGeometry().size(); i++)
	//			std::cout << " " << GetGeometry()[i].Id() << " " << GetGeometry()[i].Coordinates() << std::endl;

	//	KRATOS_WATCH("initial coordintates")
	//		for (unsigned int i = 0; i < GetGeometry().size(); i++)
	//			std::cout << " " << GetGeometry()[i].Id() << " " << GetGeometry()[i].GetInitialPosition() << std::endl;
	//}

	KRATOS_CATCH("")
}


////************************************************************************************
////************************************************************************************
//void MeshlessShellElement::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
//	std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
//{
//	//std::cout << "GetValueOnIntegrationPoints" << std::endl;
//	//if (rVariable == GREEN_LAGRANGE_STRAIN_TENSOR)
//	//{
//	//	double damage_t;
//	//	mConstitutiveLawVector->GetValue(DAMAGE_T, damage_t);
//
//	//}
//
//	//if (rVariable == PK2_STRESS_TENSOR)
//	//{
//	//	CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
//	//}
//	//// VM
//	//if (rVariable == CAUCHY_STRESS_TENSOR)
//	//{
//	//	CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
//	//}
//	// VM
//}
void MeshlessShellElement::GetValueOnIntegrationPoints(
	const Variable<double>& rVariable,
	std::vector<double>& rValues,
	const ProcessInfo& rCurrentProcessInfo
)
{
	if (rValues.size() != 1)
	{
		rValues.resize(1);
	}

	if (rVariable == DAMAGE_T)
	{
		mConstitutiveLawVector->GetValue(DAMAGE_T, rValues[0]);
	}
	else if (rVariable == DAMAGE_C)
	{
		mConstitutiveLawVector->GetValue(DAMAGE_C, rValues[0]);
	}
	else
	{
		CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
	}
}

int  MeshlessShellElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

	//verify that the variables are correctly initialized

	if (VELOCITY.Key() == 0)
		KRATOS_THROW_ERROR(std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "")

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


	mConstitutiveLawVector->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);

	ConstitutiveLaw::Features LawFeatures;
	mConstitutiveLawVector->GetLawFeatures(LawFeatures);

	if (LawFeatures.mOptions.IsNot(ConstitutiveLaw::PLANE_STRESS_LAW))
		KRATOS_THROW_ERROR(std::logic_error, "Constitutive law is compatible only with a plane stress 2D law for membrane element with Id", this->Id())

	if (LawFeatures.mOptions.IsNot(ConstitutiveLaw::INFINITESIMAL_STRAINS))
		KRATOS_THROW_ERROR(std::logic_error, "Constitutive law is compatible only with a law using infinitessimal strains for membrane element with Id", this->Id())

	if (LawFeatures.mStrainSize != 3) KRATOS_THROW_ERROR(std::logic_error, "Constitutive law expects a strain size different from 3 for membrane element with Id", this->Id())


	return 0;

	KRATOS_CATCH("");
}

} // Namespace Kratos


