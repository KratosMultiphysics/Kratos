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
	//std::cout << "start function Initialize in Meshless Shell Element" << std::endl;
	KRATOS_TRY
		//KRATOS_WATCH( "INITIALIZE ELEMENT" )
		//getting all "Actual" info from the geometry
		//const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

	// Get values of shape functions and derivatives. Derivatives are first and sevon column: first dreivatives. Third, fourth and fith are second derivatives
	double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
	Vector ShapeFunctionsN = this->GetValue(SHAPE_FUNCTION_VALUES);
	Matrix DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);



	// Initialize Variables
	mdensity = GetProperties()[DENSITY];
	mThickness0 = GetProperties()[THICKNESS];
	//WIRD DAS GEBRAUCHT??
	mThickness = 0.00;

	// Initialize Material
	mConstitutiveLawVector = GetProperties()[CONSTITUTIVE_LAW]->Clone();
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
	Hessian(H, DN_De);

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
	mTotalDomainInitialSize = mDetJ0 * integration_weight; //RAUS SCHMEI?EN

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

//***********************************************************************************
//***********************************************************************************

void MeshlessShellElement::CalculateOnIntegrationPoints(
	const Variable<Matrix>& rVariable,
	std::vector<Matrix>& Output,
	const ProcessInfo& rCurrentProcessInfo)

{
	std::cout << "CalculateOnIntegrationPoints not yet implemented" << std::endl;
}


//***********************************************************************************
//***********************************************************************************

void MeshlessShellElement::FinalizeSolutionStep(
	ProcessInfo& rCurrentProcessInfo)
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
					K(3*n + i, 3*m + j) += (SD[0] * B11(3*n + i, 3*m + j) + SD[1] * B22(3 * n + i, 3 * m + j) + SD[2] * B12(3 * n + i, 3 * m + j))*weight;
					K(3*m + j, 3*n + i) += (SD[0] * B11(3 * n + i, 3 * m + j) + SD[1] * B22(3 * n + i, 3 * m + j) + SD[2] * B12(3 * n + i, 3 * m + j))*weight;
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
	Hessian(H, DN_De);

	//basis vector g3
	array_1d<double, 3> g3;
	array_1d<double, 3> n;
	CrossProduct(g3, g1, g2);

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
		b(0, index)		= 0 - (DN_De(i, 2) * n[0] + H(0, 0)*dn(0, 0) + H(1, 0)*dn(0, 1) + H(2, 0)*dn(0, 2));
		b(0, index + 1) = 0 - (DN_De(i, 2) * n[1] + H(0, 0)*dn(1, 0) + H(1, 0)*dn(1, 1) + H(2, 0)*dn(1, 2));
		b(0, index + 2) = 0 - (DN_De(i, 2) * n[2] + H(0, 0)*dn(2, 0) + H(1, 0)*dn(2, 1) + H(2, 0)*dn(2, 2));

		//second line
		b(1, index)		= 0 - (DN_De(i, 3) * n[0] + H(0, 1)*dn(0, 0) + H(1, 1)*dn(0, 1) + H(2, 1)*dn(0, 2));
		b(1, index + 1) = 0 - (DN_De(i, 3) * n[1] + H(0, 1)*dn(1, 0) + H(1, 1)*dn(1, 1) + H(2, 1)*dn(1, 2));
		b(1, index + 2) = 0 - (DN_De(i, 3) * n[2] + H(0, 1)*dn(2, 0) + H(1, 1)*dn(2, 1) + H(2, 1)*dn(2, 2));

		//third line
		b(2, index)		= 0 - (DN_De(i, 4) * n[0] + H(0, 2)*dn(0, 0) + H(1, 2)*dn(0, 1) + H(2, 2)*dn(0, 2));
		b(2, index + 1) = 0 - (DN_De(i, 4) * n[1] + H(0, 2)*dn(1, 0) + H(1, 2)*dn(1, 1) + H(2, 2)*dn(1, 2));
		b(2, index + 2) = 0 - (DN_De(i, 4) * n[2] + H(0, 2)*dn(2, 0) + H(1, 2)*dn(2, 1) + H(2, 2)*dn(2, 2));
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
	array_1d<double, 3>& gab,
	array_1d<double, 3>& curvature_coefficient,
	array_1d<double, 3>& g1,
	array_1d<double, 3>& g2)
{
	Matrix J0;
	Jacobian(DN_De, J0);

	//auxiliary terms
	array_1d<double, 3> g3;

	//double IntegrationWeight = GetGeometry().IntegrationPoints()[0].Weight();

	g1[0] = J0(0, 0);
	g2[0] = J0(0, 1);
	g1[1] = J0(1, 0);
	g2[1] = J0(1, 1);
	g1[2] = J0(2, 0);
	g2[2] = J0(2, 1);

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
	Hessian(H, DN_De);

	curvature_coefficient[0] = H(0, 0)*n[0] + H(1, 0)*n[1] + H(2, 0)*n[2];
	curvature_coefficient[1] = H(0, 1)*n[0] + H(1, 1)*n[1] + H(2, 1)*n[2];
	curvature_coefficient[2] = H(0, 2)*n[0] + H(1, 2)*n[1] + H(2, 2)*n[2];
}

//***********************************************************************************
//***********************************************************************************
void MeshlessShellElement::CalculateSecondVariationStrainCurvature(Matrix DN_De,
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

	Matrix Hessian = ZeroMatrix(3, 3);

	this->Hessian(Hessian, DN_De);

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

					ddCurvature_curvilinear[0] = DN_De(n, 2)*dn_m(j,i) + DN_De(m, 2)*dn_n(i,j)
						+ Hessian(0, 0)*ddn[0] + Hessian(1, 0)*ddn[1] + Hessian(2, 0)*ddn[2];
					ddCurvature_curvilinear[1] = DN_De(n, 3)*dn_m(j,i) + DN_De(m, 3)*dn_n(i,j)
						+ Hessian(0, 1)*ddn[0] + Hessian(1, 1)*ddn[1] + Hessian(2, 1)*ddn[2];
					ddCurvature_curvilinear[2] = DN_De(n, 4)*dn_m(j,i) + DN_De(m, 4)*dn_n(i,j)
						+ Hessian(0, 2)*ddn[0] + Hessian(1, 2)*ddn[1] + Hessian(2, 2)*ddn[2];

					Curvature_in_Q_coordinates11(3*n + i, 3*m + j) = Q(0, 0)*ddCurvature_curvilinear[0] + Q(0, 1)*ddCurvature_curvilinear[1] + Q(0, 2)*ddCurvature_curvilinear[2];
					Curvature_in_Q_coordinates22(3*n + i, 3*m + j) = Q(1, 0)*ddCurvature_curvilinear[0] + Q(1, 1)*ddCurvature_curvilinear[1] + Q(1, 2)*ddCurvature_curvilinear[2];
					Curvature_in_Q_coordinates12(3*n + i, 3*m + j) = Q(2, 0)*ddCurvature_curvilinear[0] + Q(2, 1)*ddCurvature_curvilinear[1] + Q(2, 2)*ddCurvature_curvilinear[2];
				}
			}
		}
	}
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
	// definition of problem size
	const unsigned int number_of_nodes = GetGeometry().size();
	unsigned int MatSize = number_of_nodes * 3;

	//set up Constitutive Law
	ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

	Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRAIN, false);
	Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
	Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

	//Matrix B(3, MatSize);
	Vector StrainVector = ZeroVector(3);
	Vector StressVector = ZeroVector(3);

	Vector CurvatureVector = ZeroVector(3);

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

	// covariant metric in deformed system
	array_1d<double, 3> gab;
	// curvature covariant metric in deformed system
	array_1d<double, 3> curvature;

	// Transformation Matrix Q
	boost::numeric::ublas::bounded_matrix<double, 3, 3>  Q = ZeroMatrix(3, 3);
	Q = mQ;

	// basis vectors in deformed system
	array_1d<double, 3> g1;
	array_1d<double, 3> g2;
	CalculateMetricDeformed(DN_De, gab, curvature, g1, g2);
	CalculateStrain(StrainVector, gab, mGab0);

	Vector StrainVector_in_Q_coordinates = prod(Q, StrainVector);

	CalculateCurvature(CurvatureVector, curvature, mCurvature0);
	Vector CurvatureVector_in_Q_coordinates = prod(Q, CurvatureVector);

	//Constitive Matrices DMembrane and DCurvature
	Matrix DMembrane = ZeroMatrix(3, 3);
	Matrix DCurvature = ZeroMatrix(3, 3);

	Values.SetStrainVector(StrainVector); //this is the input parameter
	Values.SetStressVector(StressVector); //this is an ouput parameter
	Values.SetConstitutiveMatrix(DMembrane); //this is an ouput parameter

	mConstitutiveLawVector->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

	double thickness = this->GetProperties().GetValue(THICKNESS);
	DCurvature = DMembrane*(pow(thickness, 2) / 12);

	//Local Cartesian Foreces and Moments
	Vector ForceVector_in_Q_coordinates = ZeroVector(3);
	Vector MomentVector_in_Q_coordinates = ZeroVector(3);

	ForceVector_in_Q_coordinates = prod(trans(DMembrane), StrainVector_in_Q_coordinates);
	MomentVector_in_Q_coordinates = prod( trans(DCurvature), CurvatureVector_in_Q_coordinates);


	// calculate B MATRICES
	//B matrices:
	Matrix BMembrane = ZeroMatrix(3, MatSize);
	Matrix BCurvature = ZeroMatrix(3, MatSize);
	CalculateBMembrane(BMembrane, Q, DN_De, g1, g2);
	CalculateBCurvature(BCurvature, Q, DN_De, g1, g2);
	
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
	CalculateSecondVariationStrainCurvature(DN_De,
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

	//if (this->Id() == 30) //TODO: remove this! it is just for debugging purposes
	//{
	//	KRATOS_WATCH(rLeftHandSideMatrix)
	//		KRATOS_WATCH(rRightHandSideVector)
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


//************************************************************************************
//************************************************************************************
void MeshlessShellElement::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
	std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
	std::cout << "GetValueOnIntegrationPoints" << std::endl;
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

		//check constitutive law
		//for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
		//{
			mConstitutiveLawVector->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);

			ConstitutiveLaw::Features LawFeatures;
			mConstitutiveLawVector->GetLawFeatures(LawFeatures);

			if (LawFeatures.mOptions.IsNot(ConstitutiveLaw::PLANE_STRESS_LAW))
				KRATOS_THROW_ERROR(std::logic_error, "Constitutive law is compatible only with a plane stress 2D law for membrane element with Id", this->Id())

				if (LawFeatures.mOptions.IsNot(ConstitutiveLaw::INFINITESIMAL_STRAINS))
					KRATOS_THROW_ERROR(std::logic_error, "Constitutive law is compatible only with a law using infinitessimal strains for membrane element with Id", this->Id())

					if (LawFeatures.mStrainSize != 3) KRATOS_THROW_ERROR(std::logic_error, "Constitutive law expects a strain size different from 3 for membrane element with Id", this->Id())
		//}

	return 0;

	KRATOS_CATCH("");
}

} // Namespace Kratos


