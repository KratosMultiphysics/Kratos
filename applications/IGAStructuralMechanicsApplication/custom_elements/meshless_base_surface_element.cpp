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
#include "custom_elements/meshless_base_surface_element.h"

#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"

#include "utilities/math_utils.h"

#include "geometries/geometry.h"

namespace Kratos
{
	void MeshlessBaseSurfaceElement::Initialize()
	{
		KRATOS_TRY

		//Constitutive Law initialisation
		InitializeMaterial();

		MetricVariables initial_metric(3);
		CalculateMetric(initial_metric);
		mInitialMetric = initial_metric;

		KRATOS_CATCH("")
	}
	//***********************************************************************************
	//***********************************************************************************
	void MeshlessBaseSurfaceElement::FinalizeSolutionStep(
		ProcessInfo& rCurrentProcessInfo)
	{
		mConstitutiveLaw->FinalizeSolutionStep(GetProperties(),
			GetGeometry(),
			this->GetValue(SHAPE_FUNCTION_VALUES),
			rCurrentProcessInfo);
	}
	//************************************************************************************
	//************************************************************************************
	void MeshlessBaseSurfaceElement::InitializeMaterial()
	{
		KRATOS_TRY

		if (GetProperties()[CONSTITUTIVE_LAW] != NULL)
		{
			//get shape functions for evaluation of stresses inside the constitutive law
			const Vector&  N = this->GetValue(SHAPE_FUNCTION_VALUES);
			const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
			const double integration_weight = this->GetValue(INTEGRATION_WEIGHT);

			mConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
			ProcessInfo emptyProcessInfo = ProcessInfo();
			//mConstitutiveLaw->SetValue(INTEGRATION_WEIGHT, integration_weight, emptyProcessInfo);
			//mConstitutiveLaw->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, DN_De, emptyProcessInfo);
			mConstitutiveLaw->InitializeMaterial(GetProperties(), GetGeometry(), N);
		}
		else
		{
			KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
		}

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void MeshlessBaseSurfaceElement::CalculateAll(
		MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo,
		const bool CalculateStiffnessMatrixFlag,
		const bool CalculateResidualVectorFlag
	)
	{
		KRATOS_ERROR << "You have called to the CalculateAll from the base class for meshless surface elements" << std::endl;
	}

	//************************************************************************************
	//************************************************************************************
	void MeshlessBaseSurfaceElement::CalculateRightHandSide(
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo
	)
	{
		// Calculation flags
		const bool CalculateStiffnessMatrixFlag = false;
		const bool CalculateResidualVectorFlag = true;
		MatrixType temp = Matrix();

		CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
	}

	//************************************************************************************
	//************************************************************************************
	void MeshlessBaseSurfaceElement::CalculateLocalSystem(
		MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo
	)
	{
		//calculation flags
		const bool CalculateStiffnessMatrixFlag = true;
		const bool CalculateResidualVectorFlag = true;

		CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
	}

	//************************************************************************************
	//************************************************************************************
	void MeshlessBaseSurfaceElement::CalculateAndAddKm(
		MatrixType& rLeftHandSideMatrix,
		const Matrix& B,
		const Matrix& D,
		const double IntegrationWeight
	)
	{
		KRATOS_TRY
		noalias(rLeftHandSideMatrix) += IntegrationWeight * prod(trans(B), Matrix(prod(D, B)));
		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************
	void MeshlessBaseSurfaceElement::CalculateAndAddNonlinearKm(
		Matrix& rLeftHandSideMatrix,
		const Matrix& B11,
		const Matrix& B22,
		const Matrix& B12,
		const Vector& SD,
		double IntegrationWeight)

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

						rLeftHandSideMatrix(nbase, mbase) += (SD[0] * B11(nbase, mbase)
							+ SD[1] * B22(nbase, mbase)
							+ SD[2] * B12(nbase, mbase))*IntegrationWeight;
						rLeftHandSideMatrix(mbase, nbase) += (SD[0] * B11(nbase, mbase)
							+ SD[1] * B22(nbase, mbase)
							+ SD[2] * B12(nbase, mbase))*IntegrationWeight;
					}
				}
			}
		}
		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void MeshlessBaseSurfaceElement::CalculateMetric(
		MetricVariables& metric)
	{
		KRATOS_ERROR << "You have called to the CalculateMetric from the base class for meshless surface elements" << std::endl;
	}

	//************************************************************************************
	//************************************************************************************
	void MeshlessBaseSurfaceElement::CalculateConstitutiveVariables(
		MetricVariables& rActualMetric,
		ConstitutiveVariables& rThisConstitutiveVariables,
		ConstitutiveLaw::Parameters& rValues,
		const ConstitutiveLaw::StressMeasure ThisStressMeasure)
	{
		KRATOS_ERROR << "You have called to the CalculateMetric from the base class for meshless surface elements" << std::endl;
	}

	//***********************************************************************************
	//***********************************************************************************
	void MeshlessBaseSurfaceElement::CalculateStrain(
		Vector& StrainVector,
		Vector& gab,
		Vector& gab0)

	{
		KRATOS_TRY

		StrainVector[0] = 0.5 * (gab[0] - gab0[0]);
		StrainVector[1] = 0.5 * (gab[1] - gab0[1]);
		StrainVector[2] = 0.5 * (gab[2] - gab0[2]);

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************
	void MeshlessBaseSurfaceElement::CalculateCurvature(
		Vector& CurvatureVector,
		Vector& bv,
		Vector& bv_ref)

	{
		KRATOS_TRY

		CurvatureVector[0] = (bv[0] - bv_ref[0]);
		CurvatureVector[1] = (bv[1] - bv_ref[1]);
		CurvatureVector[2] = (bv[2] - bv_ref[2]);

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************
	void MeshlessBaseSurfaceElement::CalculateBMembrane(
		Matrix& rB,
		const MetricVariables& metric)
	{
	KRATOS_TRY

		if (this->Has(SHAPE_FUNCTION_LOCAL_DERIVATIVES))
		{
			const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

			const unsigned int number_of_nodes = GetGeometry().size();
			Matrix b(3, number_of_nodes * 3);

			for (unsigned int i = 0; i < number_of_nodes; i++)
			{
				unsigned int index = 3 * i;

				//first line
				b(0, index) = DN_De(i, 0) * metric.g1[0];
				b(0, index + 1) = DN_De(i, 0) * metric.g1[1];
				b(0, index + 2) = DN_De(i, 0) * metric.g1[2];

				//second line
				b(1, index) = DN_De(i, 1) * metric.g2[0];
				b(1, index + 1) = DN_De(i, 1) * metric.g2[1];
				b(1, index + 2) = DN_De(i, 1) * metric.g2[2];

				//third line
				b(2, index) = 0.5*(DN_De(i, 1) * metric.g1[0] + DN_De(i, 0) * metric.g2[0]);
				b(2, index + 1) = 0.5*(DN_De(i, 1) * metric.g1[1] + DN_De(i, 0) * metric.g2[1]);
				b(2, index + 2) = 0.5*(DN_De(i, 1) * metric.g1[2] + DN_De(i, 0) * metric.g2[2]);
			}

			rB = prod(metric.Q, b);
		}
		else
		{
			KRATOS_ERROR << "Element does not provide SHAPE_FUNCTION_LOCAL_DERIVATIVES" << std::endl;
		}
		KRATOS_CATCH("")
	}
	//***********************************************************************************
	//***********************************************************************************
	void MeshlessBaseSurfaceElement::CalculateBCurvature(
		Matrix& rB,
		const MetricVariables& metric)
	{
		KRATOS_TRY

		if (this->Has(SHAPE_FUNCTION_LOCAL_DERIVATIVES) && this->Has(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES))
		{
			const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
			const Matrix& DDN_DDe = this->GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);
			const unsigned int number_of_nodes = GetGeometry().size();

			Matrix dg3 = ZeroMatrix(3, 3);
			Matrix dn = ZeroMatrix(3, 3);
			Matrix b = ZeroMatrix(3, number_of_nodes * 3);

			//normal vector n
			Vector n = metric.g3 / metric.dA;

			double invdA = 1 / metric.dA;
			double inddA3 = 1 / pow(metric.dA, 3);


			for (unsigned int i = 0; i < number_of_nodes; i++)
			{
				unsigned int index = 3 * i;
				//first line
				dg3(0, 0) = 0;
				dg3(0, 1) = -DN_De(i, 0) * metric.g2[2] + DN_De(i, 1)*metric.g1[2];
				dg3(0, 2) = DN_De(i, 0) * metric.g2[1] - DN_De(i, 1)*metric.g1[1];

				//second line
				dg3(1, 0) = DN_De(i, 0) * metric.g2[2] - DN_De(i, 1)*metric.g1[2];
				dg3(1, 1) = 0;
				dg3(1, 2) = -DN_De(i, 0)*metric.g2[0] + DN_De(i, 1)*metric.g1[0];

				//third line
				dg3(2, 0) = -DN_De(i, 0) * metric.g2[1] + DN_De(i, 1) * metric.g1[1];
				dg3(2, 1) = DN_De(i, 0) * metric.g2[0] - DN_De(i, 1) * metric.g1[0];
				dg3(2, 2) = 0;

				//KRATOS_WATCH(dg3)

				for (unsigned int j = 0; j < 3; j++)
				{
					double g3dg3lg3 = (metric.g3[0] * dg3(j, 0) + metric.g3[1] * dg3(j, 1) + metric.g3[2] * dg3(j, 2))*inddA3;

					dn(j, 0) = dg3(j, 0)*invdA - metric.g3[0] * g3dg3lg3;
					dn(j, 1) = dg3(j, 1)*invdA - metric.g3[1] * g3dg3lg3;
					dn(j, 2) = dg3(j, 2)*invdA - metric.g3[2] * g3dg3lg3;
				}

				// curvature vector [K11,K22,K12] referred to curvilinear coordinate system
				b(0, index) = 0 - (DDN_DDe(i, 0) * n[0] + metric.H(0, 0)*dn(0, 0) + metric.H(1, 0)*dn(0, 1) + metric.H(2, 0)*dn(0, 2));
				b(0, index + 1) = 0 - (DDN_DDe(i, 0) * n[1] + metric.H(0, 0)*dn(1, 0) + metric.H(1, 0)*dn(1, 1) + metric.H(2, 0)*dn(1, 2));
				b(0, index + 2) = 0 - (DDN_DDe(i, 0) * n[2] + metric.H(0, 0)*dn(2, 0) + metric.H(1, 0)*dn(2, 1) + metric.H(2, 0)*dn(2, 2));

				//second line
				b(1, index) = 0 - (DDN_DDe(i, 1) * n[0] + metric.H(0, 1)*dn(0, 0) + metric.H(1, 1)*dn(0, 1) + metric.H(2, 1)*dn(0, 2));
				b(1, index + 1) = 0 - (DDN_DDe(i, 1) * n[1] + metric.H(0, 1)*dn(1, 0) + metric.H(1, 1)*dn(1, 1) + metric.H(2, 1)*dn(1, 2));
				b(1, index + 2) = 0 - (DDN_DDe(i, 1) * n[2] + metric.H(0, 1)*dn(2, 0) + metric.H(1, 1)*dn(2, 1) + metric.H(2, 1)*dn(2, 2));

				//third line
				b(2, index) = 0 - (DDN_DDe(i, 2) * n[0] + metric.H(0, 2)*dn(0, 0) + metric.H(1, 2)*dn(0, 1) + metric.H(2, 2)*dn(0, 2));
				b(2, index + 1) = 0 - (DDN_DDe(i, 2) * n[1] + metric.H(0, 2)*dn(1, 0) + metric.H(1, 2)*dn(1, 1) + metric.H(2, 2)*dn(1, 2));
				b(2, index + 2) = 0 - (DDN_DDe(i, 2) * n[2] + metric.H(0, 2)*dn(2, 0) + metric.H(1, 2)*dn(2, 1) + metric.H(2, 2)*dn(2, 2));
			}

			rB = prod(mInitialMetric.Q, b);
		}
		else
		{
			KRATOS_ERROR << "Element does not provide SHAPE_FUNCTION_LOCAL_DERIVATIVES and SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES" << std::endl;
		}
		KRATOS_CATCH("")
	}
	//***********************************************************************************
	//***********************************************************************************
	void MeshlessBaseSurfaceElement::CalculateSecondVariationStrainCurvature(
		Matrix& Strain_in_Q_coordinates11,
		Matrix& Strain_in_Q_coordinates22,
		Matrix& Strain_in_Q_coordinates12,
		Matrix& Curvature_in_Q_coordinates11,
		Matrix& Curvature_in_Q_coordinates22,
		Matrix& Curvature_in_Q_coordinates12,
		const MetricVariables& rMetric)
	{
		if (this->Has(SHAPE_FUNCTION_LOCAL_DERIVATIVES) && this->Has(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES))
		{
			const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
			const Matrix& DDN_DDe = this->GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

			const unsigned int number_of_nodes = GetGeometry().size();
			Matrix dg3_n = ZeroMatrix(3, 3);
			Matrix dg3_m = ZeroMatrix(3, 3);

			Vector ddStrain_curvilinear = ZeroVector(3);
			Vector ddCurvature_curvilinear = ZeroVector(3);

			Matrix dn_n = ZeroMatrix(3, 3);
			Matrix dn_m = ZeroMatrix(3, 3);

			//basis vector g3
			//array_1d<double, 3> g3;
			array_1d<double, 3> ddn;


			//CrossProduct(g3, g1, g2);

			//Matrix H = ZeroMatrix(3, 3);

			//this->Hessian(H, DDN_DDe);

			//differential area dA
			//double dA = norm_2(g3);

			//normal vector n
			array_1d<double, 3> n = rMetric.g3 / rMetric.dA;

			double invdA = 1 / rMetric.dA;
			double invdA3 = 1 / pow(rMetric.dA, 3);
			double invdA5 = 1 / pow(rMetric.dA, 5);

			for (unsigned int n = 0; n < number_of_nodes; n++)
			{
				//first line --- dg(1,0)*g2(2)-dg(2,0)*g2(1) + g1(1)*dg(2,1)-g1(2)*dg(1,1);
				dg3_n(0, 0) = 0;
				dg3_n(0, 1) = -DN_De(n, 0) * rMetric.g2[2] + DN_De(n, 1)*rMetric.g1[2];
				dg3_n(0, 2) = DN_De(n, 0) * rMetric.g2[1] - DN_De(n, 1)*rMetric.g1[1];

				//second line --- dg(2,0)*g2(0)-dg(0,0)*g2(2) + g1(2)*dg(0,1)-g1(0)*dg(2,1);
				dg3_n(1, 0) = DN_De(n, 0) * rMetric.g2[2] - DN_De(n, 1)*rMetric.g1[2];
				dg3_n(1, 1) = 0;
				dg3_n(1, 2) = -DN_De(n, 0)*rMetric.g2[0] + DN_De(n, 1)*rMetric.g1[0];

				//third line --- dg(0,0)*g2(1)-dg(1,0)*g2(0) + g1(0)*dg(1,1)-g1(1)*dg(0,1);
				dg3_n(2, 0) = -DN_De(n, 0) * rMetric.g2[1] + DN_De(n, 1) * rMetric.g1[1];
				dg3_n(2, 1) = DN_De(n, 0) * rMetric.g2[0] - DN_De(n, 1) * rMetric.g1[0];
				dg3_n(2, 2) = 0;

				//std::cout << "dg3_n: " << dg3_n << std::endl;

				for (unsigned int i = 0; i < 3; i++)
				{
					double g3dg3n = (rMetric.g3[0] * dg3_n(i, 0) + rMetric.g3[1] * dg3_n(i, 1) + rMetric.g3[2] * dg3_n(i, 2));
					double g3dg3lg3n = g3dg3n*invdA3;

					dn_n(i, 0) = dg3_n(i, 0)*invdA - rMetric.g3[0] * g3dg3lg3n;
					dn_n(i, 1) = dg3_n(i, 1)*invdA - rMetric.g3[1] * g3dg3lg3n;
					dn_n(i, 2) = dg3_n(i, 2)*invdA - rMetric.g3[2] * g3dg3lg3n;

					//std::cout << dn_n << std::endl;

					for (unsigned int m = 0; m <= n; m++)
					{
						//first line --- dg(1,0)*g2(2)-dg(2,0)*g2(1) + g1(1)*dg(2,1)-g1(2)*dg(1,1);
						dg3_m(0, 0) = 0;
						dg3_m(0, 1) = -DN_De(m, 0) * rMetric.g2[2] + DN_De(m, 1)*rMetric.g1[2];
						dg3_m(0, 2) = DN_De(m, 0) * rMetric.g2[1] - DN_De(m, 1)*rMetric.g1[1];

						//second line --- dg(2,0)*g2(0)-dg(0,0)*g2(2) + g1(2)*dg(0,1)-g1(0)*dg(2,1);
						dg3_m(1, 0) = DN_De(m, 0) * rMetric.g2[2] - DN_De(m, 1)*rMetric.g1[2];
						dg3_m(1, 1) = 0;
						dg3_m(1, 2) = -DN_De(m, 0)* rMetric.g2[0] + DN_De(m, 1)*rMetric.g1[0];

						//third line --- dg(0,0)*g2(1)-dg(1,0)*g2(0) + g1(0)*dg(1,1)-g1(1)*dg(0,1);
						dg3_m(2, 0) = -DN_De(m, 0) * rMetric.g2[1] + DN_De(m, 1) * rMetric.g1[1];
						dg3_m(2, 1) = DN_De(m, 0) * rMetric.g2[0] - DN_De(m, 1) * rMetric.g1[0];
						dg3_m(2, 2) = 0;


						//std::cout << "dg3_m: " << dg3_m << std::endl;

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

								Strain_in_Q_coordinates11(3 * n + i, 3 * m + j) = mInitialMetric.Q(0, 0)*ddStrain_curvilinear[0]
									+ mInitialMetric.Q(0, 1)*ddStrain_curvilinear[1] + mInitialMetric.Q(0, 2)*ddStrain_curvilinear[2];
								Strain_in_Q_coordinates22(3 * n + i, 3 * m + j) = mInitialMetric.Q(1, 0)*ddStrain_curvilinear[0]
									+ mInitialMetric.Q(1, 1)*ddStrain_curvilinear[1] + mInitialMetric.Q(1, 2)*ddStrain_curvilinear[2];
								Strain_in_Q_coordinates12(3 * n + i, 3 * m + j) = mInitialMetric.Q(2, 0)*ddStrain_curvilinear[0]
									+ mInitialMetric.Q(2, 1)*ddStrain_curvilinear[1] + mInitialMetric.Q(2, 2)*ddStrain_curvilinear[2];

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

							double g3dg3m = (rMetric.g3[0] * dg3_m(j, 0) + rMetric.g3[1] * dg3_m(j, 1) + rMetric.g3[2] * dg3_m(j, 2));
							double g3dg3lg3m = g3dg3m*invdA3;

							dn_m(j, 0) = dg3_m(j, 0)*invdA - rMetric.g3[0] * g3dg3lg3m;
							dn_m(j, 1) = dg3_m(j, 1)*invdA - rMetric.g3[1] * g3dg3lg3m;
							dn_m(j, 2) = dg3_m(j, 2)*invdA - rMetric.g3[2] * g3dg3lg3m;


							double c = -(ddg3[0] * rMetric.g3[0] + ddg3[1] * rMetric.g3[1] + ddg3[2] * rMetric.g3[2]
								+ dg3_n(i, 0)*dg3_m(j, 0) + dg3_n(i, 1)*dg3_m(j, 1) + dg3_n(i, 2)*dg3_m(j, 2)
								)*invdA3;


							double d = 3.0*g3dg3n * g3dg3m * invdA5;

							ddn[0] = ddg3[0] * invdA3 - g3dg3lg3m * dg3_n(i, 0) - g3dg3lg3n * dg3_m(j, 0) + (c + d)*rMetric.g3[0];
							ddn[1] = ddg3[1] * invdA3 - g3dg3lg3m * dg3_n(i, 1) - g3dg3lg3n * dg3_m(j, 1) + (c + d)*rMetric.g3[1];
							ddn[2] = ddg3[2] * invdA3 - g3dg3lg3m * dg3_n(i, 2) - g3dg3lg3n * dg3_m(j, 2) + (c + d)*rMetric.g3[2];

							ddCurvature_curvilinear[0] = DDN_DDe(n, 0)*dn_m(j, i) + DDN_DDe(m, 0)*dn_n(i, j)
								+ rMetric.H(0, 0)*ddn[0] + rMetric.H(1, 0)*ddn[1] + rMetric.H(2, 0)*ddn[2];
							ddCurvature_curvilinear[1] = DDN_DDe(n, 1)*dn_m(j, i) + DDN_DDe(m, 1)*dn_n(i, j)
								+ rMetric.H(0, 1)*ddn[0] + rMetric.H(1, 1)*ddn[1] + rMetric.H(2, 1)*ddn[2];
							ddCurvature_curvilinear[2] = DDN_DDe(n, 2)*dn_m(j, i) + DDN_DDe(m, 2)*dn_n(i, j)
								+ rMetric.H(0, 2)*ddn[0] + rMetric.H(1, 2)*ddn[1] + rMetric.H(2, 2)*ddn[2];

							Curvature_in_Q_coordinates11(3 * n + i, 3 * m + j) =
								mInitialMetric.Q(0, 0)*ddCurvature_curvilinear[0]
								+ mInitialMetric.Q(0, 1)*ddCurvature_curvilinear[1]
								+ mInitialMetric.Q(0, 2)*ddCurvature_curvilinear[2];
							Curvature_in_Q_coordinates22(3 * n + i, 3 * m + j) =
								mInitialMetric.Q(1, 0)*ddCurvature_curvilinear[0]
								+ mInitialMetric.Q(1, 1)*ddCurvature_curvilinear[1]
								+ mInitialMetric.Q(1, 2)*ddCurvature_curvilinear[2];
							Curvature_in_Q_coordinates12(3 * n + i, 3 * m + j) =
								mInitialMetric.Q(2, 0)*ddCurvature_curvilinear[0]
								+ mInitialMetric.Q(2, 1)*ddCurvature_curvilinear[1]
								+ mInitialMetric.Q(2, 2)*ddCurvature_curvilinear[2];
						}
					}
				}
			}
		}
	}
} // Namespace Kratos


