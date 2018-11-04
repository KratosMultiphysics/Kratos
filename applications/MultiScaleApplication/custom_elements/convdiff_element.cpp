//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#include "convdiff_element.hpp"
#include "multiscale_application_variables.h"
//#include "../applications/convection_diffusion_application/convection_diffusion_application.h"
//#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "includes/constitutive_law.h"

#include <string>
#include <iomanip>
#include <typeinfo>

namespace Kratos
{
	// =====================================================================================
	//
	// Class ConvDiffElement
	//
	// =====================================================================================

	ConvDiffElement::ConvDiffElement(IndexType NewId,
		GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{
	}

	ConvDiffElement::ConvDiffElement(IndexType NewId,
		GeometryType::Pointer pGeometry,
		PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
		mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
	}

	ConvDiffElement::~ConvDiffElement()
	{
	}

	Element::Pointer ConvDiffElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
	{
		GeometryType::Pointer newGeom(GetGeometry().Create(ThisNodes));
		return Element::Pointer(new ConvDiffElement(NewId, newGeom, pProperties));
	}

	ConvDiffElement::IntegrationMethod ConvDiffElement::GetIntegrationMethod() const
	{
		return mThisIntegrationMethod;
	}

	void ConvDiffElement::Initialize()
	{
		KRATOS_TRY
		InitializeMaterial();
		KRATOS_CATCH("")
	}

	void ConvDiffElement::InitializeMaterial()
	{
		KRATOS_TRY

		// NOTE:
		// This is the standard (previous) implementation:
		// If we are here, it means that no one already set up the constitutive law vector
		// through the method SetValue<CONSTITUTIVE_LAW_POINTER>

		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

		//Constitutive Law initialization
		if (mConstitutiveLawVector.size() != integration_points.size())
		{
			mConstitutiveLawVector.resize(integration_points.size());
		}
		else
		{
			// check whether the constitutive law pointers have been already set up
			bool already_set_up = true;
			for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
			{
				if (mConstitutiveLawVector[i] == NULL)
					already_set_up = false;
			}
			if (already_set_up)
			{
				for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
				{
					mConstitutiveLawVector[i]->InitializeMaterial(GetProperties(), GetGeometry(),
						row(GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod), i));
				}
				return; // if so, we are done here!
			}
		}

		// NOTE:
		// This is the standard (previous) implementation:
		// If we are here, it means that no one already set up the constitutive law vector
		// through the method SetValue<CONSTITUTIVE_LAW_POINTER>
		if (GetProperties()[CONSTITUTIVE_LAW] != NULL)
		{
			for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
			{
				mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
				mConstitutiveLawVector[i]->InitializeMaterial(GetProperties(), GetGeometry(),
					row(GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod), i));
			}
		}
		else
		{
			KRATOS_THROW_ERROR(std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id())
		}

		KRATOS_CATCH("")
	}

	void ConvDiffElement::ResetConstitutiveLaw()
	{
		KRATOS_TRY

		if (GetProperties()[CONSTITUTIVE_LAW] != NULL)
		{
			for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
				mConstitutiveLawVector[i]->ResetMaterial(GetProperties(), GetGeometry(), row(GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod), i));
		}

		KRATOS_CATCH("")
	}

	void ConvDiffElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if (rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes, false);

		for (unsigned int i = 0; i < number_of_nodes; i++)
			rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
	}

	void ConvDiffElement::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
	{
		ElementalDofList.resize(0);

		for (unsigned int i = 0; i < GetGeometry().size(); i++)
		{
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(TEMPERATURE));
		}
	}

	int ConvDiffElement::Check(const ProcessInfo& rCurrentProcessInfo)
	{
		GeometryType& geom = GetGeometry();

		SizeType ndim = geom.WorkingSpaceDimension();

		const PropertiesType & props = this->GetProperties();

		if (ndim == 2)
			if (props.Has(THICKNESS) == false)
				KRATOS_THROW_ERROR(std::invalid_argument, "THICKNESS not provided for element ", this->Id());
		return 0;
	}

	void ConvDiffElement::CleanMemory()
	{
	}

	void ConvDiffElement::GetValuesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();

		if (values.size() != number_of_nodes) values.resize(number_of_nodes, false);

		for (unsigned int i = 0; i < number_of_nodes; i++)
		{
			unsigned int index = i;
			values[index] = GetGeometry()[i].GetSolutionStepValue(TEMPERATURE, Step);
		}
	}

	void ConvDiffElement::GetFirstDerivativesVector(Vector& values, int Step)
	{
	}

	void ConvDiffElement::GetSecondDerivativesVector(Vector& values, int Step)
	{
	}

	void ConvDiffElement::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
	{
	}

	void ConvDiffElement::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
	{
	}

	void ConvDiffElement::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		const PropertiesType& props = GetProperties();
		const GeometryType& geom = GetGeometry();
		const Matrix& N = geom.ShapeFunctionsValues(mThisIntegrationMethod);
		for (SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			mConstitutiveLawVector[i]->InitializeSolutionStep(props, geom, row(N, i), CurrentProcessInfo);
		KRATOS_CATCH("");
	}

	void ConvDiffElement::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		const PropertiesType& props = GetProperties();
		const GeometryType& geom = GetGeometry();
		const Matrix& N = geom.ShapeFunctionsValues();
		for (SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			mConstitutiveLawVector[i]->FinalizeSolutionStep(props, geom, row(N, i), CurrentProcessInfo);
		KRATOS_CATCH("");
	}

	void ConvDiffElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ConvDiffElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)
	{
		CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
	}

	void ConvDiffElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)
	{
		Matrix dummy;
		CalculateAll(dummy, rRightHandSideVector, rCurrentProcessInfo, true, true);
	}

	// =====================================================================================
	//
	// Class Q4RIStabElement - Results on Gauss Points
	//
	// =====================================================================================

	void ConvDiffElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
		std::vector<double>& rOutput,
		const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

	void ConvDiffElement::CalculateOnIntegrationPoints(const Variable<Vector >& rVariable,
		std::vector< Vector >& rOutput,
		const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

	void ConvDiffElement::CalculateOnIntegrationPoints(const Variable< Matrix >& rVariable,
		std::vector< Matrix >& rOutput,
		const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

	void ConvDiffElement::SetValueOnIntegrationPoints(const Variable<double>& rVariable,
		std::vector<double>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		for (unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++)
		{
			mConstitutiveLawVector[PointNumber]->SetValue(rVariable, rValues[PointNumber], rCurrentProcessInfo);
		}
	}

	void ConvDiffElement::SetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
		std::vector<Vector>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		for (unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++)
		{
			mConstitutiveLawVector[PointNumber]->SetValue(rVariable, rValues[PointNumber], rCurrentProcessInfo);
		}
	}

	void ConvDiffElement::SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
		std::vector<Matrix>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		for (unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++)
		{
			mConstitutiveLawVector[PointNumber]->SetValue(rVariable, rValues[PointNumber], rCurrentProcessInfo);
		}
	}

	void ConvDiffElement::SetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
		std::vector<ConstitutiveLaw::Pointer>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if (mConstitutiveLawVector.size() != rValues.size())
		{
			mConstitutiveLawVector.resize(rValues.size());
			if (mConstitutiveLawVector.size() != GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod))
				KRATOS_THROW_ERROR(std::logic_error, "constitutive law not has the correct size ", mConstitutiveLawVector.size());
		}
		for (unsigned int i = 0; i<rValues.size(); i++)
			mConstitutiveLawVector[i] = rValues[i];
	}

	void ConvDiffElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
		std::vector<double>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		const unsigned int integration_points_number = GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod);

		if (rValues.size() != integration_points_number)
			rValues.resize(integration_points_number);

		for (unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++)
		{
			rValues[PointNumber] = mConstitutiveLawVector[PointNumber]->GetValue(rVariable, rValues[PointNumber]);
		}
	}

	void ConvDiffElement::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
		std::vector<Vector>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber(mThisIntegrationMethod);

		if (rValues.size() != integration_points_number)
			rValues.resize(integration_points_number);

		//if (rVariable == HEAT_FLUX_RVE)
		//{
		//	const PropertiesType& props = GetProperties();
		//	const GeometryType& geom = GetGeometry();
		//	const unsigned int number_of_points = GetGeometry().size();
		//	const unsigned int nDim = GetGeometry().WorkingSpaceDimension();
		//	//reading integration points and local gradients
		//	const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		//	const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients();
		//	const Matrix& N = GetGeometry().ShapeFunctionsValues();

		//	Element::GeometryType::JacobiansType J0;
		//	Vector mN(number_of_points);
		//	Matrix DN_DX(number_of_points, nDim);
		//	Matrix InvJ0(nDim, nDim);
		//	double IntToReferenceWeight = 0.0;

		//	GetGeometry().Jacobian(J0);
		//	double DetJ0;

		//	for (size_t PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
		//	{
		//		ConstitutiveLaw::Parameters matpar;

		//		Flags& options = matpar.GetOptions();
		//		options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
		//		options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
		//		options.Set(ConstitutiveLaw::INITIAL_CONFIGURATION);
		//		double detF = 1.0;
		//		double detF0 = 1.0;
		//		Matrix F(IdentityMatrix(2, 2));
		//		Matrix F0(IdentityMatrix(2, 2));
		//		matpar.SetDeterminantF(detF);
		//		matpar.SetDeterminantF0(detF0);
		//		matpar.SetDeformationGradientF(F);
		//		matpar.SetDeformationGradientF0(F0);

		//		matpar.SetElementGeometry(geom);
		//		matpar.SetMaterialProperties(props);
		//		matpar.SetProcessInfo(rCurrentProcessInfo);
		//		Vector q(nDim, 0.0); // physical stress vector
		//		Vector grad_T(nDim, 0.0); // (gT_x, gT_y)
		//		Matrix K(nDim, nDim, 0.0); // physical material tangent
		//		matpar.SetStrainVector(grad_T);
		//		matpar.SetStressVector(q);
		//		matpar.SetConstitutiveMatrix(K);

		//		//calculating inverse jacobian and jacobian determinant
		//		MathUtils<double>::InvertMatrix(J0[PointNumber], InvJ0, DetJ0);

		//		//Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
		//		noalias(mN) = row(N, PointNumber);
		//		noalias(DN_DX) = prod(DN_De[PointNumber], InvJ0);

		//		IntToReferenceWeight = integration_points[PointNumber].Weight() * DetJ0;
		//		if (GetGeometry().WorkingSpaceDimension() == 2)
		//			IntToReferenceWeight *= GetProperties()[THICKNESS];

		//		// E = B * T
		//		for (size_t i = 0; i < nDim; i++)
		//		{
		//			for (size_t ii = 0; ii < number_of_points; ii++)
		//			{
		//				grad_T[i] += DN_DX(ii, i) * GetGeometry()[ii].FastGetSolutionStepValue(TEMPERATURE);
		//			}
		//		}

		//		/******************** equivalent CLaw for determine Conductivity and HeatFlux *********************/

		//		// compute material response
		//		matpar.SetShapeFunctionsDerivatives(DN_DX);
		//		matpar.SetShapeFunctionsValues(mN);

		//		mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(matpar);

		//		//Vector q_aux = prod(K, grad_T);
		//		//std::stringstream ss;
		//		//ss << "K = " << K << ", " << std::endl;
		//		//ss << "grad_T = " << grad_T(1) << ", " << std::endl;
		//		//ss << "HEAT_FLUX_RVE_FROM_ELEM = " << q_aux(1) << ", " << std::endl;
		//		//std::cout << ss.str();

		//		Vector heat_flux;
		//		if (nDim == 2)
		//			heat_flux.resize(3,false);
		//		else
		//			heat_flux.resize(6,false);

		//		heat_flux(0) = q(0); //[W/mm^2]
		//		heat_flux(1) = q(1);
		//		heat_flux(2) = 0.0;
		//		if (nDim == 3)
		//		{
		//			heat_flux(2) = q(2); //[W/mm^2]
		//			heat_flux(3) = 0.0;
		//			heat_flux(4) = 0.0;
		//			heat_flux(5) = 0.0;
		//		}

		//		rValues[PointNumber] = heat_flux; // prod(K, grad_T); //

		//	}
		//}
		//else
		if (rVariable == FLUX_RVE || rVariable == HEAT_FLUX_RVE)
		{
			for (unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++)
			{
				rValues[PointNumber] = mConstitutiveLawVector[PointNumber]->GetValue(rVariable, rValues[PointNumber]);
			}
		}
		else
		{
			for (unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++)
			{
				rValues[PointNumber] = mConstitutiveLawVector[PointNumber]->GetValue(rVariable, rValues[PointNumber]);
			}
		}
	}

	void ConvDiffElement::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
		std::vector<Matrix>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		const unsigned int& integration_points_number = mConstitutiveLawVector.size();

		if (rValues.size() != integration_points_number)
			rValues.resize(integration_points_number);

		for (unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++)
		{
			rValues[PointNumber] = mConstitutiveLawVector[PointNumber]->GetValue(rVariable, rValues[PointNumber]);
		}
	}

	void ConvDiffElement::GetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
		std::vector<ConstitutiveLaw::Pointer>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if (rValues.size() != mConstitutiveLawVector.size())
			rValues.resize(mConstitutiveLawVector.size());

		for (unsigned int i = 0; i<rValues.size(); i++)
			rValues[i] = mConstitutiveLawVector[i];
	}


	// =====================================================================================
	//
	// Class ConvDiffElement - Private methods
	//
	// =====================================================================================

	void ConvDiffElement::CalculateBMatrix(double& A, Matrix& B)
	{
	}

	void ConvDiffElement::AddBodyForces(double V, VectorType& rRightHandSideVector)
	{
	}

	void ConvDiffElement::CalculateAll(MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo,
		const bool LHSrequired,
		const bool RHSrequired)
	{
		KRATOS_TRY
		const PropertiesType& props = GetProperties();
		const GeometryType& geom = GetGeometry();
		const unsigned int number_of_points = GetGeometry().size();
		const unsigned int nDim = GetGeometry().WorkingSpaceDimension();

		//resizing as needed the LHS
		if (rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points, number_of_points, false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points, number_of_points); //resetting LHS


		//resizing as needed the RHS
		if (rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points, false);
		rRightHandSideVector = ZeroVector(number_of_points); //resetting RHS

		//reading integration points and local gradients
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients();
		const Matrix& N = GetGeometry().ShapeFunctionsValues();

		Element::GeometryType::JacobiansType J0;
		Vector mN(number_of_points);
		Matrix DN_DX(number_of_points, nDim);
		Matrix InvJ0(nDim, nDim);
		double IntToReferenceWeight = 0.0;
		Vector q(nDim, 0.0); // physical stress vector
		Matrix K(nDim, nDim, 0.0); // physical material tangent

		GetGeometry().Jacobian(J0);
		double DetJ0;
		for (size_t PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
		{
			//calculating inverse jacobian and jacobian determinant
			MathUtils<double>::InvertMatrix(J0[PointNumber], InvJ0, DetJ0);

			//Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
			noalias(mN) = row(N, PointNumber);
			noalias(DN_DX) = prod(DN_De[PointNumber], InvJ0);

			IntToReferenceWeight = integration_points[PointNumber].Weight() * DetJ0;
			if (GetGeometry().WorkingSpaceDimension() == 2)
				IntToReferenceWeight *= GetProperties()[THICKNESS];

			Vector grad_T(nDim, 0.0); // (gT_x, gT_y)
			// E = B * T
			for (size_t i = 0; i < nDim; i++)
			{
				for (size_t ii = 0; ii < number_of_points; ii++)
				{
					grad_T[i] += DN_DX(ii, i) * GetGeometry()[ii].FastGetSolutionStepValue(TEMPERATURE);
				}
			}

			/******************** equivalent CLaw for determine Conductivity and HeatFlux *********************/

			// compute material response

			ConstitutiveLaw::Parameters matpar;
			matpar.SetElementGeometry(geom);
			matpar.SetMaterialProperties(props);
			matpar.SetProcessInfo(rCurrentProcessInfo);
			matpar.SetStrainVector(grad_T);
			matpar.SetStressVector(q);
			matpar.SetConstitutiveMatrix(K);
			matpar.SetShapeFunctionsDerivatives(DN_DX);
			matpar.SetShapeFunctionsValues(mN);

			Flags& options = matpar.GetOptions();
			options.Set(ConstitutiveLaw::COMPUTE_STRESS, RHSrequired);
			options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, LHSrequired);

			if (RHSrequired || LHSrequired)
			{
				mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(matpar);
				if (LHSrequired)
				{
					//CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
					// LHS = B' * D * B
					noalias(rLeftHandSideMatrix) += prod(DN_DX, Matrix(prod(K*IntToReferenceWeight, trans(DN_DX))));
				}

				// Add all contributions to the residual vector
				noalias(rRightHandSideVector) -= prod(DN_DX, q*IntToReferenceWeight);
			}
		}
		//// Add all contributions to the residual vector
		//noalias(rRightHandSideVector) -= prod(DN_DX, q);
		//rRightHandSideVector *= IntToReferenceWeight;

		//// RHS -= LHS*temperatures
		//noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, temp);
		//std::cout << "rRightHandSideVector = " << rRightHandSideVector << std::endl;
		KRATOS_CATCH("");
	}

	// =====================================================================================
	//
	// Class ConvDiffElement - Serialization
	//
	// =====================================================================================

	void ConvDiffElement::save(Serializer& rSerializer) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
		int IntMethod = int(mThisIntegrationMethod);
		rSerializer.save("IntegrationMethod", IntMethod);
		rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
	}

	void ConvDiffElement::load(Serializer& rSerializer)
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
		int IntMethod;
		rSerializer.load("IntegrationMethod", IntMethod);
		mThisIntegrationMethod = IntegrationMethod(IntMethod);
		rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
	}

}
