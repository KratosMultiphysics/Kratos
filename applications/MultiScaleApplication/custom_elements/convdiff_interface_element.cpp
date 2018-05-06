//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#include "convdiff_interface_element.hpp"
#include "multiscale_application_variables.h"

#include <string>
#include <iomanip>
#include "custom_utilities/math_helpers.h"

#define K_GET_SIGN(X) return X < 0.0 ? -1.0 : 1.0

namespace Kratos
{

	namespace Utilities
	{

	}

	// =====================================================================================
	//
	// Class ConvDiffInterfaceElement
	//
	// =====================================================================================

	ConvDiffInterfaceElement::ConvDiffInterfaceElement(IndexType NewId,
		GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
		, mInitialized(false)
		, m_delta_disp()
	{
	}

	ConvDiffInterfaceElement::ConvDiffInterfaceElement(IndexType NewId,
		GeometryType::Pointer pGeometry,
		PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
		, mInitialized(false)
		, m_delta_disp()
	{
	}

	ConvDiffInterfaceElement::~ConvDiffInterfaceElement()
	{
	}

	Element::Pointer ConvDiffInterfaceElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
	{
		GeometryType::Pointer newGeom(GetGeometry().Create(ThisNodes));
		return Element::Pointer(new ConvDiffInterfaceElement(NewId, newGeom, pProperties));
	}

	void ConvDiffInterfaceElement::Initialize()
	{
		if (mInitialized == false)
		{
			ConstitutiveLaw::Pointer& pLaw = GetProperties()[CONSTITUTIVE_LAW];
			const GeometryType::IntegrationPointsArrayType& integrationPoints = GetGeometry().IntegrationPoints();
			SizeType half_ngp = integrationPoints.size() / 2;

			mConstitutiveLawVector.clear();
			for (SizeType i = 0; i < half_ngp; i++)
			{
				ConstitutiveLaw::Pointer newCLaw = pLaw->Clone();
				newCLaw->InitializeMaterial(GetProperties(), GetGeometry(), row(GetGeometry().ShapeFunctionsValues(), i));
				mConstitutiveLawVector.push_back(newCLaw);
			}

			mInitialized = true;

			m_delta_disp = ZeroVector(GetGeometry().WorkingSpaceDimension());
		}
	}

	void ConvDiffInterfaceElement::ResetConstitutiveLaw()
    {
        for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			mConstitutiveLawVector[i]->ResetMaterial(GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues(), i ) );
	}

	void ConvDiffInterfaceElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		SizeType number_of_nodes = GetGeometry().PointsNumber();
		if (rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes, false);

		for (SizeType i = 0; i < number_of_nodes; i++)
			rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
	}

	void ConvDiffInterfaceElement::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
	{
        GeometryType & geom = this->GetGeometry();
		SizeType nnode = geom.size();
		ElementalDofList.resize(0);
        ElementalDofList.reserve(nnode);

		for (SizeType i = 0; i < nnode; i++)
		{
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(TEMPERATURE));
		}
	}

	int ConvDiffInterfaceElement::Check(const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		GeometryType& geom = GetGeometry();

		SizeType ndim = geom.WorkingSpaceDimension();

		// verify that the variables are correctly initialized
		if (TEMPERATURE.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument, "TEMPERATURE has Key zero! (check if the application is correctly registered", "");
		if (ndim == 2)
			if (THICKNESS.Key() == 0)
				KRATOS_THROW_ERROR(std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "");
		if (CONSTITUTIVE_LAW.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument, "CONSTITUTIVE_LAW has Key zero! (check if the application is correctly registered", "");

		// verify that the dofs exist
		for (unsigned int i = 0; i<geom.size(); i++)
		{
			if (geom[i].SolutionStepsDataHas(TEMPERATURE) == false)
				KRATOS_THROW_ERROR(std::invalid_argument, "missing variable TEMPERATURE on node ", geom[i].Id());
		}

		// check properties
		if (this->pGetProperties() == NULL)
			KRATOS_THROW_ERROR(std::logic_error, "Properties not provided for element ", this->Id());

		const PropertiesType & props = this->GetProperties();

		if (ndim == 2)
			if (props.Has(THICKNESS) == false)
				KRATOS_THROW_ERROR(std::invalid_argument, "THICKNESS not provided for element ", this->Id());

		if (props.Has(CONSTITUTIVE_LAW) == false)
			KRATOS_THROW_ERROR(std::invalid_argument, "CONSTITUTIVE_LAW not provided for element ", this->Id());

		const ConstitutiveLaw::Pointer& pLaw = props[CONSTITUTIVE_LAW];
		if (pLaw == NULL)
			KRATOS_THROW_ERROR(std::invalid_argument, "CONSTITUTIVE_LAW not provided for element ", this->Id());

		pLaw->Check(props, geom, rCurrentProcessInfo);

		return 0;

		KRATOS_CATCH("")
	}

	void ConvDiffInterfaceElement::CleanMemory()
	{
	}

	void ConvDiffInterfaceElement::GetValuesVector(Vector& values, int Step)
	{
		const SizeType number_of_nodes = GetGeometry().size();

		if (values.size() != number_of_nodes) values.resize(number_of_nodes, false);

		for (SizeType i = 0; i < number_of_nodes; i++)
		{
			SizeType index = i;
			values[index] = GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE, Step);
		}
	}

	void ConvDiffInterfaceElement::GetFirstDerivativesVector(Vector& values, int Step)
	{
	}

	void ConvDiffInterfaceElement::GetSecondDerivativesVector(Vector& values, int Step)
	{
	}

	void ConvDiffInterfaceElement::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
	{
	}

	void ConvDiffInterfaceElement::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
	{
	}

	void ConvDiffInterfaceElement::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		const PropertiesType& props = GetProperties();
		const GeometryType& geom = GetGeometry();
		const Matrix& N = geom.ShapeFunctionsValues();
		for (SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			mConstitutiveLawVector[i]->InitializeSolutionStep(props, geom, row(N, i), CurrentProcessInfo);
		KRATOS_CATCH("");
	}

	void ConvDiffInterfaceElement::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		const PropertiesType& props = GetProperties();
		const GeometryType& geom = GetGeometry();
		const Matrix& N = geom.ShapeFunctionsValues();
		for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			mConstitutiveLawVector[i]->FinalizeSolutionStep( props, geom, row( N, i ), CurrentProcessInfo );
		KRATOS_CATCH("");
    }

	void ConvDiffInterfaceElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        SizeType nnode = GetGeometry().size();
        if((rMassMatrix.size1() != nnode) || (rMassMatrix.size2() != nnode))
            rMassMatrix.resize(nnode, nnode, false);

        noalias( rMassMatrix ) = ZeroMatrix(nnode, nnode);
    }

	void ConvDiffInterfaceElement::CalculateDampingMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
    {
		SizeType nnode = GetGeometry().size();
        if((rDampMatrix.size1() != nnode) || (rDampMatrix.size2() != nnode))
            rDampMatrix.resize(nnode, nnode, false);

        noalias( rDampMatrix ) = ZeroMatrix(nnode, nnode);
    }

	void ConvDiffInterfaceElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)
	{
		CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
	}

	void ConvDiffInterfaceElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)
	{
		Matrix dummy;
		CalculateAll(dummy, rRightHandSideVector, rCurrentProcessInfo, true, true);
	}

	// =====================================================================================
	//
	// Class ConvDiffInterfaceElement - Results on Gauss Points
	//
	// =====================================================================================

	void ConvDiffInterfaceElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
																		 std::vector<double>& rOutput,
																		 const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

	void ConvDiffInterfaceElement::CalculateOnIntegrationPoints(const Variable< array_1d< double, 3 > >& rVariable,
																		 std::vector< array_1d<double, 3 > >& rOutput,
																		 const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}
	void ConvDiffInterfaceElement::CalculateOnIntegrationPoints(const Variable<Vector >& rVariable,
		std::vector< Vector >& rOutput,
		const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

	void ConvDiffInterfaceElement::CalculateOnIntegrationPoints(const Variable< Matrix >& rVariable,
																		 std::vector< Matrix >& rOutput,
																		 const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

	void ConvDiffInterfaceElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
		std::vector<double>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		SizeType num_gp = GetGeometry().IntegrationPoints().size();
		if (rValues.size() != num_gp)
			rValues.resize(num_gp);
		SizeType half_num_gp = num_gp / 2;
		SizeType ndim = GetGeometry().WorkingSpaceDimension();

		double res(0.0);
		if(ndim == 2)
		{
			res = 0.0;
			mConstitutiveLawVector[0]->GetValue(rVariable, res);
			rValues[0] = rValues[3] = res;
			res = 0.0;
			mConstitutiveLawVector[1]->GetValue(rVariable, res);
			rValues[1] = rValues[2] = res;
		}
		else
		{
			for (SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			{
				SizeType j = i + half_num_gp;
				res = 0.0;
				mConstitutiveLawVector[i]->GetValue(rVariable, res);
				rValues[i] = res;
				rValues[j] = res;
			}
		}
	}

	void ConvDiffInterfaceElement::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
		std::vector<Vector>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		SizeType num_gp = GetGeometry().IntegrationPoints().size();
		if (rValues.size() != num_gp)
			rValues.resize(num_gp);
		SizeType half_num_gp = num_gp / 2;
		SizeType ndim = GetGeometry().WorkingSpaceDimension();

		Vector res;
		if(ndim == 2)
		{
			res.clear();
			mConstitutiveLawVector[0]->GetValue(rVariable, res);
			rValues[0] = rValues[3] = res;
			res.clear();
			mConstitutiveLawVector[1]->GetValue(rVariable, res);
			rValues[1] = rValues[2] = res;
		}
		else
		{
			for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			{
				SizeType j = i+half_num_gp;
				res.clear();
				mConstitutiveLawVector[i]->GetValue(rVariable, res);
				rValues[i] = res;
				rValues[j] = res;
			}
		}
	}

	void ConvDiffInterfaceElement::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
		std::vector<Matrix>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
    {
		SizeType num_gp = GetGeometry().IntegrationPoints().size();
		if(rValues.size() != num_gp) rValues.resize(num_gp);
		SizeType half_num_gp = num_gp / 2;
		SizeType ndim = GetGeometry().WorkingSpaceDimension();

		Matrix res;
		if(ndim == 2)
		{
			res.clear();
			mConstitutiveLawVector[0]->GetValue(rVariable, res);
			rValues[0] = rValues[3] = res;
			res.clear();
			mConstitutiveLawVector[1]->GetValue(rVariable, res);
			rValues[1] = rValues[2] = res;
		}
		else
		{
			for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			{
				SizeType j = i+half_num_gp;
				res.clear();
				mConstitutiveLawVector[i]->GetValue(rVariable, res);
				rValues[i] = res;
				rValues[j] = res;
			}
		}
    }

	void ConvDiffInterfaceElement::GetValueOnIntegrationPoints(const Variable<array_1d<double, 3> >& rVariable,
                                                           std::vector<array_1d<double,3> >& rValues,
                                                           const ProcessInfo& rCurrentProcessInfo)
    {
		SizeType num_gp = GetGeometry().IntegrationPoints().size();
		if(rValues.size() != num_gp) rValues.resize(num_gp);
		SizeType half_num_gp = num_gp / 2;
		SizeType ndim = GetGeometry().WorkingSpaceDimension();

		array_1d<double, 3> res;
		if (ndim == 2)
		{
			res.clear();
			mConstitutiveLawVector[0]->GetValue(rVariable, res);
			rValues[0] = rValues[3] = res;
			res.clear();
			mConstitutiveLawVector[1]->GetValue(rVariable, res);
			rValues[1] = rValues[2] = res;
		}
		else
		{
			for (SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			{
				SizeType j = i + half_num_gp;
				res.clear();
				mConstitutiveLawVector[i]->GetValue(rVariable, res);
				rValues[i] = res;
				rValues[j] = res;
			}
		}
    }

	void ConvDiffInterfaceElement::GetValueOnIntegrationPoints(const Variable<array_1d<double, 6> >& rVariable,
                                                           std::vector<array_1d<double,6> >& rValues,
                                                           const ProcessInfo& rCurrentProcessInfo)
    {
		SizeType num_gp = GetGeometry().IntegrationPoints().size();
		if(rValues.size() != num_gp) rValues.resize(num_gp);
		SizeType half_num_gp = num_gp / 2;
		SizeType ndim = GetGeometry().WorkingSpaceDimension();

		array_1d<double,6> res;
		if(ndim == 2)
		{
			res.clear();
			mConstitutiveLawVector[0]->GetValue(rVariable, res);
			rValues[0] = rValues[3] = res;
			res.clear();
			mConstitutiveLawVector[1]->GetValue(rVariable, res);
			rValues[1] = rValues[2] = res;
		}
		else
		{
			for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			{
				SizeType j = i+half_num_gp;
				res.clear();
				mConstitutiveLawVector[i]->GetValue(rVariable, res);
				rValues[i] = res;
				rValues[j] = res;
			}
		}
    }

	// =====================================================================================
	//
	// Class ConvDiffInterfaceElement - Private methods
	//
	// =====================================================================================

	void ConvDiffInterfaceElement::DecimalCorrection(Vector& a)
	{
		double norm = norm_2(a);
		double tolerance = std::max(norm * 1.0E-12, 1.0E-12);
		for (SizeType i = 0; i < a.size(); i++)
			if (std::abs(a(i)) < tolerance)
				a(i) = 0.0;
	}

	void ConvDiffInterfaceElement::CalculatePermutation_disp(InterfaceIndexPermutation& p)
	{
		p.HasPermutation_disp = false;
		p.Permutation_disp.clear();
		GeometryType& geom = GetGeometry();
		if (geom.WorkingSpaceDimension() == 2)
		{
			// handle case of quad 2d 4n (this is the only one supported by now)
			if (geom.size() == 4)
			{
				p.HasPermutation_disp = true;
				// 1 2   3 4   7 8   5 6
				p.Permutation_disp.resize(8);
				p.Permutation_disp[0] = 0;
				p.Permutation_disp[1] = 1;
				p.Permutation_disp[2] = 2;
				p.Permutation_disp[3] = 3;
				p.Permutation_disp[4] = 6;
				p.Permutation_disp[5] = 7;
				p.Permutation_disp[6] = 4;
				p.Permutation_disp[7] = 5;
			}
		}
	}

	void ConvDiffInterfaceElement::CalculatePermutation_temp(InterfaceIndexPermutation& p)
	{
		p.HasPermutation_temp = false;
		p.Permutation_temp.clear();
		GeometryType& geom = GetGeometry();
		if (geom.WorkingSpaceDimension() == 2)
		{
			// handle case of quad 2d 4n (this is the only one supported by now)
			if (geom.size() == 4)
			{
				p.HasPermutation_temp = true;
				// Permutation_temp[LOCAL_NODE] = GLOBAL_NODE
				p.Permutation_temp.resize(4);
				p.Permutation_temp[0] = 0;
				p.Permutation_temp[1] = 1;
				p.Permutation_temp[2] = 3;
				p.Permutation_temp[3] = 2;
			}
		}
	}

	void ConvDiffInterfaceElement::CalculateDeltaPosition(Matrix& rDeltaPosition)
	{
		GeometryType& geom = GetGeometry();
		const unsigned int number_of_nodes = geom.PointsNumber();
		unsigned int dimension = geom.WorkingSpaceDimension();

		rDeltaPosition = zero_matrix<double>(number_of_nodes, dimension);

		for (unsigned int i = 0; i < number_of_nodes; i++)
		{
			const NodeType& iNode = geom[i];
			rDeltaPosition(i, 0) = iNode.X() - iNode.X0();
			rDeltaPosition(i, 1) = iNode.Y() - iNode.Y0();
			if (dimension == 3)
				rDeltaPosition(i, 2) = iNode.Z() - iNode.Z0();
		}
	}

	void ConvDiffInterfaceElement::CalculateJacobianAndTransformationMatrix(const SizeType pointID,
		Matrix& delta_position,
		Matrix& jacobian,
		double& J,
		Matrix& iR)
	{
		GeometryType& geom = GetGeometry();
		SizeType ndim = geom.WorkingSpaceDimension();

		geom.Jacobian(jacobian, pointID, GetIntegrationMethod(), delta_position);

		if (ndim == 2)
		{
			array_1d<double, 2> vx;
			vx[0] = jacobian(0, 0);
			vx[1] = jacobian(1, 0);
			J = std::sqrt(vx[0] * vx[0] + vx[1] * vx[1]);
			vx /= J;
			iR(0, 0) =  vx[0]; iR(0, 1) = vx[1];
			iR(1, 0) = -vx[1]; iR(1, 1) = vx[0];
		}
		else
		{
			array_1d<double, 3> vx;
			array_1d<double, 3> vy;
			vx[0] = jacobian(0, 0); vx[1] = jacobian(1, 0); vx[2] = jacobian(2, 0);
			vy[0] = jacobian(0, 1); vy[1] = jacobian(1, 1); vy[2] = jacobian(2, 1);
			array_1d<double, 3> vz;
			MathUtils<double>::CrossProduct(vz, vx, vy);
			MathUtils<double>::CrossProduct(vy, vz, vx);
			vx /= MathUtils<double>::Norm3(vx);
			vy /= MathUtils<double>::Norm3(vy);
			J = MathUtils<double>::Norm3(vz);
			vz /= J;
			for (SizeType i = 0; i < 3; i++)
			{
				iR(0, i) = vx[i];
				iR(1, i) = vy[i];
				iR(2, i) = vz[i];
			}
		}
	}

	double ConvDiffInterfaceElement::CalculateIntegrationWeight(double J,
		double iw)
	{
		double dV = J * iw * 2.0;
		if (GetGeometry().WorkingSpaceDimension() == 2)
			dV *= GetProperties()[THICKNESS];
		return dV;
	}

	void ConvDiffInterfaceElement::CalculateLocalDisplacementVector(const InterfaceIndexPermutation& P,
		const Matrix& R,
		const Vector& UG,
		Vector& UL)
	{
		GeometryType& geom = GetGeometry();
		SizeType nnodes = geom.size();
		SizeType ndim = geom.WorkingSpaceDimension();
		if (P.HasPermutation_disp)
		{
			for (SizeType inode = 0; inode < nnodes; inode++)
			{
				SizeType pos = inode*ndim;
				for (SizeType i = 0; i < ndim; i++)
				{
					double temp = 0.0;
					for (SizeType j = 0; j < ndim; j++)
					{
						temp += R(i, j) * UG(P.Permutation_disp[pos + j]);
					}
					UL(pos + i) = temp;
				}
			}
		}
		else
		{
			for (SizeType inode = 0; inode < nnodes; inode++)
			{
				SizeType pos = inode*ndim;
				for (SizeType i = 0; i < ndim; i++)
				{
					double temp = 0.0;
					for (SizeType j = 0; j < ndim; j++)
					{
						temp += R(i, j) * UG(pos + j);
					}
					UL(pos + i) = temp;
				}
			}
		}
	}

	void ConvDiffInterfaceElement::CalculateLocalTemperature(const InterfaceIndexPermutation& P,
		const Matrix& R, // TEMPERATURE HAS NOT TO BE ROTATED
		const Vector& TG,
		Vector& TL)
	{
		GeometryType& geom = GetGeometry();
		SizeType nnodes = geom.size();
		SizeType ndim = geom.WorkingSpaceDimension();
		if (P.HasPermutation_temp)
		{
			for (SizeType inode = 0; inode < nnodes; inode++)
			{
				TL(inode) = TG(P.Permutation_temp[inode]);
			}
		}
		else
		{
			for (SizeType inode = 0; inode < nnodes; inode++)
			{
				TL(inode) = TG(inode);
			}
		}
	}

	void ConvDiffInterfaceElement::TransformToGlobalAndAdd(const InterfaceIndexPermutation& P,
		const Matrix& R,
		const Matrix& LHS_local,
		const Vector& RHS_local,
		Matrix& LHS_global,
		Vector& RHS_global,
		const bool LHSrequired,
		const bool RHSrequired)
	{
		GeometryType& geom = GetGeometry();
		SizeType nnodes = geom.size();
		SizeType ndim = geom.WorkingSpaceDimension();

		Matrix RTK(ndim, ndim);

		if (P.HasPermutation_temp)
		{
			for (SizeType node_i = 0; node_i < nnodes; node_i++)
			{
				if (RHSrequired)
				{
					RHS_global(P.Permutation_temp[node_i]) += RHS_local(node_i);
				}

				if (LHSrequired)
				{
					for (SizeType node_j = 0; node_j < nnodes; node_j++)
					{
						LHS_global(P.Permutation_temp[node_i], P.Permutation_temp[node_j]) += LHS_local(node_i, node_j);
					}
				}
			}
		}
		else
		{
			for (SizeType node_i = 0; node_i < nnodes; node_i++)
			{
				if (RHSrequired)
				{
					RHS_global(node_i) += RHS_local(node_i);
				}

				if (LHSrequired)
				{
					for (SizeType node_j = 0; node_j < nnodes; node_j++)
					{
						LHS_global(node_i, node_j) += LHS_local(node_i, node_j);
					}
				}
			}
		}
	}

	void ConvDiffInterfaceElement::CalculateBMatrix(const SizeType pointID,
		Matrix& B)
	{
		/* DELTA DISPLACEMENT MATRIX
		Phi = [-I | I]
		where:
		I = identity with size = [n1 x 2*n1] where n1 =  num_dim * num_nodes / 2
		*/

		/* INTERPOLATION MATRIX
		H = [N0 | N1 | ... | Nn]
		where:
		Ni = Identity(num_dim x num_dim) * N(of node i),
		with i = 1 to num_nodes/2
		*/

		/* STRAIN DISPLACEMENT MATRIX
		B = H * Phi
		such that
		Delta_U(all nodes / 2) = Phi*U
		Delta_U(at gauss_points) = H * Delta_U = H * Phi * U = B * U
		*/

		GeometryType& geom = GetGeometry();
		const Matrix& shapes = geom.ShapeFunctionsValues();

		SizeType ndim = geom.WorkingSpaceDimension();
		SizeType nnodes = geom.size();
		SizeType half_nnodes = nnodes / 2;

		noalias(B) = ZeroMatrix(ndim,nnodes);
		for (SizeType k = 0; k < half_nnodes; k++)
		{
			double Nk = shapes(pointID, k);

			SizeType pos1 = k;
			SizeType pos2 = k + half_nnodes;
			for (SizeType i = 0; i < ndim; i++)
			{
				B(i, pos1) = -Nk;
				B(i, pos2) = Nk;
			}
		}
		//std::stringstream ss;
		//ss << "--------  B = " << B << ", " << std::endl;
		//std::cout << ss.str();
	}

	void ConvDiffInterfaceElement::CalculateBMatrixU(const SizeType pointID,
		Matrix& B)
	{
		/* DELTA DISPLACEMENT MATRIX
		Phi = [-I | I]
		where:
		I = identity with size = [n1 x 2*n1] where n1 =  num_dim * num_nodes / 2
		*/

		/* INTERPOLATION MATRIX
		H = [N0 | N1 | ... | Nn]
		where:
		Ni = Identity(num_dim x num_dim) * N(of node i),
		with i = 1 to num_nodes/2
		*/

		/* STRAIN DISPLACEMENT MATRIX
		B = H * Phi
		such that
		Delta_U(all nodes / 2) = Phi*U
		Delta_U(at gauss_points) = H * Delta_U = H * Phi * U = B * U
		*/

		GeometryType& geom = GetGeometry();
		const Matrix& shapes = geom.ShapeFunctionsValues();

		SizeType ndim = geom.WorkingSpaceDimension();
		SizeType nnodes = geom.size();
		SizeType half_nnodes = nnodes / 2;

		for (SizeType k = 0; k < half_nnodes; k++)
		{
			double Nk = shapes(pointID, k);

			SizeType pos1 = k * ndim;
			SizeType pos2 = pos1 + half_nnodes*ndim;
			for (SizeType i = 0; i < ndim; i++)
			{
				B(i, pos1 + i) = -Nk;
				B(i, pos2 + i) = Nk;
			}
		}
	}

	void ConvDiffInterfaceElement::CalculateGeneralizedTemperatureStrains(const SizeType pointID,
		const Matrix& B,
		const Vector& T,
		Vector& generalizedTemperatureStrains)
	{
		noalias(generalizedTemperatureStrains) = prod(B, T);
	}

	void ConvDiffInterfaceElement::CalculateGeneralizedDeltaDisplacement(const SizeType pointID,
		const Matrix& B,
		const Vector& U,
		Vector& deltaDispl)
	{
		noalias(deltaDispl) = prod(B, U);
	}

	void ConvDiffInterfaceElement::GetDisplacementVector(Vector& values)
	{
		const GeometryType & geom = GetGeometry();

		SizeType dim = geom.WorkingSpaceDimension();
		SizeType ndofs = dim * geom.size();
		if (values.size() != ndofs)
			values.resize(ndofs, false);

		for (SizeType i = 0; i < geom.size(); i++)
		{
			const NodeType & iNode = geom[i];
			const array_1d<double, 3>& disp = iNode.FastGetSolutionStepValue(DISPLACEMENT);

			int index = i * dim;
			values[index] = disp[0];
			values[index + 1] = disp[1];
			if (dim == 3)
				values[index + 2] = disp[2];
		}
	}

	void ConvDiffInterfaceElement::CalculateAll(MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo,
		const bool LHSrequired,
		const bool RHSrequired)
	{
		// general data
		PropertiesType props = GetProperties();
		const GeometryType& geom = GetGeometry();
		SizeType ndim = geom.WorkingSpaceDimension();
		SizeType nnodes = geom.size();
		SizeType ndofs = ndim * nnodes;

		// resize the LHS matrix
		if (LHSrequired) {
			if (rLeftHandSideMatrix.size1() != nnodes || rLeftHandSideMatrix.size2() != nnodes)
				rLeftHandSideMatrix.resize(nnodes, nnodes, false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(nnodes, nnodes);
		}

		// resize the RHS vector
		if (RHSrequired) {
			if (rRightHandSideVector.size() != nnodes)
				rRightHandSideVector.resize(nnodes, false);
			noalias(rRightHandSideVector) = ZeroVector(nnodes);
		}

		// global and local displacement vectors
		Vector globalDisplacements(ndofs);
		Vector localDisplacements(ndofs);
		GetDisplacementVector(globalDisplacements);

		// delta position for the jacobian computation
		// with respect to the reference configuration
		Matrix delta_position;
		CalculateDeltaPosition(delta_position);

		// Get TEMPERATURE
		Vector globalTemperature(nnodes,0.0);
		Vector localTemperature(nnodes,0.0);
		GetValuesVector(globalTemperature);

		// jacobian matrix and its determinant
		Matrix jacobian(ndim, ndim - 1);
		double J;

		// strain-temperature matrix
		Matrix B(ndim, nnodes, 0.0);
		Matrix B_u(ndim, ndofs, 0.0);

		// material point calculation data
		Matrix D(ndim, ndim);
		Vector generalizedTemperatureStrains(ndim);

		// transformation matrices
		Matrix iR(ndim, ndim);

		// permutation data
		InterfaceIndexPermutation permutation_disp;
		InterfaceIndexPermutation permutation_temp;
		CalculatePermutation_disp(permutation_disp);
		CalculatePermutation_temp(permutation_temp);

		// LHS and RHS in local coordinate system
		Matrix Kloc(nnodes, nnodes);
		Vector Rloc(nnodes);

		// auxiliary data to avoid extra memory allocations
		Matrix BTD(nnodes, ndim);

		// initialize material parameters
		ConstitutiveLaw::Parameters Parameters(GetGeometry(), GetProperties(),rCurrentProcessInfo);
		Flags &ConstitutiveLawOptions = Parameters.GetOptions();
		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
		if (LHSrequired) ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

		Parameters.SetStrainVector(generalizedTemperatureStrains);
		Parameters.SetStressVector(mgeneralizedTemperatureStresses);
		Parameters.SetConstitutiveMatrix(D);

		Matrix F(IdentityMatrix(ndim, ndim));
		Matrix F0(IdentityMatrix(ndim, ndim));
		double detF(1.0);
		double detF0(1.0);
		Vector N(nnodes);
		Matrix DN_DX(nnodes, ndim, 0.0);
		Parameters.SetDeformationGradientF(F);
		Parameters.SetDeterminantF(detF);
		Parameters.SetShapeFunctionsDerivatives(DN_DX); // our c.law doesn't use this!!!!!

		// END - INITIALIZE MATERIAL PARAMETERS ***************************************************************

		// loop over the integration points.
		// note: loop only the first half of them, then multiply the integration weight by 2
		const GeometryType::IntegrationPointsArrayType& integrationPoints = geom.IntegrationPoints();
		for (SizeType intp_id = 0; intp_id < integrationPoints.size() / 2; intp_id++)
		{
			// the current integration point
			const GeometryType::IntegrationPointType& ip = integrationPoints[intp_id];

			// jacobianobian and local transformation
			CalculateJacobianAndTransformationMatrix(intp_id, delta_position, jacobian, J, iR);

			// integration weight
			double dV = CalculateIntegrationWeight(J, ip.Weight());

			// calculate generalized temperature strains
			CalculateBMatrix(intp_id, B);
			CalculateLocalTemperature(permutation_temp,iR,globalTemperature,localTemperature);
			CalculateGeneralizedTemperatureStrains(intp_id, B, localTemperature, generalizedTemperatureStrains);

			// local displacement vector
			CalculateLocalDisplacementVector(permutation_disp, iR, globalDisplacements, localDisplacements);
			// calculate delta_disp
			CalculateBMatrixU(intp_id, B_u);
			CalculateGeneralizedDeltaDisplacement(intp_id, B_u, localDisplacements, m_delta_disp);
			// get the normal gap
			//std::cout << MathHelpers::VectorToString(m_delta_disp, 4, std::scientific);

			mConstitutiveLawVector[intp_id]->SetValue(GAP_INTERFACE, m_delta_disp, rCurrentProcessInfo);
			//std::cout << "+++++++++++++++++++++++++++++\n";

			// calculate material response ******************************************************************
			noalias(N) = row(geom.ShapeFunctionsValues(), intp_id);
			Parameters.SetShapeFunctionsValues(N);
			mConstitutiveLawVector[intp_id]->CalculateMaterialResponseCauchy(Parameters);
			// calculate local contributions, transform them to global system and add
			if (LHSrequired) {
				noalias(BTD) = prod(trans(B), dV*D);
				noalias(Kloc) = prod(BTD, B);
			}
			if (RHSrequired) {
				noalias(Rloc) = -prod(trans(B), dV*mgeneralizedTemperatureStresses);
			}
			TransformToGlobalAndAdd(permutation_temp, iR, Kloc, Rloc, rLeftHandSideMatrix, rRightHandSideVector, LHSrequired, RHSrequired);
		}
	}

	// =====================================================================================
	//
	// Class ConvDiffInterfaceElement - Serialization
	//
	// =====================================================================================

	void ConvDiffInterfaceElement::save(Serializer& rSerializer) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
	}

	void ConvDiffInterfaceElement::load(Serializer& rSerializer)
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
	}

}
