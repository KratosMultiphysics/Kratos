//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#include "small_displacement_interface_element.hpp"
#include "multiscale_application_variables.h"

#include <string>
#include <iomanip>
#include "custom_utilities/math_helpers.h"

namespace Kratos
{

    // =====================================================================================
    //
    // Class SmallDisplacementInterfaceElement
    //
    // =====================================================================================

    SmallDisplacementInterfaceElement::SmallDisplacementInterfaceElement(IndexType NewId,
                                               GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
		, mInitialized(false)
    {
    }

    SmallDisplacementInterfaceElement::SmallDisplacementInterfaceElement(IndexType NewId,
                                               GeometryType::Pointer pGeometry,
                                               PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
		, mInitialized(false)
    {
    }

    SmallDisplacementInterfaceElement::~SmallDisplacementInterfaceElement()
    {
    }

    Element::Pointer SmallDisplacementInterfaceElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        GeometryType::Pointer newGeom( GetGeometry().Create(ThisNodes) );
        return Element::Pointer( new SmallDisplacementInterfaceElement(NewId, newGeom, pProperties ));
    }

    void SmallDisplacementInterfaceElement::Initialize()
    {
		/*
		Note:
		here we allocate 1 constitutive law for each 2 gauss points.
		This is because GiD geometry for interfaces are the same of the corresponding continuum element.
		So it has twice the integration points that an interface needs.
		*/
        if(mInitialized == false)
		{
			ConstitutiveLaw::Pointer& pLaw = GetProperties()[CONSTITUTIVE_LAW];
			const GeometryType::IntegrationPointsArrayType& integrationPoints = GetGeometry().IntegrationPoints();
			SizeType half_ngp = integrationPoints.size() / 2;

			mConstitutiveLawVector.clear();
			for(SizeType i = 0; i < half_ngp; i++)
			{
				ConstitutiveLaw::Pointer newCLaw = pLaw->Clone();
				newCLaw->InitializeMaterial(GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues(), i ));
				mConstitutiveLawVector.push_back(newCLaw);
			}

			mInitialized = true;
		}
    }

    void SmallDisplacementInterfaceElement::ResetConstitutiveLaw()
    {
        for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			mConstitutiveLawVector[i]->ResetMaterial(GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues(), i ) );
    }

    void SmallDisplacementInterfaceElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
		GeometryType & geom = this->GetGeometry();
		SizeType ndim = geom.WorkingSpaceDimension();
		SizeType ndofs = ndim * geom.size();

        if(rResult.size() != ndofs)
            rResult.resize(ndofs, false);

        for(SizeType i = 0; i < geom.size(); i++)
        {
            int index = i * ndim;
            NodeType & iNode = geom[i];

            rResult[index]     = iNode.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = iNode.GetDof(DISPLACEMENT_Y).EquationId();
            if(ndim == 3)
				rResult[index + 2] = iNode.GetDof(DISPLACEMENT_Z).EquationId();
        }
    }

    void SmallDisplacementInterfaceElement::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        GeometryType & geom = this->GetGeometry();
		SizeType ndim = geom.WorkingSpaceDimension();
		SizeType ndofs = ndim * geom.size();

		ElementalDofList.resize(0);
        ElementalDofList.reserve(ndofs);

        for (SizeType i = 0; i < geom.size(); i++)
        {
            NodeType & iNode = geom[i];

            ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Y));
            if(ndim == 3)
				ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Z));
        }
    }

    int SmallDisplacementInterfaceElement::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& geom = GetGeometry();

		SizeType ndim = geom.WorkingSpaceDimension();

        // verify that the variables are correctly initialized
        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered","");
        if(VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered","");
        if(ACCELERATION.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered","");
        if(ndim == 2)
			if(THICKNESS.Key() == 0)
				KRATOS_THROW_ERROR(std::invalid_argument,"THICKNESS has Key zero! (check if the application is correctly registered","");
        if(CONSTITUTIVE_LAW.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"CONSTITUTIVE_LAW has Key zero! (check if the application is correctly registered","");

        // verify that the dofs exist
        for(unsigned int i=0; i<geom.size(); i++)
        {
            if(geom[i].SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing variable DISPLACEMENT on node ",geom[i].Id());
            if(geom[i].HasDofFor(DISPLACEMENT_X) == false || geom[i].HasDofFor(DISPLACEMENT_Y) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ",GetGeometry()[i].Id());
			if(ndim == 3)
				if(geom[i].HasDofFor(DISPLACEMENT_Z) == false)
					KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ",GetGeometry()[i].Id());
        }

        // check properties
        if(this->pGetProperties() == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "Properties not provided for element ", this->Id());

        const PropertiesType & props = this->GetProperties();

		if(ndim == 2)
			if(props.Has(THICKNESS) == false)
				KRATOS_THROW_ERROR(std::invalid_argument, "THICKNESS not provided for element ", this->Id());

		if(props.Has(CONSTITUTIVE_LAW) == false)
			KRATOS_THROW_ERROR(std::invalid_argument, "CONSTITUTIVE_LAW not provided for element ", this->Id());

		const ConstitutiveLaw::Pointer& pLaw = props[CONSTITUTIVE_LAW];
		if(pLaw == NULL)
			KRATOS_THROW_ERROR(std::invalid_argument, "CONSTITUTIVE_LAW not provided for element ", this->Id());

		pLaw->Check(props, geom, rCurrentProcessInfo);

        return 0;

        KRATOS_CATCH("")
    }

    void SmallDisplacementInterfaceElement::CleanMemory()
    {
    }

    void SmallDisplacementInterfaceElement::GetValuesVector(Vector& values, int Step)
    {
        const GeometryType & geom = GetGeometry();

		SizeType dim = geom.WorkingSpaceDimension();
		SizeType ndofs = dim * geom.size();
        if(values.size() != ndofs)
            values.resize(ndofs,false);

        for (SizeType i = 0; i < geom.size(); i++)
        {
            const NodeType & iNode = geom[i];
            const array_1d<double,3>& disp = iNode.FastGetSolutionStepValue(DISPLACEMENT, Step);

            int index = i * dim;
            values[index]     = disp[0];
            values[index + 1] = disp[1];
            if(dim == 3)
				values[index + 2] = disp[2];
        }
    }

    void SmallDisplacementInterfaceElement::GetFirstDerivativesVector(Vector& values, int Step)
    {
        const GeometryType & geom = GetGeometry();

		SizeType dim = geom.WorkingSpaceDimension();
		SizeType ndofs = dim * geom.size();
        if(values.size() != ndofs)
            values.resize(ndofs,false);

        for (SizeType i = 0; i < geom.size(); i++)
        {
            const NodeType & iNode = geom[i];
            const array_1d<double,3>& vel = iNode.FastGetSolutionStepValue(VELOCITY, Step);

            int index = i * dim;
            values[index]     = vel[0];
            values[index + 1] = vel[1];
            if(dim == 3)
				values[index + 2] = vel[2];
        }
    }

    void SmallDisplacementInterfaceElement::GetSecondDerivativesVector(Vector& values, int Step)
    {
		const GeometryType & geom = GetGeometry();

		SizeType dim = geom.WorkingSpaceDimension();
		SizeType ndofs = dim * geom.size();
        if(values.size() != ndofs)
            values.resize(ndofs,false);

        for (SizeType i = 0; i < geom.size(); i++)
        {
            const NodeType & iNode = geom[i];
            const array_1d<double,3>& acc = iNode.FastGetSolutionStepValue(ACCELERATION, Step);

            int index = i * dim;
            values[index]     = acc[0];
            values[index + 1] = acc[1];
            if(dim == 3)
				values[index + 2] = acc[2];
        }
    }

    void SmallDisplacementInterfaceElement::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
    }

    void SmallDisplacementInterfaceElement::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
    }

    void SmallDisplacementInterfaceElement::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
		const PropertiesType& props = GetProperties();
		const GeometryType& geom = GetGeometry();
		const Matrix& N = geom.ShapeFunctionsValues();
		for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			mConstitutiveLawVector[i]->InitializeSolutionStep( props, geom, row( N, i ), CurrentProcessInfo );
    }

    void SmallDisplacementInterfaceElement::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
		const PropertiesType& props = GetProperties();
		const GeometryType& geom = GetGeometry();
		const Matrix& N = geom.ShapeFunctionsValues();
		for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			mConstitutiveLawVector[i]->FinalizeSolutionStep( props, geom, row( N, i ), CurrentProcessInfo );
    }

    void SmallDisplacementInterfaceElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        SizeType ndofs = GetGeometry().WorkingSpaceDimension() * GetGeometry().size();
        if((rMassMatrix.size1() != ndofs) || (rMassMatrix.size2() != ndofs))
            rMassMatrix.resize(ndofs, ndofs, false);

        noalias( rMassMatrix ) = ZeroMatrix(ndofs, ndofs);
    }

    void SmallDisplacementInterfaceElement::CalculateDampingMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
    {
		SizeType ndofs = GetGeometry().WorkingSpaceDimension() * GetGeometry().size();
        if((rDampMatrix.size1() != ndofs) || (rDampMatrix.size2() != ndofs))
            rDampMatrix.resize(ndofs, ndofs, false);

        noalias( rDampMatrix ) = ZeroMatrix(ndofs, ndofs);
    }

    void SmallDisplacementInterfaceElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                    VectorType& rRightHandSideVector,
                                                    ProcessInfo& rCurrentProcessInfo)
    {
        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
    }

    void SmallDisplacementInterfaceElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                      ProcessInfo& rCurrentProcessInfo)
    {
        Matrix dummy;
        CalculateAll(dummy, rRightHandSideVector, rCurrentProcessInfo, true, true);
    }

    // =====================================================================================
    //
    // Class SmallDisplacementInterfaceElement - Results on Gauss Points
    //
    // =====================================================================================

	void SmallDisplacementInterfaceElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
																		 std::vector<double>& rOutput,
																		 const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

	void SmallDisplacementInterfaceElement::CalculateOnIntegrationPoints(const Variable< array_1d< double, 3 > >& rVariable,
																		 std::vector< array_1d<double, 3 > >& rOutput,
																		 const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

	void SmallDisplacementInterfaceElement::CalculateOnIntegrationPoints(const Variable< Vector >& rVariable,
																		 std::vector< Vector >& rOutput,
																		 const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

	void SmallDisplacementInterfaceElement::CalculateOnIntegrationPoints(const Variable< Matrix >& rVariable,
																		 std::vector< Matrix >& rOutput,
																		 const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

    void SmallDisplacementInterfaceElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                                           std::vector<double>& rValues,
                                                           const ProcessInfo& rCurrentProcessInfo)
    {
		SizeType num_gp = GetGeometry().IntegrationPoints().size();
		if(rValues.size() != num_gp) rValues.resize(num_gp);
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
			for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			{
				SizeType j = i+half_num_gp;
				res = 0.0;
				mConstitutiveLawVector[i]->GetValue(rVariable, res);
				rValues[i] = res;
				rValues[j] = res;
			}
		}
    }

    void SmallDisplacementInterfaceElement::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                           std::vector<Vector>& rValues,
                                                           const ProcessInfo& rCurrentProcessInfo)
    {
		SizeType num_gp = GetGeometry().IntegrationPoints().size();
		if(rValues.size() != num_gp) rValues.resize(num_gp);
		SizeType half_num_gp = num_gp / 2;
		SizeType ndim = GetGeometry().WorkingSpaceDimension();

		if(rVariable == INTERFACE_DISPLACEMENT_JUMP || rVariable == INTERFACE_TRACTION)
		{
			std::vector<Vector> res;
			if(rVariable == INTERFACE_DISPLACEMENT_JUMP)
				CalculateStrain(res,rCurrentProcessInfo);
			else
				CalculateStress(res,rCurrentProcessInfo);
			if(ndim == 2)
			{
				rValues[0] = rValues[3] = res[0];
				rValues[1] = rValues[2] = res[1];
			}
			else
			{
				for(SizeType i = 0; i < res.size(); i++)
				{
					SizeType j = i+half_num_gp;
					rValues[i] = res[i];
					rValues[j] = res[i];
				}
			}
		}
		else
		{
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
    }

    void SmallDisplacementInterfaceElement::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
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

    void SmallDisplacementInterfaceElement::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
                                                           std::vector<array_1d<double,3> >& rValues,
                                                           const ProcessInfo& rCurrentProcessInfo)
    {
		SizeType num_gp = GetGeometry().IntegrationPoints().size();
		if(rValues.size() != num_gp) rValues.resize(num_gp);
		SizeType half_num_gp = num_gp / 2;
		SizeType ndim = GetGeometry().WorkingSpaceDimension();

		array_1d<double,3> res;
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

    void SmallDisplacementInterfaceElement::GetValueOnIntegrationPoints(const Variable<array_1d<double,6> >& rVariable,
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
    // Class SmallDisplacementInterfaceElement - Private methods
    //
    // =====================================================================================

    void SmallDisplacementInterfaceElement::DecimalCorrection(Vector& a)
    {
        /*double norm = norm_2(a);
        double tolerance = std::max(norm * 1.0E-12, 1.0E-12);
        for(SizeType i = 0; i < a.size(); i++)
            if(std::abs(a(i)) < tolerance)
                a(i) = 0.0;*/
    }

	void SmallDisplacementInterfaceElement::CalculatePermutation(InterfaceIndexPermutation& p)
	{
		p.HasPermutation = false;
		p.Permutation.clear();
		GeometryType& geom = GetGeometry();
		if(geom.WorkingSpaceDimension() == 2)
		{
			// handle case of quad 2d 4n (this is the only one supported by now)
			if(geom.size() == 4)
			{
                p.HasPermutation = true;
				// 1 2   3 4   7 8   5 6
				p.Permutation.resize(8);
				p.Permutation[0] = 0;
				p.Permutation[1] = 1;
				p.Permutation[2] = 2;
				p.Permutation[3] = 3;
				p.Permutation[4] = 6;
				p.Permutation[5] = 7;
				p.Permutation[6] = 4;
				p.Permutation[7] = 5;
			}
		}
	}

	void SmallDisplacementInterfaceElement::CalculateDeltaPosition(Matrix& rDeltaPosition)
	{
		GeometryType& geom = GetGeometry();
		const unsigned int number_of_nodes = geom.PointsNumber();
		unsigned int dimension = geom.WorkingSpaceDimension();

		rDeltaPosition = zero_matrix<double>( number_of_nodes , dimension);

		for ( unsigned int i = 0; i < number_of_nodes; i++ )
		{
			const NodeType& iNode = geom[i];
			rDeltaPosition(i, 0) = iNode.X() - iNode.X0();
			rDeltaPosition(i, 1) = iNode.Y() - iNode.Y0();
			if(dimension == 3)
				rDeltaPosition(i, 2) = iNode.Z() - iNode.Z0();
		}
	}

	void SmallDisplacementInterfaceElement::CalculateJacobianAndTransformationMatrix(const SizeType pointID,
												                                     Matrix& delta_position,
		                                                                             Matrix& jacobian,
												                                     double& J,
												                                     Matrix& iR)
	{
		GeometryType& geom = GetGeometry();
		SizeType ndim = geom.WorkingSpaceDimension();

		geom.Jacobian(jacobian, pointID, GetIntegrationMethod(), delta_position);

		if(ndim == 2)
		{
			array_1d<double, 2> vx;
			vx[0] = jacobian(0, 0);
			vx[1] = jacobian(1, 0);
			J = std::sqrt(vx[0]*vx[0] + vx[1]*vx[1]);
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
			for(SizeType i = 0; i < 3; i++)
			{
				iR(0, i) = vx[i];
				iR(1, i) = vy[i];
				iR(2, i) = vz[i];
			}
		}
	}

	double SmallDisplacementInterfaceElement::CalculateIntegrationWeight(double J,
		                                                                 double iw)
	{
		// Note: here we multiply the integration weight by 2,
		// because we're using the geometry of a continuum element
		// but we're using only half of the integration points.
		double dV = J * iw * 2.0;
		if(GetGeometry().WorkingSpaceDimension() == 2)
			dV *= GetProperties()[THICKNESS];
		return dV;
	}

	void SmallDisplacementInterfaceElement::CalculateLocalDisplacementVector(const InterfaceIndexPermutation& P,
		                                                                     const Matrix& R,
		                                                                     const Vector& UG,
												                             Vector& UL)
	{
		GeometryType& geom = GetGeometry();
		SizeType nnodes = geom.size();
		SizeType ndim = geom.WorkingSpaceDimension();
		if(P.HasPermutation)
		{
			for(SizeType inode = 0; inode < nnodes; inode++)
			{
				SizeType pos = inode*ndim;
				for(SizeType i = 0; i < ndim; i++)
				{
					double temp = 0.0;
					for(SizeType j = 0; j < ndim; j++)
					{
						temp += R(i, j) * UG( P.Permutation[ pos+j ] );
					}
					UL(pos+i) = temp;
				}
			}
		}
		else
		{
			for(SizeType inode = 0; inode < nnodes; inode++)
			{
				SizeType pos = inode*ndim;
				for(SizeType i = 0; i < ndim; i++)
				{
					double temp = 0.0;
					for(SizeType j = 0; j < ndim; j++)
					{
						temp += R(i, j) * UG(pos+j);
					}
					UL(pos+i) = temp;
				}
			}
		}
	}

	void SmallDisplacementInterfaceElement::TransformToGlobalAndAdd(const InterfaceIndexPermutation& P,
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

		if(P.HasPermutation)
		{
			for(SizeType node_i = 0; node_i < nnodes; node_i++)
			{
				SizeType pos_i = node_i*ndim;

				if(RHSrequired)
				{
					for(SizeType i = 0; i < ndim; i++)
					{
						double temp_RHS = 0.0;
						for(SizeType j = 0; j < ndim; j++)
						{
							temp_RHS += R(j, i) * RHS_local( pos_i+j );
						}
						RHS_global( P.Permutation[ pos_i+i ] ) += temp_RHS;
					}
				}

				if(LHSrequired)
				{
					Matrix RTK(ndim, ndim);
					for(SizeType node_j = 0; node_j < nnodes; node_j++)
					{
						SizeType pos_j = node_j*ndim;

						RTK.clear();
						for(SizeType i = 0; i < ndim; i++)
						{
							for(SizeType j = 0; j < ndim; j++)
							{
								double temp_K = 0.0;
								for(SizeType k = 0; k < ndim; k++)
								{
									temp_K += R(k, i) * LHS_local( pos_i+k, pos_j+j );
								}
								RTK(i, j) = temp_K;
							}
						}
						for(SizeType i = 0; i < ndim; i++)
						{
							for(SizeType j = 0; j < ndim; j++)
							{
								double temp_K = 0.0;
								for(SizeType k = 0; k < ndim; k++)
								{
									temp_K += RTK(i, k) * R(k, j);
								}
								LHS_global( P.Permutation[ pos_i+i ], P.Permutation[ pos_j+j ] ) += temp_K;
							}
						}
					}
				}
			}
		}
		else
		{
			for(SizeType node_i = 0; node_i < nnodes; node_i++)
			{
				SizeType pos_i = node_i*ndim;

				if(RHSrequired)
				{
					for(SizeType i = 0; i < ndim; i++)
					{
						double temp_RHS = 0.0;
						for(SizeType j = 0; j < ndim; j++)
						{
							temp_RHS += R(j, i) * RHS_local( pos_i+j );
						}
						RHS_global( pos_i+i ) += temp_RHS;
					}
				}

				if(LHSrequired)
				{
					Matrix RTK(ndim, ndim);
					for(SizeType node_j = 0; node_j < nnodes; node_j++)
					{
						SizeType pos_j = node_j*ndim;

						RTK.clear();
						for(SizeType i = 0; i < ndim; i++)
						{
							for(SizeType j = 0; j < ndim; j++)
							{
								double temp_K = 0.0;
								for(SizeType k = 0; k < ndim; k++)
								{
									temp_K += R(k, i) * LHS_local( pos_i+k, pos_j+j );
								}
								RTK(i, j) = temp_K;
							}
						}
						for(SizeType i = 0; i < ndim; i++)
						{
							for(SizeType j = 0; j < ndim; j++)
							{
								double temp_K = 0.0;
								for(SizeType k = 0; k < ndim; k++)
								{
									temp_K += RTK(i, k) * R(k, j);
								}
								LHS_global( pos_i+i, pos_j+j ) += temp_K;
							}
						}
					}
				}
			}
		}
	}

	void SmallDisplacementInterfaceElement::CalculateBMatrix(const SizeType pointID,
		                                                     Matrix& B)
	{
		GeometryType& geom = GetGeometry();
		const Matrix& shapes = geom.ShapeFunctionsValues();

		SizeType ndim = geom.WorkingSpaceDimension();
		SizeType nnodes = geom.size();
		SizeType half_nnodes = nnodes / 2;

		for(SizeType k = 0; k < half_nnodes; k++)
		{
			double Nk = shapes(pointID, k);
			SizeType pos1 = k * ndim;
			SizeType pos2 = pos1 + half_nnodes*ndim;
			for(SizeType i = 0; i < ndim; i++)
			{
				B(i, pos1+i) = -Nk;
				B(i, pos2+i) =  Nk;
			}
		}
	}

	void SmallDisplacementInterfaceElement::CalculateGeneralizedStrains(const SizeType pointID,
		                                                                const Matrix& B,
		                                                                const Vector& U,
											                            Vector& generalizedStrains)
	{
		noalias( generalizedStrains ) = prod( B, U );
	}

    void SmallDisplacementInterfaceElement::CalculateAll(MatrixType& rLeftHandSideMatrix,
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

		bool use_reduced_integration = false;
		if(props.Has(INTERFACE_REDUCED_INTEGRATION))
			use_reduced_integration = (props[INTERFACE_REDUCED_INTEGRATION] != 0);

		// resize the LHS matrix
		if(LHSrequired) {
			if(rLeftHandSideMatrix.size1() != ndofs || rLeftHandSideMatrix.size2() != ndofs)
				rLeftHandSideMatrix.resize(ndofs, ndofs, false);
			noalias( rLeftHandSideMatrix ) = ZeroMatrix(ndofs, ndofs);
		}

		// resize the RHS vector
		if(RHSrequired) {
			if(rRightHandSideVector.size() != ndofs)
				rRightHandSideVector.resize(ndofs, false);
			noalias( rRightHandSideVector ) = ZeroVector(ndofs);
		}

		// global and local displacement vectors
		Vector globalDisplacements(ndofs);
		Vector localDisplacements(ndofs);
		GetValuesVector(globalDisplacements);

		// delta position for the jacobian computation with respect to the reference configuration
		Matrix delta_position;
		CalculateDeltaPosition(delta_position);

		// jacobian matrix and its determinant
		Matrix jacobian(ndim, ndim-1);
		double J;

		// strain-displacement matrix
		Matrix B(ndim, ndofs, 0.0);

		// material point calculation data
		Matrix D(ndim, ndim);
		Vector generalizedStrains(ndim);
		Vector generalizedStresses(ndim);

		// transformation matrix
		Matrix iR(ndim, ndim);

		// permutation data
		InterfaceIndexPermutation permutation;
		CalculatePermutation(permutation);

		// LHS and RHS in local coordinate system
		Matrix Kloc(ndofs, ndofs);
		Vector Rloc(ndofs);

		// auxiliary data to avoid extra memory allocations
		Matrix BTD(ndofs, ndim);

		// initialize material parameters
		ConstitutiveLaw::Parameters Parameters(GetGeometry(), GetProperties(),rCurrentProcessInfo);
		Flags &ConstitutiveLawOptions = Parameters.GetOptions();
		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
		if(LHSrequired) ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
		Parameters.SetStrainVector(generalizedStrains);
		Parameters.SetStressVector(generalizedStresses);
		Parameters.SetConstitutiveMatrix(D);
		Matrix F( IdentityMatrix(ndim, ndim) );
		double detF(1.0);
		Parameters.SetDeformationGradientF( F );
		Parameters.SetDeterminantF( detF );
		Matrix DN_DX(nnodes, ndim, 0.0);
		Vector N(nnodes);
		Parameters.SetShapeFunctionsDerivatives(DN_DX);

		// loop over the integration points.
		// note: loop only the first half of them, then multiply the integration weight by 2
		const GeometryType::IntegrationPointsArrayType& integrationPoints = geom.IntegrationPoints();

		Matrix B0(ndim, ndofs, 0.0);
		if(use_reduced_integration)
		{
			double V0 = 0.0;
			for(SizeType intp_id = 0; intp_id < integrationPoints.size()/2; intp_id++)
			{
				// the current integration point
				const GeometryType::IntegrationPointType& ip = integrationPoints[intp_id];
				// jacobianobian and local transformation
				CalculateJacobianAndTransformationMatrix(intp_id, delta_position, jacobian, J, iR);
				// integration weight
				double dV = CalculateIntegrationWeight(J, ip.Weight());
				// local displacement vector
				CalculateLocalDisplacementVector(permutation, iR, globalDisplacements, localDisplacements);
				// strain-displacement matrix
				CalculateBMatrix(intp_id, B);
				// accumulate
				B0 += B*dV;
				V0 += dV;
			}
			B0 /= V0;
		}

		for(SizeType intp_id = 0; intp_id < integrationPoints.size()/2; intp_id++)
		{
			// the current integration point
			const GeometryType::IntegrationPointType& ip = integrationPoints[intp_id];

			// jacobianobian and local transformation
			CalculateJacobianAndTransformationMatrix(intp_id, delta_position, jacobian, J, iR);

			// integration weight
			double dV = CalculateIntegrationWeight(J, ip.Weight());

			// local displacement vector
			CalculateLocalDisplacementVector(permutation, iR, globalDisplacements, localDisplacements);

			// strain-displacement matrix
			CalculateBMatrix(intp_id, B);
			if(use_reduced_integration) {
				noalias(B) = B0;
			}

			// calculate generalized strains
			CalculateGeneralizedStrains(intp_id, B, localDisplacements, generalizedStrains);
			DecimalCorrection(generalizedStrains);

			// calculate material response
			noalias( N ) = row( geom.ShapeFunctionsValues(), intp_id );
			Parameters.SetShapeFunctionsValues( N );
			mConstitutiveLawVector[intp_id]->CalculateMaterialResponseCauchy(Parameters);

			// calculate local contributions, transform them to global system and add
			if(LHSrequired) {
				noalias( BTD )  =  prod( trans( B ), dV*D );
				noalias( Kloc ) =  prod( BTD, B );
			}
			if(RHSrequired) {
				noalias( Rloc ) = -prod( trans( B ), dV*generalizedStresses );
			}
			TransformToGlobalAndAdd(permutation, iR, Kloc, Rloc,
									rLeftHandSideMatrix, rRightHandSideVector,
									LHSrequired, RHSrequired);
		}
	}

	void SmallDisplacementInterfaceElement::CalculateStrain(std::vector<Vector>& strain, const ProcessInfo& rCurrentProcessInfo)
	{
		// init the output container. it will contain ngauss/2 output variables
		strain.clear();

		// general data
		PropertiesType props = GetProperties();
		const GeometryType& geom = GetGeometry();
		SizeType ndim = geom.WorkingSpaceDimension();
		SizeType nnodes = geom.size();
		SizeType ndofs = ndim * nnodes;

		bool use_reduced_integration = false;
		if(props.Has(INTERFACE_REDUCED_INTEGRATION))
			use_reduced_integration = (props[INTERFACE_REDUCED_INTEGRATION] != 0);

		// global and local displacement vectors
		Vector globalDisplacements(ndofs);
		Vector localDisplacements(ndofs);
		GetValuesVector(globalDisplacements);

		// delta position for the jacobian computation with respect to the reference configuration
		Matrix delta_position;
		CalculateDeltaPosition(delta_position);

		// jacobian matrix and its determinant
		Matrix jacobian(ndim, ndim-1);
		double J;

		// strain-displacement matrix
		Matrix B(ndim, ndofs, 0.0);

		// material point calculation data
		Vector generalizedStrains(ndim);

		// transformation matrix
		Matrix iR(ndim, ndim);

		// permutation data
		InterfaceIndexPermutation permutation;
		CalculatePermutation(permutation);

		Vector vec3(3);

		// loop over the integration points.
		// note: loop only the first half of them, then multiply the integration weight by 2
		const GeometryType::IntegrationPointsArrayType& integrationPoints = geom.IntegrationPoints();

		Matrix B0(ndim, ndofs, 0.0);
		if(use_reduced_integration)
		{
			double V0 = 0.0;
			for(SizeType intp_id = 0; intp_id < integrationPoints.size()/2; intp_id++)
			{
				// the current integration point
				const GeometryType::IntegrationPointType& ip = integrationPoints[intp_id];
				// jacobianobian and local transformation
				CalculateJacobianAndTransformationMatrix(intp_id, delta_position, jacobian, J, iR);
				// integration weight
				double dV = CalculateIntegrationWeight(J, ip.Weight());
				// local displacement vector
				CalculateLocalDisplacementVector(permutation, iR, globalDisplacements, localDisplacements);
				// strain-displacement matrix
				CalculateBMatrix(intp_id, B);
				// accumulate
				B0 += B*dV;
				V0 += dV;
			}
			B0 /= V0;
		}

		for(SizeType intp_id = 0; intp_id < integrationPoints.size()/2; intp_id++)
		{
			// the current integration point
			const GeometryType::IntegrationPointType& ip = integrationPoints[intp_id];

			// jacobianobian and local transformation
			CalculateJacobianAndTransformationMatrix(intp_id, delta_position, jacobian, J, iR);

			// local displacement vector
			CalculateLocalDisplacementVector(permutation, iR, globalDisplacements, localDisplacements);

			// strain-displacement matrix
			CalculateBMatrix(intp_id, B);
			if(use_reduced_integration) {
				//noalias(B) = B0;
				for(unsigned int ii=0; ii < B.size2(); ii++)
					B(1,ii)=B0(1,ii);
			}

			// calculate generalized strains
			CalculateGeneralizedStrains(intp_id, B, localDisplacements, generalizedStrains);
			DecimalCorrection(generalizedStrains);

			// store output
			vec3(0) = generalizedStrains(0);
			vec3(1) = generalizedStrains(1);
			if(ndim == 3)
				vec3(2) = generalizedStrains(2);
			else
				vec3(2) = 0.0;
			strain.push_back(vec3);
		}
	}

	void SmallDisplacementInterfaceElement::CalculateStress(std::vector<Vector>& stress, const ProcessInfo& rCurrentProcessInfo)
	{
		// init the output container. it will contain ngauss/2 output variables
		stress.clear();

		// general data
		PropertiesType props = GetProperties();
		const GeometryType& geom = GetGeometry();
		SizeType ndim = geom.WorkingSpaceDimension();
		SizeType nnodes = geom.size();
		SizeType ndofs = ndim * nnodes;

		bool use_reduced_integration = false;
		if(props.Has(INTERFACE_REDUCED_INTEGRATION))
			use_reduced_integration = (props[INTERFACE_REDUCED_INTEGRATION] != 0);

		// global and local displacement vectors
		Vector globalDisplacements(ndofs);
		Vector localDisplacements(ndofs);
		GetValuesVector(globalDisplacements);

		// delta position for the jacobian computation
		// with respect to the reference configuration
		Matrix delta_position;
		CalculateDeltaPosition(delta_position);

		// jacobian matrix and its determinant
		Matrix jacobian(ndim, ndim-1);
		double J;

		// strain-displacement matrix
		Matrix B(ndim, ndofs, 0.0);

		// material point calculation data
		Matrix D(ndim, ndim);
		Vector generalizedStrains(ndim);
		Vector generalizedStresses(ndim);

		// transformation matrix
		Matrix iR(ndim, ndim);

		// permutation data
		InterfaceIndexPermutation permutation;
		CalculatePermutation(permutation);

		// initialize material parameters
		ConstitutiveLaw::Parameters Parameters(GetGeometry(), GetProperties(),rCurrentProcessInfo);
		Flags &ConstitutiveLawOptions = Parameters.GetOptions();
		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
		Parameters.SetStrainVector(generalizedStrains);
		Parameters.SetStressVector(generalizedStresses);
		Parameters.SetConstitutiveMatrix(D);
		Matrix F( IdentityMatrix(ndim, ndim) );
		double detF(1.0);
		Parameters.SetDeformationGradientF( F );
		Parameters.SetDeterminantF( detF );
		Matrix DN_DX(nnodes, ndim, 0.0);
		Vector N(nnodes);
		Parameters.SetShapeFunctionsDerivatives(DN_DX);

		Vector vec3(3);

		// loop over the integration points.
		// note: loop only the first half of them, then multiply the integration weight by 2
		const GeometryType::IntegrationPointsArrayType& integrationPoints = geom.IntegrationPoints();

		Matrix B0(ndim, ndofs, 0.0);
		if(use_reduced_integration)
		{
			double V0 = 0.0;
			for(SizeType intp_id = 0; intp_id < integrationPoints.size()/2; intp_id++)
			{
				// the current integration point
				const GeometryType::IntegrationPointType& ip = integrationPoints[intp_id];
				// jacobianobian and local transformation
				CalculateJacobianAndTransformationMatrix(intp_id, delta_position, jacobian, J, iR);
				// integration weight
				double dV = CalculateIntegrationWeight(J, ip.Weight());
				// local displacement vector
				CalculateLocalDisplacementVector(permutation, iR, globalDisplacements, localDisplacements);
				// strain-displacement matrix
				CalculateBMatrix(intp_id, B);
				// accumulate
				B0 += B*dV;
				V0 += dV;
			}
			B0 /= V0;
		}

		for(SizeType intp_id = 0; intp_id < integrationPoints.size()/2; intp_id++)
		{
			// the current integration point
			const GeometryType::IntegrationPointType& ip = integrationPoints[intp_id];

			// jacobianobian and local transformation
			CalculateJacobianAndTransformationMatrix(intp_id, delta_position, jacobian, J, iR);

			// local displacement vector
			CalculateLocalDisplacementVector(permutation, iR, globalDisplacements, localDisplacements);

			// strain-displacement matrix
			CalculateBMatrix(intp_id, B);
			if(use_reduced_integration) {
				//noalias(B) = B0;
				for(unsigned int ii=0; ii < B.size2(); ii++)
					B(1,ii)=B0(1,ii);
			}

			// calculate generalized strains
			CalculateGeneralizedStrains(intp_id, B, localDisplacements, generalizedStrains);
			DecimalCorrection(generalizedStrains);

			// calculate material response
			noalias( N ) = row( geom.ShapeFunctionsValues(), intp_id );
			Parameters.SetShapeFunctionsValues( N );
			mConstitutiveLawVector[intp_id]->CalculateMaterialResponseCauchy(Parameters);

			// store output
			vec3(0) = generalizedStresses(0);
			vec3(1) = generalizedStresses(1);
			if(ndim == 3)
				vec3(2) = generalizedStresses(2);
			else
				vec3(2) = 0.0;
			stress.push_back(vec3);
		}
	}

	// =====================================================================================
	//
	// Class SmallDisplacementInterfaceElement - Serialization
	//
	// =====================================================================================

    void SmallDisplacementInterfaceElement::save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
    }

    void SmallDisplacementInterfaceElement::load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
    }

}
