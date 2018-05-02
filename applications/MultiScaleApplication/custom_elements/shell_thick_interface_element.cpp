//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#include "shell_thick_interface_element.hpp"
#include "multiscale_application_variables.h"

#include <string>
#include <iomanip>
#include "custom_utilities/math_helpers.h"

namespace Kratos
{

    // =====================================================================================
    //
    // Class ShellThickInterfaceElement
    //
    // =====================================================================================

    ShellThickInterfaceElement::ShellThickInterfaceElement(IndexType NewId,
                                               GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
		, mInitialized(false)
    {
    }

    ShellThickInterfaceElement::ShellThickInterfaceElement(IndexType NewId,
                                               GeometryType::Pointer pGeometry,
                                               PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
		, mInitialized(false)
    {
    }

    ShellThickInterfaceElement::~ShellThickInterfaceElement()
    {
    }

    Element::Pointer ShellThickInterfaceElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        GeometryType::Pointer newGeom( GetGeometry().Create(ThisNodes) );
        return Element::Pointer( new ShellThickInterfaceElement(NewId, newGeom, pProperties ));
    }

    void ShellThickInterfaceElement::Initialize()
    {
		/*
		Note:
		here we allocate 1 constitutive law for each 2 gauss points.
		This is because GiD geometry for interfaces are the same of the corresponding continuum element.
		So it has twice the integration points that an interface needs.
		*/
        if(mInitialized == false)
		{
			/*ConstitutiveLaw::Pointer& pLaw = GetProperties()[CONSTITUTIVE_LAW];
			const GeometryType::IntegrationPointsArrayType& integrationPoints = GetGeometry().IntegrationPoints();
			SizeType half_ngp = integrationPoints.size() / 2;

			mConstitutiveLawVector.clear();
			for(SizeType i = 0; i < half_ngp; i++)
			{
				ConstitutiveLaw::Pointer newCLaw = pLaw->Clone();
				newCLaw->InitializeMaterial(GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues(), i ));
				mConstitutiveLawVector.push_back(newCLaw);
			}*/

			mInitialized = true;
		}
    }

    void ShellThickInterfaceElement::ResetConstitutiveLaw()
    {
        /*for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			mConstitutiveLawVector[i]->ResetMaterial(GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues(), i ) );*/
    }

    void ShellThickInterfaceElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
		GeometryType & geom = this->GetGeometry();
		SizeType ndim = 6;
		SizeType ndofs = ndim * geom.size();

        if(rResult.size() != ndofs)
            rResult.resize(ndofs, false);

        for(SizeType i = 0; i < geom.size(); i++)
        {
            int index = i * ndim;
            NodeType & iNode = geom[i];

            rResult[index    ] = iNode.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = iNode.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = iNode.GetDof(DISPLACEMENT_Z).EquationId();
			rResult[index + 3] = iNode.GetDof(ROTATION_X).EquationId();
            rResult[index + 4] = iNode.GetDof(ROTATION_Y).EquationId();
            rResult[index + 5] = iNode.GetDof(ROTATION_Z).EquationId();
        }
    }

    void ShellThickInterfaceElement::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        GeometryType & geom = this->GetGeometry();
		SizeType ndim = 6;
		SizeType ndofs = ndim * geom.size();

		ElementalDofList.resize(0);
        ElementalDofList.reserve(ndofs);

        for (SizeType i = 0; i < geom.size(); i++)
        {
            NodeType & iNode = geom[i];

            ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Y));
            ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Z));
			ElementalDofList.push_back(iNode.pGetDof(ROTATION_X));
            ElementalDofList.push_back(iNode.pGetDof(ROTATION_Y));
            ElementalDofList.push_back(iNode.pGetDof(ROTATION_Z));
        }
    }

    int ShellThickInterfaceElement::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& geom = GetGeometry();

		SizeType ndim = geom.WorkingSpaceDimension();

        // verify that the variables are correctly initialized
        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered","");
		if(ROTATION.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"ROTATION has Key zero! (check if the application is correctly registered","");
        if(VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered","");
        if(ACCELERATION.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered","");
        if(THICKNESS.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument,"THICKNESS has Key zero! (check if the application is correctly registered","");
        if(CONSTITUTIVE_LAW.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"CONSTITUTIVE_LAW has Key zero! (check if the application is correctly registered","");

        // verify that the dofs exist
        for(unsigned int i=0; i<geom.size(); i++)
        {
            if(geom[i].SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing variable DISPLACEMENT on node ",geom[i].Id());
            if(geom[i].HasDofFor(DISPLACEMENT_X) == false || geom[i].HasDofFor(DISPLACEMENT_Y) == false || geom[i].HasDofFor(DISPLACEMENT_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ",GetGeometry()[i].Id());
			if(geom[i].SolutionStepsDataHas(ROTATION) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing variable ROTATION on node ",geom[i].Id());
            if(geom[i].HasDofFor(ROTATION_X) == false || geom[i].HasDofFor(ROTATION_Y) == false || geom[i].HasDofFor(ROTATION_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable ROTATION on node ",GetGeometry()[i].Id());
        }

        // check properties
        if(this->pGetProperties() == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "Properties not provided for element ", this->Id());

        const PropertiesType & props = this->GetProperties();

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

    void ShellThickInterfaceElement::CleanMemory()
    {
    }

    void ShellThickInterfaceElement::GetValuesVector(Vector& values, int Step)
    {
        const GeometryType & geom = GetGeometry();

		SizeType ndim = 6;
		SizeType ndofs = ndim * geom.size();
        if(values.size() != ndofs)
            values.resize(ndofs,false);

        for (SizeType i = 0; i < geom.size(); i++)
        {
            const NodeType & iNode = geom[i];
            const array_1d<double,3>& disp = iNode.FastGetSolutionStepValue(DISPLACEMENT, Step);
			const array_1d<double,3>& roto = iNode.FastGetSolutionStepValue(ROTATION, Step);

            int index = i * ndim;
            values[index    ] = disp[0];
            values[index + 1] = disp[1];
            values[index + 2] = disp[2];
			values[index + 3] = roto[0];
            values[index + 4] = roto[1];
            values[index + 5] = roto[2];
        }
    }

    void ShellThickInterfaceElement::GetFirstDerivativesVector(Vector& values, int Step)
    {
        const GeometryType & geom = GetGeometry();

		SizeType ndim = 6;
		SizeType ndofs = ndim * geom.size();
        if(values.size() != ndofs)
            values.resize(ndofs,false);

        for (SizeType i = 0; i < geom.size(); i++)
        {
            const NodeType & iNode = geom[i];
            const array_1d<double,3>& vel = iNode.FastGetSolutionStepValue(VELOCITY, Step);

            int index = i * ndim;
            values[index    ] = vel[0];
            values[index + 1] = vel[1];
            values[index + 2] = vel[2];
			values[index + 3] = 0.0;
            values[index + 4] = 0.0;
            values[index + 5] = 0.0;
        }
    }

    void ShellThickInterfaceElement::GetSecondDerivativesVector(Vector& values, int Step)
    {
		const GeometryType & geom = GetGeometry();

		SizeType ndim = 6;
		SizeType ndofs = ndim * geom.size();
        if(values.size() != ndofs)
            values.resize(ndofs,false);

        for (SizeType i = 0; i < geom.size(); i++)
        {
            const NodeType & iNode = geom[i];
            const array_1d<double,3>& acc = iNode.FastGetSolutionStepValue(ACCELERATION, Step);

            int index = i * ndim;
            values[index    ] = acc[0];
            values[index + 1] = acc[1];
            values[index + 2] = acc[2];
			values[index + 3] = 0.0;
            values[index + 4] = 0.0;
            values[index + 5] = 0.0;
        }
    }

    void ShellThickInterfaceElement::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
    }

    void ShellThickInterfaceElement::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
    }

    void ShellThickInterfaceElement::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
		/*const PropertiesType& props = GetProperties();
		const GeometryType& geom = GetGeometry();
		const Matrix& N = geom.ShapeFunctionsValues();
		for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			mConstitutiveLawVector[i]->InitializeSolutionStep( props, geom, row( N, i ), CurrentProcessInfo );*/
    }

    void ShellThickInterfaceElement::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
		/*const PropertiesType& props = GetProperties();
		const GeometryType& geom = GetGeometry();
		const Matrix& N = geom.ShapeFunctionsValues();
		for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
			mConstitutiveLawVector[i]->FinalizeSolutionStep( props, geom, row( N, i ), CurrentProcessInfo );*/
    }

    void ShellThickInterfaceElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        SizeType ndofs = GetGeometry().WorkingSpaceDimension() * GetGeometry().size();
        if((rMassMatrix.size1() != ndofs) || (rMassMatrix.size2() != ndofs))
            rMassMatrix.resize(ndofs, ndofs, false);

        noalias( rMassMatrix ) = ZeroMatrix(ndofs, ndofs);
    }

    void ShellThickInterfaceElement::CalculateDampingMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
    {
		SizeType ndofs = GetGeometry().WorkingSpaceDimension() * GetGeometry().size();
        if((rDampMatrix.size1() != ndofs) || (rDampMatrix.size2() != ndofs))
            rDampMatrix.resize(ndofs, ndofs, false);

        noalias( rDampMatrix ) = ZeroMatrix(ndofs, ndofs);
    }

    void ShellThickInterfaceElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                    VectorType& rRightHandSideVector,
                                                    ProcessInfo& rCurrentProcessInfo)
    {
        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
    }

    void ShellThickInterfaceElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                      ProcessInfo& rCurrentProcessInfo)
    {
        Matrix dummy;
        CalculateAll(dummy, rRightHandSideVector, rCurrentProcessInfo, true, true);
    }

    // =====================================================================================
    //
    // Class ShellThickInterfaceElement - Results on Gauss Points
    //
    // =====================================================================================

	void ShellThickInterfaceElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
																		 std::vector<double>& rOutput,
																		 const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

	void ShellThickInterfaceElement::CalculateOnIntegrationPoints(const Variable< array_1d< double, 3 > >& rVariable,
																		 std::vector< array_1d<double, 3 > >& rOutput,
																		 const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

	void ShellThickInterfaceElement::CalculateOnIntegrationPoints(const Variable< Vector >& rVariable,
																		 std::vector< Vector >& rOutput,
																		 const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

	void ShellThickInterfaceElement::CalculateOnIntegrationPoints(const Variable< Matrix >& rVariable,
																		 std::vector< Matrix >& rOutput,
																		 const ProcessInfo& rCurrentProcessInfo)
	{
		GetValueOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
	}

    void ShellThickInterfaceElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                                           std::vector<double>& rValues,
                                                           const ProcessInfo& rCurrentProcessInfo)
    {
		SizeType num_gp = GetGeometry().IntegrationPoints().size();
		if(rValues.size() != num_gp) rValues.resize(num_gp);
		SizeType half_num_gp = num_gp / 2;

		double res(0.0);

		// TODO: this is a workaround for the 4-node geometry (only one supported for shell interface)
		// generalize it asap
		res = 0.0;
		//mConstitutiveLawVector[0]->GetValue(rVariable, res);
		rValues[0] = rValues[3] = res;
		res = 0.0;
		//mConstitutiveLawVector[1]->GetValue(rVariable, res);
		rValues[1] = rValues[2] = res;
    }

    void ShellThickInterfaceElement::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                           std::vector<Vector>& rValues,
                                                           const ProcessInfo& rCurrentProcessInfo)
    {
		SizeType num_gp = GetGeometry().IntegrationPoints().size();
		if(rValues.size() != num_gp) rValues.resize(num_gp);
		SizeType half_num_gp = num_gp / 2;

		if(rVariable == INTERFACE_DISPLACEMENT_JUMP || rVariable == INTERFACE_TRACTION)
		{
			std::vector<Vector> res;
			if(rVariable == INTERFACE_DISPLACEMENT_JUMP)
				CalculateStrain(res,rCurrentProcessInfo,1);
			else
				CalculateStress(res,rCurrentProcessInfo,1);

			// TODO: this is a workaround for the 4-node geometry (only one supported for shell interface)
			// generalize it asap
			rValues[0] = rValues[3] = res[0];
			rValues[1] = rValues[2] = res[1];
		}
		else if(rVariable == INTERFACE_ROTATION_JUMP || rVariable == INTERFACE_COUPLE_TRACTION)
		{
			std::vector<Vector> res;
			if(rVariable == INTERFACE_ROTATION_JUMP)
				CalculateStrain(res,rCurrentProcessInfo,2);
			else
				CalculateStress(res,rCurrentProcessInfo,2);

			// TODO: this is a workaround for the 4-node geometry (only one supported for shell interface)
			// generalize it asap
			rValues[0] = rValues[3] = res[0];
			rValues[1] = rValues[2] = res[1];
		}
		else
		{
			// TODO: this is a workaround for the 4-node geometry (only one supported for shell interface)
			// generalize it asap
			Vector res;
			res.clear();
			//mConstitutiveLawVector[0]->GetValue(rVariable, res);
			rValues[0] = rValues[3] = res;
			res.clear();
			//mConstitutiveLawVector[1]->GetValue(rVariable, res);
			rValues[1] = rValues[2] = res;
		}
    }

    void ShellThickInterfaceElement::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                           std::vector<Matrix>& rValues,
                                                           const ProcessInfo& rCurrentProcessInfo)
    {
		SizeType num_gp = GetGeometry().IntegrationPoints().size();
		if(rValues.size() != num_gp) rValues.resize(num_gp);
		SizeType half_num_gp = num_gp / 2;

		// TODO: this is a workaround for the 4-node geometry (only one supported for shell interface)
		// generalize it asap
		Matrix res;
		res.clear();
		//mConstitutiveLawVector[0]->GetValue(rVariable, res);
		rValues[0] = rValues[3] = res;
		res.clear();
		//mConstitutiveLawVector[1]->GetValue(rVariable, res);
		rValues[1] = rValues[2] = res;
    }

    void ShellThickInterfaceElement::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
                                                           std::vector<array_1d<double,3> >& rValues,
                                                           const ProcessInfo& rCurrentProcessInfo)
    {
		SizeType num_gp = GetGeometry().IntegrationPoints().size();
		if(rValues.size() != num_gp) rValues.resize(num_gp);
		SizeType half_num_gp = num_gp / 2;

		// TODO: this is a workaround for the 4-node geometry (only one supported for shell interface)
		// generalize it asap
		array_1d<double,3> res;
		res.clear();
		//mConstitutiveLawVector[0]->GetValue(rVariable, res);
		rValues[0] = rValues[3] = res;
		res.clear();
		//mConstitutiveLawVector[1]->GetValue(rVariable, res);
		rValues[1] = rValues[2] = res;
    }

    void ShellThickInterfaceElement::GetValueOnIntegrationPoints(const Variable<array_1d<double,6> >& rVariable,
                                                           std::vector<array_1d<double,6> >& rValues,
                                                           const ProcessInfo& rCurrentProcessInfo)
    {
		SizeType num_gp = GetGeometry().IntegrationPoints().size();
		if(rValues.size() != num_gp) rValues.resize(num_gp);
		SizeType half_num_gp = num_gp / 2;

		// TODO: this is a workaround for the 4-node geometry (only one supported for shell interface)
		// generalize it asap
		array_1d<double,6> res;
		res.clear();
		//mConstitutiveLawVector[0]->GetValue(rVariable, res);
		rValues[0] = rValues[3] = res;
		res.clear();
		//mConstitutiveLawVector[1]->GetValue(rVariable, res);
		rValues[1] = rValues[2] = res;
    }

    // =====================================================================================
    //
    // Class ShellThickInterfaceElement - Private methods
    //
    // =====================================================================================

    void ShellThickInterfaceElement::DecimalCorrection(Vector& a)
    {
        /*double norm = norm_2(a);
        double tolerance = std::max(norm * 1.0E-12, 1.0E-12);
        for(SizeType i = 0; i < a.size(); i++)
            if(std::abs(a(i)) < tolerance)
                a(i) = 0.0;*/
    }

	void ShellThickInterfaceElement::CalculatePermutation(InterfaceIndexPermutation& p)
	{
		// TODO: for now quad 4n is the only available geometry for this element.
		// the standard node numbering for this geometry is not the one that we want for an interface.
		// so now we set has permutation = true. (for other geometries it might change)
		p.HasPermutation = true;
		p.Permutation.clear();
		GeometryType& geom = GetGeometry();
		// handle case of quad 3d 4n (this is the only one supported by now)
		if(geom.size() == 4)
		{
            p.HasPermutation = true;
			// 1:6  7:12  19:24  13:18
			p.Permutation.resize(24);
			p.Permutation[0] = 0;
			p.Permutation[1] = 1;
			p.Permutation[2] = 2;
			p.Permutation[3] = 3;
			p.Permutation[4] = 6;
			p.Permutation[5] = 7;
			p.Permutation[6] = 4;
			p.Permutation[7] = 5;
			for(unsigned int i = 0; i < 6; i++)
			{
				p.Permutation[     i] =      i;
				p.Permutation[ 6 + i] =  6 + i;
				p.Permutation[12 + i] = 18 + i;
				p.Permutation[18 + i] = 12 + i;
			}
		}
	}

	void ShellThickInterfaceElement::CalculateDeltaPosition(Matrix& rDeltaPosition)
	{
		GeometryType& geom = GetGeometry();
		const unsigned int number_of_nodes = geom.PointsNumber();
		unsigned int dimension = 3;

		rDeltaPosition = zero_matrix<double>( number_of_nodes , dimension);

		for ( unsigned int i = 0; i < number_of_nodes; i++ )
		{
			const NodeType& iNode = geom[i];
			rDeltaPosition(i, 0) = iNode.X() - iNode.X0();
			rDeltaPosition(i, 1) = iNode.Y() - iNode.Y0();
			rDeltaPosition(i, 2) = iNode.Z() - iNode.Z0();
		}
	}

	void ShellThickInterfaceElement::CalculateJacobianAndTransformationMatrix(const SizeType pointID,
												                                     Matrix& delta_position,
		                                                                             Matrix& jacobian,
												                                     double& J,
												                                     Matrix& iR)
	{
		GeometryType& geom = GetGeometry();
		SizeType ndim = geom.WorkingSpaceDimension();

		geom.Jacobian(jacobian, pointID, GetIntegrationMethod(), delta_position);

		// this is the only vector that we can obtain from the geometry
		array_1d<double, 3> vx;
		vx[0] = jacobian(0, 0);
		vx[1] = jacobian(1, 0);
		vx[2] = jacobian(2, 0);
		J = std::sqrt(vx[0]*vx[0] + vx[1]*vx[1] + vx[2]*vx[2]);
		vx /= J;

		// check if the element/or nodes have assigned a custom normal vector.
		// otherwise use a default one (vz)
		array_1d<double, 3> vz;
		vz.clear();
		if(std::abs(vx[2]) > 0.99999)
			vz[1] = 1.0; // align local z with global Y
		else
			vz[2] = 1.0; // align local z with global Z

		array_1d<double, 3> vy;
		MathUtils<double>::CrossProduct(vy,  vz, vx);
		MathUtils<double>::CrossProduct(vz,  vx, vy);
		vy /= MathUtils<double>::Norm3(vy);
		vz /= MathUtils<double>::Norm3(vz);
		for(SizeType i = 0; i < 3; i++)
		{
			iR(0, i) = vx[i];
			iR(1, i) = vy[i];
			iR(2, i) = vz[i];
		}
	}

	double ShellThickInterfaceElement::CalculateIntegrationWeight(double J,
		                                                                 double iw)
	{
		// Note: here we multiply the integration weight by 2,
		// because we're using the geometry of a continuum element
		// but we're using only half of the integration points.
		double dV = J * iw * 2.0;
		return dV;
	}

	void ShellThickInterfaceElement::CalculateLocalDisplacementVector(const InterfaceIndexPermutation& P,
		                                                                     const Matrix& R,
		                                                                     const Vector& UG,
												                             Vector& UL)
	{
		GeometryType& geom = GetGeometry();
		SizeType nnodes = geom.size();
		if(P.HasPermutation)
		{
			for(SizeType inode = 0; inode < nnodes; inode++)
			{
				SizeType pos = inode*6;
				for(SizeType i = 0; i < 3; i++)
				{
					double temp_u = 0.0;
					double temp_r = 0.0;
					for(SizeType j = 0; j < 3; j++)
					{
						temp_u += R(i, j) * UG( P.Permutation[ pos+j  ] );
						temp_r += R(i, j) * UG( P.Permutation[ pos+j+3] );
					}
					UL(pos+i  ) = temp_u;
					UL(pos+i+3) = temp_r;
				}
			}
		}
		else
		{
			for(SizeType inode = 0; inode < nnodes; inode++)
			{
				SizeType pos = inode*6;
				for(SizeType i = 0; i < 3; i++)
				{
					double temp_u = 0.0;
					double temp_r = 0.0;
					for(SizeType j = 0; j < 3; j++)
					{
						temp_u += R(i, j) * UG(pos+j  );
						temp_r += R(i, j) * UG(pos+j+3);
					}
					UL(pos+i  ) = temp_u;
					UL(pos+i+3) = temp_r;
				}
			}
		}
	}

	void ShellThickInterfaceElement::TransformToGlobalAndAdd(const InterfaceIndexPermutation& P,
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

		Matrix RTK;
		if(LHSrequired)
			RTK.resize(6,6,false);

		if(P.HasPermutation)
		{
			for(SizeType node_i = 0; node_i < nnodes; node_i++)
			{
				SizeType pos_i = node_i*6;

				if(RHSrequired)
				{
					for(SizeType i = 0; i < 3; i++)
					{
						double temp_RHS_u = 0.0;
						double temp_RHS_r = 0.0;
						for(SizeType j = 0; j < 3; j++)
						{
							temp_RHS_u += R(j, i) * RHS_local( pos_i+j   );
							temp_RHS_r += R(j, i) * RHS_local( pos_i+j+3 );
						}
						RHS_global( P.Permutation[ pos_i+i   ] ) += temp_RHS_u;
						RHS_global( P.Permutation[ pos_i+i+3 ] ) += temp_RHS_r;
					}
				}

				if(LHSrequired)
				{
					for(SizeType node_j = 0; node_j < nnodes; node_j++)
					{
						SizeType pos_j = node_j*6;

						RTK.clear();
						for(SizeType i = 0; i < 3; i++)
						{
							for(SizeType j = 0; j < 3; j++)
							{
								double temp_K_uu = 0.0;
								double temp_K_ur = 0.0;
								double temp_K_ru = 0.0;
								double temp_K_rr = 0.0;
								for(SizeType k = 0; k < 3; k++)
								{
									temp_K_uu += R(k, i) * LHS_local( pos_i+k  , pos_j+j   );
									temp_K_ur += R(k, i) * LHS_local( pos_i+k  , pos_j+j+3 );
									temp_K_ru += R(k, i) * LHS_local( pos_i+k+3, pos_j+j   );
									temp_K_rr += R(k, i) * LHS_local( pos_i+k+3, pos_j+j+3 );
								}
								RTK(i  , j  ) = temp_K_uu;
								RTK(i  , j+3) = temp_K_ur;
								RTK(i+3, j  ) = temp_K_ru;
								RTK(i+3, j+3) = temp_K_rr;
							}
						}
						for(SizeType i = 0; i < 3; i++)
						{
							for(SizeType j = 0; j < 3; j++)
							{
								double temp_K_uu = 0.0;
								double temp_K_ur = 0.0;
								double temp_K_ru = 0.0;
								double temp_K_rr = 0.0;
								for(SizeType k = 0; k < 3; k++)
								{
									temp_K_uu += RTK(i  , k  ) * R(k, j);
									temp_K_ur += RTK(i  , k+3) * R(k, j);
									temp_K_ru += RTK(i+3, k  ) * R(k, j);
									temp_K_rr += RTK(i+3, k+3) * R(k, j);
								}
								LHS_global( P.Permutation[ pos_i+i   ], P.Permutation[ pos_j+j   ] ) += temp_K_uu;
								LHS_global( P.Permutation[ pos_i+i   ], P.Permutation[ pos_j+j+3 ] ) += temp_K_ur;
								LHS_global( P.Permutation[ pos_i+i+3 ], P.Permutation[ pos_j+j   ] ) += temp_K_ru;
								LHS_global( P.Permutation[ pos_i+i+3 ], P.Permutation[ pos_j+j+3 ] ) += temp_K_rr;
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
				SizeType pos_i = node_i*6;

				if(RHSrequired)
				{
					for(SizeType i = 0; i < 3; i++)
					{
						double temp_RHS_u = 0.0;
						double temp_RHS_r = 0.0;
						for(SizeType j = 0; j < 3; j++)
						{
							temp_RHS_u += R(j, i) * RHS_local( pos_i+j   );
							temp_RHS_r += R(j, i) * RHS_local( pos_i+j+3 );
						}
						RHS_global( pos_i+i   ) += temp_RHS_u;
						RHS_global( pos_i+i+3 ) += temp_RHS_r;
					}
				}

				if(LHSrequired)
				{
					for(SizeType node_j = 0; node_j < nnodes; node_j++)
					{
						SizeType pos_j = node_j*6;

						RTK.clear();
						for(SizeType i = 0; i < 3; i++)
						{
							for(SizeType j = 0; j < 3; j++)
							{
								double temp_K_uu = 0.0;
								double temp_K_ur = 0.0;
								double temp_K_ru = 0.0;
								double temp_K_rr = 0.0;
								for(SizeType k = 0; k < 3; k++)
								{
									temp_K_uu += R(k, i) * LHS_local( pos_i+k  , pos_j+j   );
									temp_K_ur += R(k, i) * LHS_local( pos_i+k  , pos_j+j+3 );
									temp_K_ru += R(k, i) * LHS_local( pos_i+k+3, pos_j+j   );
									temp_K_rr += R(k, i) * LHS_local( pos_i+k+3, pos_j+j+3 );
								}
								RTK(i  , j  ) = temp_K_uu;
								RTK(i  , j+3) = temp_K_ur;
								RTK(i+3, j  ) = temp_K_ru;
								RTK(i+3, j+3) = temp_K_rr;
							}
						}
						for(SizeType i = 0; i < 3; i++)
						{
							for(SizeType j = 0; j < 3; j++)
							{
								double temp_K_uu = 0.0;
								double temp_K_ur = 0.0;
								double temp_K_ru = 0.0;
								double temp_K_rr = 0.0;
								for(SizeType k = 0; k < 3; k++)
								{
									temp_K_uu += RTK(i  , k  ) * R(k, j);
									temp_K_ur += RTK(i  , k+3) * R(k, j);
									temp_K_ru += RTK(i+3, k  ) * R(k, j);
									temp_K_rr += RTK(i+3, k+3) * R(k, j);
								}
								LHS_global( pos_i+i  , pos_j+j   ) += temp_K_uu;
								LHS_global( pos_i+i  , pos_j+j+3 ) += temp_K_ur;
								LHS_global( pos_i+i+3, pos_j+j   ) += temp_K_ru;
								LHS_global( pos_i+i+3, pos_j+j+3 ) += temp_K_rr;
							}
						}
					}
				}
			}
		}
	}

	void ShellThickInterfaceElement::CalculateBMatrix(const SizeType pointID,
		                                                     Matrix& B)
	{
		GeometryType& geom = GetGeometry();
		const Matrix& shapes = geom.ShapeFunctionsValues();

		SizeType ndof = 6;
		SizeType nstrain = 4;
		SizeType nnodes = geom.size();
		SizeType half_nnodes = nnodes / 2;

		for(SizeType k = 0; k < half_nnodes; k++)
		{
			double Nk = shapes(pointID, k);
			SizeType pos1 = k * ndof;
			SizeType pos2 = pos1 + half_nnodes*ndof;
			SizeType i,j;
			// dUn
			i = 0;
			j = 1;
			B(i, pos1+j) = -Nk;
			B(i, pos2+j) =  Nk;
			// dUx
			i = 1;
			j = 0;
			B(i, pos1+j) = -Nk;
			B(i, pos2+j) =  Nk;
			// dUy
			i = 2;
			j = 2;
			B(i, pos1+j) = -Nk;
			B(i, pos2+j) =  Nk;
			// dRx
			i = 3;
			j = 3;
			B(i, pos1+j) =  Nk;
			B(i, pos2+j) = -Nk;
		}
	}

	void ShellThickInterfaceElement::CalculateGeneralizedStrains(const SizeType pointID,
		                                                                const Matrix& B,
		                                                                const Vector& U,
											                            Vector& generalizedStrains)
	{
		noalias( generalizedStrains ) = prod( B, U );
	}

    void ShellThickInterfaceElement::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                                         VectorType& rRightHandSideVector,
                                                         ProcessInfo& rCurrentProcessInfo,
                                                         const bool LHSrequired,
                                                         const bool RHSrequired)
    {
        // general data
		PropertiesType props = GetProperties();
		const GeometryType& geom = GetGeometry();
		SizeType nstrain = 4;
		SizeType nnodes = geom.size();
		SizeType ndofs = 6 * nnodes;

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
		Matrix jacobian(3, 1);
		double J;

		// strain-displacement matrix
		Matrix B(nstrain, ndofs, 0.0);

		// material point calculation data
		Matrix D(nstrain, nstrain);
		Vector generalizedStrains(nstrain);
		Vector generalizedStresses(nstrain);

		// transformation matrix
		Matrix iR(3, 3);

		// permutation data
		InterfaceIndexPermutation permutation;
		CalculatePermutation(permutation);

		// LHS and RHS in local coordinate system
		Matrix Kloc(ndofs, ndofs);
		Vector Rloc(ndofs);

		// auxiliary data to avoid extra memory allocations
		Matrix BTD(ndofs, nstrain);

		// initialize material parameters
		ConstitutiveLaw::Parameters Parameters(GetGeometry(), GetProperties(),rCurrentProcessInfo);
		Flags &ConstitutiveLawOptions = Parameters.GetOptions();
		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
		if(LHSrequired) ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
		Parameters.SetStrainVector(generalizedStrains);
		Parameters.SetStressVector(generalizedStresses);
		Parameters.SetConstitutiveMatrix(D);
		Matrix F( IdentityMatrix(3, 3) );
		double detF(1.0);
		Parameters.SetDeformationGradientF( F );
		Parameters.SetDeterminantF( detF );
		Matrix DN_DX(nnodes, 2, 0.0);
		Vector N(nnodes);
		Parameters.SetShapeFunctionsDerivatives(DN_DX);

		// elastic material test---------------------------------------------------
		double kN = props[NORMAL_STIFFNESS];
		double kT = props[TANGENTIAL_STIFFNESS];
		double th = props[THICKNESS];
		/*
		[ Kn*h,          0,          0,           0]
		[    0, (2*Kt*h)/3,          0,           0]
		[    0,          0, (2*Kt*h)/3,           0]
		[    0,          0,          0, (Kn*h^3)/12]
		*/
		D.clear();
		D(0,0) = kN*th;
		D(1,1) = (5.0*kT*th)/6.0;
		D(2,2) = D(1,1);
		D(3,3) = (kN*th*th*th)/12.0;

		// loop over the integration points.
		// note: loop only the first half of them, then multiply the integration weight by 2
		const GeometryType::IntegrationPointsArrayType& integrationPoints = geom.IntegrationPoints();

		Matrix B0(nstrain, ndofs, 0.0);
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
			//mConstitutiveLawVector[intp_id]->CalculateMaterialResponseCauchy(Parameters);
			noalias(generalizedStresses) = prod(D, generalizedStrains);

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

	void ShellThickInterfaceElement::CalculateStrain(std::vector<Vector>& strain, const ProcessInfo& rCurrentProcessInfo, int mode)
	{
		// init the output container. it will contain ngauss/2 output variables
		strain.clear();

		// general data
		PropertiesType props = GetProperties();
		const GeometryType& geom = GetGeometry();
		SizeType nstrain = 4;
		SizeType nnodes = geom.size();
		SizeType ndofs = 6 * nnodes;

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
		Matrix jacobian(3, 1);
		double J;

		// strain-displacement matrix
		Matrix B(nstrain, ndofs, 0.0);

		// material point calculation data
		Vector generalizedStrains(nstrain);

		// transformation matrix
		Matrix iR(3, 3);

		// permutation data
		InterfaceIndexPermutation permutation;
		CalculatePermutation(permutation);

		Vector vec3(3);

		// loop over the integration points.
		// note: loop only the first half of them, then multiply the integration weight by 2
		const GeometryType::IntegrationPointsArrayType& integrationPoints = geom.IntegrationPoints();

		Matrix B0(nstrain, ndofs, 0.0);
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
			if(mode == 1)
			{
				// displacement jump
				vec3(0) = generalizedStrains(0);
				vec3(1) = generalizedStrains(1);
				vec3(2) = generalizedStrains(2);
				strain.push_back(vec3);
			}
			else
			{
				// rotation jump
				vec3(0) = generalizedStrains(3);
				vec3(1) = 0.0;
				vec3(2) = 0.0;
				strain.push_back(vec3);
			}
		}
	}

	void ShellThickInterfaceElement::CalculateStress(std::vector<Vector>& stress, const ProcessInfo& rCurrentProcessInfo, int mode)
	{
		// init the output container. it will contain ngauss/2 output variables
		stress.clear();

		// general data
		PropertiesType props = GetProperties();
		const GeometryType& geom = GetGeometry();
		SizeType nstrain = 4;
		SizeType nnodes = geom.size();
		SizeType ndofs = 6 * nnodes;

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
		Matrix jacobian(3, 1);
		double J;

		// strain-displacement matrix
		Matrix B(nstrain, ndofs, 0.0);

		// material point calculation data
		Matrix D(nstrain, nstrain);
		Vector generalizedStrains(nstrain);
		Vector generalizedStresses(nstrain);

		// transformation matrix
		Matrix iR(3, 3);

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
		Matrix F( IdentityMatrix(3, 3) );
		double detF(1.0);
		Parameters.SetDeformationGradientF( F );
		Parameters.SetDeterminantF( detF );
		Matrix DN_DX(nnodes, 2, 0.0);
		Vector N(nnodes);
		Parameters.SetShapeFunctionsDerivatives(DN_DX);

		// elastic material test---------------------------------------------------
		double kN = props[NORMAL_STIFFNESS];
		double kT = props[TANGENTIAL_STIFFNESS];
		double th = props[THICKNESS];
		/*
		[ Kn*h,          0,          0,           0]
		[    0, (2*Kt*h)/3,          0,           0]
		[    0,          0, (2*Kt*h)/3,           0]
		[    0,          0,          0, (Kn*h^3)/12]
		*/
		D.clear();
		D(0,0) = kN*th;
		D(1,1) = (5.0*kT*th)/5.0;
		D(2,2) = D(1,1);
		D(3,3) = (kN*th*th*th)/12.0;

		Vector vec3(3);

		// loop over the integration points.
		// note: loop only the first half of them, then multiply the integration weight by 2
		const GeometryType::IntegrationPointsArrayType& integrationPoints = geom.IntegrationPoints();

		Matrix B0(nstrain, ndofs, 0.0);
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
			//mConstitutiveLawVector[intp_id]->CalculateMaterialResponseCauchy(Parameters);
			noalias(generalizedStresses) = prod(D, generalizedStrains);

			// store output
			if(mode == 1)
			{
				// traction vector
				vec3(0) = generalizedStresses(0);
				vec3(1) = generalizedStresses(1);
				vec3(2) = generalizedStresses(2);
				stress.push_back(vec3);
			}
			else
			{
				// traction-couple vector
				vec3(0) = generalizedStresses(3);
				vec3(1) = 0.0;
				vec3(2) = 0.0;
				stress.push_back(vec3);
			}
		}
	}

	// =====================================================================================
	//
	// Class ShellThickInterfaceElement - Serialization
	//
	// =====================================================================================

    void ShellThickInterfaceElement::save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
    }

    void ShellThickInterfaceElement::load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
    }

}
