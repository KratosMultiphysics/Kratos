//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#include "agq4_element.hpp"
#include "multiscale_application_variables.h"
#include "custom_utilities/math_helpers.h"
#include <string>
#include <iomanip>

// Number of external dofs
#ifdef AGQ4_DRILLING_DOF
#define AGQ4_NU 12
#define AGQ4_NDOF_PER_NODE 3
#else
#define AGQ4_NU 8
#define AGQ4_NDOF_PER_NODE 2
#endif

// number of internal dofs
#define AGQ4_NQ 4

// default integration method. can be 2x2 or 3x3. what's the best is still to be evaluated
#define AGQ4_DEF_INTEGRATION_METHOD Kratos::GeometryData::GI_GAUSS_3

// this flag activate the correction matrix for the incompatible modes to pass the patch test.
// this makes BQ = BQ - 1/V*int{BQ*dV}
#define AGQ4_INCOPATIBLE_MODES_CORRECTION_MTX

// this flag can be 1 or 2. it switches the formulations AGQ4I / AGQ4II.
// note that with AGQ4II the strict patch test is not fulfilled, but the weak patch test is fullfilled.
// with AGQI also the strict patch test is fulfilled, provided that the correction matrix
// for incompatible modes is active.
// so the best choice is to leave this flag = 1
#define AGQ4_FORMULATION 1

// this flag can be used to remove the incopatible modes. just for testing.
//#define AGQ4_SUPPRESS_ENHANCEMENT

namespace Kratos
{

    AGQ4Element::AGQ4Element(IndexType NewId,
                                   GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
		mThisIntegrationMethod = AGQ4_DEF_INTEGRATION_METHOD;

		mErrorCode   = 0.0;
		mInitialized = false;

		m_Q           = ZeroVector(AGQ4_NQ);
		m_Q_converged = ZeroVector(AGQ4_NQ);
		m_U           = ZeroVector(AGQ4_NU);
		m_U_converged = ZeroVector(AGQ4_NU);
		m_Q_residual  = ZeroVector(AGQ4_NQ);
		m_KQQ_inv     = ZeroMatrix(AGQ4_NQ, AGQ4_NQ);
		m_KQU         = ZeroMatrix(AGQ4_NQ, AGQ4_NU);
		m_KUQ         = ZeroMatrix(AGQ4_NU, AGQ4_NQ);

#ifdef AGQ4_DRILLING_DOF
		mG0                    = 0.0;
		m_first_step_finalized = false;
#endif // AGQ4_DRILLING_DOF

    }

    AGQ4Element::AGQ4Element(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
		mThisIntegrationMethod = AGQ4_DEF_INTEGRATION_METHOD;

		mErrorCode   = 0.0;
		mInitialized = false;

		m_Q           = ZeroVector(AGQ4_NQ);
		m_Q_converged = ZeroVector(AGQ4_NQ);
		m_U           = ZeroVector(AGQ4_NU);
		m_U_converged = ZeroVector(AGQ4_NU);
		m_Q_residual  = ZeroVector(AGQ4_NQ);
		m_KQQ_inv     = ZeroMatrix(AGQ4_NQ, AGQ4_NQ);
		m_KQU         = ZeroMatrix(AGQ4_NQ, AGQ4_NU);
		m_KUQ         = ZeroMatrix(AGQ4_NU, AGQ4_NQ);

#ifdef AGQ4_DRILLING_DOF
		mG0                    = 0.0;
		m_first_step_finalized = false;
#endif // AGQ4_DRILLING_DOF
    }

    AGQ4Element::~AGQ4Element()
    {
    }

    Element::Pointer AGQ4Element::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        GeometryType::Pointer newGeom( GetGeometry().Create(ThisNodes) );
        return Element::Pointer( new AGQ4Element(NewId, newGeom, pProperties) );
    }

    void AGQ4Element::Initialize()
    {
		KRATOS_TRY

		// NOTE:
		// This is the standard (previous) implementation:
		// If we are here, it means that no one already set up the constitutive law vector
		// through the method SetValue<CONSTITUTIVE_LAW_POINTER>

		const GeometryType & geom = GetGeometry();
		const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints( mThisIntegrationMethod );

		// initialization for internal dofs members ===================================================
		if(!mInitialized)
		{
			m_Q.clear();
			m_Q_converged.clear();
			for(unsigned int i = 0; i < geom.PointsNumber(); i++)
			{
#ifdef AGQ4_DRILLING_DOF
				unsigned int ii = i * 3;
				const array_1d<double, 3>& initialDispl = geom[i].FastGetSolutionStepValue(DISPLACEMENT);
				const array_1d<double, 3>& initialRot   = geom[i].FastGetSolutionStepValue(ROTATION);
				m_U(ii    ) = initialDispl(0);
				m_U(ii + 1) = initialDispl(1);
				m_U(ii + 2) = initialRot(2);
				m_U_converged(ii    ) = initialDispl(0);
				m_U_converged(ii + 1) = initialDispl(1);
				m_U_converged(ii + 2) = initialRot(2);
#else
				unsigned int ii = i * 2;
				const array_1d<double, 3>& initialDispl = geom[i].FastGetSolutionStepValue(DISPLACEMENT);
				m_U(ii    ) = initialDispl(0);
				m_U(ii + 1) = initialDispl(1);
				m_U_converged(ii    ) = initialDispl(0);
				m_U_converged(ii + 1) = initialDispl(1);
#endif // AGQ4_DRILLING_DOF
			}
			mInitialized = true;
		}

#ifdef AGQ4_DRILLING_DOF
		// initialization for drilling dof data =========================================================
		mG0 = 0.0;
		m_first_step_finalized = false;
#endif

		//Constitutive Law initialization ===============================================================
		if ( mConstitutiveLawVector.size() != integration_points.size() )
		{
			mConstitutiveLawVector.resize( integration_points.size() );
		}
		else
		{
			// check whether the constitutive law pointers have been already set up
			bool already_set_up = true;
			for(unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
			{
				if(mConstitutiveLawVector[i] == NULL)
					already_set_up = false;
			}
			if(already_set_up)
			{
				for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
				{
					mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), geom,
							row( geom.ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
				}

				return; // if so, we are done here!
			}
		}

		// NOTE:
		// This is the standard (previous) implementation:
		// If we are here, it means that no one already set up the constitutive law vector
		// through the method SetValue<CONSTITUTIVE_LAW_POINTER>

		if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
		{
			for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
			{
				mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
				mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), geom,
						row( geom.ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
			}
		}
		else
		{
			KRATOS_THROW_ERROR( std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id() )
		}

		KRATOS_CATCH( "" )
    }

    void AGQ4Element::ResetConstitutiveLaw()
    {
        KRATOS_TRY

        const GeometryType & geom = GetGeometry();
        const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

        const Properties& props = GetProperties();
        for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
            mConstitutiveLawVector[i]->ResetMaterial(props, geom, row(shapeFunctionsValues, i));

        KRATOS_CATCH("")
    }

	Kratos::GeometryData::IntegrationMethod AGQ4Element::GetIntegrationMethod() const
    {
        return mThisIntegrationMethod;
    }

    void AGQ4Element::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        if(rResult.size() != AGQ4_NU)
            rResult.resize(AGQ4_NU, false);

        GeometryType & geom = this->GetGeometry();

        for(int i = 0; i < geom.PointsNumber(); i++)
        {
#ifdef AGQ4_DRILLING_DOF
			int index = i * 3;
			NodeType & iNode = geom[i];
			rResult[index]     = iNode.GetDof(DISPLACEMENT_X).EquationId();
			rResult[index + 1] = iNode.GetDof(DISPLACEMENT_Y).EquationId();
			rResult[index + 2] = iNode.GetDof(ROTATION_Z).EquationId();
#else
			int index = i * 2;
			NodeType & iNode = geom[i];
			rResult[index]     = iNode.GetDof(DISPLACEMENT_X).EquationId();
			rResult[index + 1] = iNode.GetDof(DISPLACEMENT_Y).EquationId();
#endif // AGQ4_DRILLING_DOF
        }
    }

    void AGQ4Element::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        ElementalDofList.resize(0);
        ElementalDofList.reserve(AGQ4_NU);

        GeometryType & geom = this->GetGeometry();

        for (int i = 0; i < geom.PointsNumber(); i++)
        {
            NodeType & iNode = geom[i];

            ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Y));
#ifdef AGQ4_DRILLING_DOF
			ElementalDofList.push_back(iNode.pGetDof(ROTATION_Z));
#endif
        }
    }

    int AGQ4Element::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& geom = GetGeometry();

        // verify that the variables are correctly initialized
        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered","");
#ifdef AGQ4_DRILLING_DOF
		if(ROTATION.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"ROTATION has Key zero! (check if the application is correctly registered","");
#endif
        if(VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered","");
        if(ACCELERATION.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered","");
        if(DENSITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY has Key zero! (check if the application is correctly registered","");
        if(THICKNESS.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"THICKNESS has Key zero! (check if the application is correctly registered","");
        if(CONSTITUTIVE_LAW.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"CONSTITUTIVE_LAW has Key zero! (check if the application is correctly registered","");

        // verify that the dofs exist
        for(unsigned int i=0; i<geom.PointsNumber(); i++)
        {
            if(geom[i].SolutionStepsDataHas(DISPLACEMENT) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing variable DISPLACEMENT on node ",geom[i].Id());
            if(geom[i].HasDofFor(DISPLACEMENT_X) == false || geom[i].HasDofFor(DISPLACEMENT_Y) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable DISPLACEMENT on node ",GetGeometry()[i].Id());
#ifdef AGQ4_DRILLING_DOF
            if(geom[i].SolutionStepsDataHas(ROTATION) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing variable ROTATION on node ",geom[i].Id());
            if(geom[i].HasDofFor(ROTATION_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing one of the dofs for the variable ROTATION on node ",GetGeometry()[i].Id());
#endif
            if(geom[i].GetBufferSize() < 2)
                KRATOS_THROW_ERROR(std::logic_error, "This Element needs at least a buffer size = 2", "");
        }

        // check properties
        if(this->pGetProperties() == NULL)
            KRATOS_THROW_ERROR(std::logic_error, "Properties not provided for element ", this->Id());

        const PropertiesType & props = this->GetProperties();

        if(!props.Has(CONSTITUTIVE_LAW))
			KRATOS_THROW_ERROR(std::logic_error,"CONSTITUTIVE_LAW not provided for element ",this->Id());
		const ConstitutiveLaw::Pointer & claw = props[CONSTITUTIVE_LAW];
        if(claw == NULL)
            KRATOS_THROW_ERROR(std::logic_error,"CONSTITUTIVE_LAW not provided for element ",this->Id());
        claw->Check(props, geom, rCurrentProcessInfo);

        return 0;

        KRATOS_CATCH("")
    }

    void AGQ4Element::CleanMemory()
    {
    }

    void AGQ4Element::GetValuesVector(Vector& values, int Step)
    {
        if(values.size() != AGQ4_NU)
            values.resize(AGQ4_NU, false);

        const GeometryType & geom = GetGeometry();

        for (int i = 0; i < geom.PointsNumber(); i++)
        {
            const NodeType & iNode = geom[i];

#ifdef AGQ4_DRILLING_DOF
			int index = i*3;
			const array_1d<double,3>& disp = iNode.FastGetSolutionStepValue(DISPLACEMENT, Step);
			const array_1d<double,3>& roto = iNode.FastGetSolutionStepValue(ROTATION, Step);
			values[index]     = disp[0];
			values[index + 1] = disp[1];
			values[index + 2] = roto[2];
#else
			int index = i*2;
			const array_1d<double,3>& disp = iNode.FastGetSolutionStepValue(DISPLACEMENT, Step);
			values[index]     = disp[0];
			values[index + 1] = disp[1];
#endif // AGQ4_DRILLING_DOF
        }
    }

    void AGQ4Element::GetFirstDerivativesVector(Vector& values, int Step)
    {
        if(values.size() != AGQ4_NU)
            values.resize(AGQ4_NU,false);

        const GeometryType & geom = GetGeometry();

        for (int i = 0; i < geom.PointsNumber(); i++)
        {
            const NodeType & iNode = geom[i];
			const array_1d<double,3>& vel = iNode.FastGetSolutionStepValue(VELOCITY, Step);

#ifdef AGQ4_DRILLING_DOF
			int index = i * 3;
			values[index]        = vel[0];
			values[index + 1]    = vel[1];
			values[index + 2]    = 0.0;
#else
			int index = i * 2;
            values[index]        = vel[0];
            values[index + 1]    = vel[1];
#endif // AGQ4_DRILLING_DOF

        }
    }

    void AGQ4Element::GetSecondDerivativesVector(Vector& values, int Step)
    {
        if(values.size() != AGQ4_NU)
            values.resize(AGQ4_NU,false);

        const GeometryType & geom = GetGeometry();

        for (int i = 0; i < geom.PointsNumber(); i++)
        {
            const NodeType & iNode = geom[i];
			const array_1d<double,3>& acc = iNode.FastGetSolutionStepValue(ACCELERATION, Step);

#ifdef AGQ4_DRILLING_DOF
			int index = i * 3;
			values[index]        = acc[0];
			values[index + 1]    = acc[1];
			values[index + 2]    = 0.0;
#else
			int index = i * 2;
            values[index]        = acc[0];
            values[index + 1]    = acc[1];
#endif // AGQ4_DRILLING_DOF

        }
    }

    void AGQ4Element::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
    }

    void AGQ4Element::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
		Vector displacementVector(AGQ4_NU);
		GetValuesVector(displacementVector);

		Vector incrementalDispl(AGQ4_NU);
		noalias(incrementalDispl) = displacementVector - m_U;
		noalias(m_U) = displacementVector;

		array_1d<double, AGQ4_NQ> temp;
		noalias(temp)  = prod(m_KQU, incrementalDispl);
		noalias(temp) -= m_Q_residual;
		noalias(m_Q)  -= prod(m_KQQ_inv, temp);
    }

    void AGQ4Element::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
		const PropertiesType& props = GetProperties();
        const GeometryType & geom = GetGeometry();
        const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

        for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
            mConstitutiveLawVector[i]->InitializeSolutionStep(props, geom, row(shapeFunctionsValues, i), CurrentProcessInfo);

		m_U = m_U_converged;
		m_Q = m_Q_converged;
    }

    void AGQ4Element::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
		const PropertiesType& props = GetProperties();
        const GeometryType& geom = GetGeometry();
        const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

        for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
            mConstitutiveLawVector[i]->FinalizeSolutionStep(props, geom, row(shapeFunctionsValues, i), CurrentProcessInfo);

        m_U_converged = m_U;
		m_Q_converged = m_Q;

#ifdef AGQ4_DRILLING_DOF
		if(!m_first_step_finalized)
			m_first_step_finalized = true;
#endif

		/*const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(mThisIntegrationMethod);
		unsigned int ngauss = integration_points.size();

		Matrix delta_position;
		delta_position = CalculateDeltaPosition(delta_position);
		GeometryType::JacobiansType all_jacobians;
		all_jacobians = geom.Jacobian(all_jacobians, mThisIntegrationMethod, delta_position);

		double total_volume = 0.0;
		double mean_rt = 0.0;
		double mean_rc = 0.0;
		for(SizeType i = 0; i < ngauss; i++) {

			const GeometryType::IntegrationPointType & ip = integration_points[i];
			double gpw = ip.Weight();
			double det_jacobian = MathUtils<double>::Det2(all_jacobians[i]);
			double dV = det_jacobian*gpw;
			total_volume += dV;

			ConstitutiveLaw::Pointer& pclaw = mConstitutiveLawVector[i];
			double i_rt=0.0;
			double i_rc=0.0;
			pclaw->GetValue(ENH_STRAIN_PARAM_1, i_rt);
			pclaw->GetValue(ENH_STRAIN_PARAM_2, i_rc);
			mean_rt += i_rt*dV;
			mean_rc += i_rc*dV;
		}

		double tau = 0.25;
		mean_rt /= total_volume;
		mean_rc /= total_volume;
		for(SizeType i = 0; i < ngauss; i++) {
			ConstitutiveLaw::Pointer& pclaw = mConstitutiveLawVector[i];
			double i_rt=0.0;
			double i_rc=0.0;
			pclaw->GetValue(ENH_STRAIN_PARAM_1, i_rt);
			pclaw->GetValue(ENH_STRAIN_PARAM_2, i_rc);
			pclaw->SetValue(ENH_STRAIN_PARAM_1, tau*mean_rt+(1.0-tau)*i_rt, CurrentProcessInfo);
			pclaw->SetValue(ENH_STRAIN_PARAM_2, tau*mean_rc+(1.0-tau)*i_rc, CurrentProcessInfo);
		}*/
    }

    void AGQ4Element::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        if((rMassMatrix.size1() != AGQ4_NU) || (rMassMatrix.size2() != AGQ4_NU))
            rMassMatrix.resize(AGQ4_NU, AGQ4_NU, false);
        noalias(rMassMatrix) = ZeroMatrix(AGQ4_NU, AGQ4_NU);

        // Compute the total mass

		size_t num_nodes = GetGeometry().PointsNumber();
		double total_mass = GetGeometry().DomainSize() * GetProperties()[DENSITY] * GetProperties()[THICKNESS];
        double lumped_mass = total_mass / double(num_nodes);

		// loop on nodes
        for(size_t i = 0; i < num_nodes; i++)
        {
#ifdef AGQ4_DRILLING_DOF
			size_t index = i * 3;
			rMassMatrix(index, index)            = lumped_mass;
			rMassMatrix(index + 1, index + 1)    = lumped_mass;
			rMassMatrix(index + 2, index + 2)    = 0.0;
#else
			size_t index = i * 2;
            rMassMatrix(index, index)            = lumped_mass;
            rMassMatrix(index + 1, index + 1)    = lumped_mass;
#endif // AGQ4_DRILLING_DOF

        }
    }

    void AGQ4Element::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        if((rDampingMatrix.size1() != AGQ4_NU) || (rDampingMatrix.size2() != AGQ4_NU))
            rDampingMatrix.resize(AGQ4_NU, AGQ4_NU, false);

        noalias( rDampingMatrix ) = ZeroMatrix(AGQ4_NU, AGQ4_NU);
    }

    void AGQ4Element::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                     VectorType& rRightHandSideVector,
                                                     ProcessInfo& rCurrentProcessInfo)
    {
        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
    }

    void AGQ4Element::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                       ProcessInfo& rCurrentProcessInfo)
    {
        Matrix dummy;
        CalculateAll(dummy, rRightHandSideVector, rCurrentProcessInfo, true, true);
    }

    // =====================================================================================
    //
    // Class AGQ4Element - Results on Gauss Points
    //
    // =====================================================================================

	void AGQ4Element::SetValueOnIntegrationPoints(const Variable<double>& rVariable,
													   std::vector<double>& rValues,
													   const ProcessInfo& rCurrentProcessInfo)
	{
		if(rValues.size() == mConstitutiveLawVector.size())
			for(unsigned int i = 0; i < rValues.size(); i++)
				mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
	}

    void AGQ4Element::SetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
													   std::vector<Vector>& rValues,
													   const ProcessInfo& rCurrentProcessInfo)
	{
		if(rValues.size() == mConstitutiveLawVector.size())
			for(unsigned int i = 0; i < rValues.size(); i++)
				mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
	}

    void AGQ4Element::SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
													   std::vector<Matrix>& rValues,
													   const ProcessInfo& rCurrentProcessInfo)
	{
		if(rValues.size() == mConstitutiveLawVector.size())
			for(unsigned int i = 0; i < rValues.size(); i++)
				mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
	}

	void AGQ4Element::SetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
													 std::vector<ConstitutiveLaw::Pointer>& rValues,
													 const ProcessInfo& rCurrentProcessInfo )
	{
		if ( mConstitutiveLawVector.size() != rValues.size() )
		{
			mConstitutiveLawVector.resize(rValues.size());
			if( mConstitutiveLawVector.size() != GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod  ) )
				KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", mConstitutiveLawVector.size() );
		}
		for(unsigned int i=0; i<rValues.size(); i++)
			mConstitutiveLawVector[i] = rValues[i];
	}

    void AGQ4Element::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                                     std::vector<double>& rValues,
                                                     const ProcessInfo& rCurrentProcessInfo)
    {
		unsigned int numgp = mConstitutiveLawVector.size();
        if(rValues.size() != numgp)
            rValues.resize(numgp);

		std::fill(rValues.begin(), rValues.end(), 0.0);

		if(mErrorCode != 0.0)
		{
			if(rVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE)
			{
				std::fill(rValues.begin(), rValues.end(), mErrorCode);
			}
		}

        for(int i = 0; i < numgp; i++) {
			rValues[i] = 0.0;
            mConstitutiveLawVector[i]->GetValue(rVariable, rValues[i]);
		}
    }

    void AGQ4Element::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                     std::vector<Vector>& rValues,
                                                     const ProcessInfo& rCurrentProcessInfo)
    {
		unsigned int numgp = mConstitutiveLawVector.size();
		if(rValues.size() != numgp)
            rValues.resize(numgp);

        for(int i = 0; i < numgp; i++) {
            mConstitutiveLawVector[i]->GetValue(rVariable, rValues[i]);
		}
    }

    void AGQ4Element::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                     std::vector<Matrix>& rValues,
                                                     const ProcessInfo& rCurrentProcessInfo)
	{
		unsigned int numgp = mConstitutiveLawVector.size();
		if(rValues.size() != numgp)
			rValues.resize(numgp);

		if(rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ||
		   rVariable == PK2_STRESS_TENSOR ||
		   rVariable == CAUCHY_STRESS_TENSOR)
		{



			// =============================================================
			// Get some references.
			PropertiesType & props = GetProperties();
			GeometryType & geom = GetGeometry();

			const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(mThisIntegrationMethod);
			unsigned int ngauss = integration_points.size();

			const GeometryType::ShapeFunctionsGradientsType& local_gradients = geom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
			const Matrix & shapeFunctions = geom.ShapeFunctionsValues(mThisIntegrationMethod);
			Vector iN(shapeFunctions.size2());
			// =============================================================


			// =============================================================
			// Get the current displacements in global coordinate system
			Vector U(AGQ4_NU);
			GetValuesVector(U, 0);
			// =============================================================


			// =============================================================
			// Initialize parameters for the material calculation
			Matrix D(3, 3, 0.0);         // material tangent matrix.
			Vector E(3);                 // strain vector
			Vector S(3);                 // stress vector

			ConstitutiveLaw::Parameters parameters(geom, props, rCurrentProcessInfo);
			parameters.SetStrainVector( E );
			parameters.SetStressVector( S );
			parameters.SetConstitutiveMatrix( D );
			Flags& options = parameters.GetOptions();
			options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
			options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
			double detF = 1.0;
			Matrix F(IdentityMatrix(2,2));
			parameters.SetDeterminantF(detF);
			parameters.SetDeformationGradientF(F);
			// =============================================================


			// =============================================================
			// extract the x and y coordinates
			// and compute bi = yj-yk, and ci = xk-xj
			// Eq 3
			array_1d<double, 4> X,Y,b,c;
			for(unsigned int i=0; i<4; i++) {
				X(i) = geom[i].X0();
				Y(i) = geom[i].Y0();
			}
			for(unsigned int i=0; i<4; i++) {
				unsigned int j = i+1; if(j>3) j=0;
				unsigned int k = j+1; if(k>3) k=0;
				b(i) = Y(j)-Y(k);
				c(i) = X(k)-X(j);
			}
			// =============================================================


			// =============================================================
			// areas of sub-triangles
			double A1, A2, A3, A;
			A1 = 0.5*( X(1)*Y(3)+X(0)*Y(1)+X(3)*Y(0)-X(1)*Y(0)-X(3)*Y(1)-X(0)*Y(3) );
			A2 = 0.5*( X(1)*Y(2)+X(0)*Y(1)+X(2)*Y(0)-X(1)*Y(0)-X(2)*Y(1)-X(0)*Y(2) );
			A3 = 0.5*( X(2)*Y(3)+X(1)*Y(2)+X(3)*Y(1)-X(2)*Y(1)-X(3)*Y(2)-X(1)*Y(3) );
			A  = A1+A3;
			// characteristic parameters of a quadrilateral
			// Eq 4
			array_1d<double, 4> g;
			g(0) = A1/A; g(1) = A2/A; g(2)=1.0-g(0); g(3) = 1.0-g(1);
			// =============================================================


			// =============================================================
			// iso-parametric coordinates of the nodes
			array_1d<double,4> KSAI, EITA;
			KSAI(0)=-1.0; KSAI(1)=1.0; KSAI(2)=1.0; KSAI(3)=-1.0;
			EITA(0)=-1.0; EITA(1)=-1.0; EITA(2)=1.0; EITA(3)=1.0;
			// =============================================================


			// =============================================================
			// set to zero vectors and matrices for the static condensation
			// before proceding with gauss integration
			m_KQU.clear();
			m_KUQ.clear();
			m_KQQ_inv.clear();
			m_Q_residual.clear();
			// =============================================================


			// =============================================================
			// perform the gauss integration loop

			array_1d<double,4> L; // area coordinates
			Matrix BU(3, AGQ4_NU); // B matrix for external dofs
			Matrix BQ(3, AGQ4_NQ); // B matrix for internal dofs
			array_1d<double,4> NUX,NUY,NQX,NQY;
			Matrix NUXY(4,2);
			int itype = AGQ4_FORMULATION;

#ifdef AGQ4_INCOPATIBLE_MODES_CORRECTION_MTX

			double th = props[THICKNESS];
			Matrix delta_position;
			delta_position = CalculateDeltaPosition(delta_position);
			GeometryType::JacobiansType all_jacobians;
			all_jacobians = geom.Jacobian(all_jacobians, mThisIntegrationMethod, delta_position);

			Matrix BQ_mean(3, AGQ4_NQ, 0.0);
			double total_volume = 0.0;
			for(unsigned int igauss = 0; igauss < ngauss; igauss++)
			{
				// get a reference of the current integration point and shape functions
				const GeometryType::IntegrationPointType & ip = integration_points[igauss];
				double gpx = ip.X();
				double gpy = ip.Y();
				double gpw = ip.Weight();

				double det_jacobian = MathUtils<double>::Det2(all_jacobians[igauss]);
				double dV = det_jacobian*gpw*th;
				total_volume += dV;

				// area coordinates of the gauss point (Eq 7)
				L(0) = 0.25*(1.0-gpx)*(g(1)*(1.0-gpy)+g(2)*(1.0+gpy));
				L(1) = 0.25*(1.0-gpy)*(g(3)*(1.0-gpx)+g(2)*(1.0+gpx));
				L(2) = 0.25*(1.0+gpx)*(g(0)*(1.0-gpy)+g(3)*(1.0+gpy));
				L(3) = 0.25*(1.0+gpy)*(g(0)*(1.0-gpx)+g(1)*(1.0+gpx));

				// strain matrix for internal dofs
				for(unsigned int i=0; i<2; i++)
				{ // begin: LOOP 5
					unsigned int j = i+1; if(j>3) j=0;
					unsigned int k = j+1; if(k>3) k=0;
					NQX(i)=(b(i)*L(k)+b(k)*L(i))/A/2.0;
					NQY(i)=(c(i)*L(k)+c(k)*L(i))/A/2.0;
					unsigned int index1 = i*AGQ4_NDOF_PER_NODE;
					unsigned int index2 = index1+1;
					BQ_mean(0,index1) += NQX(i)*dV;
					BQ_mean(1,index2) += NQY(i)*dV;
					BQ_mean(2,index1) += NQY(i)*dV;
					BQ_mean(2,index2) += NQX(i)*dV;
				} // end: LOOP 5
			}
			BQ_mean/=total_volume;
#endif // AGQ4_INCOPATIBLE_MODES_CORRECTION_MTX

			for(unsigned int igauss = 0; igauss < ngauss; igauss++)
			{
				Matrix& output = rValues[igauss];
				if(output.size1() != 2 || output.size2() != 2)
					output.resize(2,2,false);

				// get a reference of the current integration point and shape functions
				const GeometryType::IntegrationPointType & ip = integration_points[igauss];
				double gpx = ip.X();
				double gpy = ip.Y();
				double gpw = ip.Weight();
				noalias( iN ) = row( shapeFunctions, igauss );

				// area coordinates of the gauss point (Eq 7)
				L(0) = 0.25*(1.0-gpx)*(g(1)*(1.0-gpy)+g(2)*(1.0+gpy));
				L(1) = 0.25*(1.0-gpy)*(g(3)*(1.0-gpx)+g(2)*(1.0+gpx));
				L(2) = 0.25*(1.0+gpx)*(g(0)*(1.0-gpy)+g(3)*(1.0+gpy));
				L(3) = 0.25*(1.0+gpy)*(g(0)*(1.0-gpx)+g(1)*(1.0+gpx));

				// strain matrix for external dofs
				for(unsigned int i=0; i<4; i++)
				{ // begin: LOOP 3
					unsigned int j = i+1; if(j>3) j=0;
					unsigned int k = j+1; if(k>3) k=0;
					double SX(0.0);
					double SY(0.0);
					for(unsigned int ii=0; ii<4; ii++)
					{ // begin: LOOP 4
						unsigned int jj = ii+1; if(jj>3) jj=0;
						unsigned int kk = jj+1; if(kk>3) kk=0;
						unsigned int mm = kk+1; if(mm>3) mm=0;
						// begin: select formulation
						if(itype == 1)
						{
							SX = SX+b(ii)*KSAI(ii)*EITA(ii)*(3.0*(L(jj)-L(mm))+(g(jj)-g(kk)));
							SY = SY+c(ii)*KSAI(ii)*EITA(ii)*(3.0*(L(jj)-L(mm))+(g(jj)-g(kk)));
						}
						else // itype == 2
						{
							SX = SX+b(ii)*KSAI(ii)*EITA(ii)*(L(jj)-L(mm));
							SY = SY+c(ii)*KSAI(ii)*EITA(ii)*(L(jj)-L(mm));
						}
						// end: select formulation
					} // end: LOOP 4
					// begin: select formulation
					if(itype == 1)
					{
						NUX(i)=(b(i)+b(j))/A/2.0+KSAI(i)*EITA(i)*g(k)*SX/2.0/A/(1.0+g(0)*g(2)+g(1)*g(3));
						NUY(i)=(c(i)+c(j))/A/2.0+KSAI(i)*EITA(i)*g(k)*SY/2.0/A/(1.0+g(0)*g(2)+g(1)*g(3));
					}
					else // itype == 2
					{
						NUX(i)=(b(i)+b(j))/A/2.0+KSAI(i)*EITA(i)*g(k)*SX/2.0/A/(g(0)*g(2)+g(1)*g(3));
						NUY(i)=(c(i)+c(j))/A/2.0+KSAI(i)*EITA(i)*g(k)*SY/2.0/A/(g(0)*g(2)+g(1)*g(3));
					}
					// end: select formulation
					unsigned int index1 = i*AGQ4_NDOF_PER_NODE;
					unsigned int index2 = index1+1;
					BU(0,index1) = NUX(i); BU(0,index2) = 0.0;
					BU(1,index1) = 0.0;    BU(1,index2) = NUY(i);
					BU(2,index1) = NUY(i); BU(2,index2) = NUX(i);
					NUXY(i,0)=NUX(i);
					NUXY(i,1)=NUY(i);
				} // end: LOOP 3

				// strain matrix for internal dofs
				for(unsigned int i=0; i<2; i++)
				{ // begin: LOOP 5
					unsigned int j = i+1; if(j>3) j=0;
					unsigned int k = j+1; if(k>3) k=0;
					NQX(i)=(b(i)*L(k)+b(k)*L(i))/A/2.0;
					NQY(i)=(c(i)*L(k)+c(k)*L(i))/A/2.0;
					unsigned int index1 = i*AGQ4_NDOF_PER_NODE;
					unsigned int index2 = index1+1;
					BQ(0,index1) = NQX(i); BQ(0,index2) = 0.0;
					BQ(1,index1) = 0.0;    BQ(1,index2) = NQY(i);
					BQ(2,index1) = NQY(i); BQ(2,index2) = NQX(i);
				} // end: LOOP 5
#ifdef AGQ4_INCOPATIBLE_MODES_CORRECTION_MTX
				BQ -= BQ_mean;
#endif // AGQ4_INCOPATIBLE_MODES_CORRECTION_MTX

				// compute strain vector (compatible + enhanced)
				noalias(E)  = prod(BU, U);
#ifndef AGQ4_SUPPRESS_ENHANCEMENT
				noalias(E) += prod(BQ, m_Q);
#endif // !AGQ4_SUPPRESS_ENHANCEMENT

				// strain decimal correction
				DecimalCorrection(E);

				if(rVariable == GREEN_LAGRANGE_STRAIN_TENSOR)
				{
					noalias(output) = MathUtils<double>::StrainVectorToTensor(E);
				}
				else
				{
					ConstitutiveLaw::Pointer & mat = mConstitutiveLawVector[igauss];
					parameters.SetShapeFunctionsValues( iN );
					parameters.SetShapeFunctionsDerivatives( NUXY );
					mat->CalculateMaterialResponseCauchy( parameters );
					noalias(output) = MathUtils<double>::StressVectorToTensor(S);
				}

				// compute material response
				ConstitutiveLaw::Pointer & mat = mConstitutiveLawVector[igauss];
				parameters.SetShapeFunctionsValues( iN );
				parameters.SetShapeFunctionsDerivatives( NUXY );
				mat->CalculateMaterialResponseCauchy( parameters );


			} // and for : igauss
			// =============================================================


		}
		else
		{
			for(size_t i = 0; i < mConstitutiveLawVector.size(); i++)
				mConstitutiveLawVector[i]->GetValue(rVariable, rValues[i]);
		}
	}

	void AGQ4Element::GetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
													 std::vector<ConstitutiveLaw::Pointer>& rValues,
													 const ProcessInfo& rCurrentProcessInfo)
	{
		if ( rValues.size() != mConstitutiveLawVector.size() )
			rValues.resize(mConstitutiveLawVector.size());

		for(unsigned int i=0; i<rValues.size(); i++)
			rValues[i] = mConstitutiveLawVector[i];
	}


    // =====================================================================================
    //
    // Class AGQ4Element - Private methods
    //
    // =====================================================================================

	Matrix& AGQ4Element::CalculateDeltaPosition(Matrix & rDeltaPosition)
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

		return rDeltaPosition;
	}

	void AGQ4Element::DecimalCorrection(Vector& a)
	{
		return;
		double norm = norm_2(a);
		double tolerance = std::max(norm * 1.0E-6, 1.0E-8);
		tolerance = 1.0e-9;
		for(SizeType i = 0; i < a.size(); i++)
			if(std::abs(a(i)) < tolerance)
				a(i) = 0.0;
	}

    void AGQ4Element::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo,
                                      const bool LHSrequired,
                                      const bool RHSrequired)
    {
		// Resize the Left Hand Side if necessary,
		// and initialize it to Zero
		if((rLeftHandSideMatrix.size1() != AGQ4_NU) || (rLeftHandSideMatrix.size2() != AGQ4_NU))
			rLeftHandSideMatrix.resize(AGQ4_NU, AGQ4_NU, false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(AGQ4_NU, AGQ4_NU);

		// Resize the Right Hand Side if necessary,
		// and initialize it to Zero
		if(rRightHandSideVector.size() != AGQ4_NU)
			rRightHandSideVector.resize(AGQ4_NU, false);
		noalias(rRightHandSideVector) = ZeroVector(AGQ4_NU);

		// set no error
		mErrorCode = 0.0;


		// =============================================================
		// Get some references.
		PropertiesType & props = GetProperties();
		GeometryType & geom = GetGeometry();

		const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(mThisIntegrationMethod);
		unsigned int ngauss = integration_points.size();

		const GeometryType::ShapeFunctionsGradientsType& local_gradients = geom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
        const Matrix & shapeFunctions = geom.ShapeFunctionsValues(mThisIntegrationMethod);
		Vector iN(shapeFunctions.size2());

		// thickness
		double th = props[THICKNESS];
		// =============================================================


		// =============================================================
		// Get the current displacements in global coordinate system
		Vector U(AGQ4_NU);
		GetValuesVector(U, 0);
		Matrix delta_position;
		delta_position = CalculateDeltaPosition(delta_position);
		GeometryType::JacobiansType all_jacobians;
		all_jacobians = geom.Jacobian(all_jacobians, mThisIntegrationMethod, delta_position);
		// =============================================================


		// =============================================================
		// Initialize parameters for the material calculation
		Matrix D(3, 3, 0.0);         // material tangent matrix.
		Vector E(3);                 // strain vector
		Vector S(3);                 // stress vector

		ConstitutiveLaw::Parameters parameters(geom, props, rCurrentProcessInfo);
		parameters.SetStrainVector( E );
		parameters.SetStressVector( S );
		parameters.SetConstitutiveMatrix( D );
		Flags& options = parameters.GetOptions();
		options.Set(ConstitutiveLaw::COMPUTE_STRESS, RHSrequired);
		options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, LHSrequired);
		double detF = 1.0;
		Matrix F(IdentityMatrix(2,2));
		parameters.SetDeterminantF(detF);
		parameters.SetDeformationGradientF(F);
		// =============================================================


		// =============================================================
		// extract the x and y coordinates
		// and compute bi = yj-yk, and ci = xk-xj
		// Eq 3
		array_1d<double, 4> X,Y,b,c;
		for(unsigned int i=0; i<4; i++) {
			X(i) = geom[i].X0();
			Y(i) = geom[i].Y0();
		}
		for(unsigned int i=0; i<4; i++) {
			unsigned int j = i+1; if(j>3) j=0;
			unsigned int k = j+1; if(k>3) k=0;
			b(i) = Y(j)-Y(k);
			c(i) = X(k)-X(j);
		}
		// =============================================================


		// =============================================================
		// areas of sub-triangles
		double A1, A2, A3, A;
		A1 = 0.5*( X(1)*Y(3)+X(0)*Y(1)+X(3)*Y(0)-X(1)*Y(0)-X(3)*Y(1)-X(0)*Y(3) );
		A2 = 0.5*( X(1)*Y(2)+X(0)*Y(1)+X(2)*Y(0)-X(1)*Y(0)-X(2)*Y(1)-X(0)*Y(2) );
		A3 = 0.5*( X(2)*Y(3)+X(1)*Y(2)+X(3)*Y(1)-X(2)*Y(1)-X(3)*Y(2)-X(1)*Y(3) );
		A  = A1+A3;
		// characteristic parameters of a quadrilateral
		// Eq 4
		array_1d<double, 4> g;
		g(0) = A1/A; g(1) = A2/A; g(2)=1.0-g(0); g(3) = 1.0-g(1);
		// =============================================================


		// =============================================================
		// iso-parametric coordinates of the nodes
		array_1d<double,4> KSAI, EITA;
		KSAI(0)=-1.0; KSAI(1)=1.0; KSAI(2)=1.0; KSAI(3)=-1.0;
		EITA(0)=-1.0; EITA(1)=-1.0; EITA(2)=1.0; EITA(3)=1.0;
		// =============================================================


		// =============================================================
		// set to zero vectors and matrices for the static condensation
		// before proceding with gauss integration
		m_KQU.clear();
		m_KUQ.clear();
		m_KQQ_inv.clear();
		m_Q_residual.clear();
		// =============================================================


		// =============================================================
		// perform the gauss integration loop

		array_1d<double,4> L; // area coordinates
		Matrix BU(3, AGQ4_NU); // B matrix for external dofs
		Matrix BQ(3, AGQ4_NQ); // B matrix for internal dofs
		array_1d<double,4> NUX,NUY,NQX,NQY;
		Matrix NUXY(4,2);
		int itype = AGQ4_FORMULATION;
		// auxiliary matrices
		Matrix BUT_D(AGQ4_NU, 3);
		Matrix BQT_D(AGQ4_NQ, 3);
		Matrix D_BQ(3, AGQ4_NQ);

#ifdef AGQ4_INCOPATIBLE_MODES_CORRECTION_MTX
		Matrix BQ_mean(3, AGQ4_NQ, 0.0);
		double total_volume = 0.0;
		for(unsigned int igauss = 0; igauss < ngauss; igauss++)
		{
			// get a reference of the current integration point and shape functions
			const GeometryType::IntegrationPointType & ip = integration_points[igauss];
			double gpx = ip.X();
			double gpy = ip.Y();
			double gpw = ip.Weight();

			double det_jacobian = MathUtils<double>::Det2(all_jacobians[igauss]);
			double dV = det_jacobian*gpw*th;
			total_volume += dV;

			// area coordinates of the gauss point (Eq 7)
			L(0) = 0.25*(1.0-gpx)*(g(1)*(1.0-gpy)+g(2)*(1.0+gpy));
			L(1) = 0.25*(1.0-gpy)*(g(3)*(1.0-gpx)+g(2)*(1.0+gpx));
			L(2) = 0.25*(1.0+gpx)*(g(0)*(1.0-gpy)+g(3)*(1.0+gpy));
			L(3) = 0.25*(1.0+gpy)*(g(0)*(1.0-gpx)+g(1)*(1.0+gpx));

			// strain matrix for internal dofs
			for(unsigned int i=0; i<2; i++)
			{ // begin: LOOP 5
				unsigned int j = i+1; if(j>3) j=0;
				unsigned int k = j+1; if(k>3) k=0;
				NQX(i)=(b(i)*L(k)+b(k)*L(i))/A/2.0;
				NQY(i)=(c(i)*L(k)+c(k)*L(i))/A/2.0;
				unsigned int index1 = i*AGQ4_NDOF_PER_NODE;
				unsigned int index2 = index1+1;
				BQ_mean(0,index1) += NQX(i)*dV;
				BQ_mean(1,index2) += NQY(i)*dV;
				BQ_mean(2,index1) += NQY(i)*dV;
				BQ_mean(2,index2) += NQX(i)*dV;
			} // end: LOOP 5
		}
		BQ_mean/=total_volume;
#endif // AGQ4_INCOPATIBLE_MODES_CORRECTION_MTX

		for(unsigned int igauss = 0; igauss < ngauss; igauss++)
		{
			// get a reference of the current integration point and shape functions
			const GeometryType::IntegrationPointType & ip = integration_points[igauss];
			double gpx = ip.X();
			double gpy = ip.Y();
			double gpw = ip.Weight();
			noalias( iN ) = row( shapeFunctions, igauss );

			// area coordinates of the gauss point (Eq 7)
			L(0) = 0.25*(1.0-gpx)*(g(1)*(1.0-gpy)+g(2)*(1.0+gpy));
			L(1) = 0.25*(1.0-gpy)*(g(3)*(1.0-gpx)+g(2)*(1.0+gpx));
			L(2) = 0.25*(1.0+gpx)*(g(0)*(1.0-gpy)+g(3)*(1.0+gpy));
			L(3) = 0.25*(1.0+gpy)*(g(0)*(1.0-gpx)+g(1)*(1.0+gpx));

			// strain matrix for external dofs
			for(unsigned int i=0; i<4; i++)
			{ // begin: LOOP 3
				unsigned int j = i+1; if(j>3) j=0;
				unsigned int k = j+1; if(k>3) k=0;
				double SX(0.0);
				double SY(0.0);
				for(unsigned int ii=0; ii<4; ii++)
				{ // begin: LOOP 4
					unsigned int jj = ii+1; if(jj>3) jj=0;
					unsigned int kk = jj+1; if(kk>3) kk=0;
					unsigned int mm = kk+1; if(mm>3) mm=0;
					// begin: select formulation
					if(itype == 1)
					{
						SX = SX+b(ii)*KSAI(ii)*EITA(ii)*(3.0*(L(jj)-L(mm))+(g(jj)-g(kk)));
						SY = SY+c(ii)*KSAI(ii)*EITA(ii)*(3.0*(L(jj)-L(mm))+(g(jj)-g(kk)));
					}
					else // itype == 2
					{
						SX = SX+b(ii)*KSAI(ii)*EITA(ii)*(L(jj)-L(mm));
						SY = SY+c(ii)*KSAI(ii)*EITA(ii)*(L(jj)-L(mm));
					}
					// end: select formulation
				} // end: LOOP 4
				// begin: select formulation
				if(itype == 1)
				{
					NUX(i)=(b(i)+b(j))/A/2.0+KSAI(i)*EITA(i)*g(k)*SX/2.0/A/(1.0+g(0)*g(2)+g(1)*g(3));
					NUY(i)=(c(i)+c(j))/A/2.0+KSAI(i)*EITA(i)*g(k)*SY/2.0/A/(1.0+g(0)*g(2)+g(1)*g(3));
				}
				else // itype == 2
				{
					NUX(i)=(b(i)+b(j))/A/2.0+KSAI(i)*EITA(i)*g(k)*SX/2.0/A/(g(0)*g(2)+g(1)*g(3));
					NUY(i)=(c(i)+c(j))/A/2.0+KSAI(i)*EITA(i)*g(k)*SY/2.0/A/(g(0)*g(2)+g(1)*g(3));
				}
				// end: select formulation
				unsigned int index1 = i*AGQ4_NDOF_PER_NODE;
				unsigned int index2 = index1+1;
				BU(0,index1) = NUX(i); BU(0,index2) = 0.0;
				BU(1,index1) = 0.0;    BU(1,index2) = NUY(i);
				BU(2,index1) = NUY(i); BU(2,index2) = NUX(i);
				NUXY(i,0)=NUX(i);
				NUXY(i,1)=NUY(i);
			} // end: LOOP 3

			// strain matrix for internal dofs
			for(unsigned int i=0; i<2; i++)
			{ // begin: LOOP 5
				unsigned int j = i+1; if(j>3) j=0;
				unsigned int k = j+1; if(k>3) k=0;
				NQX(i)=(b(i)*L(k)+b(k)*L(i))/A/2.0;
				NQY(i)=(c(i)*L(k)+c(k)*L(i))/A/2.0;
				unsigned int index1 = i*AGQ4_NDOF_PER_NODE;
				unsigned int index2 = index1+1;
				BQ(0,index1) = NQX(i); BQ(0,index2) = 0.0;
				BQ(1,index1) = 0.0;    BQ(1,index2) = NQY(i);
				BQ(2,index1) = NQY(i); BQ(2,index2) = NQX(i);
			} // end: LOOP 5
#ifdef AGQ4_INCOPATIBLE_MODES_CORRECTION_MTX
			BQ -= BQ_mean;
#endif // AGQ4_INCOPATIBLE_MODES_CORRECTION_MTX

			// compute strain vector (compatible + enhanced)
			noalias(E)  = prod(BU, U);
#ifndef AGQ4_SUPPRESS_ENHANCEMENT
			noalias(E) += prod(BQ, m_Q);
#endif // !AGQ4_SUPPRESS_ENHANCEMENT

			// strain decimal correction
			DecimalCorrection(E);

			// compute material response
			ConstitutiveLaw::Pointer & mat = mConstitutiveLawVector[igauss];
			parameters.SetShapeFunctionsValues( iN );
			parameters.SetShapeFunctionsDerivatives( NUXY );
			mat->CalculateMaterialResponseCauchy( parameters );

			// integration weight
			//double det_jacobian = geom.DeterminantOfJacobian(igauss, mThisIntegrationMethod);
			double det_jacobian = MathUtils<double>::Det2(all_jacobians[igauss]);
			double dV = det_jacobian*gpw*th;
			S *= dV;
			D *= dV;

			// assemble matrix and vectors for external dofs
			noalias( BUT_D )                 = prod( trans( BU ), D );
			noalias( rLeftHandSideMatrix )  += prod( BUT_D, BU );
			noalias( rRightHandSideVector ) -= prod( trans( BU ), S );

			// assemble matrix and vectors for internal dofs
			noalias( BQT_D )         = prod( trans( BQ ), D );
			noalias( D_BQ )          = prod( D, BQ );
			noalias( m_KQQ_inv )    += prod( BQT_D, BQ );
			noalias( m_Q_residual ) -= prod( trans( BQ ), S );
			noalias( m_KQU )        += prod( BQT_D, BU );
			noalias( m_KUQ )        += prod( trans( BU ), D_BQ );
		} // and for : igauss
		// =============================================================


		// =============================================================
		// perform the static condensation
		Matrix KQQ_copy( m_KQQ_inv ); // aux
		Matrix KUQ_KQQ_inv(AGQ4_NU, AGQ4_NQ); // aux
		bool eas_error = false;
		// invert KQQ
		permutation_matrix<Matrix::size_type> pm( AGQ4_NQ );
		size_t lu_res = lu_factorize( KQQ_copy, pm );
		noalias( m_KQQ_inv ) = IdentityMatrix( AGQ4_NQ );
		lu_substitute( KQQ_copy, pm, m_KQQ_inv );
		// compute L' * H^-1
		noalias( KUQ_KQQ_inv ) = prod( m_KUQ, m_KQQ_inv );
		// modify the stiffness matrix and the residual vector
		// for the static condensation of the incompatible modes
#ifndef AGQ4_SUPPRESS_ENHANCEMENT
		noalias( rRightHandSideVector ) -= prod( KUQ_KQQ_inv, m_Q_residual );   // R_mod = R - L' * H^-1 * residual
		noalias( rLeftHandSideMatrix  ) -= prod( KUQ_KQQ_inv, m_KQU );          // K_mod = K - L' * H^-1 * L
#endif // !AGQ4_SUPPRESS_ENHANCEMENT
		//eas_error = (lu_res != 0);
		// =============================================================


		// =============================================================
		// Add body forces contributions. This doesn't depend on the coordinate system
		AddBodyForces(rRightHandSideVector);
		// =============================================================


		// set error code
		//if(eas_error) mErrorCode = -1.0;
    }

	void AGQ4Element::AddBodyForces(VectorType& rRightHandSideVector)
	{
		const GeometryType& geom = GetGeometry();

		// Get shape functions
		const Matrix & N = geom.ShapeFunctionsValues(mThisIntegrationMethod);

		double rho = GetProperties()[DENSITY];
		double th = GetProperties()[THICKNESS];

		// auxiliary
		array_1d<double, 3> bf;

		// gauss loop to integrate the external force vector
		const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(mThisIntegrationMethod);
		const GeometryType::ShapeFunctionsGradientsType& local_gradients = geom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
		for(unsigned int igauss = 0; igauss < integration_points.size(); igauss++)
		{
			const GeometryType::IntegrationPointType & ip = integration_points[igauss];
			double dV = ip.Weight() * geom.DeterminantOfJacobian(igauss, mThisIntegrationMethod) * th;

			// interpolate nodal volume accelerations to this gauss point
			// and obtain the body force vector
			bf.clear();
			for(unsigned int inode = 0; inode < geom.PointsNumber(); inode++)
			{
				if( geom[inode].SolutionStepsDataHas(VOLUME_ACCELERATION) ) //temporary, will be checked once at the beginning only
					bf += N(igauss,inode) * geom[inode].FastGetSolutionStepValue(VOLUME_ACCELERATION);
			}
			bf *= (rho * dV);

			// add it to the RHS vector
			for(unsigned int inode = 0; inode < geom.PointsNumber(); inode++)
			{
#ifdef AGQ4_DRILLING_DOF
				unsigned int index = inode*3;
#else
				unsigned int index = inode*2;
#endif
				double iN = N(igauss,inode);
				rRightHandSideVector[index + 0] += iN * bf[0];
				rRightHandSideVector[index + 1] += iN * bf[1];
			}
		}
	}

	// =====================================================================================
    //
    // Class AGQ4Element - Serialization
    //
    // =====================================================================================

	void AGQ4Element::save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
        rSerializer.save("CLaw", mConstitutiveLawVector);
        rSerializer.save("IntM", (int)mThisIntegrationMethod);
		// TODO save incompatible modes data
    }

    void AGQ4Element::load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
        rSerializer.save("CLaw", mConstitutiveLawVector);
		int temp;
        rSerializer.load("IntM", temp);
		mThisIntegrationMethod = (IntegrationMethod)temp;
		// TODO restore incompatible modes data
    }

}
