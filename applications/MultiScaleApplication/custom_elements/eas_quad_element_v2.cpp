//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#include "eas_quad_element_v2.hpp"
#include "multiscale_application_variables.h"
#include "custom_utilities/math_helpers.h"
#include <string>
#include <iomanip>

/* Enhanced strain parameters */
#define EASQ4_NE 2 // EAS2 test
//#define EASQ4_NE 3 // EAS2 + directional linear
//#define EASQ4_NE 4 // Taylor
//#define EASQ4_NE 5 // Simo and Rifai
//#define EASQ4_NE 6 // Pantuso
//#define EASQ4_NE 7 // Ramm

/* Number of displacement dofs (for all nodes) */
#ifdef EASQ4_DRILLING_DOF
#define EASQ4_NU 12
#else
#define EASQ4_NU 8
#endif

/* Number of strains */
#ifdef EASQ4_DRILLING_DOF
#define EASQ4_NSTRAIN 4
#else
#define EASQ4_NSTRAIN 3
#endif

//#define SUPPRESS_EAS_ENHANCEMENT

//#define EASQ4_DRILLING_DOF_MOD
#define EASQ4_DRILLING_PENALTY_SCALE 0.1

// default integration method. can be 2x2 or 3x3. what's the best is still to be evaluated
#define EASQ4_DEF_INTEGRATION_METHOD Kratos::GeometryData::GI_GAUSS_2

namespace Kratos
{

	// =====================================================================================
    //
    // Utilties
    //
    // =====================================================================================

	namespace Utilities
	{

		template<class TVec>
		inline void ShapeFunc(double xi, double eta, TVec & N)
		{
			N(0) = 0.25 * (1.0 - xi) * (1.0 - eta); // node 1
			N(1) = 0.25 * (1.0 + xi) * (1.0 - eta); // node 2
			N(2) = 0.25 * (1.0 + xi) * (1.0 + eta); // node 3
			N(3) = 0.25 * (1.0 - xi) * (1.0 + eta); // node 4
		}

		template<class TMat>
		inline void ShapeFunc_NaturalDerivatives(double xi, double eta, TMat & dN)
		{
			dN(0, 0) = -(1.0 - eta) * 0.25;
			dN(1, 0) =  (1.0 - eta) * 0.25;
			dN(2, 0) =  (1.0 + eta) * 0.25;
			dN(3, 0) = -(1.0 + eta) * 0.25;

			dN(0, 1) = -(1.0 - xi)  * 0.25;
			dN(1, 1) = -(1.0 + xi)  * 0.25;
			dN(2, 1) =  (1.0 + xi)  * 0.25;
			dN(3, 1) =  (1.0 - xi)  * 0.25;
		}

		template<class TVec>
		inline void ShapeFunc_Drill(double xi, double eta, TVec & N)
		{
			N(0) = 0.5 * (1.0 - xi*xi) * (1.0 - eta);     // mid-side node 5
			N(1) = 0.5 * (1.0 + xi)    * (1.0 - eta*eta); // mid-side node 6
			N(2) = 0.5 * (1.0 - xi*xi) * (1.0 + eta);     // mid-side node 7
			N(3) = 0.5 * (1.0 - xi)    * (1.0 - eta*eta); // mid-side node 8
		}

		template<class TMat>
		inline void ShapeFunc_NaturalDerivatives_Drill(double xi, double eta, TMat & dN)
		{
			dN(0, 0) = xi*(eta - 1.0);
			dN(1, 0) = 0.5 - eta*eta*0.5;
			dN(2, 0) = -xi*(eta + 1.0);
			dN(3, 0) = eta*eta*0.5 - 0.5;

			dN(0, 1) = xi*xi*0.5 - 0.5;
			dN(1, 1) = -eta*(xi + 1.0);
			dN(2, 1) = 0.5 - xi*xi*0.5;
			dN(3, 1) = eta*(xi - 1.0);
		}

	}

	// =====================================================================================
    //
    // Class JacobianOperator
    //
    // =====================================================================================

	EASQuadElementV2::JacobianOperator::JacobianOperator()
        : mJac(2, 2, 0.0)
        , mInv(2, 2, 0.0)
        , mXYDeriv(4, 2, 0.0)
        , mDet(0.0)
    {
    }

    void EASQuadElementV2::JacobianOperator::Calculate(const Element::GeometryType & geom, const Matrix & dN)
    {
		const NodeType& p1 = geom[0];
		const NodeType& p2 = geom[1];
		const NodeType& p3 = geom[2];
		const NodeType& p4 = geom[3];

        mJac(0, 0) = dN(0, 0) * p1.X0() + dN(1, 0) * p2.X0() + dN(2, 0) * p3.X0() + dN(3, 0) * p4.X0();
        mJac(0, 1) = dN(0, 0) * p1.Y0() + dN(1, 0) * p2.Y0() + dN(2, 0) * p3.Y0() + dN(3, 0) * p4.Y0();
        mJac(1, 0) = dN(0, 1) * p1.X0() + dN(1, 1) * p2.X0() + dN(2, 1) * p3.X0() + dN(3, 1) * p4.X0();
        mJac(1, 1) = dN(0, 1) * p1.Y0() + dN(1, 1) * p2.Y0() + dN(2, 1) * p3.Y0() + dN(3, 1) * p4.Y0();

        mDet = mJac(0, 0) * mJac(1, 1) - mJac(1, 0) * mJac(0, 1);
        double mult = 1.0 / mDet;

        mInv(0, 0) =   mJac(1, 1) * mult;
        mInv(0, 1) = - mJac(0, 1) * mult;
        mInv(1, 0) = - mJac(1, 0) * mult;
        mInv(1, 1) =   mJac(0, 0) * mult;

        noalias( mXYDeriv ) = prod( dN, trans( mInv ) );
    }

	// =====================================================================================
    //
    // Class EASOperatorStorage
    //
    // =====================================================================================

	EASQuadElementV2::EASOperatorStorage::EASOperatorStorage()
		: alpha(EASQ4_NE, 0.0)
		, alpha_converged(EASQ4_NE, 0.0)
		, displ(EASQ4_NU, 0.0)
		, displ_converged(EASQ4_NU, 0.0)
		, residual(EASQ4_NE, 0.0)
		, Hinv(EASQ4_NE, EASQ4_NE, 0.0)
		, L(EASQ4_NE, EASQ4_NU, 0.0)
		, LT(EASQ4_NU, EASQ4_NE, 0.0)
		, mInitialized(false)
	{
	}

	void EASQuadElementV2::EASOperatorStorage::Initialize(const GeometryType& geom)
	{
		if(!mInitialized)
		{
			alpha.clear();
			alpha_converged.clear();

			for(size_t i = 0; i < geom.PointsNumber(); i++)
			{
#ifdef EASQ4_DRILLING_DOF
				size_t ii = i * 3;
				const array_1d<double, 3>& initialDispl = geom[i].FastGetSolutionStepValue(DISPLACEMENT);
				const array_1d<double, 3>& initialRot   = geom[i].FastGetSolutionStepValue(ROTATION);
				displ(ii    ) = initialDispl(0);
				displ(ii + 1) = initialDispl(1);
				displ(ii + 2) = initialRot(2);
				displ_converged(ii    ) = initialDispl(0);
				displ_converged(ii + 1) = initialDispl(1);
				displ_converged(ii + 2) = initialRot(2);
#else
				size_t ii = i * 2;
				const array_1d<double, 3>& initialDispl = geom[i].FastGetSolutionStepValue(DISPLACEMENT);
				displ(ii    ) = initialDispl(0);
				displ(ii + 1) = initialDispl(1);
				displ_converged(ii    ) = initialDispl(0);
				displ_converged(ii + 1) = initialDispl(1);
#endif // EASQ4_DRILLING_DOF
			}

			mInitialized = true;
		}
	}

	void EASQuadElementV2::EASOperatorStorage::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		displ = displ_converged;
		alpha = alpha_converged;
	}

	void EASQuadElementV2::EASOperatorStorage::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		displ_converged = displ;
		alpha_converged = alpha;
	}

	void EASQuadElementV2::EASOperatorStorage::FinalizeNonLinearIteration(const Vector& displacementVector, ProcessInfo& CurrentProcessInfo)
	{
		Vector incrementalDispl(EASQ4_NU);
		noalias(incrementalDispl) = displacementVector - displ;
		noalias(displ) = displacementVector;

		array_1d<double, EASQ4_NE> temp;
		noalias(temp) = prod(L, incrementalDispl);
		noalias(temp) -= residual;
		noalias(alpha) -= prod(Hinv, temp);
	}

	// =====================================================================================
    //
    // Class EASOperator
    //
    // =====================================================================================

	EASQuadElementV2::EASOperator::EASOperator(const GeometryType& geom, EASOperatorStorage& storage)
		: mF0inv(EASQ4_NSTRAIN, EASQ4_NSTRAIN)
		, mEnhancedStrains(EASQ4_NSTRAIN)
		, mG(EASQ4_NSTRAIN, EASQ4_NE)
	{

		// compute the jacobian at the element center

		double xi(0.0);
		double eta(0.0);

		Matrix dN(4, 2);
		Utilities::ShapeFunc_NaturalDerivatives(xi, eta, dN);

		const NodeType& p1 = geom[0];
		const NodeType& p2 = geom[1];
		const NodeType& p3 = geom[2];
		const NodeType& p4 = geom[3];
		Matrix Jac0(2, 2);
        Jac0(0, 0) = dN(0, 0) * p1.X0() + dN(1, 0) * p2.X0() + dN(2, 0) * p3.X0() + dN(3, 0) * p4.X0();
        Jac0(0, 1) = dN(0, 0) * p1.Y0() + dN(1, 0) * p2.Y0() + dN(2, 0) * p3.Y0() + dN(3, 0) * p4.Y0();
        Jac0(1, 0) = dN(0, 1) * p1.X0() + dN(1, 1) * p2.X0() + dN(2, 1) * p3.X0() + dN(3, 1) * p4.X0();
        Jac0(1, 1) = dN(0, 1) * p1.Y0() + dN(1, 1) * p2.Y0() + dN(2, 1) * p3.Y0() + dN(3, 1) * p4.Y0();

		// save the jacobian determinant at center

		mJ0 = Jac0(0, 0) * Jac0(1, 1) - Jac0(1, 0) * Jac0(0, 1);

		// compute the transformation matrix used in the implementation of the EAS method
		// which operates in the natural coordinate system

		double j11 = Jac0(0,0);
		double j22 = Jac0(1,1);
		double j12 = Jac0(0,1);
		double j21 = Jac0(1,0);

		Matrix F0(3, 3);
		F0(0,0) = j11*j11;        F0(0,1) = j21*j12;        F0(0,2) = 2.0*j11*j12;
		F0(1,0) = j12*j21;        F0(1,1) = j22*j22;        F0(1,2) = 2.0*j21*j22;
		F0(2,0) = j11*j21;        F0(2,1) = j12*j22;        F0(2,2) = j11*j22 + j12*j21;

		//noalias(F0) = IdentityMatrix(3,3);

		// invert the transformation matrix

		double F0det;
		//MathUtils<double>::InvertMatrix3(F0, mF0inv, F0det);
		mF0inv(0,0) = F0(1,1)*F0(2,2) - F0(1,2)*F0(2,1);
		mF0inv(1,0) = -F0(1,0)*F0(2,2) + F0(1,2)*F0(2,0);
		mF0inv(2,0) = F0(1,0)*F0(2,1) - F0(1,1)*F0(2,0);
		mF0inv(0,1) = -F0(0,1)*F0(2,2) + F0(0,2)*F0(2,1);
		mF0inv(1,1) = F0(0,0)*F0(2,2) - F0(0,2)*F0(2,0);
		mF0inv(2,1) = -F0(0,0)*F0(2,1) + F0(0,1)*F0(2,0);
		mF0inv(0,2) = F0(0,1)*F0(1,2) - F0(0,2)*F0(1,1);
		mF0inv(1,2) = -F0(0,0)*F0(1,2) + F0(0,2)*F0(1,0);
		mF0inv(2,2) = F0(0,0)*F0(1,1) - F0(0,1)*F0(1,0);
		F0det = F0(0,0)*mF0inv(0,0) + F0(0,1)*mF0inv(1,0) + F0(0,2)*mF0inv(2,0);
		mF0inv /= F0det;
#ifdef EASQ4_DRILLING_DOF
		mF0inv(3,0) = 0.0;
		mF0inv(3,1) = 0.0;
		mF0inv(3,2) = 0.0;
		mF0inv(0,3) = 0.0;
		mF0inv(1,3) = 0.0;
		mF0inv(2,3) = 0.0;
		mF0inv(3,3) = 1.0;
#endif

		// initialize these data to zero because they will
		// be integrated during the gauss loop
		storage.L.clear();
		storage.LT.clear();
		storage.Hinv.clear();
		storage.residual.clear();
	}

	void EASQuadElementV2::EASOperator::GaussPointComputation_Step1(double xi, double eta, const JacobianOperator& jac,
			                                                             Vector& strainVector,
			                                                             EASOperatorStorage& storage)
	{
		// construct the interpolation matrix in natural coordinate system
		Matrix E(EASQ4_NSTRAIN, EASQ4_NE, 0.0);
#if   EASQ4_NE == 2
		E(2,0) = xi;
		E(2,1) = eta;
#ifdef EASQ4_DRILLING_DOF
		/*E(3,2) = jac.Determinant()/(2.0*mJ0*mJ0)*(xi);
		E(3,3) = jac.Determinant()/(2.0*mJ0*mJ0)*(-eta);*/
#endif
#elif EASQ4_NE == 3
		E(0,0) = 0.0;
		E(1,0) = xi;
		E(2,1) = xi;
		E(2,2) = eta;
#elif EASQ4_NE == 4
		// Taylor
		E(0,0) = xi ;
		E(1,1) = eta;
		E(2,2) = xi;
		E(2,3) = eta;
#ifdef EASQ4_DRILLING_DOF
		/*E(3,2) = jac.Determinant()/(2.0*mJ0*mJ0)*(xi);
		E(3,3) = jac.Determinant()/(2.0*mJ0*mJ0)*(-eta);*/
#endif
#elif EASQ4_NE == 5
		// Simo and Rifai
		E(0,0) = xi;
		E(1,1) = eta;
		E(2,2) = xi;
		E(2,3) = eta;
		E(0,4) = xi*eta;
		E(1,4) = -xi*eta;
		E(2,4) = xi*xi - eta*eta;
#ifdef EASQ4_DRILLING_DOF
		E(3,2) = jac.Determinant()/(2.0*mJ0*mJ0)*(xi);
		E(3,3) = jac.Determinant()/(2.0*mJ0*mJ0)*(-eta);
		E(3,4) = jac.Determinant()/(2.0*mJ0*mJ0)*(xi*xi + eta*eta);
#endif
#elif EASQ4_NE == 6
		// Pantuso and Bathe
		E(0,0) = xi;
		E(1,1) = eta;
		E(2,2) = xi;
		E(2,3) = eta;
		E(0,4) = xi*eta;
		E(1,5) = xi*eta;
#elif EASQ4_NE == 7
		// Ramm
		E(0,0) = xi;
		E(1,1) = eta;
		E(2,2) = xi;
		E(2,3) = eta;
		E(0,4) = xi*eta;
		E(1,5) = xi*eta;
		E(2,6) = xi*eta;
#endif

		// construct the interpolation matrix in local coordinate system
		noalias( mG ) = mJ0 / jac.Determinant() * prod( mF0inv, E );

		// assuming that the input generalized strains has been already calculated on this
		// integration point, we just need to add the contribution of the enhanced strains
		noalias( mEnhancedStrains ) = prod( mG, storage.alpha );


#ifdef SUPPRESS_EAS_ENHANCEMENT
		return;
#endif // SUPPRESS_EAS_ENHANCEMENT



		strainVector(0) += mEnhancedStrains(0); // e.xx
		strainVector(1) += mEnhancedStrains(1); // e.yy
		strainVector(2) += mEnhancedStrains(2); // e.xy
#ifdef EASQ4_DRILLING_DOF
		strainVector(3) += mEnhancedStrains(3); // e.xy
#endif
	}

	void EASQuadElementV2::EASOperator::GaussPointComputation_Step2(const Matrix& D,
			                                                             const Matrix& B,
			                                                             const Vector& S,
			                                                             EASOperatorStorage& storage)
	{
#ifdef SUPPRESS_EAS_ENHANCEMENT
		return;
#endif // SUPPRESS_EAS_ENHANCEMENT


		Matrix GTC(EASQ4_NE, EASQ4_NSTRAIN);
		Matrix CG(EASQ4_NSTRAIN, EASQ4_NE);
		noalias( GTC )                = prod( trans( mG ), D );
		noalias( CG )                 = prod( D, mG );
		noalias( storage.Hinv )      += prod( GTC, mG );
		noalias( storage.residual )  -= prod( trans( mG ), S );
		noalias( storage.L )         += prod( GTC, B );
		noalias( storage.LT )        += prod( trans( B ), CG );
	}

	bool EASQuadElementV2::EASOperator::ComputeModfiedTangentAndResidual(Matrix& rLeftHandSideMatrix,
			                                                                  Vector& rRightHandSideVector,
			                                                                  EASOperatorStorage& storage)
	{
#ifdef SUPPRESS_EAS_ENHANCEMENT
		return true;
#endif // SUPPRESS_EAS_ENHANCEMENT



		// invert H
		/*std::cout << "H\n";
		std::cout << MathHelpers::MatrixToString(storage.Hinv, 3, std::scientific);*/

		Matrix Hcopy( storage.Hinv );
		permutation_matrix<Matrix::size_type> pm( EASQ4_NE );
		size_t lu_res = lu_factorize( Hcopy, pm );
		noalias( storage.Hinv ) = IdentityMatrix( EASQ4_NE );
		lu_substitute( Hcopy, pm, storage.Hinv );

		/*std::cout << "H^-1\n";
		std::cout << MathHelpers::MatrixToString(storage.Hinv, 3, std::scientific);*/

		// compute L' * H^-1
		Matrix LTHinv(EASQ4_NU, EASQ4_NE);
		/*noalias( LTHinv ) = prod( trans( storage.L ), storage.Hinv );*/
		noalias( LTHinv ) = prod( storage.LT, storage.Hinv );

		// modify the stiffness matrix and the residual vector
		// for the static condensation of the enhanced strain parameters
		noalias( rRightHandSideVector ) -= prod( LTHinv, storage.residual );   // R_mod = R - L' * H^-1 * residual
		noalias( rLeftHandSideMatrix  ) -= prod( LTHinv, storage.L );          // K_mod = K - L' * H^-1 * L

		//KRATOS_WATCH(storage.residual);

		return (lu_res == 0);
	}

	// =====================================================================================
    //
    // Class EASQuadElementV2
    //
    // =====================================================================================

    EASQuadElementV2::EASQuadElementV2(IndexType NewId,
                                   GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
		mThisIntegrationMethod = EASQ4_DEF_INTEGRATION_METHOD;
    }

    EASQuadElementV2::EASQuadElementV2(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
		mThisIntegrationMethod = EASQ4_DEF_INTEGRATION_METHOD;
    }

    EASQuadElementV2::~EASQuadElementV2()
    {
    }

    Element::Pointer EASQuadElementV2::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        GeometryType::Pointer newGeom( GetGeometry().Create(ThisNodes) );
        return Element::Pointer( new EASQuadElementV2(NewId, newGeom, pProperties) );
    }

    void EASQuadElementV2::Initialize()
    {
		KRATOS_TRY

		// NOTE:
		// This is the standard (previous) implementation:
		// If we are here, it means that no one already set up the constitutive law vector
		// through the method SetValue<CONSTITUTIVE_LAW_POINTER>

		const GeometryType & geom = GetGeometry();

		const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints( mThisIntegrationMethod );


		mEASStorage.Initialize(geom);

#ifdef EASQ4_DRILLING_DOF
		mG0 = 0.0;
		m_first_step_finalized = false;
#endif

		//Constitutive Law initialization
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

    void EASQuadElementV2::ResetConstitutiveLaw()
    {
        KRATOS_TRY

        const GeometryType & geom = GetGeometry();
        const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

        const Properties& props = GetProperties();
        for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
            mConstitutiveLawVector[i]->ResetMaterial(props, geom, row(shapeFunctionsValues, i));

        KRATOS_CATCH("")
    }

	Kratos::GeometryData::IntegrationMethod EASQuadElementV2::GetIntegrationMethod() const
    {
        return mThisIntegrationMethod;
    }

    void EASQuadElementV2::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        if(rResult.size() != EASQ4_NU)
            rResult.resize(EASQ4_NU, false);

        GeometryType & geom = this->GetGeometry();

        for(int i = 0; i < geom.PointsNumber(); i++)
        {
#ifdef EASQ4_DRILLING_DOF
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
#endif // EASQ4_DRILLING_DOF
        }
    }

    void EASQuadElementV2::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        ElementalDofList.resize(0);
        ElementalDofList.reserve(EASQ4_NU);

        GeometryType & geom = this->GetGeometry();

        for (int i = 0; i < geom.PointsNumber(); i++)
        {
            NodeType & iNode = geom[i];

            ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Y));
#ifdef EASQ4_DRILLING_DOF
			ElementalDofList.push_back(iNode.pGetDof(ROTATION_Z));
#endif
        }
    }

    int EASQuadElementV2::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& geom = GetGeometry();

        // verify that the variables are correctly initialized
        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered","");
#ifdef EASQ4_DRILLING_DOF
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
#ifdef EASQ4_DRILLING_DOF
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

    void EASQuadElementV2::CleanMemory()
    {
    }

    void EASQuadElementV2::GetValuesVector(Vector& values, int Step)
    {
        if(values.size() != EASQ4_NU)
            values.resize(EASQ4_NU, false);

        const GeometryType & geom = GetGeometry();

        for (int i = 0; i < geom.PointsNumber(); i++)
        {
            const NodeType & iNode = geom[i];

#ifdef EASQ4_DRILLING_DOF
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
#endif // EASQ4_DRILLING_DOF
        }
    }

    void EASQuadElementV2::GetFirstDerivativesVector(Vector& values, int Step)
    {
        if(values.size() != EASQ4_NU)
            values.resize(EASQ4_NU,false);

        const GeometryType & geom = GetGeometry();

        for (int i = 0; i < geom.PointsNumber(); i++)
        {
            const NodeType & iNode = geom[i];
			const array_1d<double,3>& vel = iNode.FastGetSolutionStepValue(VELOCITY, Step);

#ifdef EASQ4_DRILLING_DOF
			int index = i * 3;
			values[index]        = vel[0];
			values[index + 1]    = vel[1];
			values[index + 2]    = 0.0;
#else
			int index = i * 2;
            values[index]        = vel[0];
            values[index + 1]    = vel[1];
#endif // EASQ4_DRILLING_DOF

        }
    }

    void EASQuadElementV2::GetSecondDerivativesVector(Vector& values, int Step)
    {
        if(values.size() != EASQ4_NU)
            values.resize(EASQ4_NU,false);

        const GeometryType & geom = GetGeometry();

        for (int i = 0; i < geom.PointsNumber(); i++)
        {
            const NodeType & iNode = geom[i];
			const array_1d<double,3>& acc = iNode.FastGetSolutionStepValue(ACCELERATION, Step);

#ifdef EASQ4_DRILLING_DOF
			int index = i * 3;
			values[index]        = acc[0];
			values[index + 1]    = acc[1];
			values[index + 2]    = 0.0;
#else
			int index = i * 2;
            values[index]        = acc[0];
            values[index + 1]    = acc[1];
#endif // EASQ4_DRILLING_DOF

        }
    }

    void EASQuadElementV2::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
    }

    void EASQuadElementV2::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
		Vector U(EASQ4_NU);
		GetValuesVector(U);
		mEASStorage.FinalizeNonLinearIteration(U, CurrentProcessInfo);
    }

    void EASQuadElementV2::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
		const PropertiesType& props = GetProperties();
        const GeometryType & geom = GetGeometry();
        const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

        for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
            mConstitutiveLawVector[i]->InitializeSolutionStep(props, geom, row(shapeFunctionsValues, i), CurrentProcessInfo);

        mEASStorage.InitializeSolutionStep(CurrentProcessInfo);
    }

    void EASQuadElementV2::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
		 const PropertiesType& props = GetProperties();
        const GeometryType& geom = GetGeometry();
        const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

        for(SizeType i = 0; i < mConstitutiveLawVector.size(); i++)
            mConstitutiveLawVector[i]->FinalizeSolutionStep(props, geom, row(shapeFunctionsValues, i), CurrentProcessInfo);

        mEASStorage.FinalizeSolutionStep(CurrentProcessInfo);

#ifdef EASQ4_DRILLING_DOF
		if(!m_first_step_finalized)
			m_first_step_finalized = true;
#endif
    }

    void EASQuadElementV2::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        if((rMassMatrix.size1() != EASQ4_NU) || (rMassMatrix.size2() != EASQ4_NU))
            rMassMatrix.resize(EASQ4_NU, EASQ4_NU, false);
        noalias(rMassMatrix) = ZeroMatrix(EASQ4_NU, EASQ4_NU);

        // Compute the total mass

		size_t num_nodes = GetGeometry().PointsNumber();
		double total_mass = GetGeometry().DomainSize() * GetProperties()[DENSITY] * GetProperties()[THICKNESS];
        double lumped_mass = total_mass / double(num_nodes);

		// loop on nodes
        for(size_t i = 0; i < num_nodes; i++)
        {
#ifdef EASQ4_DRILLING_DOF
			size_t index = i * 3;
			rMassMatrix(index, index)            = lumped_mass;
			rMassMatrix(index + 1, index + 1)    = lumped_mass;
			rMassMatrix(index + 2, index + 2)    = 0.0;
#else
			size_t index = i * 2;
            rMassMatrix(index, index)            = lumped_mass;
            rMassMatrix(index + 1, index + 1)    = lumped_mass;
#endif // EASQ4_DRILLING_DOF

        }
    }

    void EASQuadElementV2::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        if((rDampingMatrix.size1() != EASQ4_NU) || (rDampingMatrix.size2() != EASQ4_NU))
            rDampingMatrix.resize(EASQ4_NU, EASQ4_NU, false);

        noalias( rDampingMatrix ) = ZeroMatrix(EASQ4_NU, EASQ4_NU);
    }

    void EASQuadElementV2::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                     VectorType& rRightHandSideVector,
                                                     ProcessInfo& rCurrentProcessInfo)
    {
        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
    }

    void EASQuadElementV2::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                       ProcessInfo& rCurrentProcessInfo)
    {
        Matrix dummy;
        CalculateAll(dummy, rRightHandSideVector, rCurrentProcessInfo, true, true);
    }

    // =====================================================================================
    //
    // Class EASQuadElementV2 - Results on Gauss Points
    //
    // =====================================================================================

	void EASQuadElementV2::SetValueOnIntegrationPoints(const Variable<double>& rVariable,
													   std::vector<double>& rValues,
													   const ProcessInfo& rCurrentProcessInfo)
	{
#if   EASQ4_NE == 2
		if(rVariable == ENH_STRAIN_PARAM_1) {
			this->mEASStorage.alpha_converged(0) = rValues[0];
			return;
		}
		else if(rVariable == ENH_STRAIN_PARAM_2) {
			this->mEASStorage.alpha_converged(0) = rValues[1];
			return;
		}
#elif EASQ4_NE == 4
		if(rVariable == ENH_STRAIN_PARAM_1) {
			this->mEASStorage.alpha_converged(0) = rValues[0];
			return;
		}
		else if(rVariable == ENH_STRAIN_PARAM_2) {
			this->mEASStorage.alpha_converged(0) = rValues[1];
			return;
		}
		else if(rVariable == ENH_STRAIN_PARAM_3) {
			this->mEASStorage.alpha_converged(0) = rValues[2];
			return;
		}
		else if(rVariable == ENH_STRAIN_PARAM_4) {
			this->mEASStorage.alpha_converged(0) = rValues[3];
			return;
		}
#endif

		if(rValues.size() == mConstitutiveLawVector.size())
			for(unsigned int i = 0; i < rValues.size(); i++)
				mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
	}

    void EASQuadElementV2::SetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
													   std::vector<Vector>& rValues,
													   const ProcessInfo& rCurrentProcessInfo)
	{
		if(rValues.size() == mConstitutiveLawVector.size())
			for(unsigned int i = 0; i < rValues.size(); i++)
				mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
	}

    void EASQuadElementV2::SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
													   std::vector<Matrix>& rValues,
													   const ProcessInfo& rCurrentProcessInfo)
	{
		if(rValues.size() == mConstitutiveLawVector.size())
			for(unsigned int i = 0; i < rValues.size(); i++)
				mConstitutiveLawVector[i]->SetValue(rVariable, rValues[i], rCurrentProcessInfo);
	}

	void EASQuadElementV2::SetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
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

    void EASQuadElementV2::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                                     std::vector<double>& rValues,
                                                     const ProcessInfo& rCurrentProcessInfo)
    {
		unsigned int numgp = mConstitutiveLawVector.size();
        if(rValues.size() != numgp)
            rValues.resize(numgp);

		std::fill(rValues.begin(), rValues.end(), 0.0);

		if(rVariable == TEMPERATURE)
		{
			Vector U(EASQ4_NU);
			GetValuesVector(U, 0);
			PropertiesType & props = GetProperties();
			GeometryType & geom = GetGeometry();
			const Matrix & shapeFunctions = geom.ShapeFunctionsValues(mThisIntegrationMethod);
			Vector iN(shapeFunctions.size2());
			JacobianOperator jacOp;      // jacobian operator
			Matrix B(EASQ4_NSTRAIN, EASQ4_NU, 0.0);  // strain-displacement matrix.
			Matrix D(EASQ4_NSTRAIN, EASQ4_NSTRAIN, 0.0);         // material tangent matrix.
			Vector E(EASQ4_NSTRAIN);                 // strain vector
			Vector S(EASQ4_NSTRAIN);                 // stress vector
			Matrix Du(3, 3, 0.0);         // material tangent matrix.
			Vector Eu(3);                 // strain vector
			Vector Su(3);                 // stress vector
			EASOperator EASOp(geom, mEASStorage);
			ConstitutiveLaw::Parameters parameters(geom, props, rCurrentProcessInfo);
			parameters.SetStrainVector( Eu );
			parameters.SetStressVector( Su );
			parameters.SetConstitutiveMatrix( Du );
			Flags& options = parameters.GetOptions();
			options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
			options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
			double detF = 1.0;
			Matrix F(IdentityMatrix(2,2));
			parameters.SetDeterminantF(detF);
			parameters.SetDeformationGradientF(F);
			const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(mThisIntegrationMethod);
			const GeometryType::ShapeFunctionsGradientsType& local_gradients = geom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
			for(size_t i = 0; i < numgp; i++)
			{
				const GeometryType::IntegrationPointType & ip = integration_points[i];
				noalias( iN ) = row( shapeFunctions, i );
				jacOp.Calculate(geom, local_gradients[i]);
				CalculateBMatrix( ip.X(), ip.Y(), jacOp, iN, B );
				noalias( E ) = prod( B, U );
				rValues[i] = E(3); // drilling strain = skew-grad u
			}
			return;
		}

		if(mErrorCode != 0.0)
		{
			if(rVariable == CONSTITUTIVE_INTEGRATION_ERROR_CODE)
			{
				std::fill(rValues.begin(), rValues.end(), mErrorCode);
			}
		}

#if   EASQ4_NE == 2
		if(rVariable == ENH_STRAIN_PARAM_1) {
			std::fill(rValues.begin(), rValues.end(), this->mEASStorage.alpha_converged(0));
			return;
		}
		else if(rVariable == ENH_STRAIN_PARAM_2) {
			std::fill(rValues.begin(), rValues.end(), this->mEASStorage.alpha_converged(1));
			return;
		}
#elif EASQ4_NE == 4
		if(rVariable == ENH_STRAIN_PARAM_1) {
			std::fill(rValues.begin(), rValues.end(), this->mEASStorage.alpha_converged(0));
			return;
		}
		else if(rVariable == ENH_STRAIN_PARAM_2) {
			std::fill(rValues.begin(), rValues.end(), this->mEASStorage.alpha_converged(1));
			return;
		}
		else if(rVariable == ENH_STRAIN_PARAM_3) {
			std::fill(rValues.begin(), rValues.end(), this->mEASStorage.alpha_converged(2));
			return;
		}
		else if(rVariable == ENH_STRAIN_PARAM_4) {
			std::fill(rValues.begin(), rValues.end(), this->mEASStorage.alpha_converged(3));
			return;
		}
#endif

        for(int i = 0; i < numgp; i++) {
			rValues[i] = 0.0;
            mConstitutiveLawVector[i]->GetValue(rVariable, rValues[i]);
		}
    }

    void EASQuadElementV2::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
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

    void EASQuadElementV2::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
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
			// Get some references.
			PropertiesType & props = GetProperties();
			GeometryType & geom = GetGeometry();
			const Matrix & shapeFunctions = geom.ShapeFunctionsValues(mThisIntegrationMethod);
			Vector iN(shapeFunctions.size2());

			// some data
			JacobianOperator jacOp;      // jacobian operator
			Matrix B(EASQ4_NSTRAIN, EASQ4_NU, 0.0);  // strain-displacement matrix.
			Matrix D(EASQ4_NSTRAIN, EASQ4_NSTRAIN, 0.0);         // material tangent matrix.
			Vector E(EASQ4_NSTRAIN);                 // strain vector
			Vector S(EASQ4_NSTRAIN);                 // stress vector
#ifdef EASQ4_DRILLING_DOF
			Matrix Du(3, 3, 0.0);         // material tangent matrix.
			Vector Eu(3);                 // strain vector
			Vector Su(3);                 // stress vector
#endif

			// Get the current displacements in global coordinate system
			Vector U(EASQ4_NU);
			GetValuesVector(U, 0);

			// Instantiate the EAS Operator.
			EASOperator EASOp(geom, mEASStorage);

			// Initialize parameters for the material calculation
			ConstitutiveLaw::Parameters parameters(geom, props, rCurrentProcessInfo);
#ifdef EASQ4_DRILLING_DOF
			parameters.SetStrainVector( Eu );
			parameters.SetStressVector( Su );
			parameters.SetConstitutiveMatrix( Du );
#else
			parameters.SetStrainVector( E );
			parameters.SetStressVector( S );
			parameters.SetConstitutiveMatrix( D );
#endif // EASQ4_DRILLING_DOF
			Flags& options = parameters.GetOptions();
			options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
			options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
			double detF = 1.0;
			Matrix F(IdentityMatrix(2,2));
			parameters.SetDeterminantF(detF);
			parameters.SetDeformationGradientF(F);

			// gauss loop
			const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(mThisIntegrationMethod);
			const GeometryType::ShapeFunctionsGradientsType& local_gradients = geom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
			for(size_t i = 0; i < numgp; i++)
			{
				Matrix& output = rValues[i];
				if(output.size1() != 2 || output.size2() != 2)
					output.resize(2,2,false);

				// get a reference of the current integration point and shape functions
				const GeometryType::IntegrationPointType & ip = integration_points[i];
				noalias( iN ) = row( shapeFunctions, i );

				// Compute Jacobian, Inverse of Jacobian, Determinant of Jacobian
				// and Shape functions derivatives in the local coordinate system
				jacOp.Calculate(geom, local_gradients[i]);

				// Compute all strain-displacement matrices
				CalculateBMatrix( ip.X(), ip.Y(), jacOp, iN, B );

				// Calculate strain vector
				noalias( E ) = prod( B, U );

				// Apply the EAS method.
				EASOp.GaussPointComputation_Step1(ip.X(), ip.Y(), jacOp, E, mEASStorage);

				if(rVariable == GREEN_LAGRANGE_STRAIN_TENSOR)
				{
#ifdef EASQ4_DRILLING_DOF
					Eu(0) = E(0);
					Eu(1) = E(1);
					Eu(2) = E(2);
					noalias(output) = MathUtils<double>::StrainVectorToTensor(Eu);
#else
					noalias(output) = MathUtils<double>::StrainVectorToTensor(E);
#endif
				}
				else
				{
#ifdef EASQ4_DRILLING_DOF
					Eu(0) = E(0);
					Eu(1) = E(1);
					Eu(2) = E(2);
#endif
					ConstitutiveLaw::Pointer & mat = mConstitutiveLawVector[i];
					parameters.SetShapeFunctionsValues( iN );
					parameters.SetShapeFunctionsDerivatives( jacOp.XYDerivatives() );
					mat->CalculateMaterialResponseCauchy( parameters );
#ifdef EASQ4_DRILLING_DOF
					noalias(output) = MathUtils<double>::StressVectorToTensor(Su);
#else
					noalias(output) = MathUtils<double>::StressVectorToTensor(S);
#endif
				}
			}
		}
		else
		{
			for(size_t i = 0; i < mConstitutiveLawVector.size(); i++)
				mConstitutiveLawVector[i]->GetValue(rVariable, rValues[i]);
		}
	}

	void EASQuadElementV2::GetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
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
    // Class EASQuadElementV2 - Private methods
    //
    // =====================================================================================

	void EASQuadElementV2::DecimalCorrection(Vector& a)
	{
		double norm = norm_2(a);
		double tolerance = std::max(norm * 1.0E-12, 1.0E-12);
		for(SizeType i = 0; i < a.size(); i++)
			if(std::abs(a(i)) < tolerance)
				a(i) = 0.0;
	}

	inline void drutils__calculate_indices(
		unsigned int i,
		unsigned int &m, unsigned int &l, unsigned int &j, unsigned int &k)
	{
		m = i+4;
		l = m-1+4*std::floor(1/(i+1));
		k = std::fmod(m+1,4);
		j = l-4;
		m-=4;
		l-=4;
	}

	inline void drutils__calculate_edges(
		const Element::GeometryType& geom,
		unsigned int i, unsigned int j, unsigned int k,
		double &xij, double &yij, double &xik, double &yik)
	{
		xij = geom[j].X0() - geom[i].X0();
		yij = geom[j].Y0() - geom[i].Y0();
		xik = geom[k].X0() - geom[i].X0();
		yik = geom[k].Y0() - geom[i].Y0();
	}

	void EASQuadElementV2::CalculateBMatrix(double xi, double eta, const JacobianOperator& Jac, const Vector& N, Matrix& B)
	{
		const Matrix& dNxy = Jac.XYDerivatives();

#ifdef EASQ4_DRILLING_DOF

		B(0,  0) =  dNxy(0, 0);   B(0, 3) =  dNxy(1, 0);   B(0, 6) =  dNxy(2, 0);   B(0, 9)  =  dNxy(3, 0);
		B(1,  1) =  dNxy(0, 1);   B(1, 4) =  dNxy(1, 1);   B(1, 7) =  dNxy(2, 1);   B(1, 10) =  dNxy(3, 1);
		B(2,  0) =  dNxy(0, 1);   B(2, 3) =  dNxy(1, 1);   B(2, 6) =  dNxy(2, 1);   B(2, 9)  =  dNxy(3, 1);
		B(2,  1) =  dNxy(0, 0);   B(2, 4) =  dNxy(1, 0);   B(2, 7) =  dNxy(2, 0);   B(2, 10) =  dNxy(3, 0);

#ifndef EASQ4_DRILLING_DOF_MOD

		/*const GeometryType& geom = GetGeometry();
		JacobianOperator jac_0;
		const GeometryType::IntegrationPointsArrayType& integration_points_0 = geom.IntegrationPoints(IntegrationMethod::GI_GAUSS_1);
		const GeometryType::ShapeFunctionsGradientsType& local_gradients_0 = geom.ShapeFunctionsLocalGradients(IntegrationMethod::GI_GAUSS_1);
		const Matrix & shape_functions_0 = geom.ShapeFunctionsValues(IntegrationMethod::GI_GAUSS_1);
		jac_0.Calculate(geom, local_gradients_0[0]);
		const Matrix & dNxy_0 = jac_0.XYDerivatives();

		B(3,  0) = -0.5*dNxy_0(0, 1);
		B(3,  1) =  0.5*dNxy_0(0, 0);
		B(3,  2) = -N(0);
		B(3,  3) = -0.5*dNxy_0(1, 1);
		B(3,  4) =  0.5*dNxy_0(1, 0);
		B(3,  5) = -N(1);
		B(3,  6) = -0.5*dNxy_0(2, 1);
		B(3,  7) =  0.5*dNxy_0(2, 0);
		B(3,  8) = -N(2);
		B(3,  9) = -0.5*dNxy_0(3, 1);
		B(3, 10) =  0.5*dNxy_0(3, 0);
		B(3, 11) = -N(3);*/

		B(3,  0) = -0.5*dNxy(0, 1);
		B(3,  1) =  0.5*dNxy(0, 0);
		B(3,  2) = -N(0);
		B(3,  3) = -0.5*dNxy(1, 1);
		B(3,  4) =  0.5*dNxy(1, 0);
		B(3,  5) = -N(1);
		B(3,  6) = -0.5*dNxy(2, 1);
		B(3,  7) =  0.5*dNxy(2, 0);
		B(3,  8) = -N(2);
		B(3,  9) = -0.5*dNxy(3, 1);
		B(3, 10) =  0.5*dNxy(3, 0);
		B(3, 11) = -N(3);

#else

		const GeometryType& geom = GetGeometry();

		unsigned int i,m,l,k,j;
		double xij,yij,xik,yik;

		// get shape function and their gradients for the drilling part
		Matrix dN_d(4,2);
		Matrix dNxy_d(4,2);
		Utilities::ShapeFunc_NaturalDerivatives_Drill(xi,eta,dN_d);
		noalias( dNxy_d ) = prod( dN_d, trans( Jac.Inverse() ) );

		// compute the B-matrix entries for the drilling part
		for(i = 0; i < 4; i++)
		{
			// drilling utilities
			drutils__calculate_indices(i,l,m,k,j);
			drutils__calculate_edges(geom,i,j,k,xij,yij,xik,yik);
			// Allman-type shape functions gradients
			double NXi_x = (yij*dNxy_d(l,0) - yik*dNxy_d(m,0))/8.0;
			double NXi_y = (yij*dNxy_d(l,1) - yik*dNxy_d(m,1))/8.0;
			double NYi_x = (xij*dNxy_d(l,0) - xik*dNxy_d(m,0))/8.0;
			double NYi_y = (xij*dNxy_d(l,1) - xik*dNxy_d(m,1))/8.0;
			if(i==0 || i==2) {
				NXi_x *= -1.0;
				NXi_y *= -1.0;
				NYi_x *= -1.0;
				NYi_y *= -1.0;
			}
			// Fill the membrane B-matrix with the drilling contribution
			unsigned int index = i*3+2;
			B(0, index) = NXi_x;
			B(1, index) = NYi_y;
			B(2, index) = NXi_y + NYi_x;

			/*B(3, index-2) = -0.5*dNxy(i, 1);
			B(3, index-1) =  0.5*dNxy(i, 0);
			B(3, index)   = -N(i) + (-yij*dNxy_d(l,1)+yik*dNxy_d(m,1)+xij*dNxy_d(l,0)-xik*dNxy_d(m,0))/16.0;*/
		}

		// Compute the correction to pass the constant stress patch test
		// [B_drilling_bar] = [B_drilling] - 1/V*int{ [B_drilling] dV }
		JacobianOperator i_jac;
		const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(mThisIntegrationMethod);
		const GeometryType::ShapeFunctionsGradientsType& local_gradients = geom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
		Matrix Bd_mean(3,4,0.0);
		Vector bi_mean(12,0.0);
		double V_tot(0.0);
		for(int igp = 0; igp < integration_points.size(); igp++)
		{
			// compute jacobian for this integration point
			const GeometryType::IntegrationPointType & ip = integration_points[igp];
			i_jac.Calculate(geom, local_gradients[igp]);
			double dV = ip.Weight() * i_jac.Determinant();
			V_tot += dV;
			// re-calculate the drilling shape functions and local gradients for this integration point
			Utilities::ShapeFunc_NaturalDerivatives_Drill(ip.X(),ip.Y(),dN_d);
			noalias( dNxy_d ) = prod( dN_d, trans( i_jac.Inverse() ) );
			for(i = 0; i < 4; i++)
			{
				// drilling utilities
				drutils__calculate_indices(i,l,m,k,j);
				drutils__calculate_edges(geom,i,j,k,xij,yij,xik,yik);
				// Allman-type shape functions gradients
				double NXi_x = (yij*dNxy_d(l,0) - yik*dNxy_d(m,0))/8.0;
				double NXi_y = (yij*dNxy_d(l,1) - yik*dNxy_d(m,1))/8.0;
				double NYi_x = (xij*dNxy_d(l,0) - xik*dNxy_d(m,0))/8.0;
				double NYi_y = (xij*dNxy_d(l,1) - xik*dNxy_d(m,1))/8.0;
				if(i==0 || i==2) {
					NXi_x *= -1.0;
					NXi_y *= -1.0;
					NYi_x *= -1.0;
					NYi_y *= -1.0;
				}
				// Fill the membrane B-matrix with the drilling contribution
				Bd_mean(0, i) += dV*(NXi_x);
				Bd_mean(1, i) += dV*(NYi_y);
				Bd_mean(2, i) += dV*(NXi_y + NYi_x);
			}
		}
		Bd_mean /= V_tot;
		for(i = 0; i < 4; i++)
		{
			// Fill the membrane B-matrix with the drilling contribution
			unsigned int index = i*3+2;
			B(0, index) -= Bd_mean(0, i);
			B(1, index) -= Bd_mean(1, i);
			B(2, index) -= Bd_mean(2, i);
		}

		// Compute the vector b rank-one update (Vol*gamma*[b]*[b]^T)
		JacobianOperator jac_0;
		const GeometryType::IntegrationPointsArrayType& integration_points_0 = geom.IntegrationPoints(IntegrationMethod::GI_GAUSS_1);
		const GeometryType::ShapeFunctionsGradientsType& local_gradients_0 = geom.ShapeFunctionsLocalGradients(IntegrationMethod::GI_GAUSS_1);
		const Matrix & shape_functions_0 = geom.ShapeFunctionsValues(IntegrationMethod::GI_GAUSS_1);
		jac_0.Calculate(geom, local_gradients_0[0]);
		const Matrix & dNxy_0 = jac_0.XYDerivatives();
		Matrix dN_d_0(4,2);
		Matrix dNxy_d_0(4,2);
		Utilities::ShapeFunc_NaturalDerivatives_Drill(0.0,0.0,dN_d_0);
		noalias( dNxy_d_0 ) = prod( dN_d_0, trans( jac_0.Inverse() ) );
		for(i = 0; i < 4; i++)
		{
			// drilling utilities
			drutils__calculate_indices(i,l,m,k,j);
			drutils__calculate_edges(geom,i,j,k,xij,yij,xik,yik);
			// Fill the 4th row of the B-matrix for the rank-one update to remove the hourglass drilling mode.
			unsigned int index = i*3+2;
			B(3, index-2) = -0.5*dNxy(i, 1);
			B(3, index-1) =  0.5*dNxy(i, 0);
			B(3, index)   = -N(i) + (-yij*dNxy_d_0(l,1)+yik*dNxy_d_0(m,1)+xij*dNxy_d_0(l,0)-xik*dNxy_d_0(m,0))/16.0;
		}

#endif // !EASQ4_DRILLING_DOF_MOD

#else

		B(0,  0) =  dNxy(0, 0);   B(0, 2) =  dNxy(1, 0);   B(0, 4) =  dNxy(2, 0);   B(0, 6) =  dNxy(3, 0);
		B(1,  1) =  dNxy(0, 1);   B(1, 3) =  dNxy(1, 1);   B(1, 5) =  dNxy(2, 1);   B(1, 7) =  dNxy(3, 1);
		B(2,  0) =  dNxy(0, 1);   B(2, 2) =  dNxy(1, 1);   B(2, 4) =  dNxy(2, 1);   B(2, 6) =  dNxy(3, 1);
		B(2,  1) =  dNxy(0, 0);   B(2, 3) =  dNxy(1, 0);   B(2, 5) =  dNxy(2, 0);   B(2, 7) =  dNxy(3, 0);

#endif // EASQ4_DRILLING_DOF

	}

    void EASQuadElementV2::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo,
                                      const bool LHSrequired,
                                      const bool RHSrequired)
    {
        // Resize the Left Hand Side if necessary,
        // and initialize it to Zero
        if((rLeftHandSideMatrix.size1() != EASQ4_NU) || (rLeftHandSideMatrix.size2() != EASQ4_NU))
            rLeftHandSideMatrix.resize(EASQ4_NU, EASQ4_NU, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(EASQ4_NU, EASQ4_NU);

        // Resize the Right Hand Side if necessary,
        // and initialize it to Zero
        if(rRightHandSideVector.size() != EASQ4_NU)
            rRightHandSideVector.resize(EASQ4_NU, false);
        noalias(rRightHandSideVector) = ZeroVector(EASQ4_NU);

		// set no error
		mErrorCode = 0.0;

        // Get some references.
        PropertiesType & props = GetProperties();
        GeometryType & geom = GetGeometry();
        const Matrix & shapeFunctions = geom.ShapeFunctionsValues(mThisIntegrationMethod);
		Vector iN(shapeFunctions.size2());

		// thickness
		double th = props[THICKNESS];

        // some data
        JacobianOperator jacOp;      // jacobian operator
        Matrix B(EASQ4_NSTRAIN, EASQ4_NU, 0.0);  // strain-displacement matrix.
		Matrix BTD(EASQ4_NU, EASQ4_NSTRAIN);     // auxiliary matrix to store the product B'*D
        Matrix D(EASQ4_NSTRAIN, EASQ4_NSTRAIN, 0.0);         // material tangent matrix.
        Vector E(EASQ4_NSTRAIN);                 // strain vector
        Vector S(EASQ4_NSTRAIN);                 // stress vector
#ifdef EASQ4_DRILLING_DOF
		Matrix Du(3, 3, 0.0);         // material tangent matrix.
        Vector Eu(3);                 // strain vector
        Vector Su(3);                 // stress vector
#endif

        // Get the current displacements in global coordinate system
        Vector U(EASQ4_NU);
        GetValuesVector(U, 0);

        // Instantiate the EAS Operator.
		EASOperator EASOp(geom, mEASStorage);

		// Initialize parameters for the material calculation
		ConstitutiveLaw::Parameters parameters(geom, props, rCurrentProcessInfo);
#ifdef EASQ4_DRILLING_DOF
		parameters.SetStrainVector( Eu );
		parameters.SetStressVector( Su );
		parameters.SetConstitutiveMatrix( Du );
#else
		parameters.SetStrainVector( E );
		parameters.SetStressVector( S );
		parameters.SetConstitutiveMatrix( D );
#endif // EASQ4_DRILLING_DOF
		Flags& options = parameters.GetOptions();
		options.Set(ConstitutiveLaw::COMPUTE_STRESS, RHSrequired);
		options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, LHSrequired);
		double detF = 1.0;
		Matrix F(IdentityMatrix(2,2));
		parameters.SetDeterminantF(detF);
		parameters.SetDeformationGradientF(F);

		// Gauss Loop.
		const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(mThisIntegrationMethod);
		const GeometryType::ShapeFunctionsGradientsType& local_gradients = geom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
        for(int i = 0; i < integration_points.size(); i++)
        {
            // get a reference of the current integration point and shape functions
            const GeometryType::IntegrationPointType & ip = integration_points[i];
			noalias( iN ) = row( shapeFunctions, i );

            // Compute Jacobian, Inverse of Jacobian, Determinant of Jacobian
            // and Shape functions derivatives in the local coordinate system
			jacOp.Calculate(geom, local_gradients[i]);

            // compute the 'volume' of the current integration point
            double dV = ip.Weight() * jacOp.Determinant() * th;

            // Compute all strain-displacement matrices
			CalculateBMatrix( ip.X(), ip.Y(), jacOp, iN, B );

            // Calculate strain vector
            noalias( E ) = prod( B, U );

			// Apply the EAS method.
            EASOp.GaussPointComputation_Step1(ip.X(), ip.Y(), jacOp, E, mEASStorage);

            // Calculate material response
#ifdef EASQ4_DRILLING_DOF
			Eu(0) = E(0);
			Eu(1) = E(1);
			Eu(2) = E(2);
#endif
			ConstitutiveLaw::Pointer & mat = mConstitutiveLawVector[i];
			parameters.SetShapeFunctionsValues( iN );
			parameters.SetShapeFunctionsDerivatives( jacOp.XYDerivatives() );
			mat->CalculateMaterialResponseCauchy( parameters );
#ifdef EASQ4_DRILLING_DOF
			if(LHSrequired)
			{
				if(!m_first_step_finalized)
					mG0 = Du(2,2);
			}
			double drilling_modulus = mG0*EASQ4_DRILLING_PENALTY_SCALE;
			D(0,0) = Du(0,0); D(0,1) = Du(0,1); D(0,2) = Du(0,2);
			D(1,0) = Du(1,0); D(1,1) = Du(1,1); D(1,2) = Du(1,2);
			D(2,0) = Du(2,0); D(2,1) = Du(2,1); D(2,2) = Du(2,2);
			D(3,3) = drilling_modulus;
			S(0) = Su(0);
			S(1) = Su(1);
			S(2) = Su(2);
			S(3) = D(3,3)*E(3);
#endif
			// integration weight
			S *= dV;
			D *= dV;

			// Stiffness Matrix
			noalias( BTD ) = prod( trans( B ), D );
			noalias( rLeftHandSideMatrix ) += prod( BTD, B );

			// Residual vector
			noalias( rRightHandSideVector ) -= prod( trans( B ), S );

			// Continue the calculation of the EAS method now that the contitutive response
			// has been computed
			EASOp.GaussPointComputation_Step2(D, B, S, mEASStorage);
        }

        // Now that the gauss integration is over, let the EAS operator modify
        // the stiffness matrix and residual vector.
        // It will perform a static condensation to remove the enhanced strain parameters
        // at the element level
        bool eas_error = EASOp.ComputeModfiedTangentAndResidual(rLeftHandSideMatrix, rRightHandSideVector, mEASStorage);

        // Add body forces contributions. This doesn't depend on the coordinate system
		AddBodyForces(rRightHandSideVector);

		// set error code
		if(eas_error) mErrorCode = -1.0;
    }

	void EASQuadElementV2::AddBodyForces(VectorType& rRightHandSideVector)
	{
		const GeometryType& geom = GetGeometry();

		// Get shape functions
		const Matrix & N = geom.ShapeFunctionsValues(mThisIntegrationMethod);

		double rho = GetProperties()[DENSITY];
		double th = GetProperties()[THICKNESS];

		// auxiliary
		array_1d<double, 3> bf;

		// gauss loop to integrate the external force vector
		JacobianOperator jacOp;
		const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(mThisIntegrationMethod);
		const GeometryType::ShapeFunctionsGradientsType& local_gradients = geom.ShapeFunctionsLocalGradients(mThisIntegrationMethod);
		for(unsigned int igauss = 0; igauss < integration_points.size(); igauss++)
		{
			const GeometryType::IntegrationPointType & ip = integration_points[igauss];
			jacOp.Calculate(geom, local_gradients[igauss]);
			double dV = ip.Weight() * jacOp.Determinant() * th;

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
#ifdef EASQ4_DRILLING_DOF
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
    // Class EASQuadElementV2 - Serialization
    //
    // =====================================================================================

	void EASQuadElementV2::save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
        rSerializer.save("CLaw", mConstitutiveLawVector);
        rSerializer.save("IntM", (int)mThisIntegrationMethod);
		rSerializer.save("EAS", mEASStorage);
    }

    void EASQuadElementV2::load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
        rSerializer.save("CLaw", mConstitutiveLawVector);
		int temp;
        rSerializer.load("IntM", temp);
		rSerializer.load("EAS", mEASStorage);
    }

}
