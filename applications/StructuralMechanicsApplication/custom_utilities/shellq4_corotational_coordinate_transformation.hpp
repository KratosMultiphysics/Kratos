//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-12-29 23:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(KRATOS_SHELLQ4_COROTATIONAL_COORDINATE_TRANSFORMATION_H_INCLUDED )
#define  KRATOS_SHELLQ4_COROTATIONAL_COORDINATE_TRANSFORMATION_H_INCLUDED

#include "EICR.hpp"
#include "shellq4_coordinate_transformation.hpp"

namespace Kratos
{

	/** \brief EICR ShellQ4_CorotationalCoordinateTransformation
	*
	* This class represents a corotational (nonlinear) coordinate transformation
	* that can be used by any element whose geometry is a QUAD 4 in 3D space,
	* with 6 D.O.F.s per node.
	* It's main aim is to:
	* 1) Create the local coordinate system
	* 2) Transform the incoming global displacements in local coordinate system
	*    removing rigid body displacements and rotations.
	* 3) Transform the outgoing matrices and vectors in global coordinate system
	*    with rigid body displacements and rotations.
	*
	* Updated version:
	* - Makes use of Quaternions (Euler Parameters) to parametrize finite rotations
	*   in an efficient and robust way.
	*
	* References:
	*
	* - C.A.Felippa,B.Haugen, "Unified formulation of small-strain corotational
	*   finite elements: I. Theory",
	*   CU-CAS-05-02, January 2005
	* - C.A.Felippa, AFEM.Ch.38, "Quadrilateral Shell Elements",
	*   Chapter 5 of B.Haugen's Thesis.
	*   link: http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/
	*/
	class ShellQ4_CorotationalCoordinateTransformation : public ShellQ4_CoordinateTransformation
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION( ShellQ4_CorotationalCoordinateTransformation );

		typedef double RealType;

		typedef BoundedMatrix<RealType, 3, 3> TransformationMatrixType;

		typedef array_1d<RealType, 3> Vector3Type;

		typedef Vector VectorType;

		typedef Matrix MatrixType;

		typedef Quaternion<RealType> QuaternionType;

	public:

		ShellQ4_CorotationalCoordinateTransformation(const GeometryType::Pointer & pGeometry)
			: ShellQ4_CoordinateTransformation(pGeometry)
			, mInitialized(false)
		{
		}

		~ShellQ4_CorotationalCoordinateTransformation() override
		{
		}

	public:

		ShellQ4_CoordinateTransformation::Pointer Create(GeometryType::Pointer pGeometry)const override
		{
			return ShellQ4_CorotationalCoordinateTransformation::Pointer( new ShellQ4_CorotationalCoordinateTransformation( pGeometry ) );
		}

		void Initialize() override
		{
			KRATOS_TRY

			if(!mInitialized)
			{
				const GeometryType & geom = GetGeometry();

				ShellQ4_LocalCoordinateSystem LCS( CreateReferenceCoordinateSystem() );

				mQ0 = QuaternionType::FromRotationMatrix( LCS.Orientation() );
				mC0 = LCS.Center();

				for(int i = 0; i < 4; i++)
				{
					mRV[i] = geom[i].FastGetSolutionStepValue( ROTATION );
					mQN[i] = QuaternionType::FromRotationVector( mRV[i] );

					mRV_converged[i] = mRV[i];
					mQN_converged[i] = mQN[i];
				}

				mInitialized = true;
			}

			KRATOS_CATCH("")
		}

		void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override
		{
			for(int i = 0; i < 4; i++)
			{
				mRV[i] = mRV_converged[i];
				mQN[i] = mQN_converged[i];
			}
		}

		void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo) override
		{
			for(int i = 0; i < 4; i++)
			{
				mRV_converged[i] = mRV[i];
				mQN_converged[i] = mQN[i];
			}
		}

		void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo) override
		{
		}

		void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo) override
		{
			const GeometryType & geom = GetGeometry();
			Vector3Type incrementalRotation;

			for(int i = 0; i < 4; i++)
			{
				const Vector3Type & currentRotVec = geom[i].FastGetSolutionStepValue( ROTATION );
				noalias( incrementalRotation ) = currentRotVec - mRV[i];
				noalias( mRV[i] ) = currentRotVec;

				QuaternionType incrementalQuaternion = QuaternionType::FromRotationVector( incrementalRotation );

				mQN[i] = incrementalQuaternion * mQN[i];
			}
		}

		ShellQ4_LocalCoordinateSystem CreateLocalCoordinateSystem()const override
		{
			const GeometryType & geom = GetGeometry();

			//return ShellQ4_LocalCoordinateSystem(geom[0], geom[1], geom[2], geom[3]);

			// reference coordinate system
			ShellQ4_LocalCoordinateSystem a( CreateReferenceCoordinateSystem() );

			// current coordinate system using the 1-2 side alignement
			ShellQ4_LocalCoordinateSystem b( geom[0], geom[1], geom[2], geom[3] );

			double aX1 = a.X1(); double aY1 = a.Y1();
			double bX1 = b.X1(); double bY1 = b.Y1();
			double aX2 = a.X2(); double aY2 = a.Y2();
			double bX2 = b.X2(); double bY2 = b.Y2();
			double aX3 = a.X3(); double aY3 = a.Y3();
			double bX3 = b.X3(); double bY3 = b.Y3();
			double aX4 = a.X4(); double aY4 = a.Y4();
			double bX4 = b.X4(); double bY4 = b.Y4();

			// now we are in the local coordinate systems (reference and current), i.e. we are looking in the local Z direction
			// which is the same for both coordinate systems.
			// now we can compute the 2D deformation gradient between the 2 configurations, at the element center.

			double C1 = 1/(aX1*aY2 - aX2*aY1 - aX1*aY4 + aX2*aY3 - aX3*aY2 + aX4*aY1 + aX3*aY4 - aX4*aY3);
			double C2 = bY1/4 + bY2/4 - bY3/4 - bY4/4;
			double C3 = bY1/4 - bY2/4 - bY3/4 + bY4/4;
			double C4 = bX1/4 + bX2/4 - bX3/4 - bX4/4;
			double C5 = bX1/4 - bX2/4 - bX3/4 + bX4/4;
			double C6 = aX1 + aX2 - aX3 - aX4;
			double C7 = aX1 - aX2 - aX3 + aX4;
			double C8 = aY1 + aY2 - aY3 - aY4;
			double C9 = aY1 - aY2 - aY3 + aY4;
			double f11 = 2*C1*C5*C8 - 2*C1*C4*C9;
			double f12 = 2*C1*C4*C7 - 2*C1*C5*C6;
			double f21 = 2*C1*C3*C8 - 2*C1*C2*C9;
			double f22 = 2*C1*C2*C7 - 2*C1*C3*C6;

			// now we can extrapolate the rotation angle that makes this deformation gradient symmetric.
			// F = R*U -> find R such that R'*F = U
			double alpha = std::atan2( f21 - f12, f11 + f22 );

			// this final coordinate system is the one in which
			// the deformation gradient is equal to the stretch tensor
			return ShellQ4_LocalCoordinateSystem( geom[0], geom[1], geom[2], geom[3], alpha );
		}

		VectorType CalculateLocalDisplacements(const ShellQ4_LocalCoordinateSystem & LCS,
													   const VectorType & globalDisplacements) override
		{
			const GeometryType & geom = GetGeometry();

			Vector3Type deformationalDisplacements;

			QuaternionType Q = QuaternionType::FromRotationMatrix( LCS.Orientation() );

			const Vector3Type & C = LCS.Center();

			VectorType localDisplacements(24);

			TransformationMatrixType T, T0;
			Q.ToRotationMatrix( T );
			mQ0.ToRotationMatrix( T0 );

			for(unsigned int i = 0; i < 4; i++)
			{
				unsigned int index = i * 6;

				// get deformational displacements

				noalias( deformationalDisplacements )  = prod( T , geom[i] - C );
				noalias( deformationalDisplacements ) -= prod( T0, geom[i].GetInitialPosition() - mC0 );

				localDisplacements[index]     = deformationalDisplacements[0];
				localDisplacements[index + 1] = deformationalDisplacements[1];
				localDisplacements[index + 2] = deformationalDisplacements[2];

				// get deformational rotations

				QuaternionType Qd = Q * mQN[i] * mQ0.conjugate();

				Qd.ToRotationVector( localDisplacements[index + 3],
									 localDisplacements[index + 4],
									 localDisplacements[index + 5] );
			}

			return localDisplacements;
		}

		void FinalizeCalculations(const ShellQ4_LocalCoordinateSystem & LCS,
										  const VectorType & globalDisplacements,
										  const VectorType & localDisplacements,
										  MatrixType & rLeftHandSideMatrix,
										  VectorType & rRightHandSideVector,
										  const bool RHSrequired,
										  const bool LHSrequired) override
		{
			// Get the total rotation matrix (local - to - global)
			// Note: do NOT include the warpage correction matrix!
			// Explanation:
			// The Warpage correction matrix computed by the LocalCoordinateSystem is a Linear Projector.
			// It should be used in a LinearCoordinateTransformation.
			// Here instead we already calculate a nonlinear Projector (P = Pu - S * G)!

			MatrixType T(24, 24);
			LCS.ComputeTotalRotationMatrix( T );

			// Form all matrices:
			// S: Spin-Fitter matrix
			// G: Spin-Lever matrix
			// P: Projector (Translational & Rotational)

			MatrixType P( EICR::Compute_Pt(4) );
			MatrixType S( EICR::Compute_S(LCS.Nodes()) );
			//MatrixType G( Compute_G(LCS.P1(), LCS.P2(), LCS.P3(), LCS.P4(), LCS.Area()) );
			MatrixType G( RotationGradient() );
			noalias( P ) -= prod( S, G );

			// Compute the projected local forces ( pe = P' * RHS ).
			// Note: here the RHS is already given as a residual vector -> - internalForces -> (pe = - Ke * U)
			// so projectedLocalForces = - P' * Ke * U

			VectorType projectedLocalForces( prod( trans( P ), rRightHandSideVector ) );

			// Compute the Right-Hand-Side vector in global coordinate system (- T' * P' * Km * U).
			// At this point the computation of the Right-Hand-Side is complete.

			noalias( rRightHandSideVector ) = prod( trans( T ), projectedLocalForces );

			// Begin the computation of the Left-Hand-Side Matrix :

			if(!LHSrequired) return; // avoid useless calculations!

			// This is a temporary matrix to store intermediate values
			// to avoid extra dynamic memory allocations!

			MatrixType temp(24, 24);

			// H: Axial Vector Jacobian
			MatrixType H( EICR::Compute_H(localDisplacements) );

			// Step 1: ( K.M : Material Stiffness Matrix )
			// Apply the projector to the Material Stiffness Matrix (Ke = P' * Km * H * P)
			// At this point 'LHS' contains the 'projected' Material Stiffness matrix
			// in local corotational coordinate system

			noalias( temp ) = prod( rLeftHandSideMatrix, H );
			noalias( rLeftHandSideMatrix ) = prod( temp, P );
			noalias( temp ) = prod( trans( P ), rLeftHandSideMatrix );
			rLeftHandSideMatrix.swap( temp );

			// Step 2: ( K.GP: Equilibrium Projection Geometric Stiffness Matrix )
			// First assemble the 'Fnm' matrix with the Spins of the nodal forces.
			// Actually at this point the 'Fnm' Matrix is the 'Fn' Matrix,
			// because it only contains the spins of the 'translational' forces.
			// At this point 'LHS' contains also this term of the Geometric stiffness
			// (Ke = (P' * Km * H * P) - (G' * Fn' * P))

			MatrixType Fnm(24, 3, 0.0);

			EICR::Spin_AtRow( projectedLocalForces, Fnm, 0  );
			EICR::Spin_AtRow( projectedLocalForces, Fnm, 6  );
			EICR::Spin_AtRow( projectedLocalForces, Fnm, 12 );
			EICR::Spin_AtRow( projectedLocalForces, Fnm, 18 );

			noalias( temp ) = prod( trans( G ), trans( Fnm ) );
			noalias( rLeftHandSideMatrix ) += prod( temp, P ); // note: '+' not '-' because the RHS already has the negative sign

			// Step 3: ( K.GR: Rotational Geometric Stiffness Matrix )
			// Add the Spins of the nodal moments to 'Fnm'.
			// At this point 'LHS' contains also this term of the Geometric stiffness
			// (Ke = (P' * Km * H * P) - (G' * Fn' * P) - (Fnm * G))

			EICR::Spin_AtRow( projectedLocalForces, Fnm, 3  );
			EICR::Spin_AtRow( projectedLocalForces, Fnm, 9  );
			EICR::Spin_AtRow( projectedLocalForces, Fnm, 15 );
			EICR::Spin_AtRow( projectedLocalForces, Fnm, 21 );

			noalias( rLeftHandSideMatrix ) += prod( Fnm, G ); // note: '+' not '-' because the RHS already has the negative sign

			// Step 4: (Global Stiffness Matrix)
			// Transform the LHS to the Global coordinate system.
			// T' * [(P' * Km * H * P) - (G' * Fn' * P) - (Fnm * G)] * T
			noalias( temp ) = prod( rLeftHandSideMatrix, T );
			noalias( rLeftHandSideMatrix ) = prod( trans( T ), temp );
		}

		MatrixType GetNodalDeformationalRotationTensor(const ShellQ4_LocalCoordinateSystem & LCS,
			                                                   const Vector& globalDisplacements,
															   size_t nodeid) override
		{
			if(nodeid>3) return IdentityMatrix(3,3);

			QuaternionType Q = QuaternionType::FromRotationMatrix( LCS.Orientation() );

			QuaternionType Qd = Q * mQN[nodeid] * mQ0.conjugate();

			MatrixType R(3,3);
			Qd.ToRotationMatrix(R);
			return R;
		}

		MatrixType GetNodalDeformationalRotationTensor(const ShellQ4_LocalCoordinateSystem & LCS,
			                                                   const Vector& globalDisplacements,
															   const Vector& N) override
		{
			QuaternionType Q = QuaternionType::FromRotationMatrix( LCS.Orientation() );

			RealType qx(0.0);
			RealType qy(0.0);
			RealType qz(0.0);
			RealType qw(0.0);

			for(int i = 0; i < 4; i++)
			{
				QuaternionType iQd = Q * mQN[i] * mQ0.conjugate();
				iQd.normalize();
				RealType iN = N[i];
				qx += iN * iQd.x();
				qy += iN * iQd.y();
				qz += iN * iQd.z();
				qw += iN * iQd.w();
			}

			MatrixType R(3, 3);

			QuaternionType Qd(qw, qx, qy, qz);
			Qd.normalize();
			Qd.ToRotationMatrix(R);

			return R;
		}

	private:

		/**
		* Computes the Spin Fitter Matrix.
		* This is the only matrix not included in the EICR, because it depends on how
		* the corotational frame follows the element.
		* In this implementation the local x axis of the corotational frame is kept
		* parallel to the 1-2 side.
		* @param P1 the 1st node
		* @param P2 the 2nd node
		* @param P3 the 3rd node
		* @param P4 the 4th node
		* @param area the element area
		* @return the Spin Fitter Matrix
		*/
		inline MatrixType Compute_G(const Vector3Type & P1,
									const Vector3Type & P2,
									const Vector3Type & P3,
									const Vector3Type & P4,
									RealType area)
		{
			RealType Ap = 2.0 * area;
			RealType m = 1.0 / Ap;

			Vector3Type D12( P2 - P1 );
			Vector3Type D24( P4 - P2 );
			Vector3Type D13( P3 - P1 );

			RealType x42 = D24(0);
			RealType x24 = - x42;
			RealType y42 = D24(1);
			RealType y24 = - y42;
			RealType x31 = D13(0);
			RealType x13 = - x31;
			RealType y31 = D13(1);
			RealType y13 = - y31;

			// Note, assuming the input vectors are in local CR,
			// l12 is the length of the side 1-2 projected onto the xy plane.
			RealType l12 = std::sqrt( D12(0)*D12(0) + D12(1)*D12(1) );

			MatrixType G(3, 24, 0.0);

			// G1

			G(0,  2) = x42 * m;
			G(1,  2) = y42 * m;
			G(2,  1) = - 1.0 / l12;

			// G2

			G(0,  8) = x13 * m;
			G(1,  8) = y13 * m;
			G(2,  7) = 1.0 / l12;

			// G3

			G(0, 14) = x24 * m;
			G(1, 14) = y24 * m;

			// G4

			G(0, 20) = x31 * m;
			G(1, 20) = y31 * m;

			return G;
		}

		/**
		* Computes the Spin Fitter Matrix.
		* This is the only matrix not included in the EICR, because it depends on how
		* the corotational frame follows the element.
		* This implementation works for the fitting based on polar decomposition.
		* For the moment the gradient is calculated numerically (perturbation).
		* TODO: Find a closed form of the derivative.
		* @return the Spin Fitter Matrix
		*/
		inline MatrixType RotationGradient()
		{
			Matrix G(3, 24, 0.0);

			ShellQ4_LocalCoordinateSystem a( CreateReferenceCoordinateSystem() ); // current coordinate system

			ShellQ4_LocalCoordinateSystem::Vector3ContainerType nodes = a.Nodes();

			double aX1 = a.X1(); double aY1 = a.Y1();
			double aX2 = a.X2(); double aY2 = a.Y2();
			double aX3 = a.X3(); double aY3 = a.Y3();
			double aX4 = a.X4(); double aY4 = a.Y4();

			double pert = std::sqrt(a.Area()) * 1.0E-2;

			Vector3Type e3;
			Vector3Type e1;

			for(int i = 0; i < 4; i++) // for each node
			{
				Vector3Type& iNode = nodes[i];

				int index = i * 6;

				for(int j = 0; j < 3; j++) // for each component in [x, y, z]
				{
					double saved_coord = iNode[j]; // save the current coordinate

					iNode[j] += pert; // apply perturbation

					ShellQ4_LocalCoordinateSystem b( nodes[0], nodes[1], nodes[2], nodes[3] ); // perturbed coordinate system (1-2 side alignement)

					double bX1 = b.X1(); double bY1 = b.Y1();
					double bX2 = b.X2(); double bY2 = b.Y2();
					double bX3 = b.X3(); double bY3 = b.Y3();
					double bX4 = b.X4(); double bY4 = b.Y4();

					double C1 = 1/(aX1*aY2 - aX2*aY1 - aX1*aY4 + aX2*aY3 - aX3*aY2 + aX4*aY1 + aX3*aY4 - aX4*aY3);
					double C2 = bY1/4 + bY2/4 - bY3/4 - bY4/4;
					double C3 = bY1/4 - bY2/4 - bY3/4 + bY4/4;
					double C4 = bX1/4 + bX2/4 - bX3/4 - bX4/4;
					double C5 = bX1/4 - bX2/4 - bX3/4 + bX4/4;
					double C6 = aX1 + aX2 - aX3 - aX4;
					double C7 = aX1 - aX2 - aX3 + aX4;
					double C8 = aY1 + aY2 - aY3 - aY4;
					double C9 = aY1 - aY2 - aY3 + aY4;
					double f11 = 2*C1*C5*C8 - 2*C1*C4*C9;
					double f12 = 2*C1*C4*C7 - 2*C1*C5*C6;
					double f21 = 2*C1*C3*C8 - 2*C1*C2*C9;
					double f22 = 2*C1*C2*C7 - 2*C1*C3*C6;

					double alpha = std::atan2( f21 - f12, f11 + f22 );

					ShellQ4_LocalCoordinateSystem c( nodes[0], nodes[1], nodes[2], nodes[3], alpha ); // perturbed coordinate system (polar decomposition)

					// save the (numerical) rotation gradient

					noalias( e3 ) = c.Vz();
					noalias( e1 ) = c.Vx();

					G(0, index + j) = -e3(1) / pert; // - d(e3.y)/dxj , where xj is x,y,z for j=0,1,2
					G(1, index + j) =  e3(0) / pert; // + d(e3.x)/dxj , where xj is x,y,z for j=0,1,2
					G(2, index + j) =  e1(1) / pert; // + d(e1.y)/dxj , where xj is x,y,z for j=0,1,2

					iNode[j] = saved_coord; // restore the current coordinate
				}
			}

			return G;
		}

	protected:

		ShellQ4_CorotationalCoordinateTransformation()
			: ShellQ4_CoordinateTransformation()
		{
		}

	private:

		bool mInitialized;
		QuaternionType mQ0;
		Vector3Type mC0;
		array_1d< QuaternionType, 4 > mQN;
		array_1d< Vector3Type, 4 > mRV;
		array_1d< QuaternionType, 4 > mQN_converged;
		array_1d< Vector3Type, 4 > mRV_converged;

	private:

		friend class Serializer;

		void save(Serializer& rSerializer) const override
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  ShellQ4_CoordinateTransformation );
			rSerializer.save("init", mInitialized);
			rSerializer.save("Q0", mQ0);
			rSerializer.save("C0", mC0);
			rSerializer.save("QN", mQN);
			rSerializer.save("RV", mRV);
			rSerializer.save("QN_conv", mQN_converged);
			rSerializer.save("RV_conv", mRV_converged);
		}

		void load(Serializer& rSerializer) override
		{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  ShellQ4_CoordinateTransformation );
			rSerializer.load("init", mInitialized);
			rSerializer.load("Q0", mQ0);
			rSerializer.load("C0", mC0);
			rSerializer.load("QN", mQN);
			rSerializer.load("RV", mRV);
			rSerializer.load("QN_conv", mQN_converged);
			rSerializer.load("RV_conv", mRV_converged);

		}

	};

}


#endif // KRATOS_SHELLQ4_COROTATIONAL_COORDINATE_TRANSFORMATION_H_INCLUDED
