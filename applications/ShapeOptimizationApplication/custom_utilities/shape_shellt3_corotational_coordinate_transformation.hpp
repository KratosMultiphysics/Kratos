// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Massimo Petracca
//

#if !defined(KRATOS_SHAPE_SHELLT3_COROTATIONAL_COORDINATE_TRANSFORMATION_H_INCLUDED )
#define  KRATOS_SHAPE_SHELLT3_COROTATIONAL_COORDINATE_TRANSFORMATION_H_INCLUDED

#include "shape_EICR.hpp"
#include "shape_shellt3_coordinate_transformation.hpp"

#define COROT_POLAR_DECOMPOSITION

namespace Kratos
{
/** \brief SHAPE_EICR ShapeShellT3_CorotationalCoordinateTransformation
*
* This class represents a corotational (nonlinear) coordinate transformation
* that can be used by any element whose geometry is a TRIANGLE 3 in 3D space,
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
* - C.A.Felippa, AFEM.Ch.37, "Triangular ShapeShell Elements",
*   link: http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/
*/
class ShapeShellT3_CorotationalCoordinateTransformation : public ShapeShellT3_CoordinateTransformation
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ShapeShellT3_CorotationalCoordinateTransformation);

    typedef double RealType;

    typedef BoundedMatrix<RealType, 3, 3> TransformationMatrixType;

    typedef array_1d<RealType, 3> Vector3Type;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef Quaternion<RealType> QuaternionType;

public:

    ShapeShellT3_CorotationalCoordinateTransformation(const GeometryType::Pointer& pGeometry)
        : ShapeShellT3_CoordinateTransformation(pGeometry)
        , mInitialized(false)
    {
    }

    ~ShapeShellT3_CorotationalCoordinateTransformation() override
    {
    }

public:

    ShapeShellT3_CoordinateTransformation::Pointer Create(GeometryType::Pointer pGeometry)const override
    {
        return ShapeShellT3_CorotationalCoordinateTransformation::Pointer(new ShapeShellT3_CorotationalCoordinateTransformation(pGeometry));
    }

    void Initialize() override
    {
        KRATOS_TRY

        if (!mInitialized) {
            const GeometryType& geom = GetGeometry();

            ShapeShellT3_LocalCoordinateSystem LCS(CreateReferenceCoordinateSystem());

            mQ0 = QuaternionType::FromRotationMatrix(LCS.Orientation());
            mC0 = LCS.Center();

            for (int i = 0; i < 3; i++) {
                mRV[i] = geom[i].FastGetSolutionStepValue(HELMHOLTZ_ROTS);
                mQN[i] = QuaternionType::FromRotationVector(mRV[i]);

                mRV_converged[i] = mRV[i];
                mQN_converged[i] = mQN[i];
            }

            mInitialized = true;
        }

        KRATOS_CATCH("")
    }

    void InitializeSolutionStep() override
    {
        for (int i = 0; i < 3; i++) {
            mRV[i] = mRV_converged[i];
            mQN[i] = mQN_converged[i];
        }
    }

    void FinalizeSolutionStep() override
    {
        for (int i = 0; i < 3; i++) {
            mRV_converged[i] = mRV[i];
            mQN_converged[i] = mQN[i];
        }
    }

    void FinalizeNonLinearIteration() override
    {
        const GeometryType& geom = GetGeometry();
        Vector3Type incrementalRotation;

        for (int i = 0; i < 3; i++) {
            const Vector3Type& currentRotVec = geom[i].FastGetSolutionStepValue(HELMHOLTZ_ROTS);
            noalias(incrementalRotation) = currentRotVec - mRV[i];
            noalias(mRV[i]) = currentRotVec;

            QuaternionType incrementalQuaternion = QuaternionType::FromRotationVector(incrementalRotation);

            mQN[i] = incrementalQuaternion * mQN[i];
        }
    }

    ShapeShellT3_LocalCoordinateSystem CreateLocalCoordinateSystem()const override
    {
        const GeometryType& geom = GetGeometry();

#ifdef COROT_POLAR_DECOMPOSITION

        // reference coordinate system
        ShapeShellT3_LocalCoordinateSystem a(CreateReferenceCoordinateSystem());

        // current coordinate system using the 1-2 side alignement
        ShapeShellT3_LocalCoordinateSystem b(geom[0], geom[1], geom[2]);

        double aX1 = a.X1();
        double aY1 = a.Y1();
        double bX1 = b.X1();
        double bY1 = b.Y1();
        double aX2 = a.X2();
        double aY2 = a.Y2();
        double bX2 = b.X2();
        double bY2 = b.Y2();
        double aX3 = a.X3();
        double aY3 = a.Y3();
        double bX3 = b.X3();
        double bY3 = b.Y3();

        // now we are in the local coordinate systems (reference and current), i.e. we are looking in the local Z direction
        // which is the same for both coordinate systems.
        // now we can compute the 2D deformation gradient between the 2 configurations, at the element center.

        double C1 = 1/(aX1*aY2 - aX2*aY1 - aX1*aY3 + aX3*aY1 + aX2*aY3 - aX3*aY2);
        double f11 = C1*(aY1 - aY3)*(bX1 - bX2) - C1*(aY1 - aY2)*(bX1 - bX3);
        double f12 = C1*(aX1 - aX2)*(bX1 - bX3) - C1*(aX1 - aX3)*(bX1 - bX2);
        double f21 = C1*(aY1 - aY3)*(bY1 - bY2) - C1*(aY1 - aY2)*(bY1 - bY3);
        double f22 = C1*(aX1 - aX2)*(bY1 - bY3) - C1*(aX1 - aX3)*(bY1 - bY2);

        // now we can extrapolate the rotation angle that makes this deformation gradient symmetric.
        // F = R*U -> find R such that R'*F = U
        double alpha = std::atan2(f21 - f12, f11 + f22);

        // this final coordinate system is the one in which
        // the deformation gradient is equal to the stretch tensor
        return ShapeShellT3_LocalCoordinateSystem(geom[0], geom[1], geom[2], alpha);

#else

        return ShapeShellT3_LocalCoordinateSystem(geom[0], geom[1], geom[2]);

#endif // COROT_POLAR_DECOMPOSITION

    }

    VectorType CalculateLocalDisplacements(const ShapeShellT3_LocalCoordinateSystem& LCS,
                                           const VectorType& globalDisplacements) override
    {
        const GeometryType& geom = GetGeometry();

        Vector3Type deformationalDisplacements;

        QuaternionType Q = QuaternionType::FromRotationMatrix(LCS.Orientation());

        const Vector3Type& C = LCS.Center();

        VectorType localDisplacements(18);

        TransformationMatrixType T, T0;
        Q.ToRotationMatrix(T);
        mQ0.ToRotationMatrix(T0);

        for (unsigned int i = 0; i < 3; i++) {
            unsigned int index = i * 6;

            // get deformational displacements

            noalias(deformationalDisplacements)  = prod(T , geom[i] - C);
            noalias(deformationalDisplacements) -= prod(T0, geom[i].GetInitialPosition() - mC0);

            localDisplacements[index]     = deformationalDisplacements[0];
            localDisplacements[index + 1] = deformationalDisplacements[1];
            localDisplacements[index + 2] = deformationalDisplacements[2];

            // get deformational rotations

            QuaternionType Qd = Q * mQN[i] * mQ0.conjugate();

            Qd.ToRotationVector(localDisplacements[index + 3],
                                localDisplacements[index + 4],
                                localDisplacements[index + 5]);
        }

        return localDisplacements;
    }

    void FinalizeCalculations(const ShapeShellT3_LocalCoordinateSystem& LCS,
                              const VectorType& globalDisplacements,
                              const VectorType& localDisplacements,
                              MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const bool RHSrequired,
                              const bool LHSrequired) override
    {
        // Get the total rotation matrix (local - to - global)
        // Note: do NOT include the warpage correction matrix!
        // Explanation:
        // The Warpage correction matrix computed by the LocalCoordinateSystem is a Linear Projector.
        // It should be used in a LinearCoordinateTransformation.
        // Here instead we already calculate a nonlinear Projector (P = Pu - S * G)!

        MatrixType T(18, 18);
        LCS.ComputeTotalRotationMatrix(T);

        // Form all matrices:
        // S: Spin-Fitter matrix
        // G: Spin-Lever matrix
        // P: Projector (Translational & Rotational)

        MatrixType P(SHAPE_EICR::Compute_Pt(3));
        MatrixType S(SHAPE_EICR::Compute_S(LCS.Nodes()));
        MatrixType G(RotationGradient(LCS));
        noalias(P) -= prod(S, G);

        // Compute the projected local forces ( pe = P' * RHS ).
        // Note: here the RHS is already given as a residual vector -> - internalForces -> (pe = - Ke * U)
        // so projectedLocalForces = - P' * Ke * U

        VectorType projectedLocalForces(prod(trans(P), rRightHandSideVector));

        // Compute the Right-Hand-Side vector in global coordinate system (- T' * P' * Km * U).
        // At this point the computation of the Right-Hand-Side is complete.

        noalias(rRightHandSideVector) = prod(trans(T), projectedLocalForces);

        // Begin the computation of the Left-Hand-Side Matrix :

        if (!LHSrequired) {
            return;    // avoid useless calculations!
        }

        // This is a temporary matrix to store intermediate values
        // to avoid extra dynamic memory allocations!

        MatrixType temp(18, 18);

        // H: Axial Vector Jacobian
        MatrixType H(SHAPE_EICR::Compute_H(localDisplacements));

        // Step 1: ( K.M : Material Stiffness Matrix )
        // Apply the projector to the Material Stiffness Matrix (Ke = P' * Km * H * P)
        // At this point 'LHS' contains the 'projected' Material Stiffness matrix
        // in local corotational coordinate system

        noalias(temp) = prod(rLeftHandSideMatrix, H);
        noalias(rLeftHandSideMatrix) = prod(temp, P);
        noalias(temp) = prod(trans(P), rLeftHandSideMatrix);
        rLeftHandSideMatrix.swap(temp);

        // Step 2: ( K.GP: Equilibrium Projection Geometric Stiffness Matrix )
        // First assemble the 'Fnm' matrix with the Spins of the nodal forces.
        // Actually at this point the 'Fnm' Matrix is the 'Fn' Matrix,
        // because it only contains the spins of the 'translational' forces.
        // At this point 'LHS' contains also this term of the Geometric stiffness
        // (Ke = (P' * Km * H * P) - (G' * Fn' * P))

        MatrixType Fnm(18, 3, 0.0);

        SHAPE_EICR::Spin_AtRow(projectedLocalForces, Fnm, 0);
        SHAPE_EICR::Spin_AtRow(projectedLocalForces, Fnm, 6);
        SHAPE_EICR::Spin_AtRow(projectedLocalForces, Fnm, 12);

        noalias(temp) = prod(trans(G), trans(Fnm));
        noalias(rLeftHandSideMatrix) += prod(temp, P);     // note: '+' not '-' because the RHS already has the negative sign

        // Step 3: ( K.GR: Rotational Geometric Stiffness Matrix )
        // Add the Spins of the nodal moments to 'Fnm'.
        // At this point 'LHS' contains also this term of the Geometric stiffness
        // (Ke = (P' * Km * H * P) - (G' * Fn' * P) - (Fnm * G))

        SHAPE_EICR::Spin_AtRow(projectedLocalForces, Fnm, 3);
        SHAPE_EICR::Spin_AtRow(projectedLocalForces, Fnm, 9);
        SHAPE_EICR::Spin_AtRow(projectedLocalForces, Fnm, 15);

        noalias(rLeftHandSideMatrix) += prod(Fnm, G);     // note: '+' not '-' because the RHS already has the negative sign

        // Step 4: (Global Stiffness Matrix)
        // Transform the LHS to the Global coordinate system.
        // T' * [(P' * Km * H * P) - (G' * Fn' * P) - (Fnm * G)] * T
        noalias(temp) = prod(rLeftHandSideMatrix, T);
        noalias(rLeftHandSideMatrix) = prod(trans(T), temp);
    }

    MatrixType GetNodalDeformationalRotationTensor(const ShapeShellT3_LocalCoordinateSystem& LCS,
            const Vector& globalDisplacements,
            size_t nodeid) override
    {
        if (nodeid>2) {
            return IdentityMatrix(3);
        }

        QuaternionType Q = QuaternionType::FromRotationMatrix(LCS.Orientation());

        QuaternionType Qd = Q * mQN[nodeid] * mQ0.conjugate();

        MatrixType R(3,3);
        Qd.ToRotationMatrix(R);
        return R;
    }

    MatrixType GetNodalDeformationalRotationTensor(const ShapeShellT3_LocalCoordinateSystem& LCS,
            const Vector& globalDisplacements,
            const Vector& N) override
    {
        QuaternionType Q = QuaternionType::FromRotationMatrix(LCS.Orientation());

        RealType qx(0.0);
        RealType qy(0.0);
        RealType qz(0.0);
        RealType qw(0.0);

        for (int i = 0; i < 3; i++) {
            QuaternionType iQd = Q * mQN[i] * mQ0.conjugate();
            iQd.normalize();
            RealType iN = N[i];
            qx += iN * iQd.X();
            qy += iN * iQd.Y();
            qz += iN * iQd.Z();
            qw += iN * iQd.W();
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
    * This is the only matrix not included in the SHAPE_EICR, because it depends on how
    * the corotational frame follows the element.
    * This implementation works for the fitting based on polar decomposition.
    * For the moment the gradient is calculated numerically (perturbation).
    * TODO: Find a closed form of the derivative.
    * @return the Spin Fitter Matrix
    */
#ifdef COROT_POLAR_DECOMPOSITION
    inline MatrixType RotationGradient(const ShapeShellT3_LocalCoordinateSystem& LCS)
    {
        Matrix G(3, 18, 0.0);
        //const GeometryType & geom = GetGeometry();

        ShapeShellT3_LocalCoordinateSystem a(CreateReferenceCoordinateSystem());   // REF coordinate system

        ShapeShellT3_LocalCoordinateSystem::Vector3ContainerType nodes = a.Nodes();

        double aX1 = a.X1();
        double aY1 = a.Y1();
        double aX2 = a.X2();
        double aY2 = a.Y2();
        double aX3 = a.X3();
        double aY3 = a.Y3();

        double pert = std::sqrt(2.0*a.Area()) * 1.0E-3;

        Vector3Type e3;
        Vector3Type e1;

        for (int i = 0; i < 3; i++) { // for each node
            Vector3Type& iNode = nodes[i];

            int index = i * 6;

            for (int j = 0; j < 3; j++) { // for each component in [x, y, z]
                double saved_coord = iNode[j]; // save the current coordinate

                iNode[j] += pert; // apply perturbation

                ShapeShellT3_LocalCoordinateSystem b(nodes[0], nodes[1], nodes[2]);   // perturbed coordinate system (1-2 side alignement)

                double bX1 = b.X1();
                double bY1 = b.Y1();
                double bX2 = b.X2();
                double bY2 = b.Y2();
                double bX3 = b.X3();
                double bY3 = b.Y3();

                // now we are in the local coordinate systems (reference and current), i.e. we are looking in the local Z direction
                // which is the same for both coordinate systems.
                // now we can compute the 2D deformation gradient between the 2 configurations, at the element center.

                double C1 = 1/(aX1*aY2 - aX2*aY1 - aX1*aY3 + aX3*aY1 + aX2*aY3 - aX3*aY2);
                double f11 = C1*(aY1 - aY3)*(bX1 - bX2) - C1*(aY1 - aY2)*(bX1 - bX3);
                double f12 = C1*(aX1 - aX2)*(bX1 - bX3) - C1*(aX1 - aX3)*(bX1 - bX2);
                double f21 = C1*(aY1 - aY3)*(bY1 - bY2) - C1*(aY1 - aY2)*(bY1 - bY3);
                double f22 = C1*(aX1 - aX2)*(bY1 - bY3) - C1*(aX1 - aX3)*(bY1 - bY2);

                double alpha = std::atan2(f21 - f12, f11 + f22);

                ShapeShellT3_LocalCoordinateSystem c(nodes[0], nodes[1], nodes[2], alpha);  // perturbed coordinate system (polar decomposition)

                // save the (numerical) rotation gradient

                noalias(e3) = c.Vz();
                noalias(e1) = c.Vx();

                G(0, index + j) = -e3(1) / pert; // - d(e3.y)/dxj , where xj is x,y,z for j=0,1,2
                G(1, index + j) =  e3(0) / pert; // + d(e3.x)/dxj , where xj is x,y,z for j=0,1,2
                G(2, index + j) =  e1(1) / pert; // + d(e1.y)/dxj , where xj is x,y,z for j=0,1,2

                iNode[j] = saved_coord; // restore the current coordinate
            }
        }

        return G;
    }
#else
    inline MatrixType RotationGradient(const ShapeShellT3_LocalCoordinateSystem& LCS)
    {
        //Matrix G(3, 18, 0.0);
        //const GeometryType & geom = GetGeometry();

        //double m = 1.0 / 2.0 * LCS.Area();

        //double x32 = LCS.X3() - LCS.X2();
        //double y32 = LCS.Y3() - LCS.Y2();

        //double x13 = LCS.X1() - LCS.X3();
        //double y13 = LCS.Y1() - LCS.Y3();

        //double x21 = LCS.X2() - LCS.X1();
        //double y21 = LCS.Y2() - LCS.Y1();

        //double l12 = norm_2( LCS.P2() - LCS.P1() );

        //// G1

        //G(0,  2) = m * x32;
        //G(1,  2) = m * y32;
        //G(2,  1) = -1.0 / l12;

        //// G2

        //G(0,  8) = m * x13;
        //G(1,  8) = m * y13;
        //G(2,  7) = 1.0 / l12;

        //// G3

        //G(0, 14) = m * x21;
        //G(1, 14) = m * y21;

        //return G;

        Matrix G(3, 18, 0.0);
        const GeometryType& geom = GetGeometry();

        ShapeShellT3_LocalCoordinateSystem a(CreateLocalCoordinateSystem());   // REF coordinate system

        ShapeShellT3_LocalCoordinateSystem::Vector3ContainerType nodes = a.Nodes();

        double aX1 = a.X1();
        double aY1 = a.Y1();
        double aX2 = a.X2();
        double aY2 = a.Y2();
        double aX3 = a.X3();
        double aY3 = a.Y3();

        double pert = std::sqrt(2.0*a.Area()) * 1.0E-4;

        Vector3Type e3;
        Vector3Type e1;

        for (int i = 0; i < 3; i++) { // for each node
            Vector3Type& iNode = nodes[i];

            int index = i * 6;

            for (int j = 0; j < 3; j++) { // for each component in [x, y, z]
                double saved_coord = iNode[j]; // save the current coordinate

                iNode[j] += pert; // apply perturbation

                ShapeShellT3_LocalCoordinateSystem b(nodes[0], nodes[1], nodes[2]);   // perturbed coordinate system (1-2 side alignement)

                // save the (numerical) rotation gradient

                noalias(e3) = b.Vz();
                noalias(e1) = b.Vx();

                G(0, index + j) = -e3(1) / pert; // - d(e3.y)/dxj , where xj is x,y,z for j=0,1,2
                G(1, index + j) =  e3(0) / pert; // + d(e3.x)/dxj , where xj is x,y,z for j=0,1,2
                G(2, index + j) =  e1(1) / pert; // + d(e1.y)/dxj , where xj is x,y,z for j=0,1,2

                iNode[j] = saved_coord; // restore the current coordinate
            }
        }

        return G;
    }
#endif // COROT_NUM_ROT_GRAD


protected:

    ShapeShellT3_CorotationalCoordinateTransformation()
        : ShapeShellT3_CoordinateTransformation()
    {
    }

private:

    bool mInitialized;
    QuaternionType mQ0;
    Vector3Type mC0;
    array_1d< QuaternionType, 3 > mQN;
    array_1d< Vector3Type, 3 > mRV;
    array_1d< QuaternionType, 3 > mQN_converged;
    array_1d< Vector3Type, 3 > mRV_converged;

private:

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  ShapeShellT3_CoordinateTransformation);
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
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,  ShapeShellT3_CoordinateTransformation);
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


#endif // KRATOS_SHAPE_SHELLT3_COROTATIONAL_COORDINATE_TRANSFORMATION_H_INCLUDED
