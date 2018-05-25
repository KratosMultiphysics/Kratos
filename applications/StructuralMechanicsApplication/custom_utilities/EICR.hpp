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

#if !defined(EICR_H_INCLUDED )
#define  EICR_H_INCLUDED

namespace Kratos
{
/** \brief EICR Element Independent CoRotational formulation
*
* E.I.C.R. is a utility class containing static methods related to
* the Element Independent Corotational Formulation.
* This class implements methods that do not depend on the element type,
* and so they can be used by any implementation of a corotational coordinate transformation.
*/
class EICR
{

public:

    typedef double RealType;

    typedef BoundedMatrix<RealType, 3, 3> BoundedMatrixType3x3;

    typedef array_1d<RealType, 3> Vector3Type;

    typedef std::vector<Vector3Type> Vector3ContainerType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t SizeType;

public:

    /**
    * Computes the Spin of the input vector rV, and saves the result into the output matrix rS.
    * Note: no check is made on the size of the input-output arguments.
    * @param rV the input vector (assumed size: >= 3)
    * @param rS the output matrix (assumed size: >= 3x3)
    */
    template< class TVec, class TMat>
    inline static void Spin(const TVec& rV, TMat& rS)
    {
        rS(0, 0) =	 0.00;
        rS(0, 1) = - rV(2);
        rS(0, 2) =   rV(1);
        rS(1, 0) =	 rV(2);
        rS(1, 1) =   0.00;
        rS(1, 2) = - rV(0);
        rS(2, 0) = - rV(1);
        rS(2, 1) =   rV(0);
        rS(2, 2) =   0.00;
    }

    /**
    * Computes the Spin of the input vector rV, and saves the result into the output matrix rS,
    * at the specified row index.
    * Note: no check is made on the size of the input-output arguments.
    * @param rV the input vector (assumed size: >= 3)
    * @param rS the output matrix (assumed size: >= 3x3)
    * @param RowIndex the index of the first row in the output matrix where the spin has to be saved
    */
    template< class TVec, class TMat>
    inline static void Spin_AtRow(const TVec& rV, TMat& rS, const SizeType RowIndex)
    {
        const SizeType i0 = RowIndex;
        const SizeType i1 = RowIndex + 1;
        const SizeType i2 = RowIndex + 2;
        const double v0 = rV(i0);
        const double v1 = rV(i1);
        const double v2 = rV(i2);
        rS(i0, 0) =	  0.00;
        rS(i0, 1) = - v2;
        rS(i0, 2) =   v1;
        rS(i1, 0) =	  v2;
        rS(i1, 1) =   0.00;
        rS(i1, 2) = - v0;
        rS(i2, 0) = - v1;
        rS(i2, 1) =   v0;
        rS(i2, 2) =   0.00;
    }

    /**
    * Computes the Spin of the input vector rV, from the specified index, and saves the result into the output matrix rS,
    * at the specified row index.
    * Note: no check is made on the size of the input-output arguments.
    * @param rV the input vector (assumed size: >= 3)
    * @param rS the output matrix (assumed size: >= 3x3)
    * @param VectorIndex the index of the first component of the input vector to be used to compute the spin
    * @param MatrixRowIndex the index of the first row in the output matrix where the spin has to be saved
    */
    template< class TVec, class TMat>
    inline static void Spin_AtRow(const TVec& rV, TMat& rS,
                                  const SizeType VectorIndex,
                                  const SizeType MatrixRowIndex)
    {
        const SizeType i0 = MatrixRowIndex;
        const SizeType i1 = MatrixRowIndex + 1;
        const SizeType i2 = MatrixRowIndex + 2;
        const double v0 = rV(VectorIndex);
        const double v1 = rV(VectorIndex + 1);
        const double v2 = rV(VectorIndex + 2);
        rS(i0, 0) =	  0.00;
        rS(i0, 1) = - v2;
        rS(i0, 2) =   v1;
        rS(i1, 0) =   v2;
        rS(i1, 1) =   0.00;
        rS(i1, 2) = - v0;
        rS(i2, 0) = - v1;
        rS(i2, 1) =   v0;
        rS(i2, 2) =   0.00;
    }

    /**
    * Computes the Spin of the input vector rV, and saves the result into the output matrix rS.
    * This version uses a multiplier for the output values.
    * Note: no check is made on the size of the input-output arguments.
    * @param rV the input vector (assumed size: >= 3)
    * @param rS the output matrix (assumed size: >= 3x3)
    * @param Multiplier the multiplier for the output values
    */
    template< class TVec, class TMat>
    inline static void Spin(const TVec& rV, TMat& rS, double Multiplier)
    {
        rS(0, 0) =	 0.00;
        rS(0, 1) = - Multiplier * rV(2);
        rS(0, 2) =   Multiplier * rV(1);
        rS(1, 0) =	 Multiplier * rV(2);
        rS(1, 1) =   0.00;
        rS(1, 2) = - Multiplier * rV(0);
        rS(2, 0) = - Multiplier * rV(1);
        rS(2, 1) =   Multiplier * rV(0);
        rS(2, 2) =   0.00;
    }

    /**
    * Computes the Spin of the input vector rV, and saves the result into the output matrix rS,
    * at the specified row index.
    * This version uses a multiplier for the output values.
    * Note: no check is made on the size of the input-output arguments.
    * @param rV the input vector (assumed size: >= 3)
    * @param rS the output matrix (assumed size: >= 3x3)
    * @param Multiplier the multiplier for the output values
    * @param RowIndex the index of the first row in the output matrix where the spin has to be saved
    */
    template< class TVec, class TMat>
    inline static void Spin_AtRow(const TVec& rV, TMat& rS,
                                  const double Multiplier,
                                  const SizeType RowIndex)
    {
        const SizeType i0 = RowIndex;
        const SizeType i1 = RowIndex + 1;
        const SizeType i2 = RowIndex + 2;
        const double v0 = Multiplier * rV(i0);
        const double v1 = Multiplier * rV(i1);
        const double v2 = Multiplier * rV(i2);
        rS(i0, 0) =	  0.00;
        rS(i0, 1) = - v2;
        rS(i0, 2) =   v1;
        rS(i1, 0) =	  v2;
        rS(i1, 1) =   0.00;
        rS(i1, 2) = - v0;
        rS(i2, 0) = - v1;
        rS(i2, 1) =   v0;
        rS(i2, 2) =   0.00;
    }

    /**
    * Computes the Spin of the input vector rV, from the specified index, and saves the result into the output matrix rS,
    * at the specified row index.
    * This version uses a multiplier for the output values.
    * Note: no check is made on the size of the input-output arguments.
    * @param rV the input vector (assumed size: >= 3)
    * @param rS the output matrix (assumed size: >= 3x3)
    * @param Multiplier the multiplier for the output values
    * @param VectorIndex the index of the first component of the input vector to be used to compute the spin
    * @param MatrixRowIndex the index of the first row in the output matrix where the spin has to be saved
    */
    template< class TVec, class TMat>
    inline static void Spin_AtRow(const TVec& rV, TMat& rS,
                                  const double Multiplier,
                                  const SizeType VectorIndex,
                                  const SizeType MatrixRowIndex)
    {
        const SizeType i0 = MatrixRowIndex;
        const SizeType i1 = MatrixRowIndex + 1;
        const SizeType i2 = MatrixRowIndex + 2;
        const double v0 = Multiplier * rV(VectorIndex);
        const double v1 = Multiplier * rV(VectorIndex + 1);
        const double v2 = Multiplier * rV(VectorIndex + 2);
        rS(i0, 0) =	  0.00;
        rS(i0, 1) = - v2;
        rS(i0, 2) =   v1;
        rS(i1, 0) =	  v2;
        rS(i1, 1) =   0.00;
        rS(i1, 2) = - v0;
        rS(i2, 0) = - v1;
        rS(i2, 1) =   v0;
        rS(i2, 2) =   0.00;
    }

    /**
    * Computes the Translational Projector Matrix.
    * The output is a square matrix of size NumNodes*6.
    * Note that 6 Degrees Of Freedom are assumed for each node.
    * @param NumNodes the number of nodes
    * @return the Translational Projector Matrix
    */
    inline static MatrixType Compute_Pt(const SizeType NumNodes)
    {
        const RealType a = RealType(NumNodes - 1) / RealType(NumNodes);
        const RealType b = -1.0 / RealType(NumNodes);

        const SizeType num_dofs = NumNodes * 6;

        MatrixType P( IdentityMatrix(num_dofs, num_dofs) );

        for(SizeType i = 0; i < NumNodes; i++)
        {
            const SizeType j = i * 6;

            // diagonal block
            P(j    ,     j) = a;
            P(j + 1, j + 1) = a;
            P(j + 2, j + 2) = a;

            // out-of-diagonal block
            for(SizeType k = i + 1; k < NumNodes; k++)
            {
                const SizeType w = k * 6;

                P(j    , w    ) = b;
                P(j + 1, w + 1) = b;
                P(j + 2, w + 2) = b;

                P(w    , j    ) = b;
                P(w + 1, j + 1) = b;
                P(w + 2, j + 2) = b;
            }
        }

        return P;
    }

    /**
    * Computes the Spin Lever Matrix.
    * The output is a rectangular matrix of 3 columns and rNodes.size()*6 rows.
    * Note that 6 Degrees Of Freedom are assumed for each node.
    * @param rNodes the input nodes
    * @return the Spin Lever Matrix
    */
    inline static MatrixType Compute_S(const Vector3ContainerType& rNodes)
    {
        const SizeType num_nodes = rNodes.size();
        const SizeType num_dofs = num_nodes * 6;

        MatrixType S(num_dofs, 3, 0.0);

        for(SizeType i = 0; i < num_nodes; i++)
        {
            SizeType j = i * 6;

            Spin_AtRow( rNodes[i], S, -1.0, 0, j );

            S(j + 3, 0) = 1.0;
            S(j + 4, 1) = 1.0;
            S(j + 5, 2) = 1.0;
        }

        return S;
    }

    /**
    * Computes the Axial Vector Jacobian.
    * The output is a square matrix of size rDisplacements.size() (which is num_nodes * 6).
    * Note that 6 Degrees Of Freedom are assumed for each node.
    * @param rDisplacements the vector of nodal displacements and rotations in the local corotational coordinate system. (assumed size = num_nodes*6)
    * @return the H matrix
    */
    inline static MatrixType Compute_H(const VectorType& rDisplacements)
    {
        const SizeType num_dofs = rDisplacements.size();
        const SizeType num_nodes = num_dofs / 6;

        MatrixType H( IdentityMatrix(num_dofs, num_dofs) );

        BoundedMatrixType3x3 omega(3, 3);
        BoundedMatrixType3x3 Hi(3, 3);

        for(SizeType i = 0; i < num_nodes; i++)
        {
            const SizeType index = i * 6;
            Vector3Type rv = project( rDisplacements, range(index + 3, index + 6) );

            double angle = norm_2(rv);

            if(angle >= 2.0 * Globals::Pi)
                angle = std::fmod(angle, 2.0 * Globals::Pi);

            double eta;
            if(angle < 0.05)
            {
                double angle2 = angle * angle;
                double angle4 = angle2 * angle2;
                double angle6 = angle4 * angle2;
                eta = 1.0 / 12.0 + 1.0 / 270.0 * angle2 + 1.0 / 30240.0 * angle4 + 1.0 / 1209600.0 * angle6;
            }
            else
            {
                eta = ( 1.0 - 0.5 * angle * std::tan( 0.5 * Globals::Pi - 0.5 * angle ) ) / (angle * angle);
            }

            Spin( rv, omega );

            noalias( Hi ) = IdentityMatrix(3, 3);
            noalias( Hi ) -= 0.5 * omega;
            noalias( Hi ) += eta * prod( omega, omega );

            range i_range(index + 3, index + 6);
            project( H, i_range, i_range ) = Hi;
        }

        return H;
    }

    /**
    * Computes the Spin derivative of (Axial Vector Jacobian)^T contracted with the nodal moment vector.
    * The output is a square matrix of size rDisplacements.size() (which is num_nodes * 6).
    * Note that 6 Degrees Of Freedom are assumed for each node.
    * @param rDisplacements the vector of nodal displacements and rotations in the local corotational coordinate system. (assumed size = num_nodes*6)
    * @param rForces the vector of nodal forces and moments in the local corotational coordinate system. (assumed size = num_nodes*6)
    * @param rH the Axial Vector Jacobian Matrix computed with a previous call to EICR::Compute_H(rDisplacements)
    * @return the L matrix
    */
    inline static MatrixType Compute_L(const VectorType& rDisplacements,
                                       const VectorType& rForces,
                                       const MatrixType& rH)
    {
        const SizeType num_dofs = rDisplacements.size();
        const SizeType num_nodes = num_dofs / 6;

        MatrixType L(num_dofs, num_dofs, 0.0);

        Vector3Type rotation_vector;
        Vector3Type moment_vector;
        BoundedMatrixType3x3 omega(3, 3);
        BoundedMatrixType3x3 omega_2(3, 3);
        BoundedMatrixType3x3 Li(3, 3);
        BoundedMatrixType3x3 LiTemp1(3, 3);
        BoundedMatrixType3x3 LiTemp2(3, 3);

        for(SizeType i = 0; i < num_nodes; i++)
        {
            const SizeType index = i * 6;
            range i_range(index + 3, index + 6);
            noalias( rotation_vector ) = project( rDisplacements, i_range );
            noalias( moment_vector ) = project( rForces, i_range );

            double angle = norm_2(rotation_vector);

            if(angle >= 2.0 * Globals::Pi)
                angle = std::fmod(angle, 2.0 * Globals::Pi);

            const double angle2 = angle * angle;
            const double angle4 = angle2 * angle2;
            const double angle6 = angle4 * angle2;

            double eta;
            double mu;
            if(angle < 0.05)
            {
                eta = 1.0 / 12.0 + angle2 / 270.0 + angle4 / 30240.0 + angle6 / 1209600.0;
                mu  = 1.0 / 360.0 + angle2 / 7560.0 + angle4 / 201600.0 + angle6 / 5987520.0;
            }
            else
            {
                eta = ( 1.0 - 0.5 * angle * std::tan( 0.5 * Globals::Pi - 0.5 * angle ) ) / (angle * angle);
                double sin_h_angle = std::sin(0.5 * angle);
                mu  = ( angle2 + 4.0 * std::cos(angle) + angle * std::sin(angle) - 4.0 ) / ( 4.0 * angle4 * sin_h_angle * sin_h_angle );
            }

            Spin( rotation_vector, omega );
            noalias( omega_2 ) = prod( omega, omega );

            noalias( LiTemp2 ) = outer_prod( moment_vector, rotation_vector );

            noalias( Li ) = inner_prod( rotation_vector, moment_vector ) * IdentityMatrix(3, 3);
            noalias( Li ) += outer_prod( rotation_vector, moment_vector );
            noalias( Li ) -= LiTemp2;

            noalias( LiTemp1 ) = mu * prod( omega_2, LiTemp2 );
            Spin( moment_vector, LiTemp2, 0.5 );
            noalias( LiTemp1 ) -= LiTemp2;

            noalias( LiTemp1 ) += eta * Li;

            noalias( Li ) = prod( LiTemp1, project( rH, i_range, i_range ) );

            project( L, i_range, i_range ) = Li;
        }

        return L;
    }

};

}

#endif // EICR_H_INCLUDED
