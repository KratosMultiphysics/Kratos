//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:       Massimo Petracca $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            December 2013 $
//   Revision:            $Revision:                  0.0 $
//
//


#if !defined(EICR_H_INCLUDED)
#define  EICR_H_INCLUDED

#include "utilities/quaternion.h"

namespace Kratos
{
/** \brief EICR Element Independent CoRotational formulation
 *
 * E.I.C.R. is a utility class containing static methods related to
 * the Element Independent Corotational Formulation.
 * This class implements methods that do not depend on the element type,
 * and so they can be used by any implementation of a corotational coordinate transformation.
 */
class KRATOS_API(SOLID_MECHANICS_APPLICATION) EICR
{

 public:

  typedef double RealType;

  typedef BoundedMatrix<RealType, 3, 3> TransformationMatrixType;

  typedef array_1d<RealType, 3> Vector3Type;

  typedef std::vector<Vector3Type> Vector3ContainerType;

  typedef Vector VectorType;

  typedef Matrix MatrixType;

  typedef Quaternion<RealType> QuaternionType;

 public:

  /**
   * Computes the Spin of the input vector V, and saves the result into the output matrix S.
   * Note: no check is made on the size of the input-output arguments.
   * @param V the input vector (assumed size: >= 3)
   * @param S the output matrix (assumed size: >= 3x3)
   */
  template< class TVec, class TMat>
  inline static void Spin(const TVec & V, TMat & S)
  {
    S(0, 0) =	0.00;		S(0, 1) = - V(2);		S(0, 2) =   V(1);
    S(1, 0) =	V(2);		S(1, 1) =   0.00;		S(1, 2) = - V(0);
    S(2, 0) = - V(1);		S(2, 1) =   V(0);		S(2, 2) =   0.00;
  }

  /**
   * Computes the Spin of the input vector V, and saves the result into the output matrix S,
   * at the specified row index.
   * Note: no check is made on the size of the input-output arguments.
   * @param V the input vector (assumed size: >= 3)
   * @param S the output matrix (assumed size: >= 3x3)
   * @param row_index the index of the first row in the output matrix where the spin has to be saved
   */
  template< class TVec, class TMat>
  inline static void Spin_AtRow(const TVec & V, TMat & S, size_t row_index)
  {
    size_t i0 = row_index;
    size_t i1 = 1 + row_index;
    size_t i2 = 2 + row_index;
    double v0 = V(i0);
    double v1 = V(i1);
    double v2 = V(i2);
    S(i0, 0) =	0.00;		S(i0, 1) = - v2;		S(i0, 2) =   v1;
    S(i1, 0) =	v2;		S(i1, 1) =   0.00;		S(i1, 2) = - v0;
    S(i2, 0) = - v1;		S(i2, 1) =   v0;		S(i2, 2) =   0.00;
  }

  /**
   * Computes the Spin of the input vector V, from the specified index, and saves the result into the output matrix S,
   * at the specified row index.
   * Note: no check is made on the size of the input-output arguments.
   * @param V the input vector (assumed size: >= 3)
   * @param S the output matrix (assumed size: >= 3x3)
   * @param vector_index the index of the first component of the input vector to be used to compute the spin
   * @param row_index the index of the first row in the output matrix where the spin has to be saved
   */
  template< class TVec, class TMat>
  inline static void Spin_AtRow(const TVec & V, TMat & S, size_t vector_index, size_t matrix_row_index)
  {
    size_t i0 = matrix_row_index;
    size_t i1 = 1 + matrix_row_index;
    size_t i2 = 2 + matrix_row_index;
    double v0 = V(vector_index);
    double v1 = V(vector_index + 1);
    double v2 = V(vector_index + 2);
    S(i0, 0) =	0.00;		S(i0, 1) = - v2;		S(i0, 2) =   v1;
    S(i1, 0) =	v2;			S(i1, 1) =   0.00;		S(i1, 2) = - v0;
    S(i2, 0) = - v1;		S(i2, 1) =   v0;		S(i2, 2) =   0.00;
  }

  /**
   * Computes the Spin of the input vector V, and saves the result into the output matrix S.
   * This version uses a multiplier for the output values.
   * Note: no check is made on the size of the input-output arguments.
   * @param V the input vector (assumed size: >= 3)
   * @param S the output matrix (assumed size: >= 3x3)
   * @param mult the multiplier for the output values
   */
  template< class TVec, class TMat>
  inline static void Spin(const TVec & V, TMat & S, double mult)
  {
    S(0, 0) =	0.00;			S(0, 1) = - mult * V(2);	S(0, 2) =   mult * V(1);
    S(1, 0) =	mult * V(2);	S(1, 1) =   0.00;			S(1, 2) = - mult * V(0);
    S(2, 0) = - mult * V(1);	S(2, 1) =   mult * V(0);	S(2, 2) =   0.00;
  }

  /**
   * Computes the Spin of the input vector V, and saves the result into the output matrix S,
   * at the specified row index.
   * This version uses a multiplier for the output values.
   * Note: no check is made on the size of the input-output arguments.
   * @param V the input vector (assumed size: >= 3)
   * @param S the output matrix (assumed size: >= 3x3)
   * @param mult the multiplier for the output values
   * @param row_index the index of the first row in the output matrix where the spin has to be saved
   */
  template< class TVec, class TMat>
  inline static void Spin_AtRow(const TVec & V, TMat & S, double mult, size_t row_index)
  {
    size_t i0 = row_index;
    size_t i1 = 1 + row_index;
    size_t i2 = 2 + row_index;
    double v0 = mult * V(i0);
    double v1 = mult * V(i1);
    double v2 = mult * V(i2);
    S(i0, 0) =	0.00;		S(i0, 1) = - v2;		S(i0, 2) =   v1;
    S(i1, 0) =	v2;			S(i1, 1) =   0.00;		S(i1, 2) = - v0;
    S(i2, 0) = - v1;		S(i2, 1) =   v0;		S(i2, 2) =   0.00;
  }

  /**
   * Computes the Spin of the input vector V, from the specified index, and saves the result into the output matrix S,
   * at the specified row index.
   * This version uses a multiplier for the output values.
   * Note: no check is made on the size of the input-output arguments.
   * @param V the input vector (assumed size: >= 3)
   * @param S the output matrix (assumed size: >= 3x3)
   * @param mult the multiplier for the output values
   * @param vector_index the index of the first component of the input vector to be used to compute the spin
   * @param row_index the index of the first row in the output matrix where the spin has to be saved
   */
  template< class TVec, class TMat>
  inline static void Spin_AtRow(const TVec & V, TMat & S, double mult, size_t vector_index, size_t matrix_row_index)
  {
    size_t i0 = matrix_row_index;
    size_t i1 = 1 + matrix_row_index;
    size_t i2 = 2 + matrix_row_index;
    double v0 = mult * V(vector_index);
    double v1 = mult * V(vector_index + 1);
    double v2 = mult * V(vector_index + 2);
    S(i0, 0) =	0.00;		S(i0, 1) = - v2;		S(i0, 2) =   v1;
    S(i1, 0) =	v2;			S(i1, 1) =   0.00;		S(i1, 2) = - v0;
    S(i2, 0) = - v1;		S(i2, 1) =   v0;		S(i2, 2) =   0.00;
  }

 public:

  /**
   * Computes the Translational Projector Matrix.
   * The output is a square matrix of size num_nodes*6.
   * Note that 6 Degrees Of Freedom are assumed for each node.
   * @param num_nodes the number of nodes
   * @return the Translational Projector Matrix
   */
  inline static MatrixType Compute_Pt(size_t num_nodes)
  {
    RealType a = RealType(num_nodes - 1) / RealType(num_nodes);
    RealType b = -1.0 / RealType(num_nodes);

    size_t num_dofs = num_nodes * 6;

    MatrixType P( IdentityMatrix(num_dofs, num_dofs) );

    for(size_t i = 0; i < num_nodes; i++)
    {
      size_t j = i * 6;

      // diagonal block
      P(j    ,     j) = a;
      P(j + 1, j + 1) = a;
      P(j + 2, j + 2) = a;

      // out-of-diagonal block
      for(size_t k = i + 1; k < num_nodes; k++)
      {
        size_t w = k * 6;

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
   * The output is a rectangular matrix of 3 columns and nodes.size()*6 rows.
   * Note that 6 Degrees Of Freedom are assumed for each node.
   * @param nodes the input nodes
   * @return the Spin Lever Matrix
   */
  inline static MatrixType Compute_S(const Vector3ContainerType& nodes)
  {
    size_t num_nodes = nodes.size();
    size_t num_dofs = num_nodes * 6;

    MatrixType S(num_dofs, 3, 0.0);

    for(size_t i = 0; i < num_nodes; i++)
    {
      size_t j = i * 6;

      Spin_AtRow( nodes[i], S, -1.0, 0, j );

      S(j + 3, 0) = 1.0;
      S(j + 4, 1) = 1.0;
      S(j + 5, 2) = 1.0;
    }

    return S;
  }

  /**
   * Computes the Axial Vector Jacobian.
   * The output is a square matrix of size displacements.size() (which is num_nodes * 6).
   * Note that 6 Degrees Of Freedom are assumed for each node.
   * @param displacements the vector of nodal displacements and rotations in the local corotational coordinate system. (assumed size = num_nodes*6)
   * @return the H matrix
   */
  inline static MatrixType Compute_H(const VectorType & displacements)
  {
    size_t num_dofs = displacements.size();
    size_t num_nodes = num_dofs / 6;

    MatrixType H( IdentityMatrix(num_dofs, num_dofs) );

    MatrixType Omega(3, 3);
    MatrixType Hi(3, 3);

    for(size_t i = 0; i < num_nodes; i++)
    {
      size_t index = i * 6;
      Vector3Type rv = project( displacements, range(index + 3, index + 6) );

      double angle = norm_2(rv);

      if(angle >= 2.0 * Globals::Pi)
        angle = std::fmod(angle, 2.0 * Globals::Pi);

      double eta;
      if(angle < 0.05) {
        double angle2 = angle * angle;
        double angle4 = angle2 * angle2;
        double angle6 = angle4 * angle2;
        eta = 1.0 / 12.0 + 1.0 / 270.0 * angle2 + 1.0 / 30240.0 * angle4 + 1.0 / 1209600.0 * angle6;
      }
      else {
        eta = ( 1.0 - 0.5 * angle * std::tan( 0.5 * Globals::Pi - 0.5 * angle ) ) / (angle * angle);
      }

      Spin( rv, Omega );

      noalias( Hi ) = IdentityMatrix(3, 3);
      noalias( Hi ) -= 0.5 * Omega;
      noalias( Hi ) += eta * prod( Omega, Omega );

      range iRange(index + 3, index + 6);
      project( H, iRange, iRange ) = Hi;
    }

    return H;
  }

  /**
   * Computes the Spin derivative of (Axial Vector Jacobian)^T contracted with the nodal moment vector.
   * The output is a square matrix of size displacements.size() (which is num_nodes * 6).
   * Note that 6 Degrees Of Freedom are assumed for each node.
   * @param displacements the vector of nodal displacements and rotations in the local corotational coordinate system. (assumed size = num_nodes*6)
   * @param forces the vector of nodal forces and moments in the local corotational coordinate system. (assumed size = num_nodes*6)
   * @param H the Axial Vector Jacobian Matrix computed with a previous call to EICR::Compute_H(displacements)
   * @return the L matrix
   */
  inline static MatrixType Compute_L(const VectorType & displacements, const VectorType & forces, const MatrixType & H)
  {
    size_t num_dofs = displacements.size();
    size_t num_nodes = num_dofs / 6;

    MatrixType L(num_dofs, num_dofs, 0.0);

    Vector3Type rotationVector;
    Vector3Type momentVector;
    MatrixType Omega(3, 3);
    MatrixType Omega2(3, 3);
    MatrixType Li(3, 3);
    MatrixType LiTemp1(3, 3);
    MatrixType LiTemp2(3, 3);

    for(size_t i = 0; i < num_nodes; i++)
    {
      size_t index = i * 6;
      range iRange(index + 3, index + 6);
      noalias( rotationVector ) = project( displacements, iRange );
      noalias( momentVector ) = project( forces, iRange );

      double angle = norm_2(rotationVector);

      if(angle >= 2.0 * Globals::Pi)
        angle = std::fmod(angle, 2.0 * Globals::Pi);

      double angle2 = angle * angle;
      double angle4 = angle2 * angle2;
      double angle6 = angle4 * angle2;

      double eta;
      double mu;
      if(angle < 0.05) {
        eta = 1.0 / 12.0 + angle2 / 270.0 + angle4 / 30240.0 + angle6 / 1209600.0;
        mu  = 1.0 / 360.0 + angle2 / 7560.0 + angle4 / 201600.0 + angle6 / 5987520.0;
      }
      else {
        eta = ( 1.0 - 0.5 * angle * std::tan( 0.5 * Globals::Pi - 0.5 * angle ) ) / (angle * angle);
        double sin_h_angle = std::sin(0.5 * angle);
        mu  = ( angle2 + 4.0 * std::cos(angle) + angle * std::sin(angle) - 4.0 ) / ( 4.0 * angle4 * sin_h_angle * sin_h_angle );
      }

      Spin( rotationVector, Omega );
      noalias( Omega2 ) = prod( Omega, Omega );

      noalias( LiTemp2 ) = outer_prod( momentVector, rotationVector );

      noalias( Li ) = inner_prod( rotationVector, momentVector ) * IdentityMatrix(3, 3);
      noalias( Li ) += outer_prod( rotationVector, momentVector );
      noalias( Li ) -= LiTemp2;

      noalias( LiTemp1 ) = mu * prod( Omega2, LiTemp2 );
      Spin( momentVector, LiTemp2, 0.5 );
      noalias( LiTemp1 ) -= LiTemp2;

      noalias( LiTemp1 ) += eta * Li;

      noalias( Li ) = prod( LiTemp1, project( H, iRange, iRange ) );

      project( L, iRange, iRange ) = Li;
    }

    return L;
  }


};

}


#endif // EICR_H_INCLUDED
