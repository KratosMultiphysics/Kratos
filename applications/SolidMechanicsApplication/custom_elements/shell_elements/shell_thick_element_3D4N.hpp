//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:       Massimo Petracca $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:           September 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(SHELL_THICK_ELEMENT_3D4N_H_INCLUDED)
#define  SHELL_THICK_ELEMENT_3D4N_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "custom_utilities/shell_cross_section.hpp"
#include "utilities/quaternion.h"

namespace Kratos
{


///@name Kratos Globals
///@{
///@}

///@name Type Definitions
///@{
///@}

class ShellQ4_CoordinateTransformation;

class ShellQ4_LocalCoordinateSystem;

///@name  Enum's
///@{
///@}

///@name  Functions
///@{
///@}

///@name Kratos Classes
///@{

/** \brief ShellThickElement3D4N
 *
 * This element represents a 4-node bilinear Shell element
 * based on the Enhanced Assumed Strain Method (E.A.S.) for the membrane part
 * and on the Mixed Interpolation of Tensorial Components (M.I.T.C.)
 * for the trasverse shear part.
 * This element is formulated for small strains,
 * but can be used in Geometrically nonlinear problems
 * involving large displacements and rotations
 * using a Corotational Coordinate Transformation.
 * Material nonlinearity is handled by means of the cross section object.
 */
class KRATOS_API(SOLID_MECHANICS_APPLICATION) ShellThickElement3D4N : public Element
{
 public:

  ///@name Type Definitions
  ///@{

  KRATOS_CLASS_POINTER_DEFINITION(ShellThickElement3D4N);

  typedef std::vector< ShellCrossSection::Pointer > CrossSectionContainerType;

  typedef ShellQ4_CoordinateTransformation CoordinateTransformationBaseType;

  typedef Kratos::shared_ptr<CoordinateTransformationBaseType> CoordinateTransformationBasePointerType;

  typedef array_1d<double, 3> Vector3Type;

  typedef Quaternion<double> QuaternionType;

  ///@}

  ///@name Classes
  ///@{

  /** \brief JacobianOperator
   *
   * This class is a utility to compute at a given integration point,
   * the Jacobian, its inverse, its determinant
   * and the derivatives of the shape functions in the local
   * cartesian coordinate system.
   */
  class JacobianOperator
  {
   public:

    JacobianOperator();

    void Calculate(const ShellQ4_LocalCoordinateSystem & CS, const Matrix & dN);

    inline const Matrix & Jacobian()const { return mJac; }

    inline const Matrix & Inverse()const { return mInv; }

    inline const Matrix & XYDerivatives()const { return mXYDeriv; }

    inline const double Determinant()const { return mDet; }

   private:

    Matrix mJac;     /*!< Jacobian matrix */
    Matrix mInv;     /*!< Inverse of the Jacobian matrix */
    Matrix mXYDeriv; /*!< Shape function derivatives in cartesian coordinates */
    double mDet;     /*!< Determinant of the Jacobian matrix */
  };

  /** \brief MITC4Params
   *
   * This class performs some operations and stores some data to compute
   * the transverse shear contribution of the stiffness matrix using the
   * M.I.T.C. formulation.
   *
   * References:
   * -   Dvorkin,Bathe, "A continuum mechanics based four node shell
   *     element for general nonlinear analysis",
   *     Eng.Comput.,vol. 1, 77-88, 1984
   * -   Bathe, Dvorkin, "Short communication A four-node plate bending element
   *     based on Mindlin/Reissner plate theory and a Mixed Interpolation",
   *     International Journal for Numerical Methods in Eng.,
   *     vol. 21, 367-383, 1985
   */
  struct MITC4Params
  {

    double Ax;
    double Ay;
    double Bx;
    double By;
    double Cx;
    double Cy;
    Matrix Transformation;
    Matrix ShearStrains;

    MITC4Params(const ShellQ4_LocalCoordinateSystem & LCS);

  };

  class EASOperator; // forward declaration

  /** \brief EASOperatorStorage
   *
   * This class is meant to store persistent data for the EAS calculations.
   * This class is stored in the element and used by the EASOperator.
   */
  class EASOperatorStorage
  {

   public:

    friend class EASOperator;

    typedef Element::GeometryType GeometryType;

   public:

    EASOperatorStorage();

    inline void Initialize(const GeometryType& geom);

    inline void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

    inline void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);

    inline void FinalizeNonLinearIteration(const Vector& displacementVector, ProcessInfo& CurrentProcessInfo);

   private:

    array_1d<double, 5> alpha;              /*!< (trial) vector containing the 5 enhanced strain parameters */
    array_1d<double, 5> alpha_converged;    /*!< (converged) vector containing the 5 enhanced strain parameters */

    array_1d<double, 24> displ;             /*!< (trial) vector containing the displacement vector */
    array_1d<double, 24> displ_converged;   /*!< (converged) vector containing the displacement vector */

    array_1d<double, 5>           residual; /*!< vector containing the 5 residuals for the 5 enhanced strain parameters */
    BoundedMatrix<double, 5, 5>  Hinv;     /*!< 5x5 matrix that stores H^-1 */
    BoundedMatrix<double, 5, 24> L;        /*!< 5x24 coupling matrix */

    bool mInitialized;                      /*!< Initialization flag */

   private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

  };

  /** \brief EASOperator
   *
   * This class performs some operations and stores some data to compute
   * the membrane contribution of the stiffness matrix using the
   * Enhanced Assumed Strain formulation.
   *
   * References:
   * -   J.C.Simo,M.S.Rifai, "A class of mixed assumed strain methods
   *     and the method of incompatible modes",
   *     International Journal for Numerical Methods in Eng.,
   *     vol. 29, 1595-1638, 1990
   */
  class EASOperator
  {

   public:

    /**
     * Constructor
     */
    EASOperator(const ShellQ4_LocalCoordinateSystem& LCS, EASOperatorStorage& storage);

   public:

    /**
     * this method should be called in the Gauss Loop
     * after the standard strain-displacement matrix has been computed, as well as the standard
     * generalized strains, but before the computation of the constitutive laws.
     */
    inline void GaussPointComputation_Step1(double xi, double eta, const JacobianOperator& jac,
                                            Vector& generalizedStrains,
                                            EASOperatorStorage& storage);

    /**
     * this method should be called in the Gauss Loop
     * after the standard computation of the constitutive laws.
     * note:
     * the input algorithmic tangent and generalized stress vector are assumed to be already multiplied
     * by the integration weight, and their size is assumed to be those of a standard shell element
     * (i.e. algorithmicTangent(8x8), generalizedStresses(8x1))
     */
    inline void GaussPointComputation_Step2(const Matrix& D,
                                            const Matrix& B,
                                            const Vector& S,
                                            EASOperatorStorage& storage);

    /**
     * this method should be called at the end of the Gauss Loop,
     * when the integration is terminated, but before transforming everything
     * to the global system: Here we are still operating in the element local
     * coordinate system.
     */
    inline void ComputeModfiedTangentAndResidual(Matrix& rLeftHandSideMatrix,
                                                 Vector& rRightHandSideVector,
                                                 EASOperatorStorage& storage);

   private:

    Matrix mF0inv;           /*!< 3x3 inverse deformation matrix at the element center */
    double mJ0;              /*!< determinant of the jacobian at the element center */
    Vector mEnhancedStrains; /*!< vector of 3 enhanced strains [e.xx, e.yy, 2e.xy] */
    Matrix mG;               /*!< 3x5 interpolation matrix in cartesian coordinates */
  };

  ///@}

  ///@name Life Cycle
  ///@{

  ShellThickElement3D4N(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        bool NLGeom = false);

  ShellThickElement3D4N(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties,
                        bool NLGeom = false);

  ShellThickElement3D4N(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties,
                        CoordinateTransformationBasePointerType pCoordinateTransformation);

  ~ShellThickElement3D4N() override;

  ///@}

  ///@name Operations
  ///@{

  // Basic

  Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

  void Initialize() override;

  void ResetConstitutiveLaw() override;

  void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

  void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) override;

  int Check(const ProcessInfo& rCurrentProcessInfo) override;

  void CleanMemory() override;

  void GetValuesVector(Vector& values, int Step = 0) override;

  void GetFirstDerivativesVector(Vector& values, int Step = 0) override;

  void GetSecondDerivativesVector(Vector& values, int Step = 0) override;

  void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo) override;

  void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo) override;

  void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

  void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

  void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;

  void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override;

  void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                            VectorType& rRightHandSideVector,
                            ProcessInfo& rCurrentProcessInfo) override;

  void CalculateRightHandSide(VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override;

  // Results calculation on integration points

  void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  void GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable, std::vector<array_1d<double,3> >& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  void GetValueOnIntegrationPoints(const Variable<array_1d<double,6> >& rVariable, std::vector<array_1d<double,6> >& rValues, const ProcessInfo& rCurrentProcessInfo) override;

  ///@}

  ///@name Public specialized Access - Temporary
  ///@{

  void SetCrossSectionsOnIntegrationPoints(std::vector< ShellCrossSection::Pointer >& crossSections);

  ///@}

 protected:

  ///@name Protected Lyfe Cycle
  ///@{

  /**
   * Protected empty constructor
   */
  ShellThickElement3D4N() : Element()
  {
  }

  ///@}

 private:

  ///@name Private Operations
  ///@{

  void DecimalCorrection(Vector& a);

  void SetupOrientationAngles();

  void CalculateBMatrix(double xi, double eta,
                        const JacobianOperator& Jac, const MITC4Params& params,
                        const Vector& N,
                        Matrix& B, Vector& Bdrill);

  void CalculateAll(MatrixType& rLeftHandSideMatrix,
                    VectorType& rRightHandSideVector,
                    ProcessInfo& rCurrentProcessInfo,
                    const bool LHSrequired,
                    const bool RHSrequired);

  void AddBodyForces(const array_1d<double,4> & dA, VectorType& rRightHandSideVector);

  bool TryGetValueOnIntegrationPoints_MaterialOrientation(const Variable<array_1d<double,3> >& rVariable,
                                                          std::vector<array_1d<double,3> >& rValues,
                                                          const ProcessInfo& rCurrentProcessInfo);

  bool TryGetValueOnIntegrationPoints_GeneralizedStrainsOrStresses(const Variable<Matrix>& rVariable,
                                                                   std::vector<Matrix>& rValues,
                                                                   const ProcessInfo& rCurrentProcessInfo);

  ///@}

  ///@name Static Member Variables
  ///@{
  ///@}

  ///@name Member Variables
  ///@{

  CoordinateTransformationBasePointerType mpCoordinateTransformation; /*!< The Coordinate Transformation */

  CrossSectionContainerType mSections; /*!< Container for cross section associated to each integration point */

  EASOperatorStorage mEASStorage; /*!< The storage instance for the EAS Operator */

  ///@}

  ///@name Serialization
  ///@{

  friend class Serializer;

  void save(Serializer& rSerializer) const override;

  void load(Serializer& rSerializer) override;

  ///@}

  ///@name Private  Access
  ///@{
  ///@}

  ///@name Private Inquiry
  ///@{
  ///@}

  ///@name Un accessible methods
  ///@{
  ///@}

};

}
#endif // SHELL_THICK_ELEMENT_3D4N_H_INCLUDED
