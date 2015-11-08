//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//


#if !defined(EAS_QUAD_ELEMENT_V2_H_INCLUDED )
#define  EAS_QUAD_ELEMENT_V2_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/element.h"

/* Drilling dof flag */
//#define EASQ4_DRILLING_DOF

namespace Kratos
{


  ///@name Kratos Globals
  ///@{
  ///@}

  ///@name Type Definitions
  ///@{
  ///@}

  ///@name  Enum's
  ///@{
  ///@}

  ///@name  Functions
  ///@{
  ///@}

  ///@name Kratos Classes
  ///@{

  /** \brief EASQuadElementV2
   */
  class EASQuadElementV2 : public Element
  {
  public:

    ///@name Type Definitions
    ///@{
    
    KRATOS_CLASS_POINTER_DEFINITION(EASQuadElementV2);

    typedef array_1d<double, 3> Vector3Type;

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

      void Calculate(const Element::GeometryType & geom, const Matrix & dN);

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

    class EASOperator; // forward declaration

    /** \brief EASOperatorStorage
     *
     * This class is meant to store persistent data for the EAS calculations.
     * This class is stored in the element and used by the EASOperator.
     */
    class EASOperatorStorage
    {

    public:

	  friend class EASQuadElementV2;
      friend class EASOperator;

      typedef Element::GeometryType GeometryType;

    public:

      EASOperatorStorage();

      void Initialize(const GeometryType& geom);

      void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

      void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);

      void FinalizeNonLinearIteration(const Vector& displacementVector, ProcessInfo& CurrentProcessInfo);

	public:

		inline const Matrix& H()const { return Hinv; }

    private:

      Vector alpha;
      Vector alpha_converged;

      Vector displ;
      Vector displ_converged; 

      Vector residual;
      Matrix Hinv;
      Matrix L;
	  Matrix LT;

      bool   mInitialized;

    private:

      friend class Serializer;

      virtual void save(Serializer& rSerializer) const
      {
		rSerializer.save("A0", alpha);
		rSerializer.save("A1", alpha_converged);
		rSerializer.save("U0", displ);
		rSerializer.save("U1", displ_converged);
		rSerializer.save("init", mInitialized);
      }

      virtual void load(Serializer& rSerializer)
      {
	rSerializer.load("A0", alpha);
		rSerializer.load("A1", alpha_converged);
		rSerializer.load("U0", displ);
		rSerializer.load("U1", displ_converged);
		rSerializer.load("init", mInitialized);
      }
    };

    /** \brief EASOperator
     *
     * This class performs some operations and stores some data for the
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
      EASOperator(const GeometryType& geom, EASOperatorStorage& storage);

    public:

      /**
       * this method should be called in the Gauss Loop
       * after the standard strain-displacement matrix has been computed, as well as the strains, 
	   * but before the computation of the constitutive laws.
       */
      void GaussPointComputation_Step1(double xi, double eta, const JacobianOperator& jac, 
					      Vector& strainVector,
					      EASOperatorStorage& storage);

      /**
       * this method should be called in the Gauss Loop
       * after the standard computation of the constitutive laws.
       * note:
       * the input algorithmic tangent and stress vector are assumed to be already multiplied
       * by the integration weight, and their size is assumed to be those of a standard 2D planestress element
       * (i.e. algorithmicTangent(3x3), generalizedStresses(3x1))
       */
      void GaussPointComputation_Step2(const Matrix& D, 
					      const Matrix& B,
					      const Vector& S,
					      EASOperatorStorage& storage);

      /**
       * this method should be called at the end of the Gauss Loop,
       * when the integration is terminated.
       */
      bool ComputeModfiedTangentAndResidual(Matrix& rLeftHandSideMatrix,
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

    EASQuadElementV2(IndexType NewId, 
                   GeometryType::Pointer pGeometry);
    
    EASQuadElementV2(IndexType NewId, 
                   GeometryType::Pointer pGeometry, 
                   PropertiesType::Pointer pProperties);

    virtual ~EASQuadElementV2();

    ///@}
    
    ///@name Operations
    ///@{

    // Basic

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    void Initialize();

    void ResetConstitutiveLaw();

	IntegrationMethod GetIntegrationMethod() const;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo);

    int Check(const ProcessInfo& rCurrentProcessInfo);

    void CleanMemory();

    void GetValuesVector(Vector& values, int Step = 0);

    void GetFirstDerivativesVector(Vector& values, int Step = 0);
    
    void GetSecondDerivativesVector(Vector& values, int Step = 0);

    void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo);

    void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo);

    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

    void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);

    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo);

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo);

    // Results calculation on integration points

	void SetValueOnIntegrationPoints(const Variable<double>& rVariable, 
									 std::vector<double>& rValues, 
									 const ProcessInfo& rCurrentProcessInfo);

    void SetValueOnIntegrationPoints(const Variable<Vector>& rVariable, 
									 std::vector<Vector>& rValues, 
									 const ProcessInfo& rCurrentProcessInfo);

    void SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, 
									 std::vector<Matrix>& rValues, 
									 const ProcessInfo& rCurrentProcessInfo);

    void SetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
									 std::vector<ConstitutiveLaw::Pointer>& rValues,
									 const ProcessInfo& rCurrentProcessInfo );

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, 
									 std::vector<double>& rValues, 
									 const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, 
									 std::vector<Vector>& rValues, 
									 const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, 
									 std::vector<Matrix>& rValues, 
									 const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
									 std::vector<ConstitutiveLaw::Pointer>& rValues,
									 const ProcessInfo& rCurrentProcessInfo);

    ///@}

  protected:
    
    ///@name Protected Lyfe Cycle
    ///@{
    
    /**
     * Protected empty constructor
     */
    EASQuadElementV2() : Element()
    {
    }
    
    ///@}

  private:

    ///@name Private Operations
    ///@{

    void DecimalCorrection(Vector& a);

    void CalculateBMatrix(double xi, double eta, const JacobianOperator& Jac, const Vector& N, Matrix& B);

    void CalculateAll(MatrixType& rLeftHandSideMatrix,
                      VectorType& rRightHandSideVector,
                      ProcessInfo& rCurrentProcessInfo,
                      const bool LHSrequired,
                      const bool RHSrequired);

	void AddBodyForces(VectorType& rRightHandSideVector);

    ///@}

    ///@name Static Member Variables
    ///@{
    ///@}
    
    ///@name Member Variables
    ///@{
    
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; 

    IntegrationMethod mThisIntegrationMethod;
    
    EASOperatorStorage mEASStorage;

	double mErrorCode;

#ifdef EASQ4_DRILLING_DOF
	double mG0;
	bool   m_first_step_finalized;
#endif

    ///@}
    
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);
    
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
#endif // EAS_QUAD_ELEMENT_V2_H_INCLUDED
