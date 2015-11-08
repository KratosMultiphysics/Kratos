//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//


#if !defined(AGQ4_ELEMENT_V2_H_INCLUDED )
#define  AGQ4_ELEMENT_V2_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/element.h"

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

  /** \brief AGQ4Element
   */
  class AGQ4Element : public Element
  {
  public:

    ///@name Type Definitions
    ///@{
    
    KRATOS_CLASS_POINTER_DEFINITION(AGQ4Element);

    typedef array_1d<double, 3> Vector3Type;

    ///@}

    ///@name Life Cycle
    ///@{

    AGQ4Element(IndexType NewId, 
                   GeometryType::Pointer pGeometry);
    
    AGQ4Element(IndexType NewId, 
                   GeometryType::Pointer pGeometry, 
                   PropertiesType::Pointer pProperties);

    virtual ~AGQ4Element();

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
    AGQ4Element() : Element()
    {
    }
    
    ///@}

  private:

    ///@name Private Operations
    ///@{

	Matrix& CalculateDeltaPosition(Matrix & rDeltaPosition);

    void DecimalCorrection(Vector& a);

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
	
	// standard element members
	std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; 
	IntegrationMethod mThisIntegrationMethod;
	double mErrorCode;
	bool   mInitialized;
	
	// members for non-linear treatement of internal dofs
	Vector m_Q;
	Vector m_Q_converged;
	Vector m_U;
	Vector m_U_converged; 
	Vector m_Q_residual;
	Matrix m_KQQ_inv;
	Matrix m_KQU; // L = G'*C*B
	Matrix m_KUQ; // L^T = B'*C'*G
	
#ifdef EASQ4_DRILLING_DOF
	// members for drilling dofs
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
#endif // AGQ4_ELEMENT_V2_H_INCLUDED
