//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(Q4RISTAB_ELEMENT_H_INCLUDED )
#define  Q4RISTAB_ELEMENT_H_INCLUDED


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

  /** \brief Q4RIStabElement
   *
   * This element represents a 4-node plane stress element
   * with reduced integration plus hourglass stabilization.
   */
  class Q4RIStabElement : public Element
  {
  public:

    ///@name Type Definitions
    ///@{
    
    KRATOS_CLASS_POINTER_DEFINITION(Q4RIStabElement);

    ///@}

    ///@name Classes
    ///@{

    // TODO: Add Calulation Data
	
    ///@}

    ///@name Life Cycle
    ///@{

    Q4RIStabElement(IndexType NewId, 
                         GeometryType::Pointer pGeometry);
    
    Q4RIStabElement(IndexType NewId, 
                         GeometryType::Pointer pGeometry, 
                         PropertiesType::Pointer pProperties);

    virtual ~Q4RIStabElement();

    ///@}
    
    ///@name Operations
    ///@{

    // Basic

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    IntegrationMethod GetIntegrationMethod() const;
    
    void Initialize();

    void ResetConstitutiveLaw();

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

    // Get/Set values on/from integration points

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

	void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
									  std::vector<double>& rOutput,
									  const ProcessInfo& rCurrentProcessInfo);

	void CalculateOnIntegrationPoints(const Variable< array_1d< double, 3 > >& rVariable,
									  std::vector< array_1d<double, 3 > >& rOutput,
									  const ProcessInfo& rCurrentProcessInfo);

	void CalculateOnIntegrationPoints(const Variable< Vector >& rVariable,
									  std::vector< Vector >& rOutput,
									  const ProcessInfo& rCurrentProcessInfo);

	void CalculateOnIntegrationPoints(const Variable< Matrix >& rVariable,
									  std::vector< Matrix >& rOutput,
									  const ProcessInfo& rCurrentProcessInfo);

    ///@}

  protected:
    
    ///@name Protected Lyfe Cycle
    ///@{
    
    /**
     * Protected empty constructor
     */
    Q4RIStabElement() : Element()
    {
    }
    
    ///@}

  private:

    ///@name Private Classes
    ///@{

    ///@}

    ///@name Private Operations
    ///@{
	
	void Precompute(ProcessInfo& CurrentProcessInfo);

    void AddBodyForces(double V, VectorType& rRightHandSideVector);

    void CalculateAll(MatrixType& rLeftHandSideMatrix,
                      VectorType& rRightHandSideVector,
                      ProcessInfo& rCurrentProcessInfo,
                      const bool LHSrequired,
                      const bool RHSrequired);
    ///@}

    ///@name Static Member Variables
    ///@{
    ///@}
    
    ///@name Member Variables
    ///@{

    ConstitutiveLaw::Pointer mConstitutiveLaw;
	Matrix m_H;
	double m_dQ0;
	double m_dQn;
	double m_dQ_trial;
	bool   m_H_calculated;

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
#endif // Q4RISTAB_ELEMENT_H_INCLUDED
