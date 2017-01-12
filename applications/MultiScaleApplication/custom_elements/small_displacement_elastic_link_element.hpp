//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//


#if !defined(SMALL_DISPLACEMENT_ELASTIC_LINK_ELEMENT_H_INCLUDED )
#define  SMALL_DISPLACEMENT_ELASTIC_LINK_ELEMENT_H_INCLUDED


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

  /** \brief SmallDisplacementElasticLinkElement
   */
  class SmallDisplacementElasticLinkElement : public Element
  {
  public:

    ///@name Type Definitions
    ///@{
    
    KRATOS_CLASS_POINTER_DEFINITION(SmallDisplacementElasticLinkElement);

    ///@}

    ///@name Life Cycle
    ///@{

    SmallDisplacementElasticLinkElement(IndexType NewId, 
                   GeometryType::Pointer pGeometry);
    
    SmallDisplacementElasticLinkElement(IndexType NewId, 
                   GeometryType::Pointer pGeometry, 
                   PropertiesType::Pointer pProperties);

    virtual ~SmallDisplacementElasticLinkElement();

    ///@}
    
    ///@name Operations
    ///@{

    // Basic

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

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
    SmallDisplacementElasticLinkElement() : Element()
    {
    }
    
    ///@}

  private:

    ///@name Private Operations
    ///@{

	void CalculateLocalAxes(array_1d<double,3>& vx, array_1d<double,3>& vy, array_1d<double,3>& vz, double& L);
	
    ///@}

    ///@name Static Member Variables
    ///@{
    ///@}
    
    ///@name Member Variables
    ///@{
	
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
#endif // SMALL_DISPLACEMENT_ELASTIC_LINK_ELEMENT_H_INCLUDED
