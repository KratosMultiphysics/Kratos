//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(SMALL_DISPLACEMENT_INTERFACE_ELEMENT_H_INCLUDED )
#define  SMALL_DISPLACEMENT_INTERFACE_ELEMENT_H_INCLUDED


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

  /** \brief SmallDisplacementInterfaceElement
   *
   * Da aggiungere
   */
  class SmallDisplacementInterfaceElement : public Element
  {
  public:

    ///@name Type Definitions
    ///@{
    
    KRATOS_CLASS_POINTER_DEFINITION(SmallDisplacementInterfaceElement);

    ///@}

    ///@name Classes
    ///@{

    struct InterfaceIndexPermutation
	{
		bool HasPermutation;
		std::vector<SizeType> Permutation;
	};

    ///@}

    ///@name Life Cycle
    ///@{

    SmallDisplacementInterfaceElement(IndexType NewId, 
                         GeometryType::Pointer pGeometry);
    
    SmallDisplacementInterfaceElement(IndexType NewId, 
                         GeometryType::Pointer pGeometry, 
                         PropertiesType::Pointer pProperties);

    virtual ~SmallDisplacementInterfaceElement();

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

    void CalculateDampingMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo);

    // Results calculation on integration points

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

	void GetValueOnIntegrationPoints(const Variable<double>& rVariable, 
									 std::vector<double>& rValues, 
									 const ProcessInfo& rCurrentProcessInfo);
	
	void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, 
									 std::vector<Vector>& rValues, 
									 const ProcessInfo& rCurrentProcessInfo);
	
	void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, 
									 std::vector<Matrix>& rValues, 
									 const ProcessInfo& rCurrentProcessInfo);
	
	void GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable, 
									 std::vector<array_1d<double,3> >& rValues, 
									 const ProcessInfo& rCurrentProcessInfo);
	
	void GetValueOnIntegrationPoints(const Variable<array_1d<double,6> >& rVariable, 
									 std::vector<array_1d<double,6> >& rValues, 
									 const ProcessInfo& rCurrentProcessInfo);
	
	///@}

  protected:
    
    ///@name Protected Lyfe Cycle
    ///@{
    
    /**
     * Protected empty constructor
     */
    SmallDisplacementInterfaceElement() : Element()
    {
    }
    
    ///@}

  private:

    ///@name Private Operations
    ///@{

    void DecimalCorrection(Vector& a);

	void InitializeContactData();

	virtual void CalculatePermutation(InterfaceIndexPermutation& p);

	virtual void CalculateDeltaPosition(Matrix& rDeltaPosition);

	virtual void CalculateJacobianAndTransformationMatrix(const SizeType pointID,
												          Matrix& delta_position,
		                                                  Matrix& jacobian, 
												          double& J, 
												          Matrix& iR);

	virtual double CalculateIntegrationWeight(double J, 
		                                      double iw);

	virtual void CalculateLocalDisplacementVector(const InterfaceIndexPermutation& P, 
		                                          const Matrix& R,
		                                          const Vector& UG, 
												  Vector& UL);

	virtual void TransformToGlobalAndAdd(const InterfaceIndexPermutation& P,
		                                 const Matrix& R,
										 const Matrix& LHS_local,
										 const Vector& RHS_local,
										 Matrix& LHS_global,
										 Vector& RHS_global,
										 const bool LHSrequired,
										 const bool RHSrequired);

	virtual void CalculateBMatrix(const SizeType pointID, 
		                          Matrix& B);

	virtual void CalculateGeneralizedStrains(const SizeType pointID,
		                                     const Matrix& B, 
		                                     const Vector& U,
											 Vector& generalizedStrains);

    void CalculateAll(MatrixType& rLeftHandSideMatrix,
                      VectorType& rRightHandSideVector,
                      ProcessInfo& rCurrentProcessInfo,
                      const bool LHSrequired,
                      const bool RHSrequired);

	void CalculateStrain(std::vector<Vector>& strain, const ProcessInfo& rCurrentProcessInfo);
	
	void CalculateStress(std::vector<Vector>& stress, const ProcessInfo& rCurrentProcessInfo);

    ///@}

    ///@name Static Member Variables
    ///@{
    ///@}
    
    ///@name Member Variables
    ///@{

	bool mInitialized;
	std::vector< ConstitutiveLaw::Pointer > mConstitutiveLawVector;



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
#endif // SMALL_DISPLACEMENT_INTERFACE_ELEMENT_H_INCLUDED
