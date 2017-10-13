//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(OPT_TRIANGLE_ELEMENT_H_INCLUDED )
#define  OPT_TRIANGLE_ELEMENT_H_INCLUDED


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

  /** \brief OptTriangleElement
   *
   * This element represents a 3-node plane stress element
   * based on the Assumed Natural DEviatoric Strain (ANDES) by Felippa.
   */
  class OptTriangleElement : public Element
  {
  public:

    ///@name Type Definitions
    ///@{
    
    KRATOS_CLASS_POINTER_DEFINITION(OptTriangleElement);

    ///@}

    ///@name Classes
    ///@{

    // TODO: Add Calulation Data
	
    ///@}

    ///@name Life Cycle
    ///@{

    OptTriangleElement(IndexType NewId, 
                         GeometryType::Pointer pGeometry);
    
    OptTriangleElement(IndexType NewId, 
                         GeometryType::Pointer pGeometry, 
                         PropertiesType::Pointer pProperties);

    virtual ~OptTriangleElement();

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
    OptTriangleElement() : Element()
    {
    }
    
    ///@}

  private:

    ///@name Private Classes
    ///@{

    class CalculationData
    {

    public:

      // ---------------------------------------
      // calculation-constant data
      // ----------------------------------------
      // these data are allocated and constructed
      // at the beginning of the calculation

      MatrixType L;

      MatrixType Q1;
      MatrixType Q2;
      MatrixType Q3;

      MatrixType Te;
      MatrixType TTu;

      double dV;
      std::vector< array_1d<double,3> > gpLocations;

      MatrixType dNxy; /*!< shape function cartesian derivatives */

      VectorType U; /*!< global displacement vector */

      bool CalculateRHS; /*!< flag for the calculation of the right-hand-side vector */
      bool CalculateLHS; /*!< flag for the calculation of the left-hand-side vector */

      // ---------------------------------------
      // calculation-variable data
      // ---------------------------------------
      // these data are updated during the
      // calculations

      double beta0;
      size_t gpIndex;

      // ---------------------------------------
      // calculation-variable data
      // ---------------------------------------
      // these data are updated during the
      // calculations, but they are allocated
      // only once(the first time they are used)
      // to avoid useless re-allocations

      MatrixType B;   /*!< total strain-displacement matrix at the current integration point */
      MatrixType D;   /*!< section constitutive matrix at the current integration point */
      MatrixType BTD; /*!< auxiliary matrix to store the product B'*D */

      VectorType E;  /*!< generalized strain vector at the current integration point */
      VectorType S; /*!< generalized stress vector at the current integration point */

      VectorType N; /*!< shape function vector at the current integration point */

      MatrixType Q; /*!< 3x3 - stores the weighted sum of Q1, Q2 and Q3 */
      MatrixType Qh; /*!< 3x9 - the higher order B matrix */
      MatrixType TeQ; /*!< 3x3 - stores the product Te*Q */

      ConstitutiveLaw::Parameters MaterialParameters; /*!< parameters for material calculations */

	  double detF ;
	  Matrix F;

    public:

      const ProcessInfo& CurrentProcessInfo;

    public:

      CalculationData(const ProcessInfo& rCurrentProcessInfo);

    };

    ///@}

    ///@name Private Operations
    ///@{

    void DecimalCorrection(Vector& a);

    void InitializeCalculationData(CalculationData& data);

    void CalculateBMatrix(CalculationData& data);

    void CalculateBeta0(CalculationData& data);

	void CalculateConstitutiveLawResponse(CalculationData& data);

    void CalculateGaussPointContribution(CalculationData& data, MatrixType& LHS, VectorType& RHS);

    void AddBodyForces(CalculationData& data, VectorType& rRightHandSideVector);

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

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; 

    IntegrationMethod mThisIntegrationMethod; 
    
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
#endif // OPT_TRIANGLE_ELEMENT_H_INCLUDED
