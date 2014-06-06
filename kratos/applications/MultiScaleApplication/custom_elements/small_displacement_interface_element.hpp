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

    struct ContactData
	{
	public:
		double GapT1;
		double GapT2;
		double LastGapT1;
		double LastGapT2;
		bool IsSliding;
		double Tau1;
		double Tau2;
		double LastTau1;
		double LastTau2;

	private:
		bool WasSliding;

	public:
		ContactData()
			: GapT1(0.0)
			, GapT2(0.0)
			, LastGapT1(0.0)
			, LastGapT2(0.0)
            , IsSliding(false)
			, Tau1(0.0)
			, Tau2(0.0)
			, LastTau1(0.0)
			, LastTau2(0.0)
            , WasSliding(false)
		{
		}

		inline void FinalizeSolutionStep()
		{
			/*if(WasSliding)
			{
				if(!IsSliding)
				{
					LastGapT1 = GapT1;
					LastGapT2 = GapT2;

					LastTau1 = Tau1;
					LastTau2 = Tau2;
				}
			}
			WasSliding = IsSliding;*/

			if(IsSliding)
			{
				LastGapT1 = GapT1;
				LastGapT2 = GapT2;

				LastTau1 = Tau1;
				LastTau2 = Tau2;
			}

			/*LastGapT1 = GapT1;
			LastGapT2 = GapT2;

			LastTau1 = Tau1;
			LastTau2 = Tau2;*/
		}

		friend class Serializer;

		void save(Serializer& rSerializer) const
		{
			rSerializer.save("sp1",  GapT1);
			rSerializer.save("sp2",  GapT2);
			rSerializer.save("spc1", LastGapT1);
			rSerializer.save("spc2", LastGapT2);
			rSerializer.save("bs",  IsSliding);
		}

		void load(Serializer& rSerializer)
		{
			rSerializer.load("sp1",  GapT1);
			rSerializer.load("sp2",  GapT2);
			rSerializer.load("spc1", LastGapT1);
			rSerializer.load("spc2", LastGapT2);
			rSerializer.load("bs",  IsSliding);
		}
	};

    struct InterfaceIndexPermutation
	{
		bool HasPermutation;
		std::vector<SizeType> Permutation;
	};

	typedef std::vector<ContactData> ContactDataCollectionType;

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

	void CalculateOnIntegrationPoints(const Variable<Vector >& rVariable, std::vector< Vector >& rOutput, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable, std::vector<array_1d<double,3> >& rValues, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<array_1d<double,6> >& rVariable, std::vector<array_1d<double,6> >& rValues, const ProcessInfo& rCurrentProcessInfo);

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
												          Matrix& iR, 
												          Matrix& R);

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
										 Vector& RHS_global);

	virtual void CalculateBMatrix(const SizeType pointID, 
		                          Matrix& B);

	virtual void CalculateGeneralizedStrains(const SizeType pointID,
		                                     const Matrix& B, 
		                                     const Vector& U,
											 Vector& generalizedStrains);

	// OUTDATED
	void CalcMat(const Vector& dU, 
		         Vector& sigma, 
				 Matrix& D);

	// OUTADATED
	void CalculateCoulombFrictionLaw(const Vector& sigma,
		                             double& phi,
									 double& norm_tau,
									 double& Fs);

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

	ContactDataCollectionType mContactData;
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
