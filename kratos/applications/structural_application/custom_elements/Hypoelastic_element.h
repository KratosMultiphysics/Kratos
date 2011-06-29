//   Project Name:        Kratos       
//   Last Modified by:    $Author: kazem $
//   Date:                $Date: 2008-07-24 16:23:50 $
//   Revision:            $Revision: 1.1 $
//

#if !defined(KRATOS_HYPOELASTIC_ELEMENT_H_INCLUDED )

#define  KRATOS_HYPOELASTIC_ELEMENT_H_INCLUDED
// System includes 
// External includes 
// Project includes

#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h" //for vec & mat manipulation
#include "includes/variables.h" //insert kratos varibles
#include "includes/constitutive_law.h"

namespace Kratos

{



  class HypoelasticElement
	  : public Element
    {
    public:

      typedef ConstitutiveLaw ConstitutiveLawType;
      typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

      KRATOS_CLASS_POINTER_DEFINITION(HypoelasticElement);
      HypoelasticElement(IndexType NewId, GeometryType::Pointer pGeometry);
      HypoelasticElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
      /// Destructor.
      virtual ~HypoelasticElement();

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;
      void Initialize();
      void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);
      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
      void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);
  //    void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo);
		// order of DOF in numbers
	  void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
		//order of DOF in x,y
	 // void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo);
		// calculate stiffnes and other in element
      void GetSecondDerivativesVector(Vector& values, int Step = 0);
      void GetFirstDerivativesVector(Vector& values,int Step=0);
      void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo);
      //void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo);
      void DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);
      void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
      void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo);
     
  
	//  void GetValuesVector(Vector& values, int Step = 0);
	 

    protected:


    private:

      ///@name Static Member Variables 

      ///@{ 

		static Matrix msB;
		static Matrix msF;
		static Matrix msD;
		//static Matrix msC;
		static Vector msStrainVector;
		static Vector msStressVector;
		static Matrix msDN_DX;


	  double mTotalDomainInitialSize;

		std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

	  std::vector< Matrix > mInvJ0;
		Vector mDetJ0;
	std::vector< Matrix > mInvJ;
		Vector mDetJ;
		        

        

      ///@} 

      ///@name Private Operators

      ///@{ 

		/** K += weight*Btrans*D*B */

		void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										ProcessInfo& rCurrentProcessInfo,
										bool CalculateStiffnessMatrixFlag,
										bool CalculateResidualVectorFlag);






		void CalculateBodyForces(
								Vector& BodyForce,
								const ProcessInfo& CurrentProcessInfo
								);


		void InitializeVariables();
		
 		virtual void InitializeMaterial();

		double CalculateIntegrationWeight
								(double GaussPointWeight,
							     double DetJ0);



		void CalculateAndAdd_ExtForceContribution(
					const Vector& N,
					const ProcessInfo& CurrentProcessInfo,
					Vector& BodyForce,
					VectorType& mResidualVector,
					double weight
					);



		void CalculateStrain(
					const Matrix& B,Vector& StrainVector);


		void CalculateB(		  Matrix& B,Matrix& DN_DX
								  );
      		void Hypoconstitutive(Matrix& D,Matrix& F,const ProcessInfo& rCurrentProcessInfo);
		void CalculateStress(Vector& StrainVector, Vector& StressVector, Matrix& F, Matrix& D, const ProcessInfo& rCurrentProcessInfo);
		void CalculateRotatedStress(Vector& StressVector,Matrix& F,const ProcessInfo& rCurrentProcessInfo);
		void ResizeAndInitializeAuxiliaries();

		
	///@}
	///@name Serialization
	///@{	

	friend class Serializer; 
	// A private default constructor necessary for serialization 
	HypoelasticElement(){}

	virtual void save(Serializer& rSerializer) const
	{
	   rSerializer.save("Name","HypoelasticElement");
	   KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
	}

	virtual void load(Serializer& rSerializer)
	{
	   KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,  Element );
	}
		
		

	}; // Class HypoelasticElement 

 } // namespace Kratos.



#endif // KRATOS_HYPOELASTIC_ELEMENT_H_INCLUDED  defined 





