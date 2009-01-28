//   Project Name:        Kratos       
//   Last Modified by:    $Author: kazem $
//   Date:                $Date: 2007-12-14 08:49:51 $
//   Revision:            $Revision: 1.1 $
//

#if !defined(KRATOS_LINEAR_ELEMENT_H_INCLUDED )

#define  KRATOS_LINEAR_ELEMENT_H_INCLUDED
// System includes 
// External includes 
// Project includes

#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h" //for vec & mat manipulation
#include "includes/variables.h" //insert kratos varibles
#include "includes/constitutive_law.h"

namespace Kratos

{



  class LinearElement
	  : public Element
    {
    public:

      typedef ConstitutiveLaw<Node<3> > ConstitutiveLawType;
      typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

      KRATOS_CLASS_POINTER_DEFINITION(LinearElement);
      LinearElement(IndexType NewId, GeometryType::Pointer pGeometry);
      LinearElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
      /// Destructor.
      virtual ~LinearElement();

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;
      void Initialize();
      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
      void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
		// order of DOF in numbers
	  void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
		//order of DOF in x,y
	  void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo);
		// calculate stiffnes and other in element


	//  void GetValuesVector(Vector& values, int Step = 0);
	 

    protected:


    private:

      ///@name Static Member Variables 

      ///@{ 

		static Matrix msB;
		//static Matrix msF;
		static Matrix msD;
		//static Matrix msC;
		static Vector msStrainVector;
		static Vector msStressVector;
		static Matrix msDN_DX;


	  double mTotalDomainInitialSize;

		std::vector<ConstitutiveLaw<Node<3> >::Pointer> mConstitutiveLawVector;

	  std::vector< Matrix > mInvJ0;
		Vector mDetJ0;
		        

        

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
      

        

	}; // Class LinearElement 

 } // namespace Kratos.



#endif // KRATOS_LINEAR_ELEMENT_H_INCLUDED  defined 





