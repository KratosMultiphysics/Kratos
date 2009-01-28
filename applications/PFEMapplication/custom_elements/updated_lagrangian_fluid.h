//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-01-23 15:38:00 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_UPDATED_LAGRANGIAN_FLUID_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_FLUID_H_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
//#include "constitutive_laws/constitutive_laws.h"

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
  
  /// Short class definition.
  /** Detail class definition.
  */
  class UpdatedLagrangianFluid
	  : public Element
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of UpdatedLagrangianFluid
      KRATOS_CLASS_POINTER_DEFINITION(UpdatedLagrangianFluid);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
	  UpdatedLagrangianFluid(IndexType NewId, GeometryType::Pointer pGeometry);
      UpdatedLagrangianFluid(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~UpdatedLagrangianFluid();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

//	  void Initialize();

  	  void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
	  //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

	  void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
	  
	  void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);

	  void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);

	  void DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);
	  
	  void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo);
	  //void Calculate(array_1d<double,3>& rVariable, array_1d<double,3>& Output, const ProcessInfo& rCurrentProcessInfo);
	  void Calculate(const Variable< array_1d<double,3> >& rVariable, array_1d<double,3>& Output, const ProcessInfo& rCurrentProcessInfo);
  	  //void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);

	  void GetValuesVector(Vector& values, int Step = 0);
	  void GetFirstDerivativesVector(Vector& values, int Step = 0);
	  void GetSecondDerivativesVector(Vector& values, int Step = 0);


      ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
//      virtual String Info() const;
      
      /// Print information about this object.
//      virtual void PrintInfo(std::ostream& rOStream) const;

      /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;
      
            
      ///@}      
      ///@name Friends
      ///@{

            
      ///@}
      
    protected:
      ///@name Protected static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operators
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operations
      ///@{ 
        
        
      ///@} 
      ///@name Protected  Access 
      ///@{ 
        
        
      ///@}      
      ///@name Protected Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Protected LifeCycle 
      ///@{ 
      
            
      ///@}
      
    private:
      ///@name Static Member Variables 
      ///@{ 
		static Matrix msB;
		static Matrix msF;
		static Matrix msD;
		static Matrix msC;
		//approximate constitutive tensor
		static Matrix msCapx;
		static Vector msStrainVector;
		static Vector msStressVector;

		static Matrix msDN_DX;
        
      ///@} 
      ///@name Member Variables 
      ///@{ 

	//  std::vector<ConstitutiveLaws::Pointer> mConstitutiveLawVector;
		
	  Geometry< Point<3,double> >::Pointer  mpReferenceGeometry; 

	  double mTotalDomainInitialSize;
	  std::vector< Matrix > mInvJ0;
	  Vector mDetJ0;
		        
        
      ///@} 
      ///@name Private Operators
      ///@{ 
		/** K += weight*Btrans*D*B */
		
		//this function computes the Cauchy stress in the fluid which will
		//consist of the "shear part" (velocity grads) and the volumetric part
		//and then transforms the Cauchy tensor into PK2 stress tensor
		//then using Voigt rule we store PK2 stress as a vector
		void CalculateIncPK2Stress(double& BulkModulus, Vector& StressVector, Matrix& F);
		
		void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										ProcessInfo& rCurrentProcessInfo,
										bool CalculateStiffnessMatrixFlag,
										bool CalculateResidualVectorFlag);

		void CalculateAndAddKm(
									MatrixType& K,
									Matrix& B,
									Matrix& D,
									double weight);
		
		/** Calculation of the Geometric Stiffness Matrix. Kg = dB * S*/
		void CalculateAndAddKg(
					MatrixType& K,
					Matrix& DN_DX,
					Vector& StressVector,
					double weight
					);
		
		
		

		void CalculateBodyForces(
								Vector& BodyForce,
								const ProcessInfo& CurrentProcessInfo
								);


		void InitializeVariables();
		
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

		void CalculateStrain(const Matrix& C,
							Vector& StrainVector);

		void CalculateB(		  Matrix& B,
								  Matrix& F,
								  Matrix& DN_DX,
								  unsigned int StrainSize);
      ///@} 
      ///@name Private Operations
      ///@{ 
        
        
      ///@} 
      ///@name Private  Access 
      ///@{ 
        
        
      ///@}    
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      //UpdatedLagrangianFluid& operator=(const UpdatedLagrangianFluid& rOther);

      /// Copy constructor.
      //UpdatedLagrangianFluid(const UpdatedLagrangianFluid& rOther);

        
      ///@}    
        
    }; // Class UpdatedLagrangianFluid 

  ///@} 
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream, 
				    UpdatedLagrangianFluid& rThis);
*/
  /// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream, 
				    const UpdatedLagrangianFluid& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  ///@} 

}  // namespace Kratos.

#endif // KRATOS_UPDATED_LAGRANGIAN_FLUID_H_INCLUDED  defined 


