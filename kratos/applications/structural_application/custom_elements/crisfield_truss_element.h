/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
//   Project Name:        Kratos       
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2009-01-14 17:14:42 $
//   Revision:            $Revision: 1.2 $
//

#if !defined(KRATOS_CRISFIELD_TRUSS_ELEMENT_H_INCLUDED )
#define  KRATOS_CRISFIELD_TRUSS_ELEMENT_H_INCLUDED

// System includes 

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

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
  

  class CrisfieldTrussElement
	  : public Element
    {
    public:
      ///@name Type Definitions
      ///@{
		      
      /// Counted pointer of CrisfieldTrussElement
      KRATOS_CLASS_POINTER_DEFINITION(CrisfieldTrussElement);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
	/**
	* Constructor.
	* This deals without DOFs
	* @param NewId element ID
	* @param pGeometry geometry pointer
	* @see CrisfieldTrussElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
	*/
	CrisfieldTrussElement(IndexType NewId, GeometryType::Pointer pGeometry);
	
	/**
	* Constructor.
	* This deals with DOFs
	* @param NewId element ID
	* @param pGeometry geometry pointer
	* @param pProperties properties pointer
	* @see CrisfieldTrussElement(IndexType NewId, GeometryType::Pointer pGeometry)
	*/
	CrisfieldTrussElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

      /// Destructor.
      	virtual ~CrisfieldTrussElement();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      	
	/**
	* Create a new Crisfield truss element.
	* @return crisfield truss elment
	* @param NewId element ID
	* @param ThisNodes array of nodes
	* @param pProperties properties pointer
	*/
	Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;
	
	/**
	* Initialization of the Crisfield truss element. 
	* This initializes the cross-section, length, position vector and matrix A for the element
	*/
	void Initialize();
	
	/**
	* Calculation of the local system.
	* This calculates both the elemental stiffness matrix and the elemental residual vector
	* @param rLeftHandSideMatrix elemental stiffness matrix
	* @param rRightHandSideVector elemental residual vetor
	* @param rCurrentProcessInfo process info
	* @see CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
	* @see CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	*/
	void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
	
	/**
	* Calculation of the left hand side.
	* This calculates only the elemental stiffness matrix
	* @param rLeftHandSideMatrix elemental stiffness matrix
	* @param rCurrentProcessInfo process info
	* @see CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	* @see CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	*/
	void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
	
	/**
	* Calculation of the right hand side.
	* This calculates only the elemental residual vector
	* @param rRightHandSideVector elemental residual vetor
	* @param rCurrentProcessInfo process info
	* @see CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	* @see CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
	*/
	void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
	
	/**
	* Get the equation ID vector of the element.
	* @param rResult equation ID vector
	* @param rCurrentProcessInfo process info
	*/
	void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
	
	/**
	* Get the DOF list of the element.
	* @param ElementalDofList elemental DOF vector
	* @param rCurrentProcessInfo process info
	*/
	void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
	
	/**
	* Get the mass matrix of the element.
	* @param rMassMatrix mass matrix
	* @param rCurrentProcessInfo process info
	* TODO: assign the mass matrix
	*/
	void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
	
	/**
	* Get the damping matrix of the element.
	* @param rDampMatrix damping matrix
	* @param rCurrentProcessInfo process info
	* TODO: assign the damping matrix
	*/
	void DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);
	
	/**
	* Get the displacement vector of the element
	* @param values displacement vector
	* @param Step solution step
	* @see GetFirstDerivativesVector(Vector& values, int Step)
	* @see GetSecondDerivativesVector(Vector& values, int Step)
	*/
	void GetValuesVector(Vector& values, int Step = 0); 			
	
	/**
	* Get the velocity vector of the element
	* @param values velocity vector
	* @param Step solution step
	* @see GetValuesVector(Vector& values, int Step)
	* @see GetSecondDerivativesVector(Vector& values, int Step)
	*/
	void GetFirstDerivativesVector(Vector& values, int Step = 0); 	
	
	/**
	* Get the acceleration vector of the element
	* @param values acceleration vector
	* @param Step solution step
	* @see GetValuesVector(Vector& values, int Step)
	* @see GetFirstDerivativesVector(Vector& values, int Step)
	*/
	void GetSecondDerivativesVector(Vector& values, int Step = 0); 	
	  
// 	  void CalculateOnIntegrationPoints(const Variable<double>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo);
	
// 	  void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);

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
		Matrix msA;
		Vector msX;
		Vector msU;
		//double msStrain;
		//double msStress;
		//double msKappa; //history parameter for Scalar Damage Model
		double mLength;
		double mArea;
		int dimension;
		unsigned int number_of_nodes;

        
      ///@} 
      ///@name Member Variables 
      ///@{ 
		
//	  Geometry< Point<3,double> >::Pointer  mpReferenceGeometry; 

	  //double mTotalDomainInitialSize;
	  
		        
        
      ///@} 
      ///@name Private Operators
      ///@{ 
		
	/**
	* Auxiliary function.
	* This calculates the elemental stiffness matrix and the elemental residual vector, when required
	* @param rLeftHandSideMatrix elemental stiffness matrix
	* @param rRightHandSideVector elemental residual vector
	* @param rCurrentProcessInfo process info
	* @param CalculateStiffnessMatrixFlag flag for elemental stiffness matrix
	* @param CalculateResidualVectorFlag flag for elemental residual vector
	*/
	void CalculateAll(	MatrixType& rLeftHandSideMatrix, 
				VectorType& rRightHandSideVector, 
				ProcessInfo& rCurrentProcessInfo,
				bool CalculateStiffnessMatrixFlag,
				bool CalculateResidualVectorFlag);
	
	/**
	* Auxiliary function.
	* This adds the external force vector to the elemental residual vector
	* @param rRightHandSideVector elemental residual vector
	* @param CurrentProcessInfo process info
	*/						
	void CalculateAndAdd_ExtForce(
				VectorType& rRightHandSideVector,
				const ProcessInfo& CurrentProcessInfo);
	
	/**
	* Auxiliary function.
	* This substracts the internal force vector to the elemental residual vector
	* @param rRightHandSideVector elemental residual vector
	* @param CurrentProcessInfo process info
	* @param X position vector
	* @param U displacement vector
	* @param weight weighting factor
	*/								
	void CalculateAndMinus_IntForce(
				VectorType& rRightHandSideVector,
				const ProcessInfo& CurrentProcessInfo,
				Vector& X,
				Vector& U,
				double weight);
	
	/**
	* Auxiliary function.
	* This add the material element stiffness matrix to the elemental stiffness matrix
	* @param rLeftHandSideMatrix elemental stiffness matrix
	* @param A matrix A
	* @param X position vector
	* @param U displacement vector
	* @param weight weighting factor
	*/
	void CalculateAndAddKm(	MatrixType& rLeftHandSideMatrix,
				const Matrix& A,
				Vector& X,
				Vector& U,
				double weight);
	
	/**
	* Auxiliary function.
	* This add the geometrial element stiffness matrix to the elemental stiffness matrix
	* @param rLeftHandSideMatrix elemental stiffness matrix
	* @param A matrix A
	* @param weight weighting factor
	*/
	void CalculateAndAddKg(	MatrixType& rLeftHandSideMatrix,
				const Matrix& A,
				double weight);
	
	/**
	* Auxiliary function.
	* This calculates the GREEN-LAGRANGE strain.
	* @return GREEN-LAGRANGE strain
	* @param A matrix A
	* @param X position vector
	* @param U displacement vector
	* @param weight weighting factor
	*/						
	double CalculateStrain(	const Matrix& A,
				const Vector& X,
				const Vector& U,
				double weight );

        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
        
      ///@} 
      ///@name Private  Access 
      ///@{ 
        
        
      ///@}
      ///@name Serialization
      ///@{	
      friend class Serializer;
      
       // A private default constructor necessary for serialization 
      CrisfieldTrussElement(){};

      virtual void save(Serializer& rSerializer)
      {
        rSerializer.save("Name","CrisfieldTrussElement");
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
      }

      virtual void load(Serializer& rSerializer)
      {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
      }

      ///@}    
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      //CrisfieldTrussElement& operator=(const CrisfieldTrussElement& rOther);

      /// Copy constructor.
      //CrisfieldTrussElement(const CrisfieldTrussElement& rOther);

        
      ///@}    
        
    }; // Class CrisfieldTrussElement 

  ///@} 
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream, 
				    CrisfieldTrussElement& rThis);
*/
  /// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream, 
				    const CrisfieldTrussElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  ///@} 

}  // namespace Kratos.

#endif // KRATOS_TOTAL_CRISFIELD_ELEMENT_H_INCLUDED  defined 


