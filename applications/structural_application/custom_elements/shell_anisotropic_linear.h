/*
==============================================================================
KratosIncompressibleFluidApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


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

//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-11-27 16:13:28 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_SHELL_ANISOTROPIC_LINEAR_H_INCLUDED )
#define  KRATOS_SHELL_ANISOTROPIC_LINEAR_H_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"

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
	class ShellAnisotropicLinear
		: public Element
	{
	public:
		///@name Type Definitions
		///@{
		
		/// Counted pointer of ShellAnisotropicLinear
		KRATOS_CLASS_POINTER_DEFINITION(ShellAnisotropicLinear);
		
		///@}
		///@name Life Cycle 
		///@{ 
		
		/// Default constructor.
		ShellAnisotropicLinear(IndexType NewId, GeometryType::Pointer pGeometry);
		ShellAnisotropicLinear(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
		
		/// Destructor.
		virtual ~ShellAnisotropicLinear();
		
		
		///@}
		///@name Operators 
		///@{
		
		
		///@}
		///@name Operations
		///@{
		
		Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;
		
		void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
		
//		void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
			
		void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
		
		void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
		
//		void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);
		
		void GetValuesVector(Vector& values, int Step);
		
		void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo);
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
			
		///@} 
		///@name Member Variables 
		///@{ 
			
		///@} 
		///@name Private Operators
		///@{ 
		void CalculateLocalGlobalTransformation(
					double& x12,
					double& x23,
					double& x31,
					double& y12,
					double& y23,
					double& y31,
					array_1d<double,3>& v1,
					array_1d<double,3>& v2,
					array_1d<double,3>& v3,
					double& area
					);
		
		void CalculateMembraneB(
					boost::numeric::ublas::bounded_matrix<double,9,3>& B,
					const double&  beta0,
					const double& loc1,
					const double& loc2,
					const double& loc3,
					const double& x12,
					const double& x23,
					const double& x31,
					const double& y12,
					const double& y23,
					const double& y31
					);
		
		
		void CalculateBendingB( 
					boost::numeric::ublas::bounded_matrix<double,9,3>& Bb,
					const double& loc2,
					const double& loc3,
					const double& x12,
					const double& x23,
					const double& x31,
					const double& y12,
					const double& y23,
					const double& y31
					);
		
		void CalculateMembraneContribution( 
					const boost::numeric::ublas::bounded_matrix<double,9,3>& Bm,
					const boost::numeric::ublas::bounded_matrix<double,3,3>& Em,
					boost::numeric::ublas::bounded_matrix<double,9,9>& Km
					);
		
		void AssembleMembraneContribution( 
					const boost::numeric::ublas::bounded_matrix<double,9,9>& Km,
					const double& coeff, 	
					boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system 
					);
				
		void CalculateBendingContribution(
					const boost::numeric::ublas::bounded_matrix<double,9,3>& Bb, 
					const boost::numeric::ublas::bounded_matrix<double,3,3>& Eb, 
					boost::numeric::ublas::bounded_matrix<double,9,9>& Kb
					);
		
		void AssembleBendingContribution( 
					const boost::numeric::ublas::bounded_matrix<double,9,9>& Kb, 
					const double& coeff, 
					boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system 
					);
		
		void CalculateMixedContribution(
					const boost::numeric::ublas::bounded_matrix<double,9,3>& Bm, 
					const boost::numeric::ublas::bounded_matrix<double,9,3>& Bb, 
					const boost::numeric::ublas::bounded_matrix<double,3,3>& B, 
					boost::numeric::ublas::bounded_matrix<double,9,9>& Kmix_top
					);
		
		void AssembleMixedContribution( 
					const boost::numeric::ublas::bounded_matrix<double,9,9>& Kmix_top, 
					const double& coeff, 
					boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system 
					);	
			
		void CalculateGaussPointContribution(
					boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system ,
					const boost::numeric::ublas::bounded_matrix<double,3,3>& Em,
					const boost::numeric::ublas::bounded_matrix<double,3,3>& B,	
					const boost::numeric::ublas::bounded_matrix<double,3,3>& Eb,	
					const double& weight,
					const double& h, /*thickness*/
					const double& loc1, /*local coords*/
					const double& loc2,
					const double& loc3,
					const double& x12,
					const double& x23,
					const double& x31,
					const double& y12,
					const double& y23,
					const double& y31
					);
								
		double CalculateBeta(	
					const boost::numeric::ublas::bounded_matrix<double,3,3>& Em 
				);
							
		void CalculateAllMatrices(	
					MatrixType& rLeftHandSideMatrix,
					VectorType& rRightHandSideVector,
					ProcessInfo& rCurrentProcessInfo
					);
	
		void ApplyBodyForce(	
				    	const double& h,
	 				const double& Area,
 					VectorType& rRightHandSideVector
					);
		
		void RotateToGlobal(
					const array_1d<double,3>& v1,
					const array_1d<double,3>& v2,
					const array_1d<double,3>& v3,
					const boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system,
					Matrix& rLeftHandSideMatrix
				);
		
		void NicePrint(const Matrix& A);
		
		void  AddVoigtTensorComponents(
					const double local_component,
					array_1d<double,6>& v,
					const array_1d<double,3>& a,
					const array_1d<double,3>& b
					);
	
	
		  
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
		ShellAnisotropicLinear(){}

		virtual void save(Serializer& rSerializer) const
		{
		    rSerializer.save("Name","ShellAnisotropicLinear");
		    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
		}

		virtual void load(Serializer& rSerializer)
		{
		   KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer,  Element );
		}
			
		///@}    
		///@name Private Inquiry 
		///@{ 
			
			
		///@}    
		///@name Un accessible methods 
		///@{ 
		
		/// Assignment operator.
		//ShellAnisotropicLinear& operator=(const ShellAnisotropicLinear& rOther);
		
		/// Copy constructor.
		//ShellAnisotropicLinear(const ShellAnisotropicLinear& rOther);
		
			
		///@}    
			
	}; // Class ShellAnisotropicLinear 
		
	/// input stream function
	/*  inline std::istream& operator >> (std::istream& rIStream, 
					ShellAnisotropicLinear& rThis);
	*/
	/// output stream function
	/*  inline std::ostream& operator << (std::ostream& rOStream, 
					const ShellAnisotropicLinear& rThis)
	{
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);
	
	return rOStream;
	}*/
	///@} 
	
}  // namespace Kratos.
#endif // KRATOS_SHELL_ANISOTROPIC_LINEAR_H_INCLUDED  defined 
	
	
