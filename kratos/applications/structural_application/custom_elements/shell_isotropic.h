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
//   Date:                $Date: 2008-02-21 09:36:17 $
//   Revision:            $Revision: 1.7 $
//
//


#if !defined(KRATOS_SHELL_ISOTROPIC_H_INCLUDED )
#define  KRATOS_SHELL_ISOTROPIC_H_INCLUDED



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
	class ShellIsotropic
		: public Element
	{
	public:
		///@name Type Definitions
		///@{
		
		/// Counted pointer of ShellIsotropic
		KRATOS_CLASS_POINTER_DEFINITION(ShellIsotropic);
		
		///@}
		///@name Life Cycle 
		///@{ 
		
		/// Default constructor.
		ShellIsotropic(IndexType NewId, GeometryType::Pointer pGeometry);
		ShellIsotropic(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
		
		/// Destructor.
		virtual ~ShellIsotropic();
		
		
		///@}
		///@name Operators 
		///@{
		
		
		///@}
		///@name Operations
		///@{
		
		Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;
		
		void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
		
		void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
			
		void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
		
		void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
		
//		void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);
		
		void GetValuesVector(Vector& values, int Step);
		void GetFirstDerivativesVector(Vector& values, int Step = 0);
		void GetSecondDerivativesVector(Vector& values, int Step = 0);
		
		void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo);

		void CalculateOnIntegrationPoints(const Variable<double >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo);
		
                void GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, 
                                           std::vector<Matrix>& rValues, 
                                           const ProcessInfo& rCurrentProcessInfo);
					   
		void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo);
		
		void Initialize();
		
		void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo);

		void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);

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
		array_1d< boost::numeric::ublas::bounded_matrix<double,3,3> , 3 > mTs;
		boost::numeric::ublas::bounded_matrix<double,3,3> mTE0;
		
		array_1d< array_1d<double,3>, 3> rot_oldit;
			
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
		
		void CalculateGaussPointContribution(
					boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system ,
					const boost::numeric::ublas::bounded_matrix<double,3,3>& Em,
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
		
		void CalculateMembraneElasticityTensor( 
					boost::numeric::ublas::bounded_matrix<double,3,3>& Em, 
					const double& h 
					);
					
		void CalculateBendingElasticityTensor( 
					boost::numeric::ublas::bounded_matrix<double,3,3>& Eb, 
     					const double& h );
					
		void CalculateAllMatrices(	
					MatrixType& rLeftHandSideMatrix,
					VectorType& rRightHandSideVector,
					ProcessInfo& rCurrentProcessInfo
					);
	
		void AddBodyForce(	
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
		
		void RotateToGlobal(
				    	const array_1d<double,3>& v1,
					const array_1d<double,3>& v2,
					const array_1d<double,3>& v3,
					const boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system,
					Matrix& rLeftHandSideMatrix,
					VectorType& rRightHandSideVector
				   );

		
		
		void NicePrint(const Matrix& A);
	
		void  AddVoigtTensorComponents(
					const double local_component,
					array_1d<double,6>& v,
					const array_1d<double,3>& a,
					const array_1d<double,3>& b
					);
		
		void CalculateAndAddKg(
				        MatrixType& LHS, 
				  	boost::numeric::ublas::bounded_matrix<double,18,18>& rWorkMatrix,
					const double& x12,
					const double& x23,
					const double& x31,
					const double& y12,
					const double& y23,
					const double& y31,
					const array_1d<double,3>& v1,
					const array_1d<double,3>& v2,
					const array_1d<double,3>& v3,
					const double& A
				);		
		
		void CalculateKg_GaussPointContribution(
					boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system ,
					const boost::numeric::ublas::bounded_matrix<double,3,3>& Em,
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
					const double& y31,
					const array_1d<double,9>& membrane_disp
					);
		
		void CalculateLocalShapeDerivatives(
				double alpha,
				boost::numeric::ublas::bounded_matrix<double,2,9>& DNu_loc ,
				boost::numeric::ublas::bounded_matrix<double,2,9>& DNv_loc ,
				boost::numeric::ublas::bounded_matrix<double,2,9>& DNw_loc ,
				const double& a, /*local coords*/ //loc1
				const double& b, //loc2
				const double& c, //loc3
				const double& x12,
				const double& x23,
				const double& x31,
				const double& y12,
				const double& y23,
				const double& y31
				);

		void CalculateProjectionOperator(
				boost::numeric::ublas::bounded_matrix<double,18,18>& rProjOperator, 
				const double& x12,
				const double& x23,
				const double& x31,
				const double& y12,
				const double& y23,
				const double& y31
				);
				
		void ApplyProjection(
				boost::numeric::ublas::bounded_matrix<double,18,18>& rLeftHandSideMatrix, 
	 			VectorType& rRightHandSideVector,
				boost::numeric::ublas::bounded_matrix<double,18,18>& rWorkMatrix, 
				array_1d<double,18>& rWorkArray,
				const boost::numeric::ublas::bounded_matrix<double,18,18>& rProjOperator
				);
		
		void UpdateNodalReferenceSystem(
				const double& x12,
				const double& x23,
				const double& x31,
				const double& y12,
				const double& y23,
				const double& y31
				);	
		
		void SaveOriginalReference(
				const array_1d<double,3>& v1,
				const array_1d<double,3>& v2,
				const array_1d<double,3>& v3
				);
		
		void CalculatePureDisplacement(
				Vector& values,
				const array_1d<double,3>& v1,
				const array_1d<double,3>& v2,
				const array_1d<double,3>& v3
				);
		
		void CalculatePureMembraneDisplacement(
				array_1d<double,9>& values,
				const array_1d<double,3>& v1,
				const array_1d<double,3>& v2,
				const array_1d<double,3>& v3
				);

		void CalculatePureBendingDisplacement(
				array_1d<double,9>& values,
				const array_1d<double,3>& v1,
				const array_1d<double,3>& v2,
				const array_1d<double,3>& v3
				);
		
		void InvertMatrix(
				const boost::numeric::ublas::bounded_matrix<double,3,3>& InputMatrix, 
				boost::numeric::ublas::bounded_matrix<double,3,3>& InvertedMatrix,
				double& InputMatrixDet
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
		ShellIsotropic(){}

		virtual void save(Serializer& rSerializer)
		{
		rSerializer.save("Name","ShellIsotropic");
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
		//ShellIsotropic& operator=(const ShellIsotropic& rOther);
		
		/// Copy constructor.
		//ShellIsotropic(const ShellIsotropic& rOther);
		
			
		///@}    
			
	}; // Class ShellIsotropic 
		
	/// input stream function
	/*  inline std::istream& operator >> (std::istream& rIStream, 
					ShellIsotropic& rThis);
	*/
	/// output stream function
	/*  inline std::ostream& operator << (std::ostream& rOStream, 
					const ShellIsotropic& rThis)
	{
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);
	
	return rOStream;
	}*/
	///@} 
	
}  // namespace Kratos.
#endif // KRATOS_SHELL_ISOTROPIC_H_INCLUDED  defined 
	
	
