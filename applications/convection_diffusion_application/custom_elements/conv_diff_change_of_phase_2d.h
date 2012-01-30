/*
==============================================================================
KratosConvectionDiffusionApplication 
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
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-12-15 15:42:14 $
//   Revision:            $Revision: 1.1 $
//
//
 

#if !defined(KRATOS_TRIANGULAR_CONVDIFF_2_ELEM_H_INCLUDED )
#define  KRATOS_TRIANGULAR_CONVDIFF_2_ELEM_H_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

#include "includes/serializer.h"
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
  class ConvDiffChangeOfPhase2D
	  : public Element
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of ConvDiffChangeOfPhase2D
      KRATOS_CLASS_POINTER_DEFINITION(ConvDiffChangeOfPhase2D);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
	  ConvDiffChangeOfPhase2D(IndexType NewId, GeometryType::Pointer pGeometry);
      ConvDiffChangeOfPhase2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~ConvDiffChangeOfPhase2D();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      //void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

	void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
	
	void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

	void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

	double AA(double, double, double, double, double, double);
	
	Vector Z(Matrix, double, double, double, double, double, double, double, double, double);
	
	void CC(double, double, double, double, double, double, double, MatrixType&, VectorType&, ProcessInfo&);

	void DD(double, double, double, double, double, double, double, double, MatrixType&, VectorType& , ProcessInfo&);
	
	
	void A2(double ,double, double, double, double,double ,double,double, MatrixType&, VectorType& , ProcessInfo&);

	void A1(double,double, double, double, double,double,double, MatrixType& , VectorType& , ProcessInfo&);

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
		static boost::numeric::ublas::bounded_matrix<double,3,3> msMassFactors;
		static boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
  		static array_1d<double,3> msN; //dimension = number of nodes
		//static Matrix msDN_DX;
		//static Matrix msMassFactors;
		static array_1d<double,2> ms_vel_gauss; //dimesion coincides with space dimension
  		static array_1d<double,3> ms_temp_vec_np; //dimension = number of nodes
		static array_1d<double,3> ms_u_DN;

	double f(double T); 

	double k(double a);
        
	double g(double T); 
	
	///@}
        ///@name Serialization
        ///@{
	friend class Serializer;

        // A private default constructor necessary for serialization
        ConvDiffChangeOfPhase2D() : Element()
	{
	}
	
        virtual void save(Serializer& rSerializer) const
	{
	KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
	}
	
	virtual void load(Serializer& rSerializer)
	{
	KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
	}


      ///@} 
      ///@name Member Variables 
      ///@{ 
		
        
        
      ///@} 
      ///@name Private Operators
      ///@{ 
          
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
      //ConvDiffChangeOfPhase2D& operator=(const ConvDiffChangeOfPhase2D& rOther);

      /// Copy constructor.
      //ConvDiffChangeOfPhase2D(const ConvDiffChangeOfPhase2D& rOther);

        
      ///@}    
        
    }; // Class ConvDiffChangeOfPhase2D 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream, 
				    ConvDiffChangeOfPhase2D& rThis);
*/
  /// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream, 
				    const ConvDiffChangeOfPhase2D& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  ///@} 

}  // namespace Kratos.

#endif // KRATOS_TRIANGULAR_CONVDIFF_2_ELEM_H_INCLUDED  defined 


