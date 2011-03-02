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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:32 $
//   Revision:            $Revision: 1.2 $
//
//
 

#if !defined(KRATOS_TRIANGULAR_CONVDIFF_ELEM_H_INCLUDED )
#define  KRATOS_TRIANGULAR_CONVDIFF_ELEM_H_INCLUDED



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
  
  /// A stabilized element for solving the convection diffusion equations in 2D.
  
   /** 
	The problem can be written in matrix form as 
     \f$ \rho C M \frac{\partial T}{\partial t} + \rho C S T  = - \kappa L T \f$


  */


  class ConvDiff2D
	  : public Element
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of ConvDiff2D
     // typedef GeometryData::IntegrationMethod IntegrationMethod;
      KRATOS_CLASS_POINTER_DEFINITION(ConvDiff2D);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
	  ConvDiff2D(IndexType NewId, GeometryType::Pointer pGeometry);
      ConvDiff2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~ConvDiff2D();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
	//IntegrationMethod GetIntegrationMethod();
      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;
      
      ///Evaluates  \f$ L h s = \frac{\rho C}{\Delta t} (W, N) + (W, v. \nabla N) + \kappa (\nabla W, \nabla N)  \f$ and \f$R h s = \rho (W, Q) + \frac{\rho C}{\Delta t} (W, T^n) - L h s \ast T \f$



      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

      ///Provides the global indices for each one of this element's local rows. @see NoNewtonianASGS2D 
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

      ///Returns a list of the element's Dofs. @see NoNewtonianASGS2D 
	  void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
	/// Calculates the temperature convective projection 
	  void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

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
	//IntegrationMethod mThisIntegrationMethod;
        //std::vector< Matrix > mInvJ0;
        //Vector mDetJ0;
      ///@name Static Member Variables 
      ///@{ 
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
		
	///@}
        ///@name Serialization
        ///@{
	friend class Serializer;

        // A private default constructor necessary for serialization
        ConvDiff2D() : Element()
	{
	}
	
        virtual void save(Serializer& rSerializer)
	{
	rSerializer.save("Name", "ConvDiff2D");
	KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
	}
	
	virtual void load(Serializer& rSerializer)
	{
	KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
	}
	


        
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
      //ConvDiff2D& operator=(const ConvDiff2D& rOther);

      /// Copy constructor.
      //ConvDiff2D(const ConvDiff2D& rOther);

        
      ///@}    
        
    }; // Class ConvDiff2D 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream, 
				    ConvDiff2D& rThis);
*/
  /// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream, 
				    const ConvDiff2D& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  ///@} 

}  // namespace Kratos.

#endif // KRATOS_TRIANGULAR_CONVDIFF_ELEM_H_INCLUDED  defined 


