/*
==============================================================================
KratosR1PoissonApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2008
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
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_POISSON_2D_ELEM_H_INCLUDED)
#define  KRATOS_POISSON_2D_ELEM_H_INCLUDED

// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

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

  class Poisson2D
 	  : public Element
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of Poisson2D
      KRATOS_CLASS_POINTER_DEFINITION(Poisson2D);
 
 
      ///@}
      ///@name Life Cycle 
      ///@{ 

     /// Default constructor.
      Poisson2D(IndexType NewId, GeometryType::Pointer pGeometry);
      Poisson2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
 
      /// Destructor.
      virtual ~ Poisson2D();
 
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
// 		static boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
// 		static boost::numeric::ublas::bounded_matrix<double,2,2> Poisson2D::msD;
//   		static array_1d<double,3> msN; //dimension = number of nodes
// 		static array_1d<double,3> ms_temp; //dimension = number of nodes
// 		static array_1d<double,3> Poisson2D::point_sources; //dimension = number of nodes
// 		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
// 		boost::numeric::ublas::bounded_matrix<double,2,2> Poisson2D::msD;
//   		array_1d<double,3> msN; //dimension = number of nodes
// 		array_1d<double,3> ms_temp; //dimension = number of nodes
// 		array_1d<double,3> Poisson2D::point_sources; //dimension = number of nodes
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
      // Poisson2D& operator=(const Poisson2D& rOther);
 
      /// Copy constructor.
      // Poisson2D(const Poisson2D& rOther);
        
      ///@}    
        
    }; // Class Poisson2D

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function

  /// output stream function
 
  ///@} 

}  // namespace Kratos.

#endif // KRATOS_POISSON_2D_ELEM_H_INCLUDED  defined
