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
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2010-10-26 14:14:49 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_COUETTE_NONNEWTONIAN_ASGS_2D_H_INCLUDED )
#define  KRATOS_COUETTE_NONNEWTONIAN_ASGS_2D_H_INCLUDED


// System includes  


// External includes 
#include "boost/smart_ptr.hpp"
 

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "custom_elements/nonewtonian_asgs_2d.h"

#include "includes/serializer.h"

namespace Kratos
{
  ///@addtogroup IncompressibleFluidApplication
  ///@{
  
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
  
  /// This class allows the calculation of a Non-Newtonian  fluid using a BINGHAM constitutive model.
  /** @author  Antonia Larese De Tetto <antoldt@cimne.upc.edu>
  * 
  * This class implements a 2D linear triangular element. A non-newtonian constituve law is developed using a BINGHAM model.
  * It is a derived class @see NoNewtonianASGS2D
  * The only difference between the two approaches is in the calculation of the variable viscosity. In the present element a 
  * Bingham plastic is used with the help of an exponencial variation of the viscosity in function of the strain rate
  * Reference Papanastasiou, T. C. Flows of materials with yield. Journal of Rheology, 1987, 31, 385-404

  */
  class CouetteNonNewtonianASGS2D
	  : public NoNewtonianASGS2D
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of CouetteNonNewtonianASGS2D
      KRATOS_CLASS_POINTER_DEFINITION(CouetteNonNewtonianASGS2D);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
	  CouetteNonNewtonianASGS2D(IndexType NewId, GeometryType::Pointer pGeometry);
      CouetteNonNewtonianASGS2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~CouetteNonNewtonianASGS2D();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

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
		virtual std::string Info() const
		{
			return "CouetteNonNewtonianASGS2D #" ;
		}

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << Info() << Id();
	}

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
       virtual void CalculateApparentViscosity(double & ApparentViscosity, double & ApparentViscosityDerivative , array_1d<double,3> & grad_sym_vel, double & gamma_dot, const boost::numeric::ublas::bounded_matrix<double, 3, 6 > & B, const double & mu, const double & m_coef);
//        virtual void CalculateApparentViscosityStbl(double & ApparentViscosity, double & ApparentViscosityDerivative , array_1d<double,3> & grad_sym_vel, double & gamma_dot, const boost::numeric::ublas::bounded_matrix<double, 3, 6 > & B, const double & mu);
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
      ///@name Serialization
      ///@{

      friend class Serializer;

      // A private default constructor necessary for serialization
      CouetteNonNewtonianASGS2D() : NoNewtonianASGS2D()
      {
      }

      virtual void save(Serializer& rSerializer) const
      {
	  rSerializer.save("Name", "CouetteNonNewtonianASGS2D");
	  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, NoNewtonianASGS2D);
      }

      virtual void load(Serializer& rSerializer)
      {
	  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, NoNewtonianASGS2D);
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

      /// Copy constructor.

        
      ///@}    
        
    }; // Class Fluid2DASGS 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function


  /// output stream function

  ///@} 

  ///@} IncompressibleFluidApplication group

}  // namespace Kratos.

#endif // KRATOS_COUETTE_NONNEWTONIAN_ASGS_2D_H_INCLUDED  defined 


