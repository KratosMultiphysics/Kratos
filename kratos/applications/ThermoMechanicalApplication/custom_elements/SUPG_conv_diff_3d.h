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
//   Last Modified by:    $Author: kazem $
//   Date:                $Date: 2009-01-21 14:14:49 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_SUPG_CONV_DIFF_3D_INCLUDED )
#define  KRATOS_SUPG_CONV_DIFF_3D_INCLUDED


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
  
  /// ASGS, Incompressible fluid, Variational multi scale method, Quasi-static subscales, Implicit method.
  /** 
  ASGS is an abriviation for Algebraic Sub-Grid Scale element. It is implemented to solve
  Implicitly the NS equations in a variotionally consistant sub-grid scale methid. It also has the OSS swith
  to use Orthogonal Sub Scales to use impose explicity the orthogonality condition on subscales´ estimation.
  The "Dynamic_Tau" swith allows the use of "Dt", time step, in calculation of Tau.
  This element just work with Monolithic schemes like "monolithic_solver_eulerian" or "monolithic_solver_lagranigan".
  The detailed description of the formulation could be fined in
     "Stabilized finite element approximation of transient incompressible flows using orthogonal subscales, Comput. Methods Appl. Mech. Engrg. 191 (2002) 4295?4321"
  
  */
  class SUPGConvDiff3D
	  : public Element
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of Fluid2DASGS
      KRATOS_CLASS_POINTER_DEFINITION(SUPGConvDiff3D);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.

	  SUPGConvDiff3D(IndexType NewId, GeometryType::Pointer pGeometry);
      SUPGConvDiff3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~SUPGConvDiff3D();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

      	        /// The DOF´s are VELOCITY_X, VELOCITY_Y and PRESSURE
	        /**
	         * @param ElementalDofList: the list of DOFs
	         * @param rCurrentProcessInfo: the current process info instance
		 */
	         void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
	  


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
			return "SUPGConvDiff3D #" ;
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
 		 
      ///@}
      ///@name Protected Operators
      ///@{ 
       SUPGConvDiff3D() : Element(){}
        
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
      
   // private:
       	         
 	         /// To Calculate tau
	        /**
	         * @return tau: multiplied by Residual of the momentum equation		    
		 */		 
	         virtual void CalculateTau(array_1d<double, 3 >& ms_adv_vel, double& tau,const double& K, const double time, const double area, const ProcessInfo& rCurrentProcessInfo);
                              
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
//         ASGS2D() : Element()
//         {
//         }
        
        virtual void save(Serializer& rSerializer) const
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
        }

        virtual void load(Serializer& rSerializer)
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
        }
        
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
      //Fluid2DASGS& operator=(const Fluid2DASGS& rOther);

      /// Copy constructor.
      //Fluid2DASGS(const Fluid2DASGS& rOther);

        
      ///@}    
        
    }; // Class Fluid2DASGS 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream, 
				    Fluid2DASGS& rThis);
*/
  /// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Fluid2DASGS& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  ///@} 

}  // namespace Kratos.

#endif // KRATOS_ASGS_2D_H_INCLUDED  defined 


