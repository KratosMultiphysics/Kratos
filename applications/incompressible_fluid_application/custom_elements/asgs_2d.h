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


#if !defined(KRATOS_ASGS_2D_H_INCLUDED )
#define  KRATOS_ASGS_2D_H_INCLUDED


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
  class ASGS2D
	  : public Element
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of Fluid2DASGS
      KRATOS_CLASS_POINTER_DEFINITION(ASGS2D);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
	  ASGS2D(IndexType NewId, GeometryType::Pointer pGeometry);
      ASGS2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~ASGS2D();
      

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
	  
	         /// To print a scalar value on Gausse points
	        /**
	         * @param rVariable: The "rVariable" must be either "TAUONE" or "TAUTWO"
	         * @param rCurrentProcessInfo: the current process info instance
		 */
	         void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);
//	  void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

	         /// To calculate NODAL_MASS or The OSS projections
	        /**
	         * @param rVariable: If the variable is "NODAL_MASS" it returns the nodal mass vector which is "lumped mass" for first component and the rest zero,
		                     else it calulates the Orthogonal projections
		 * @return Output:  Has value just in the case of NODAL_MASS   
	         * @param rCurrentProcessInfo: the current process info instance
		 */
                 void Calculate( const Variable<array_1d<double,3> >& rVariable, 
                       array_1d<double,3>& Output, 
                       const ProcessInfo& rCurrentProcessInfo);
	         /// To calculate minimum length inside element 
	        /**
	         * @param rVariable: Is not used
		 * @return Output: returns the min value   
	         * @param rCurrentProcessInfo: the current process info instance
		 */		       
                 void Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);

		 /// Returns vx, vy, p for each node
                 void GetFirstDerivativesVector(Vector& values, int Step = 0);

		 /// Returns ax, ay, 0 for each node
                 void GetSecondDerivativesVector(Vector& values, int Step = 0);

 //      void DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);
       void CalculateLocalVelocityContribution(MatrixType& rDampMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo);


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
			return "ASGS2D #" ;
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
       virtual void calculatedensity(Geometry< Node<3> > geom, double& density, double& viscosity);
       virtual void CalculateResidual(const MatrixType& K, VectorType& F);

 	         /// To compute projections for OSS 
	        /**
	         * @return adv_proj: projection due to advection
		 * @return adv_proj: projection due to divergence 
	         * @param rCurrentProcessInfo: the current process info instance
		 */      
                 virtual void ComputeProjections(array_1d<double,6>& adv_proj , array_1d<double,3>& div_proj, const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const double thawone,const double thawtwo,const array_1d<double,3>& N,const double area, const double time); 

 	         /// To Calculate stabilization of the form (a.grad(U) , U_dot) 
	        /**
		    It is assembeled directly to LHS and RHS
		 */   		 
		 virtual void CalculateAdvMassStblTerms(MatrixType& M,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N, const double thawone,const double area);

 	         /// To Calculate stabilization of the form (grad(q) , U_dot) 
	        /**
		    It is assembeled directly to LHS and RHS
		 */  		 
		 virtual void CalculateGradMassStblTerms(MatrixType& M,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N,const double thawone,const double area);

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
      
   // private:
      ///@name Static Member Variables 
      ///@{ 
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
      ///@} 
      ///@name Private Operators
      ///@{ 
       	         /// To Calculate Lumped mass matrix 
	        /**
		    It is assembeled directly to LHS
		 */
                virtual void CalculateMassContribution(MatrixType& K,const double time,const double area);
		
       	         /// To Calculate viscouse term  
	        /**
		    It is assembeled directly to LHS
		 */		
	         virtual void CalculateViscousTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double area);
		 
       	         /// To Calculate advective term  
	        /**
		    It is assembeled directly to LHS
		 */			 
	         virtual void CalculateAdvectiveTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N, const double thawone, const double thawtwo, const double time,const double area);


       	         /// To Calculate Pressure term, (Div(V), p) and (q, Div(U))  
	        /**
		    It is assembeled directly to LHS
		 */			 
		 virtual void CalculatePressureTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N,const double time ,const double area);



 	         /// To Calculate stabilization of the form ( Div(u) , Div(v) ) 
	        /**
		    It is assembeled directly to LHS and is scaled by Tau2
		 */  			 
	         virtual void CalculateDivStblTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double thawtwo,const double area);

 	         /// To Calculate stabilization of the form ( a.grad(V) , a.grad(U) ) 
	        /**
		    It is assembeled directly to LHS and is scaled by Tau1
		 */ 		 
		 virtual void CalculateAdvStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const array_1d<double,3>& N, const double thawone,const double time,const double area);

 	         /// To Calculate stabilization of the form ( grad(q),grad(p) ) 
	        /**
		    It is assembeled directly to LHS and is scaled by Tau1
		 */ 		 
		 virtual void CalculateGradStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX,const array_1d<double,3>& N, const double time,const double thawone,const double area);

 	         /// To add body force 
	        /**
		    It is assembeled directly to RHS and is scaled by Tau1
		 */
		 virtual void AddBodyForceAndMomentum(VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX,const array_1d<double,3>& N, const double time,const double area,const double thawone,const double thawtwo);

 	         /// To Calculate tau1 & tau2 
	        /**
	         * @return tau1: multiplied by Residual of the momentum equation
	         * @return tau2: multiplied by Residual of the constrain equation		    
		 */		 
	         virtual void CalculateTau(const array_1d<double,3>& N,double& thawone, double& thawtwo, const double time,const double area,const ProcessInfo& rCurrentProcessInfo);

 	         /// To add projection forces 
	        /**
                  This function is called by Calculate and assemble explicitly the projection forces to the RHS
		 */  		 
	         virtual void AddProjectionForces(VectorType& F, const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX, const double area,const double thawone,const double thawtwo);

		 virtual void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
	    private:
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


