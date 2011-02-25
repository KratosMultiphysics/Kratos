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
//   Date:                $Date: 2009-01-21 14:14:49 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_NONEWTONIAN_ASGS_2D_H_INCLUDED )
#define  KRATOS_NONEWTONIAN_ASGS_2D_H_INCLUDED


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
  
  /// This class allows the calculation of a Non-Newtonian fluid using a visco rigid constitutive model.
  /** @author  Antonia Larese De Tetto <antoldt@cimne.upc.edu>
  * 
  * This class implements a 2D linear triangular element. A non-newtonian constituve law is developed using a visco-rigid model.
  * The solutions system is solved via a monolithic approach. 
  * The dofs per node: velocity and pressure. 
  * Time integration scheme: Bossak residualbased predictor-corrector @see residualbased_predictorcorrector_velocity_bossak_scheme
  * Strategy: @see residualbased_newton_raphson_strategy
  * Stabilization technique: ASGS (optionally OSS can be used activating the OSS_SWITCH parameter)
  * Python solvers that can use this elelement: @see monolithic_solver_lagrangian_nonnewtonian,  monolithic_solver_eulerian
  * 
  */
  class NoNewtonianASGS2D
	  : public Element
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of NoNewtonianASGS2D
      KRATOS_CLASS_POINTER_DEFINITION(NoNewtonianASGS2D);

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      NoNewtonianASGS2D(IndexType NewId, GeometryType::Pointer pGeometry);
      NoNewtonianASGS2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~NoNewtonianASGS2D();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      /// Create a new element of this type
      /**
	* Returns a pointer to a new NoNewtonianASGS2D element, created using given input
	* @param NewId: the ID of the new element
	* @param ThisNodes: the nodes of the new element
	* @param pProperties: the properties assigned to the new element
	* @return a Pointer to the new element
	*/
      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;
      
      ///Calculate the local external force vector and resize local sistem matrices
      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      ///Calulate the residual (RHS) of the solution system
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      
      /// Provides the global indices for each one of this element's local rows
	/**
	* this determines the elemental equation ID vector for all elemental
	* DOFs
	* @param rResult: A vector containing the global Id of each row
	* @param rCurrentProcessInfo: the current process info object
	*/
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
      
	/// Returns a list of the element's Dofs
	/**
	* @param ElementalDofList: the list of DOFs
	* @param rCurrentProcessInfo: the current process info instance
	*/
	  void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
	  
	/// Returns the values on the integration points
	/**
	* @param rVariable: Kratos vector variable to compute (implemented for the variable viscosity variable and for the rate of strain variable)
	* @param Output: Values of variable on integration points
	* @param rCurrentProcessInfo: Process info instance
	*/
	void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);
//	  void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

      void Calculate( const Variable<array_1d<double,3> >& rVariable, 
		      array_1d<double,3>& Output, 
		      const ProcessInfo& rCurrentProcessInfo);

      void GetFirstDerivativesVector(Vector& values, int Step = 0);
      void GetSecondDerivativesVector(Vector& values, int Step = 0);

//      void DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);
      ///Calculate all the lhs contribution multiplied by velocity: i.e. the convective, pressure, viscous, darcy contributions
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
			return "NoNewtonianASGS2D #" ;
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
      ///Evaluates the elemental density, and (fluid) viscosity
      virtual void calculatedensity(Geometry< Node<3> > geom, double& density, double& viscosity);
      ///Evaluates the residual of the solution system including the viscous contribution \f$ rhs = -lhs  u - B^{T} \tau  \f$
      virtual void CalculateResidual(const MatrixType& K, VectorType& F, const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double area);
      ///Compute the projection in case OSS stabilization thechnique is chosen (OSS_SWITCH should be set = 1.0);
      virtual void ComputeProjections(array_1d<double,6>& adv_proj , array_1d<double,3>& div_proj, const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const double thawone,const double thawtwo,const array_1d<double,3>& N,const double area, const double time); 
      ///Evaluates the following stabilization terms:  \f$  (a \cdot \nabla w, \partial_{t} u) \f$ 
      virtual void CalculateAdvMassStblTerms(MatrixType& M,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N, const double thawone,const double area);
      ///Evaluates the following stabilization terms:  \f$  (\nabla q, \partial_{t} u) \f$ 
      virtual void CalculateGradMassStblTerms(MatrixType& M,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N, const double thawone,const double area);
      ///Calculate the mass contribution to the K global matrix	
      virtual void CalculateMassContribution(MatrixType& K,const double time,const double area); 
      ///Calculate the linearized viscous contribution ONLY to the LHS @todo Make linearization works quadratically
      virtual void CalculateViscousTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,  const double area, const int it_num);
      ///Calculate the advective contribution to the lhs
      virtual void CalculateAdvectiveTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double, 3 >& N, const double thawone, const double thawtwo, const double time,const double area);
	///Calculate the pressure contribution to the lhs and divergence term as well
      virtual void CalculatePressureTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N,const double time ,const double area);

	///Calculate the following stabilization terms:  \f$  (a \nabla w, \nabla \cdot u) \f$ 
      virtual void CalculateDivStblTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double thawtwo,const double area);
	///Calculate the following stabilization terms:  \f$  (a \nabla w, a \nabla w \nabla p + f) \f$  and \f$  (\nabla q, a \nabla u ) \f$
      virtual void CalculateAdvStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const array_1d<double,3>& N, const double thawone,const double time,const double area);
	///Calculate the following stabilization terms:  \f$  (\nabla q, \nabla p + f) \f$ 
      virtual void CalculateGradStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX,const array_1d<double,3>& N, const double time,const double thawone,const double area);
	///Add body forces to the lhs
      virtual void AddBodyForceAndMomentum(VectorType& F, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double,3>& N, const double time,const double area,const double thawone,const double thawtwo);
	///Calculate stabilization parameter
       virtual void CalculateTau(const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX, const array_1d<double,3>& N, double& thawone, double& thawtwo, const double time,const double area,const ProcessInfo& rCurrentProcessInfo);
        ///Add the projection in case OSS stabilization thechnique is chosen (OSS_SWITCH should be set = 1.0);
	virtual void AddProjectionForces(VectorType& F, const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX, const double area,const double thawone,const double thawtwo);
	///Calcualte the mass contributions
	virtual void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
//         virtual void CalculateInternalForces(VectorType& rRightHandSideVector, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double,3>& N, const double area); 
	///Calculate the shape function derivatives matrix
	virtual void CalculateB(	 boost::numeric::ublas::bounded_matrix<double, 3, 6 > & B,const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX);
	///Calculate the symmetric gradient of velocity
	virtual void CalculateGradSymVel(array_1d<double,3> & grad_sym_vel, double & gamma_dot,const boost::numeric::ublas::bounded_matrix<double, 3, 6 > & B);///@} 
 ///Evaluates the viscosity of the nodes variable in function of the rate of strain
	/**
	* @param ApparentViscosity: \f$ \tilde{\nu} \f$ non-newtonian variable viscosity
	* @param ApparentViscosityDerivative: \f$ \frac{\partial \tilde{\nu}}{\partial \gamma } \f$
	* @param grad_sym_vel: Simmetric gradient of velocity \f$ \nabla^{s} u = \varepsilon \f$
	* @param gamma_dot: rate of strain \f$ \gamma = \sqrt{2 \varepsilon:\varepsilon }    \f$
	* @param B: Matrix 3x6 of the shape function derivatives
	* @param mu: fluid minimum viscosity possible
	*/
      virtual void CalculateApparentViscosity(double & ApparentViscosity, double & ApparentViscosityDerivative, array_1d<double,3> & grad_sym_vel, double & gamma_dot, const boost::numeric::ublas::bounded_matrix<double, 3, 6 > & B, const double & mu);
      virtual void CalculateApparentViscosityStbl(double & ApparentViscosity, double & ApparentViscosityDerivative, array_1d<double,3> & grad_sym_vel, double & gamma_dot, const boost::numeric::ublas::bounded_matrix<double, 3, 6 > & B, const double & mu);

      ///@name Protected Operators
      ///@{ 
	
      ///@}
      ///@name Serialization
      ///@{

      // A private default constructor necessary for serialization
      NoNewtonianASGS2D() : Element()
      {
      }

	
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
      
	  MatrixType mlhs0;
	  MatrixType mKvisc0;	
          double mDevStress;
	  
      ///@}
      ///@name Serialization
      ///@{

      friend class Serializer;

      virtual void save(Serializer& rSerializer)
      {
	  rSerializer.save("Name", "NoNewtonianASGS2D");
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

      /// Copy constructor.

	
      ///@}    
	
    }; // Class NoNewtonianASGS3D

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
  
  ///@} IncompressibleFluidApplication group

}  // namespace Kratos.

#endif // KRATOS_NONEWTONIAN_ASGS_2D_H_INCLUDED  defined 


