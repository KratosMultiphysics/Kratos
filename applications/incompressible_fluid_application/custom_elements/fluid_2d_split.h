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


#if !defined(KRATOS_FLUID_2D_SPLIT_H_INCLUDED )
#define  KRATOS_FLUID_2D_SPLIT_H_INCLUDED


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
  
  /// This class allows the calculation of a Newtonian fluid in presence of a porous medium.
  /** @author  Antonia Larese De Tetto <antoldt@cimne.upc.edu>
  * 
  * This class implements a 2D linear triangular element. A non-linear Darcy law is taken into account.
  * The fluid can be calculated both outside and inside a porous medium at once.
  * The solutions system is solved via a monolithic approach. 
  * The dofs per node: velocity and pressure. 
  * Time integration scheme: Bossak residualbased predictor-corrector @see residualbased_predictorcorrector_velocity_bossak_scheme
  * Strategy: @see residualbased_newton_raphson_strategy
  * Stabilization technique: ASGS (optionally OSS can be used activating the OSS_SWITCH parameter)
  * Python solvers that can use this elelement: @see monolithic_solver_eulerian
  * 
  */
  class Fluid2DSplit
	  : public Element
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of Fluid2DSplit
      KRATOS_CLASS_POINTER_DEFINITION(Fluid2DSplit);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Fluid2DSplit(IndexType NewId, GeometryType::Pointer pGeometry);
      Fluid2DSplit(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~Fluid2DSplit();
      

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

//        void Calculate( const Variable<array_1d<double,3> >& rVariable, 
//                        array_1d<double,3>& Output, 
//                        const ProcessInfo& rCurrentProcessInfo);

       void GetFirstDerivativesVector(Vector& values, int Step = 0);
       void GetSecondDerivativesVector(Vector& values, int Step = 0);

 //      void DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);
      ///Calculate all the lhs contributions  multiplied by velocity: i.e. the convective, pressure, viscous, darcy contributions
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
			return "Fluid2DSplit #" ;
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
      ///Evaluates the elemental density, and (fluid) viscosity, porosity and average diameter of teh porous material;
       virtual void CalculateDensity(Geometry< Node<3> > geom, double& _density, double& viscosity, double& porosity, double& diameter);
      ///Evaluates the residual of the solution system including the viscous contribution \f$ rhs = -lhs  u  \tau  \f$
       virtual void CalculateResidual(const MatrixType& K, VectorType& F);
//        virtual void ComputeProjections(array_1d<double,6>& adv_proj , array_1d<double,3>& div_proj, const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const double thawone,const double thawtwo,const array_1d<double,3>& N,const double area, const double time); 
      ///Evaluates the following stabilization terms:  \f$  (a \cdot \nabla w, \partial_{t} u) \f$ 
       virtual void CalculateAdvMassStblTerms(MatrixType& M,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N, const double thawone,const double area);
      ///Evaluates the following stabilization terms:  \f$  (\nabla q, \partial_{t} u) \f$ 
       virtual void CalculateGradMassStblTerms(MatrixType& M,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N, const double thawone,const double area);
//        virtual void CalculateDarcyMassStblTerms(MatrixType& M,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const double thawone,const double area);
      ///Calculate the mass contribution to the K global matrix	
        virtual void CalculateMassContribution(MatrixType& K,const double time,const double area); 
      ///Calculate the viscous contribution to the lhs
	virtual void CalculateViscousTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double area);
      ///Calculate the advective contribution to the lhs
	virtual void CalculateAdvectiveTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>&  N, const double thawone, const double thawtwo, const double time,const double area);
	//Evaluate the Darcy non linear term for the element cut by the free surface. Non linearity is treated using a Picard method.
	virtual void CalculateDarcyTerm_SubElem(MatrixType& K, const array_1d<double,3>&  N, const double area);
	//Evaluate the Darcy non linear term. Non linearity is treated using a Picard method.
	virtual void CalculateDarcyTerm(MatrixType& K, const double area);
	///Calculate the pressure contribution to the lhs and divergence term as well
	virtual void CalculatePressureTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N,const double time ,const double area);
	///Calculate the following stabilization terms:  \f$  (a \nabla w, \nabla \cdot u) \f$ 
	virtual void CalculateDivStblTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double thawtwo,const double area);
	///Calculate the following stabilization terms:  \f$  (a \nabla w, a \nabla w \nabla p + f) \f$  and \f$  (\nabla q, a \nabla u ) \f$
	virtual void CalculateAdvStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const array_1d<double,3>& N, const double thawone,const double time,const double area);
	///Calculate the following stabilization terms:  \f$  (\nabla q, \nabla p + f) \f$ 
	virtual void CalculateGradStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX,const array_1d<double,3>& N, const double time,const double thawone,const double area);
// 	virtual void CalculateDarcyStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX, const double time,const double thawone,const double area);
	///Add body forces to the lhs
        virtual void AddBodyForceAndMomentum(VectorType& F,const array_1d<double,3>& N, const double time,const double area,const double thawone,const double thawtwo);
        virtual void AddVolumeCorrection(VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N, const double time,const double area);
	///Calculate stabilization parameter
	virtual void CalculateTau(const array_1d<double,3>& N, double& thawone, double& thawtwo, const double time,const double area,const ProcessInfo& rCurrentProcessInfo);
	
// 	virtual void AddProjectionForces(VectorType& F, const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX, const double area,const double thawone,const double thawtwo);
 	///Calcualte the mass contributions
       virtual void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
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
      Fluid2DSplit() : Element()
      {
      }

      virtual void save(Serializer& rSerializer)
      {
	  rSerializer.save("Name", "Fluid2DSplit");
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
      //Fluid2DSplit& operator=(const Fluid2DSplit& rOther);

      /// Copy constructor.
      //Fluid2DSplit(const Fluid2DSplit& rOther);

        
      ///@}    
        
    }; // Class Fluid2DSplit 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream, 
				    Fluid2DSplit& rThis);
*/
  /// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Fluid2DSplit& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  ///@} 

  ///@} IncompressibleFluidApplication group

}  // namespace Kratos.

#endif // KRATOS_FLUID_2D_SPLIT_H_INCLUDED  defined 


