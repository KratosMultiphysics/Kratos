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

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      //virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

	  void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
	void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);
//	  void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

       void Calculate( const Variable<array_1d<double,3> >& rVariable, 
                       array_1d<double,3>& Output, 
                       const ProcessInfo& rCurrentProcessInfo);

       void GetFirstDerivativesVector(Vector& values, int Step = 0);
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
       virtual void calculatedensity(Geometry< Node<3> > geom, double& density, double& viscosity);
       virtual void CalculateResidual(const MatrixType& K, VectorType& F, const double area);
       virtual void ComputeProjections(array_1d<double,6>& adv_proj , array_1d<double,3>& div_proj, const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const double thawone,const double thawtwo,const array_1d<double,3>& N,const double area, const double time); 
       virtual void CalculateAdvMassStblTerms(MatrixType& M,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N, const double thawone,const double area);
       virtual void CalculateGradMassStblTerms(MatrixType& M,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const double thawone,const double area);
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
      
      ///@name Static Member Variables 
      ///@{ 
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        virtual void CalculateMassContribution(MatrixType& K,const double time,const double area); 
	virtual void CalculateViscousTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,  const double area);
	virtual void CalculateAdvectiveTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double thawone, const double thawtwo, const double time,const double area);
	virtual void CalculatePressureTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N,const double time ,const double area);

	virtual void CalculateDivStblTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double thawtwo,const double area);
	virtual void CalculateAdvStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const array_1d<double,3>& N, const double thawone,const double time,const double area);
	virtual void CalculateGradStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX, const double time,const double thawone,const double area);
        virtual void AddBodyForceAndMomentum(VectorType& F,const array_1d<double,3>& N, const double time,const double area,const double thawone,const double thawtwo);

	virtual void CalculateTau(double& thawone, double& thawtwo, const double time,const double area,const ProcessInfo& rCurrentProcessInfo);
	
	virtual void AddProjectionForces(VectorType& F, const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX, const double area,const double thawone,const double thawtwo);
        virtual void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
//         virtual void CalculateInternalForces(VectorType& rRightHandSideVector, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double,3>& N, const double area); 
	virtual void CalculateB(	 boost::numeric::ublas::bounded_matrix<double, 3, 6 > & B,const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX);

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

#endif // KRATOS_NONEWTONIAN_ASGS_2D_H_INCLUDED  defined 


