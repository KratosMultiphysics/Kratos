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


#if !defined(KRATOS_EXPLICIT_ASGS_COMP_PR_DC_2D_H_INCLUDED )
#define  KRATOS_EXPLICIT_ASGS_COMP_PR_DC_2D_H_INCLUDED


// System includes  


// External includes 
#include "boost/smart_ptr.hpp"
 

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "custom_elements/explicit_asgs_compressible_2d.h"
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
  class ExplicitASGSCOMPPRDC2D
	  : public ExplicitASGSCompressible2D
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of ASGSCOMPPRDC2D
      KRATOS_CLASS_POINTER_DEFINITION(ExplicitASGSCOMPPRDC2D);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      ExplicitASGSCOMPPRDC2D(IndexType NewId, GeometryType::Pointer pGeometry);
      ExplicitASGSCOMPPRDC2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~ExplicitASGSCOMPPRDC2D();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;
       void GetSecondDerivativesVector(Vector& values, int Step = 0);
       void GetFirstDerivativesVector(Vector& values, int Step = 0);
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
      void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
//       void Calculate( const Variable<array_1d<double,3> >& rVariable, 
// 		      array_1d<double,3>& Output, 
// 		      const ProcessInfo& rCurrentProcessInfo);

	void Calculate( const Variable<double>& rVariable, 
			      double& Output, 
			      const ProcessInfo& rCurrentProcessInfo);
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
			return "ExplicitASGSCOMPPRDC2D #" ;
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

	virtual void CalculateSoundVelocity(Geometry< Node<3> > geom, double& vc);
       virtual void calculatedensity(Geometry< Node<3> > geom, double& density, double& viscosity);
       virtual void CalculateResidual(const MatrixType& K, VectorType& F);
	virtual void CalculateDivPdotStblTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX,const array_1d<double,3>& N, const double time,const double thawone,const double volume);
	virtual void CalculateCharectristicLength(double& ch_length, const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,double& norm_grad );
 	virtual void CalculateArtifitialViscosity(double& art_visc,double& Pr_art_visc ,const boost::numeric::ublas::bounded_matrix<double,3,2>&DN_DX);

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
	/*	
        double m_thawone;
        double m_thawtwo ;
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        void CalculateMassContribution(MatrixType& K,const double time,const double area); 
	void CalculateViscousTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double area);
	void CalculateAdvectiveTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, double& thawone, double& thawtwo, const double time,const double area);
	void CalculatePressureTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const array_1d<double,3>& N,const double time ,const double area);

	void CalculateDivStblTerm(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double thawtwo,const double area);
	void CalculateAdvStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const array_1d<double,3>& N, const double thawone,const double time,const double area);
	void CalculateGradStblAllTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX, const double time,const double thawone,const double area);
	void AddBodyForceAndMomentum(VectorType& F, const array_1d<double,3>& N, const double time,const double area);


	
	void AddProjectionForces(VectorType& F, const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX, const double area);
	*/
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
        ExplicitASGSCOMPPRDC2D() : ExplicitASGSCompressible2D()
        {
        }
        
        virtual void save(Serializer& rSerializer)
        {
            rSerializer.save("Name", "ExplicitASGSCOMPPRDC2D");
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ExplicitASGSCompressible2D);
        }

        virtual void load(Serializer& rSerializer)
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ExplicitASGSCompressible2D);
        }
        
      ///@}       
      
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      //ASGSCOMPPRDC2D& operator=(const ASGSCOMPPRDC2D& rOther);

      /// Copy constructor.
      //ASGSCOMPPRDC2D(const ASGSCOMPPRDC2D& rOther);

        
      ///@}    
        
    }; // Class ASGSCOMPPRDC2D 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream, 
				    ASGSCOMPPRDC2D& rThis);
*/
  /// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream, 
				    const ASGSCOMPPRDC2D& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
  ///@} 

}  // namespace Kratos.

#endif // KRATOS_ASGS_PR_DC_H_INCLUDED  defined 


