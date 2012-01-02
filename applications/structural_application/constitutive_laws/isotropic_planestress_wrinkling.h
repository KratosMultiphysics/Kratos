/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
//   Date:                $Date: 2007-11-22 11:49:23 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_ISOTROPIC_PLANESTRESS_WRINKLING_H_INCLUDED )
#define  KRATOS_ISOTROPIC_PLANESTRESS_WRINKLING_H_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


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
	class IsotropicPlaneStressWrinkling : public ConstitutiveLaw
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of IsotropicPlaneStressWrinkling
      //typedef boost::shared_ptr<IsotropicPlaneStressWrinkling> Pointer;
      KRATOS_CLASS_POINTER_DEFINITION( IsotropicPlaneStressWrinkling );
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      IsotropicPlaneStressWrinkling();
      IsotropicPlaneStressWrinkling(const double E, const double NU);

      /// Destructor.
      virtual ~IsotropicPlaneStressWrinkling();
      

      ///@}
      ///@name Operators 
      ///@{
	boost::shared_ptr<ConstitutiveLaw> Clone() const
	{
		boost::shared_ptr<ConstitutiveLaw> p_clone(new IsotropicPlaneStressWrinkling());
		return p_clone;
	}
      
	  void InitializeMaterial( 	const Properties& props,
					const GeometryType& geom,
					const Vector& ShapeFunctionsValues );

	  void CalculateConstitutiveMatrix( const Vector& StrainVector, 
					    Matrix& ElasticityTensor );

	  void CalculateStress( const Vector& StrainVector, 
				Vector& StressVector);
     		
	  void UpdateMaterial( const Vector& StrainVector,
				const Properties& props,
				const GeometryType& geom,
				const Vector& ShapeFunctionsValues ,
				const ProcessInfo& CurrentProcessInfo);
	  
	  void Calculate( const Variable<Matrix >& rVariable, 
			  Matrix& Output, 
    			  const ProcessInfo& rCurrentProcessInfo);
		
		
	  void CalculateCauchyStresses(
			Vector& Cauchy_StressVector,
			const Matrix& F,
			const Vector& PK2_StressVector,
			const Vector& GreenLagrangeStrainVector);
	  
      void CalculateMaterialResponse( const Vector& StrainVector,
                                      const Matrix& DeformationGradient,
                                      Vector& StressVector,
                                      Matrix& AlgorithmicTangent,
                                      const ProcessInfo& CurrentProcessInfo,
                                      const Properties& props, 
                                      const GeometryType& geom,
                                      const Vector& ShapeFunctionsValues,
                                      bool CalculateStresses = true,
                                      int CalculateTangent = true,
                                      bool SaveInternalVariables = true );

	  double CalculateThicknessRatio(const Vector& GreenLagrangeStrainVector);
	  
	  int Check(const Properties& props, const GeometryType& geom, const ProcessInfo& CurrentProcessInfo);

      ///@}
      ///@name Operations
      ///@{
      
      
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
      //virtual String Info() const;
      
      /// Print information about this object.
      //virtual void PrintInfo(std::ostream& rOStream) const;

      /// Print object's data.
      //virtual void PrintData(std::ostream& rOStream) const;
      
            
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
		double mE;
		double mNU;
		Matrix mCtangent;
		Vector mEw;
       
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
		void CalculateStressEigNonNormalized(const Vector& StressVector,double& smin, Vector& mineigenvect,double& smax, Vector& maxeigenvect);
		void CalculateStressEig(const Vector& StressVector,double& smin, Vector& mineigenvect,double& smax, Vector& maxeigenvect);
		void PrincipSTRAIN(const Vector& StrainVector,double& eps1, double& eps2);
		void PrincipSTRESS(const Vector& StressVector,double& str1, double& str2);

        unsigned int AssessState(const Vector& ElasticStress,const Vector& StrainVector);
		void CalculateTangentMatrix(const Vector& StrainVector);

		Matrix ConstructUnidirectionalConstitutiveMatrix(const Vector& v);
	
		Matrix CalculateElasticMatrix(const double E, const double NU);
	
		Matrix CalculateDirectionVariationTerm(const Matrix& Celastic, const Vector& ElasticStrain);

        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
        
      ///@} 
      ///@name Private  Access 
      ///@{ 
        
        
      private:

      ///@}
      ///@name Serialization
      ///@{	
      friend class Serializer;

      virtual void save(Serializer& rSerializer) const
      {
         KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
      }

      virtual void load(Serializer& rSerializer)
      {
         KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
      }
		
	
      ///@}    
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      //IsotropicPlaneStressWrinkling& operator=(const IsotropicPlaneStressWrinkling& rOther);

      /// Copy constructor.
      //IsotropicPlaneStressWrinkling(const IsotropicPlaneStressWrinkling& rOther);

        
      ///@}    
        
    }; // Class IsotropicPlaneStressWrinkling 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  ///// input stream function
  //inline std::istream& operator >> (std::istream& rIStream, 
		//		    IsotropicPlaneStressWrinkling& rThis);

  ///// output stream function
  //inline std::ostream& operator << (std::ostream& rOStream, 
		//		    const IsotropicPlaneStressWrinkling& rThis)
  //  {
  //    rThis.PrintInfo(rOStream);
  //    rOStream << std::endl;
  //    rThis.PrintData(rOStream);

  //    return rOStream;
  //  }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_ISOTROPIC_PLANESTRESS_WRINKLING_H_INCLUDED  defined 


