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
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2007-09-20 11:54:06 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_TRUSS_ELEMENT_INCLUDED )
#define  KRATOS_TRUSS_ELEMENT_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
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
    class TrussElement
    : public Element
           {
               public:
      ///@name Type Definitions
      ///@{
				typedef GeometryData::IntegrationMethod IntegrationMethod;

                typedef ConstitutiveLaw ConstitutiveLawType;

               	typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
      /// Counted pointer of TrussElement

                       KRATOS_CLASS_POINTER_DEFINITION(TrussElement);
 
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
                       TrussElement(IndexType NewId, GeometryType::Pointer pGeometry);
               TrussElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

      /// Destructor.
               virtual ~TrussElement();
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
				IntegrationMethod GetIntegrationMethod();

               Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

               void Initialize();

               void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
               void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
      
               void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

               void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);
	  
               void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);
               void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);
 			   void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);
               	//************************************************************************************

                void GetValuesVector(Vector& values, int Step);
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
      ///@name Static Member Variables 
      ///@{ 
				Matrix mA_Matrix;
				Vector mb_Vector;
				double mArea;
				double mYoungs;
				double mReference_length;
				Vector mCurrentDisplacement;
				double mCurrentStrain;
				Vector mAtimesU;
				double mPrescribedStrain;

              void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
                                         ProcessInfo& rCurrentProcessInfo,
                                         bool CalculateStiffnessMatrixFlag,
                                         bool CalculateResidualVectorFlag);
				void GreenStrain();
				double VectorMatrixVector(Vector& vector);
				Vector MatrixVector(Vector& vector);
				Vector MatrixVector();
				void CalcAtimesU();
				void CalculateRHS(Vector& rRightHandSideVector);

				void CalculateLHS(Matrix& rLeftHandSideMatrix);
      ///@} 
      ///@name Member Variables 
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
      //TrussElement& operator=(const TrussElement& rOther);

      /// Copy constructor.
      //TrussElement(const TrussElement& rOther);

        
      ///@}    
        
           }; // Class TrussElement 

  ///@} 
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream, 
                   TrussElement& rThis);
*/
  /// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream, 
                   const TrussElement& rThis)
           {
                   rThis.PrintInfo(rOStream);
                   rOStream << std::endl;
                   rThis.PrintData(rOStream);

                   return rOStream;
}*/
  ///@} 

}  // namespace Kratos.

#endif // KRATOS_UNSATURATED_SOILS_ELEMENT_INCLUDED defined 


