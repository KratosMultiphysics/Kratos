//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_DISCRETE_ELEMENT_H_INCLUDED )
#define  KRATOS_DISCRETE_ELEMENT_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cmath>

// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "includes/process_info.h"
#include "utilities/indexed_object.h"
#include "containers/weak_pointer_vector.h"
#include "includes/constitutive_law.h"

//Cfeng,RigidFace
#include "includes/condition.h"

namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
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
  
  /// Short class definition.
  /** Detail class definition.
  */
  class DiscreteElement 
       : public Element
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of DiscreteElement
      KRATOS_CLASS_POINTER_DEFINITION(DiscreteElement);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
        /// Default constructor.
        //DiscreteElement() : Element() {}
     
        DiscreteElement(IndexType NewId = 0) : Element(NewId){ } 
     
        /**
         * Constructor using an array of nodes
         */
        DiscreteElement(IndexType NewId, const NodesArrayType& ThisNodes) :
        Element(NewId, ThisNodes) 
        {
        }

        /**
         * Constructor using Geometry
         */
        DiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry) :
        Element(NewId, pGeometry)
        {
        }

        /**
         * Constructor using Properties
         */
        DiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : 
        Element(NewId, pGeometry , pProperties)
        {
        }

        /// Copy constructor.
        DiscreteElement (DiscreteElement const& rOther) : Element(rOther)
        {
        }
      
   
      /// Destructor.
      virtual ~DiscreteElement(){}
      
      /*
      virtual Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
      {
        KRATOS_TRY
           return Element::Pointer(new Element(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
      }
      */
      
      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      /**
         * is called to initialize the element.
         * Must be called before any calculation is done!
         */
        virtual void Initialize()
        {
        }

        /**
         * this is called during the assembling process in order
         * to calculate all elemental contributions to the global system
         * matrix and the right hand side
         * @param rLeftHandSideMatrix: the elemental left hand side matrix
         * @param rRightHandSideVector: the elemental right hand side
         * @param r_process_info: the current process info instance
         */

        virtual void ResetConstitutiveLaw()
        {
        }

       
        /**
         * this is called during the assembling process in order
         * to calculate the elemental right hand side vector only
         * @param rRightHandSideVector: the elemental right hand side vector
         * @param r_process_info: the current process info instance
         */
        virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
                ProcessInfo& r_process_info)
        {
            if (rRightHandSideVector.size() != 0)
                rRightHandSideVector.resize(0);
        }

       
        virtual void EquationIdVector(EquationIdVectorType& rResult,
                ProcessInfo& r_process_info)
        {
            if (rResult.size() != 0)
                rResult.resize(0);
        }

        /**
         * this is called during the assembling process in order
         * to calculate the elemental mass matrix
         * @param rMassMatrix: the elemental mass matrix
         * @param r_process_info: the current process info instance
         */
        virtual void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info)
        {
            if (rMassMatrix.size1() != 0)
                rMassMatrix.resize(0, 0);
        }

        
        virtual void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info)
        {
            if (rDampingMatrix.size1() != 0)
                rDampingMatrix.resize(0, 0);
        }


        virtual void GetDofList(DofsVectorType& ElementalDofList,
                ProcessInfo& r_process_info)
        {
            if (ElementalDofList.size() != 0)
                ElementalDofList.resize(0);
        }
      
        /**
         * this is called in the beginning of each solution step
         */
        virtual void InitializeSolutionStep(const ProcessInfo& r_process_info)
        {
        }

        /**
         * this is called at the end of each solution step
         */
        virtual void FinalizeSolutionStep(const ProcessInfo& r_process_info)
        {
        }

        /**
         * deletes all obsolete data from memory
         */
        virtual void CleanMemory()
        {
        }

        virtual void GetValuesVector(Vector& values, int Step = 0)
        {
        }

        virtual void GetFirstDerivativesVector(Vector& values, int Step = 0)
        {
        }

        virtual void GetSecondDerivativesVector(Vector& values, int Step = 0)
        {
        }

        //output On integration points
        //calculate on element

        virtual void Calculate(const Variable<double>& rVariable,
                double& Output,
                const ProcessInfo& r_process_info)
        {
        }

        virtual void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                array_1d<double, 3 > & Output,
                const ProcessInfo& r_process_info)
        {
        }

        virtual void Calculate(const Variable<Vector >& rVariable,
                Vector& Output,
                const ProcessInfo& r_process_info)
        {
        }

        virtual void Calculate(const Variable<Matrix >& rVariable,
                Matrix& Output,
                const ProcessInfo& r_process_info)
        {
        }








      ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      virtual std::string Info() const
      {
	  std::stringstream buffer;
	  buffer << "Discrete Element #" << Id();
	  return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
	  rOStream << "Discrete Element #" << Id();
      }
      
      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
	  //mpGeometry->PrintData(rOStream);
      }
      
            
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
      ///@name Serialization
      ///@{

      friend class Serializer;

      virtual void save(Serializer& rSerializer) const
      {
          KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
      }

      virtual void load(Serializer& rSerializer)
      {
          KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
      }

        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{    
      ///@}    
        
    }; // Class DiscreteElement 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    DiscreteElement& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const DiscreteElement& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DISCRETE_ELEMENT_H_INCLUDED  defined 
