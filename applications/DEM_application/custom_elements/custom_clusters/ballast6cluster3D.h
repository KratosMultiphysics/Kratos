//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Joaquín $
//   Date:                $Date: 2016-01-26 15:00:00 $
//   Revision:            $Revision: 1.1.1.1 $
//
//

#if !defined(KRATOS_BALLAST6CLUSTER3D_H_INCLUDED )
#define  KRATOS_BALLAST6CLUSTER3D_H_INCLUDED

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
#include "custom_utilities/create_and_destroy.h"

#include "includes/condition.h"
#include "custom_elements/cluster3D.h"

namespace Kratos
{

  /// Short class definition.
  /** Detail class definition.
  */
    class Ballast6Cluster3D : public Cluster3D
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Cluster3D
        KRATOS_CLASS_POINTER_DEFINITION(Ballast6Cluster3D);
  
      ///@}
      ///@name Life Cycle
      ///@{ 
      
        /// Default constructor.
        //Cluster3D() : Element() {}
    
        /**
         * Constructor using an array of nodes
         */
//         Cluster3D(IndexType NewId, const NodesArrayType& ThisNodes) :
//         Element(NewId, ThisNodes) 
//         {
//         }
//         
        Ballast6Cluster3D( );
        Ballast6Cluster3D( IndexType NewId, GeometryType::Pointer pGeometry );
        Ballast6Cluster3D( IndexType NewId, NodesArrayType const& ThisNodes);
        Ballast6Cluster3D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;      
        

        /**
         * Constructor using Geometry
         */
//         Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry) :
//         Element(NewId, pGeometry)
//         {
//         }

        /**
         * Constructor using Properties
         */
//         Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
//         Element(NewId, pGeometry , pProperties)
//         {
//         }

        /// Copy constructor.
//         Cluster3D (Cluster3D const& rOther) : Element(rOther)
//         {
//         }
      
   
      /// Destructor.
        virtual ~Ballast6Cluster3D();
      
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
        virtual void CustomInitialize();
        
        /**
         * this is called during the assembling process in order
         * to calculate all elemental contributions to the global system
         * matrix and the right hand side
         * @param rLeftHandSideMatrix: the elemental left hand side matrix
         * @param rRightHandSideVector: the elemental right hand side
         * @param rCurrentProcessInfo: the current process info instance
         */


        /**
         * this is called during the assembling process in order
         * to calculate the elemental right hand side vector only
         * @param rRightHandSideVector: the elemental right hand side vector
         * @param rCurrentProcessInfo: the current process info instance
         */
        virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
                ProcessInfo& rCurrentProcessInfo);

       
        virtual void EquationIdVector(EquationIdVectorType& rResult,
                ProcessInfo& rCurrentProcessInfo);
        

        /**
         * this is called during the assembling process in order
         * to calculate the elemental mass matrix
         * @param rMassMatrix: the elemental mass matrix
         * @param rCurrentProcessInfo: the current process info instance
         */
        virtual void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
        
        
        virtual void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo);
     

        virtual void GetDofList(DofsVectorType& ElementalDofList,
                ProcessInfo& CurrentProcessInfo);
        
      
        /**
         * this is called in the beginning of each solution step
         */
        virtual void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);


        /**
         * this is called at the end of each solution step
         */
        virtual void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);
        

        /**
         * deletes all obsolete data from memory
         */
   


        //output On integration points
        //calculate on element

        virtual void Calculate(const Variable<double>& rVariable,
                double& Output,
                const ProcessInfo& rCurrentProcessInfo);
    

        virtual void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                array_1d<double, 3 > & Output,
                const ProcessInfo& rCurrentProcessInfo);
    
        virtual void Calculate(const Variable<Vector >& rVariable,
                Vector& Output,
                const ProcessInfo& rCurrentProcessInfo);
     

        virtual void Calculate(const Variable<Matrix >& rVariable,
                Matrix& Output,
                const ProcessInfo& rCurrentProcessInfo);
   
        double SlowGetDensity();    


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
//
//        friend class Serializer;
//
//        virtual void save(Serializer& rSerializer) const
//        {
//            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Cluster3D );
//        }
//
//        virtual void load(Serializer& rSerializer)
//        {
//            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Cluster3D );
//        }

        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{    
      ///@}    
        
    }; // Class Ballast6Cluster3D

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, 
                    Ballast6Cluster3D& rThis);

  /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, 
                    const Ballast6Cluster3D& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BALLAST6CLUSTER3D_INCLUDED  defined
