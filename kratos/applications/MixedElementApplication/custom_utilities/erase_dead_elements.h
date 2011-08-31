//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_ERASE_DEAD_ELEMENTS_UTILITY_H_INCLUDED )
#define  KRATOS_ERASE_DEAD_ELEMENTS_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"


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
  class EraseDeadElementsUtility
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of EraseDeadElementsUtility
      KRATOS_CLASS_POINTER_DEFINITION(EraseDeadElementsUtility);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      EraseDeadElementsUtility(){}

      /// Destructor.
      virtual ~EraseDeadElementsUtility(){}
      

      ///@}
      ///@name Operators 
      ///@{
      void Erase(double damage_limit,ModelPart& rmodel_part)
      {
          KRATOS_TRY

          //first of all mark all nodes so that they are NOT preserved
          for(ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin(); it!=rmodel_part.NodesEnd(); it++)
          {
              it->SetValue(IS_INACTIVE,true);
          }

          //now loop on all elements and find min damage, if min_damage is below threshold keep them, and mark their nodes to be kept
          for(ModelPart::ElementsContainerType::iterator it=rmodel_part.ElementsBegin(); it!=rmodel_part.ElementsEnd(); it++)
          {
              double min_damage = 1.0;
              Geometry<Node<3>>& geom = it->GetGeometry(); it++)

              for(unsigned int i=0; i<geom.size(); i++)
              {
                  double damage = GetGeometry()[i].FastGetSolutionStepValue(DAMAGE);
                  if(damage < min_damage)
                      min_damage = damage;
              }

              if(min_damage > damage_limit)
              {
                  it->SetValue(ERASE_FLAG,true);
              }
              else
              {
                  it->SetValue(ERASE_FLAG,false);
                  for(unsigned int i=0; i<geom.size(); i++)
                      geom[i].SetValue(IS_INACTIVE,false);
              }

          }

          //remove marked elements

          //remove marked nodes
          for(ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin(); it!=rmodel_part.NodesEnd(); it++)
          {
              if( it->GetValue(IS_INACTIVE) == true)
                  remove element
          }

          //remove conditions

          std::vector<double> values;
          ProcessInfo rCurrentProcessInfo = rmodel_part.GetProcessInfo();
          
          //loop on all elements, gather the constitutive law and use it to replace nodal values
          for(ModelPart::ElementsContainerType::iterator it=rmodel_part.ElementsBegin(); it!=rmodel_part.ElementsEnd(); it++)
          {
              //get constitutive law data
              it->GetValueOnIntegrationPoints(rVariable, values, rCurrentProcessInfo);


              //get geometry
              Geometry<Node<3> >& geom = it->GetGeometry();

              if(values.size() >= geom.size())
              {
                  for(unsigned int i=0; i<geom.size(); i++)
                  {
//                      std::cout << values[i] << std::endl;
                      geom[i].FastGetSolutionStepValue(rVariable) = values[i];
                  }
              }
          }

          KRATOS_CATCH("")
      }
      
      
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
      virtual std::string Info() const
      {
	std::stringstream buffer;
        buffer << "EraseDeadElementsUtility" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "EraseDeadElementsUtility";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const {}
      
            
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
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
//      EraseDeadElementsUtility& operator=(EraseDeadElementsUtility const& rOther){}

      /// Copy constructor.
      EraseDeadElementsUtility(EraseDeadElementsUtility const& rOther){}

        
      ///@}    
        
    }; // Class EraseDeadElementsUtility

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    EraseDeadElementsUtility& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const EraseDeadElementsUtility& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_ERASE_DEAD_ELEMENTS_UTILITY_H_INCLUDED  defined


