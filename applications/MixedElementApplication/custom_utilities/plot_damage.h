//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_PLOT_NODAL_CONSTITUTIVELAW_VAR_UTILITY_H_INCLUDED )
#define  KRATOS_PLOT_NODAL_CONSTITUTIVELAW_VAR_UTILITY_H_INCLUDED



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
  class PlotUtility
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of PlotUtility
      KRATOS_CLASS_POINTER_DEFINITION(PlotUtility);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      PlotUtility(){}

      /// Destructor.
      virtual ~PlotUtility(){}
      

      ///@}
      ///@name Operators 
      ///@{
      void PlotVariable(const Variable<double>& rVariable, ModelPart& rmodel_part)
      {
          KRATOS_TRY

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
        buffer << "PlotUtility" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "PlotUtility";}

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
//      PlotUtility& operator=(PlotUtility const& rOther){}

      /// Copy constructor.
      PlotUtility(PlotUtility const& rOther){}

        
      ///@}    
        
    }; // Class PlotUtility

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    PlotUtility& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const PlotUtility& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_PLOT_NODAL_CONSTITUTIVELAW_VAR_UTILITY_H_INCLUDED  defined


