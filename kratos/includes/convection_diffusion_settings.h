//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_CONVECTION_DIFFUSION_SETTINGS_INCLUDED )
#define  KRATOS_CONVECTION_DIFFUSION_SETTINGS_INCLUDED


 
// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"


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
  class ConvectionDiffusionSettings
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of ConvectionDiffusionSettings
      KRATOS_CLASS_POINTER_DEFINITION(ConvectionDiffusionSettings);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      ConvectionDiffusionSettings(){};
      ConvectionDiffusionSettings(const ConvectionDiffusionSettings& rOther):
	    mpDensityVar(rOther.mpDensityVar), 
	    mpDiffusionVar(rOther.mpDensityVar),
	    mpUnknownVar(rOther.mpUnknownVar),
	    mpVolumeSourceVar(rOther.mpVolumeSourceVar),
	    mpSurfaceSourceVar(rOther.mpSurfaceSourceVar),
	    mpProjectionVar(rOther. mpProjectionVar),
	    mpConvectionVar(rOther.mpConvectionVar),
	    mpMeshVelocityVar(rOther.mpMeshVelocityVar)
	    {
	    }

      /// Destructor.
      virtual ~ConvectionDiffusionSettings(){};
      
      ///@}
      ///@name Operators 
      ///@{
      void SetDensityVariable(const Variable<double>& rvar){mpDensityVar = &rvar;}
      const Variable<double>& GetDensityVariable(){return *mpDensityVar;}
      
      void SetDiffusionVariable(const Variable<double>& rvar){mpDiffusionVar = &rvar;}
      const Variable<double>& GetDiffusionVariable(){return *mpDiffusionVar;}
      
      void SetUnknownVariable(const Variable<double>& rvar){mpUnknownVar = &rvar;}
      const Variable<double>& GetUnknownVariable(){return *mpUnknownVar;}
      
      void SetVolumeSourceVariable(const Variable<double>& rvar){mpVolumeSourceVar = &rvar;}
      const Variable<double>& GetVolumeSourceVariable(){return *mpVolumeSourceVar;}
      
      void SetSurfaceSourceVariable(const Variable<double>& rvar){mpSurfaceSourceVar = &rvar;}
      const Variable<double>& GetSurfaceSourceVariable(){return *mpSurfaceSourceVar;}
      
      void SetProjectionVariable(const Variable<double>& rvar){mpProjectionVar = &rvar;}
      const Variable<double>& GetProjectionVariable(){return *mpProjectionVar;}
      
      void SetConvectionVariable(const Variable<array_1d<double,3> >& rvar){mpConvectionVar = &rvar;}
      const Variable<array_1d<double,3> >& GetConvectionVariable(){return *mpConvectionVar;}
      
      void SetMeshVelocityVariable(const Variable<array_1d<double,3> >& rvar){mpMeshVelocityVar = &rvar;}
      const Variable<array_1d<double,3> >& GetMeshVelocityVariable(){return *mpMeshVelocityVar;}
     
      ///@}
      ///@name Operations
      ///@{
      
      
      ///@}
      ///@name Access
      ///@{ 
      /// Assignment operator.
      ConvectionDiffusionSettings& operator=(ConvectionDiffusionSettings const& rOther)
      {
	    mpDensityVar = rOther.mpDensityVar;
	    mpDiffusionVar = rOther.mpDensityVar;
	    mpUnknownVar = rOther.mpUnknownVar;
	    mpVolumeSourceVar = rOther.mpVolumeSourceVar;
	    mpSurfaceSourceVar = rOther.mpSurfaceSourceVar;
	    mpProjectionVar = rOther. mpProjectionVar;
	    mpConvectionVar = rOther.mpConvectionVar;
	    mpMeshVelocityVar = rOther.mpMeshVelocityVar;
	    return *this;
      }
      
      
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
	  buffer << "ConvectionDiffusionSettings #" ;
	  return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
	  rOStream << "ConvectionDiffusionSettings #";
      }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
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
      const Variable<double>* mpDensityVar;
      const Variable<double>* mpDiffusionVar;
      const Variable<double>* mpUnknownVar;
      const Variable<double>* mpVolumeSourceVar;
      const Variable<double>* mpSurfaceSourceVar;
      const Variable<double>* mpProjectionVar;
      const Variable<array_1d<double,3> >* mpConvectionVar;
      const Variable<array_1d<double,3> >* mpMeshVelocityVar;
        
        
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
      
      


        
      ///@}    
        
    }; // Class ConvectionDiffusionSettings 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    ConvectionDiffusionSettings& rThis)
    {return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const ConvectionDiffusionSettings& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
  
  
  KRATOS_DEFINE_VARIABLE(ConvectionDiffusionSettings::Pointer, CONVECTION_DIFFUSION_SETTINGS)
  
}  // namespace Kratos.

#endif // KRATOS_CONVECTION_DIFFUSION_SETTINGS_INCLUDED  defined 


