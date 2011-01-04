//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_WAVEGENERATOR_H_INCLUDED )
#define  KRATOS_WAVEGENERATOR_H_INCLUDED



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
  class WaveGenerator
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of WaveGenerator
      KRATOS_CLASS_POINTER_DEFINITION(WaveGenerator);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      WaveGenerator(){}

      /// Destructor.
      virtual ~WaveGenerator(){}
      
      void GenerateWaveXYPlane(ModelPart::NodesContainerType& rNodes, const double d, const double  H, const double  T, const double  Z0, const double t, const double g)
      {
	double C0=g*T/( 6.2831853071795862 );
	double L0=C0*T;
	double teta = (-6.2831853071795862*t / T) - (3.1415926535897931 / 2.0);
	double aux = cosh(6.2831853071795862*d/L0);
	double ux= 0.5 * H * (g * T / L0) * cos(teta) / aux;
	double uy= 0.5 * H * (g * T / L0) * sin(teta) / aux;
	double h=0.5 * H * cos(teta);
	
	const PointerVector< Node<3> >::iterator it_begin = rNodes.begin();
	array_1d<double,3> temp;
	
	#pragma omp parallel for private(temp)
	for(int i=0; i<static_cast<int>(rNodes.size()); i++)
	{
	    
	    PointerVector< Node<3> >::iterator it = it_begin + i;
	    const double Z = it->Z();
	    temp[0]  = ux * cosh(6.2831853071795862*(Z - Z0 + d)/L0);
	    temp[1]  = 0.0;
	    temp[2]  = uy * sinh(6.2831853071795862*(Z - Z0 + d)/L0);
	    
	    noalias(it->FastGetSolutionStepValue(VELOCITY))=temp;
	    it->FastGetSolutionStepValue(DISTANCE,1)=-h+Z - Z0;
	}
      }

      ///@}
      ///@name Operators 
      ///@{
      
      
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
        buffer << "WaveGenerator" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "WaveGenerator";}

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
      WaveGenerator& operator=(WaveGenerator const& rOther){return *this;}

      /// Copy constructor.
      WaveGenerator(WaveGenerator const& rOther){};

        
      ///@}    
        
    }; // Class WaveGenerator 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    WaveGenerator& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const WaveGenerator& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_WAVEGENERATOR_H_INCLUDED  defined 


