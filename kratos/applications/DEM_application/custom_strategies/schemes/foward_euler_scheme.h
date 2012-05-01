//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson Lafontaine $
//   Date:                $Date: 2012-04-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//



#if !defined(KRATOS_FOWARD_EULER_SCHEME_H_INCLUDED )
#define  KRATOS_FOWARD_EULER_SCHEME_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

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
  
  class FowardEulerScheme  
    {
    public:
      ///@name Type Definitions
      ///@{
      
      typedef ModelPart::NodesContainerType NodesArrayType;
	
	
      /// Pointer definition of FowardEulerScheme
      KRATOS_CLASS_POINTER_DEFINITION(FowardEulerScheme);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      FowardEulerScheme(){}

      /// Destructor.
      virtual ~FowardEulerScheme(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      

     /// Its the same to do a loop`in nodes or element??? Need to be compared.  
     /// Need to check if the velocity or the dispalcement are the degree of freedon. Talk to M. Celigueta
     void Calculate(ModelPart& model_part)
     {
        KRATOS_TRY
        
	
	ProcessInfo& CurrentProcessInfo  = model_part.GetProcessInfo();
	NodesArrayType& pNodes           = model_part.Nodes(); 
        
	double aux     = 0;
	double delta_t =  CurrentProcessInfo[DELTA_TIME];  
        vector<unsigned int> node_partition;
	NodesArrayType::iterator it_begin = pNodes.ptr_begin();
	NodesArrayType::iterator it_end   = pNodes.ptr_end();
	int number_of_threads             = OpenMPUtils::GetNumThreads();
	OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);
	
	
	#pragma omp parallel for firstprivate(aux) shared(delta_t) 
	for(int k=0; k<number_of_threads; k++)
	{
	  NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	  NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
	  for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)      
	  {
	     array_1d<double, 3 > & vel          = i->FastGetSolutionStepValue(VELOCITY);
	     array_1d<double, 3 > & displ        = i->FastGetSolutionStepValue(DISPLACEMENT);
	     array_1d<double, 3 > & coor         = i->Coordinates();
  	     array_1d<double, 3 > & initial_coor = i->GetInitialPosition();
  	     array_1d<double, 3 > & force        = i->FastGetSolutionStepValue(FORCE);
	     const double& mass                  = i->FastGetSolutionStepValue(NODAL_MASS); 
	     aux = delta_t / mass;
	     
	     //Evolution of position (u(n+1) = u(n) + v(n+0.5)*delta_t):
	     if(i->pGetDof(DISPLACEMENT_X)->IsFixed() == false){
	         vel[0]    += aux * force[0]; 
	         displ[0]  += delta_t * vel[0];  
	         coor[0]   += initial_coor[0] + displ[0];
	     }
	     
	     if( i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false)
	     {
	         vel[1]    += aux * force[1];
	         displ[1]  += delta_t * vel[1];  
	         coor[1]   += initial_coor[1] + displ[1];
	     }
	     
	     if( i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false)
	     {
	         vel[2]    += aux * force[2];
	         displ[2]  += delta_t * vel[2];  
	         coor[2]   += initial_coor[2] + displ[2];
	     }
	   }
	}
	KRATOS_CATCH(" ")
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

      /// Turn back information as a string.
      virtual std::string Info() const
      {
	std::stringstream buffer;
        buffer << "FowardEulerScheme" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "FowardEulerScheme";}

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
     FowardEulerScheme& operator=(FowardEulerScheme const& rOther)
     {
       return *this;
     }

      /// Copy constructor.
      FowardEulerScheme(FowardEulerScheme const& rOther)
      {
	*this = rOther;
      }

        
      ///@}    
        
    }; // Class FowardEulerScheme 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    FowardEulerScheme& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const FowardEulerScheme& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_FOWARD_EULER_SCHEME_H_INCLUDED  defined 

