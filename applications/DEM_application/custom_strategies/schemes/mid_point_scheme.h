//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson Lafontaine $
//   Date:                $Date: 2012-04-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//



#if !defined(KRATOS_MID_POINT__SCHEME_H_INCLUDED )
#define  KRATOS_MID_POINT__SCHEME_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "integration_scheme.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

#include "DEM_application.h"

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
  
  class MidPointScheme :  public IntegrationScheme
    {
    public:
      ///@name Type Definitions
      ///@{
      
      typedef ModelPart::NodesContainerType    NodesArrayType;
      typedef ModelPart::ElementsContainerType ElementsArrayType;
       
	
      /// Pointer definition of MidPointScheme
      KRATOS_CLASS_POINTER_DEFINITION(MidPointScheme);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      MidPointScheme(){}

      /// Destructor.
      virtual ~MidPointScheme(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      

     /// Its the same to do a loop`in nodes or element??? Need to be compared.  
     /// Need to check if the velocity or the dispalcement are the degree of freedon. Talk to M. Celigueta
     /// WARNING = TO BE CHECK IT
     void Calculate(ModelPart& model_part)
     {
        KRATOS_TRY
        
        typedef std::vector<array_1d<double, 3 > > ComponentVectorType;
        typedef std::vector<array_1d<double, 3 > >::iterator ComponentIteratorType;
         
	ProcessInfo& CurrentProcessInfo  = model_part.GetProcessInfo();
        NodesArrayType& pNodes           = model_part.Nodes();

 	double delta_t      =  CurrentProcessInfo[DELTA_TIME];
        double half_delta_t = 0.5 * delta_t;
	array_1d<double, 3 > aux, vel_copy, displ_copy;
	ComponentVectorType vel_old, displ_new, kf, kv;

        vector<unsigned int> node_partition;
	NodesArrayType::iterator it_begin = pNodes.ptr_begin();
	NodesArrayType::iterator it_end   = pNodes.ptr_end();
	int number_of_threads             = 1; //OpenMPUtils::GetNumThreads();
	OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);
	
	#pragma omp parallel for private(aux) shared(delta_t) 
	for(int k=0; k<number_of_threads; k++)
	{
	  NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	  NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
	  for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
	  {
	    //Element::GeometryType& geom           = i->GetGeometry();
	    array_1d<double, 3 > & vel            = i->FastGetSolutionStepValue(VELOCITY);
	    array_1d<double, 3 > & vel_old        = i->FastGetSolutionStepValue(VELOCITY,1);
	    array_1d<double, 3 > & displ          = i->FastGetSolutionStepValue(DISPLACEMENT);
	    array_1d<double, 3 > & displ_old      = i->FastGetSolutionStepValue(DISPLACEMENT,1);
	    array_1d<double, 3 > & force          = i->FastGetSolutionStepValue(RHS);
	    array_1d<double, 3 > & coor           = i->Coordinates();
	    array_1d<double, 3 > & initial_coor   = i->GetInitialPosition();

	    const double& mass                    = i->FastGetSolutionStepValue(NODAL_MASS);
	    noalias(aux)                          = (half_delta_t/ mass) * force;
	    
	    if( i->pGetDof(VELOCITY_X)->IsFixed() == false  )
            {
	      vel[0]    += aux[0]; 
	      displ[0]  += half_delta_t * vel[0];
	    }

	    if( i->pGetDof(VELOCITY_Y)->IsFixed() == false  )
            {
	      vel[1]    += aux[1]; 
	      displ[1]  += half_delta_t * vel[1];
	    }

	    if( i->pGetDof(VELOCITY_Z)->IsFixed() == false  )
            {
	    vel[2]    += aux[2];
	    displ[2]  += half_delta_t * vel[2];  
	    }
	    
	    /// TALK TO M. Angel
	    //Calculate_Forces(type_id, damp_id, delta_t, gravity);
	    //i->Calculate(FORCE, force, CurrentProcessInfo);
	    
	    if( i->pGetDof(VELOCITY_X)->IsFixed() == false  )
            {
	        vel[0]    = vel_old[0]      + (delta_t/mass) * force[0];
	        displ[0]  = displ_old[0]    + delta_t * vel_old[0] * (1 + half_delta_t);
		coor[0]   = initial_coor[0] + displ[0];
	    }

	    if(  i->pGetDof(VELOCITY_Y)->IsFixed() == false  )
            {
	        vel[1]    = vel_old[1]      + (delta_t/mass) * force[1];
	        displ[1]  = displ_old[1]    + delta_t * vel_old[1] * (1 + half_delta_t);
		coor[1]   = initial_coor[1] + displ[1];
	    }

	    if(  i->pGetDof(VELOCITY_Z)->IsFixed() == false  )
            {
	        vel[2]    = vel_old[2]      + (delta_t/mass) * force[2];
	        displ[2]  = displ_old[2]    + delta_t * vel_old[2] * (1 + half_delta_t);
		coor[2]   = initial_coor[2] + displ[2]; 
	      }
	    }
	  }
	  
	  KRATOS_CATCH("")
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
        buffer << "MidPointScheme" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "MidPointScheme";}

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
     MidPointScheme& operator=(MidPointScheme const& rOther)
     {
       return *this;
     }

      /// Copy constructor.
      MidPointScheme(MidPointScheme const& rOther)
      {
	*this = rOther;
      }

        
      ///@}    
        
    }; // Class MidPointScheme

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    MidPointScheme& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const MidPointScheme& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_MID_POINT__SCHEME_H_INCLUDED  defined 

