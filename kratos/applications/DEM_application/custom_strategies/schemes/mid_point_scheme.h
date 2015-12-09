//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//


#if !defined(KRATOS_MID_POINT__SCHEME_H_INCLUDED )
#define  KRATOS_MID_POINT__SCHEME_H_INCLUDED

// System includes
#include <string>
#include <iostream> 

// External includes 

// Project includes
#include "dem_integration_scheme.h"
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
  
  class MidPointScheme :  public DEMIntegrationScheme
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
      

     //NOTE: THIS IS NOT EXACTLY CENTRAL DIFFERENCES
     
     void Calculate(ModelPart& model_part)
      {
          ProcessInfo& rCurrentProcessInfo  = model_part.GetProcessInfo();

          CalculateTranslationalMotion(model_part);

          if(rCurrentProcessInfo[ROTATION_OPTION]!=0)
          {
             CalculateRotationalMotion(model_part);
          }

      }
    
    void CalculateTranslationalMotion(ModelPart& model_part)
    {  
      
    KRATOS_TRY
        
    //typedef std::vector<array_1d<double, 3 > > ComponentVectorType;
    //typedef std::vector<array_1d<double, 3 > >::iterator ComponentIteratorType;
         
	ProcessInfo& CurrentProcessInfo  = model_part.GetProcessInfo();
        NodesArrayType& pNodes           = model_part.Nodes();

 	double delta_t      =  CurrentProcessInfo[DELTA_TIME];
    double half_delta_t = 0.5 * delta_t;
	array_1d<double, 3 > aux, vel_copy;
    
    vector<unsigned int> node_partition;

	#ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
    #else
    int number_of_threads = 1;
    #endif
	OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);
	
	#pragma omp parallel for private(aux) shared(delta_t) 
	for(int k=0; k<number_of_threads; k++)
	{
	  NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	  NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
	  for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
	  {

        array_1d<double, 3 > & vel             = i->FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3 > & displ           = i->FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3 > & delta_displ     = i->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        array_1d<double, 3 > & coor            = i->Coordinates();
        array_1d<double, 3 > & initial_coor    = i->GetInitialPosition();
        array_1d<double, 3 > & force           = i->FastGetSolutionStepValue(TOTAL_FORCES);

        double mass                           = i->FastGetSolutionStepValue(NODAL_MASS);

	    noalias(aux)                          = (half_delta_t/ mass) * force;
	    
	   
	    if( i->pGetDof(VELOCITY_X)->IsFixed() == false  )
        {
	         vel_copy[0]    = vel [0] + aux[0];
              
	         delta_displ[0]  = delta_t * vel_copy[0];
              
                 displ[0]   += delta_displ[0];
              
                 coor[0]   = initial_coor[0] + displ[0];
              
                 vel[0] += (delta_t/ mass) * force[0];
	    }

            else
            {
                                
                displ[0]  += delta_displ[0];
                
                delta_displ[0] = delta_t * vel[0];

                coor[0]   = initial_coor[0] + displ[0];
            }

	    if( i->pGetDof(VELOCITY_Y)->IsFixed() == false  )
            {
                vel_copy[1]    = vel [1] + aux[1];
              
                delta_displ[1]  = delta_t * vel_copy[1];
              
                 displ[1]   +=delta_displ[1];
              
                 coor[1]   = initial_coor[1] + displ[1];
              
                 vel[1] += (delta_t/ mass) * force[1];
	    }

            else
            {
                 
                 displ[1]  += delta_displ[1];
                 
                 delta_displ[1] = delta_t * vel[1];
                  
                 coor[1]   = initial_coor[1] + displ[1];
            }


	    if( i->pGetDof(VELOCITY_Z)->IsFixed() == false  )
            {
                vel_copy[2]    = vel [2] + aux[2];
              
                delta_displ[2]  = delta_t * vel_copy[2];
              
                 displ[2]   +=delta_displ[2];
              
                 coor[2]   = initial_coor[2] + displ[2];
              
                 vel[2] += (delta_t/ mass) * force[2];
	    }    

            else
            {
                displ[2]  += delta_displ[2]; 
              
                delta_displ[2] = delta_t * vel[2];    

                coor[2]   = initial_coor[2] + displ[2];
            }            
	
	    } //End nodes loop
	  } //End openMP loop
	  
	  KRATOS_CATCH("")
	}


      void CalculateRotationalMotion(ModelPart& model_part)
      {
          KRATOS_TRY   

                
        ProcessInfo& rCurrentProcessInfo  = model_part.GetProcessInfo();
        NodesArrayType& pNodes           = model_part.Nodes();

    
    double delta_t =  rCurrentProcessInfo[DELTA_TIME];
    double half_delta_t = 0.5 * delta_t;
    array_1d<double, 3 > AngularVelAux;
    
     vector<unsigned int> node_partition;

    #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
    #else
        int number_of_threads = 1;
    #endif
    OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

    #pragma omp parallel for shared(delta_t)
    for(int k=0; k<number_of_threads; k++)
    {
            NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
            NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
            for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
            {

                double PMomentOfInertia = i->FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
                double coeff            = rCurrentProcessInfo[NODAL_MASS_COEFF];

                array_1d<double, 3 > & AngularVel             = i->FastGetSolutionStepValue(ANGULAR_VELOCITY);
                array_1d<double, 3 > & RotaMoment             = i->FastGetSolutionStepValue(PARTICLE_MOMENT);
                //array_1d<double, 3 > & delta_rotation_displ   = i->FastGetSolutionStepValue(DELTA_ROTA_DISPLACEMENT);
                array_1d<double, 3 > & Rota_Displace          = i->FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
                array_1d<double, 3 > delta_rotation_displ;                   
                //double & Orientation_real                     = i->FastGetSolutionStepValue(ORIENTATION_REAL); 
                //array_1d<double, 3 > & Orientation_imag       = i->FastGetSolutionStepValue(ORIENTATION_IMAG);                
                
                bool If_Fix_Rotation[3] = {false, false, false};
                If_Fix_Rotation[0] = i->pGetDof(ANGULAR_VELOCITY_X)->IsFixed();
                If_Fix_Rotation[1] = i->pGetDof(ANGULAR_VELOCITY_Y)->IsFixed();
                If_Fix_Rotation[2] = i->pGetDof(ANGULAR_VELOCITY_Z)->IsFixed();
                
                for(std::size_t iterator = 0 ; iterator < 3; iterator++)
                {
                    if(If_Fix_Rotation[iterator] == false)
                    {
                         double RotaAcc = 0.0;
                         
                         RotaAcc = (RotaMoment[iterator]) / (PMomentOfInertia);

                         if(rCurrentProcessInfo[VIRTUAL_MASS_OPTION])
                         {
                                  RotaAcc = RotaAcc * ( 1 - coeff );
                         }
                       
                       
                         AngularVelAux[iterator] = AngularVel[iterator] + RotaAcc * half_delta_t;
                         
                         delta_rotation_displ[iterator] = AngularVelAux[iterator] * delta_t;
                         
                         Rota_Displace[iterator] +=  delta_rotation_displ[iterator];  
                         
                         AngularVel[iterator] += RotaAcc * delta_t;
      
       
                    }                   

                    else
                    {

                        delta_rotation_displ[iterator]= 0.0;
                                           
                        /*
                       *
                       *
                       *  implementation of fixed rotational motion.
                       *
                       *
                       */

                    }
                    
                    //RotaMoment[iterator] = 0.0;
                }
                
                //double RotationAngle;
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

