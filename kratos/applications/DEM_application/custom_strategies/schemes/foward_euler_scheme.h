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
#include "integration_scheme.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

#include "DEM_application.h"

namespace Kratos
{
  
  
  class FowardEulerScheme :  public IntegrationScheme
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
      


     /// Its the same to do a loop`in nodes or element??? Need to be compared.  
     /// Need to check if the velocity or the dispalcement are the degree of freedon. Talk to M. Celigueta
     void Calculate(ModelPart& model_part)
     {
        CalculateTranslationalMotion(model_part);
        CalculateRotationalMotion(model_part);
     }

     void CalculateTranslationalMotion(ModelPart& model_part)
     {
        KRATOS_TRY

	ProcessInfo& rCurrentProcessInfo  = model_part.GetProcessInfo();
	NodesArrayType& pNodes           = model_part.Nodes(); 
        
	double aux     = 0;
	double delta_t =  rCurrentProcessInfo[DELTA_TIME];

        vector<unsigned int> node_partition;
	NodesArrayType::iterator it_begin = pNodes.ptr_begin();
	NodesArrayType::iterator it_end   = pNodes.ptr_end();
	int number_of_threads             = 1; //OpenMPUtils::GetNumThreads();
	OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);
	
	
	#pragma omp parallel for firstprivate(aux) shared(delta_t) 
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
  	     array_1d<double, 3 > & force           = i->FastGetSolutionStepValue(RHS);

	     const double mass                      = i->FastGetSolutionStepValue(NODAL_MASS);
            
	     aux = delta_t / mass;
          
	   
	     if( i->pGetDof(VELOCITY_X)->IsFixed() == false ) // equivalently:  i->IsFixed(VELOCITY_X) == false
             {
	         vel[0]    += aux * force[0];

                 delta_displ[0] = delta_t * vel[0];
	         
                 displ[0]  +=  delta_displ[0];

	         coor[0]   = initial_coor[0] + displ[0];
                 
             }
             else
             {
                 delta_displ[0] = delta_t * vel[0];

                 displ[0]  += delta_displ[0];

	         coor[0]   = initial_coor[0] + displ[0];

             }
	     
	     if(  i->pGetDof(VELOCITY_Y)->IsFixed() == false  )
             {
	         vel[1]    += aux * force[1];

                 delta_displ[1] = delta_t * vel[1];

                 displ[1]  +=  delta_displ[1];

	         coor[1]   = initial_coor[1] + displ[1];

             }
             else
             {
                 delta_displ[1] = delta_t * vel[1];

                 displ[1]  += delta_displ[1];

	         coor[1]   = initial_coor[1] + displ[1];

             }
	     
             if(  i->pGetDof(VELOCITY_Z)->IsFixed() == false  )
	     {
	         vel[2]    += aux * force[2];

                 delta_displ[2] = delta_t * vel[2];

                 displ[2]  +=  delta_displ[2];

	         coor[2]   = initial_coor[2] + displ[2];
               
             }
             else
             {
                 delta_displ[2] = delta_t * vel[2];

                 displ[2]  += delta_displ[2];

	         coor[2]   = initial_coor[2] + displ[2];

             }
	   }
	}
	KRATOS_CATCH(" ")
     }


     void CalculateRotationalMotion(ModelPart& model_part)
     {
        KRATOS_TRY   
     
	ProcessInfo& rCurrentProcessInfo  = model_part.GetProcessInfo();
	NodesArrayType& pNodes           = model_part.Nodes();

	
	double delta_t =  rCurrentProcessInfo[DELTA_TIME];

        //KRATOS_WATCH(delta_t)

        vector<unsigned int> node_partition;
	NodesArrayType::iterator it_begin = pNodes.ptr_begin();
	NodesArrayType::iterator it_end   = pNodes.ptr_end();
	int number_of_threads             = 1; //OpenMPUtils::GetNumThreads();
	OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);


	#pragma omp parallel for shared(delta_t)
	for(int k=0; k<number_of_threads; k++)
	{
            NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
            NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
            for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
            {

                double PMass            = i->FastGetSolutionStepValue(NODAL_MASS);
                double PMomentOfInertia = i->FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);

                array_1d<double, 3 > & AngularVel      = i->FastGetSolutionStepValue(ANGULAR_VELOCITY);
                array_1d<double, 3 > & RotaMoment      = i->FastGetSolutionStepValue(PARTICLE_MOMENT);
                array_1d<double, 3 > & Rota_Displace   = i->FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);

                bool If_Fix_Rotation[3] = {false, false, false};
                If_Fix_Rotation[0] = i->pGetDof(VELOCITY_X)->IsFixed();
                If_Fix_Rotation[1] = i->pGetDof(VELOCITY_Y)->IsFixed();
                If_Fix_Rotation[2] = i->pGetDof(VELOCITY_Z)->IsFixed();

                

                for(std::size_t iterator = 0 ; iterator < 3; iterator++)
                {
                    if(If_Fix_Rotation[iterator] == false)
                    {
                         double RotaAcc = 0.0;
                         RotaAcc = (RotaMoment[iterator]) / PMass / PMomentOfInertia;

                         double RotaVelOld = AngularVel[iterator];
                         double RotaVelNew = RotaVelOld + RotaAcc * delta_t;

                         AngularVel[iterator]  = 0.5 * (RotaVelOld + RotaVelNew);
                     

                         Rota_Displace[iterator] += AngularVel[iterator] * delta_t / M_PI * 180.0;
                     
                    }
                   

                    else
                    {
                       /*
                       *
                       *
                       *  implementation of fixed rotational motion.
                       *
                       *
                       */

                    }
                     RotaMoment[iterator] = 0.0;


                }
                    
            }

        }
        KRATOS_CATCH(" ")
     }
	
   


     

     
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

