
#if !defined(KRATOS_CONSTANT_ACERAGE_ACCELERATION_SCHEME_H_INCLUDED )
#define  KRATOS_CONSTANT_ACERAGE_ACCELERATION_SCHEME_H_INCLUDED

// System includes
#include <string>
#include <iostream> 

// External includes 

// Project includes
#include "dem_integration_scheme.h"
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"

namespace Kratos
{
  
  
  class ConstAverageAccelerationScheme :  public DEMIntegrationScheme
    {
    public:
      ///@name Type Definitions
      ///@{
      
      typedef ModelPart::NodesContainerType NodesArrayType;
	
	
      /// Pointer definition of ConstAverageAccelerationScheme
      KRATOS_CLASS_POINTER_DEFINITION(ConstAverageAccelerationScheme);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      ConstAverageAccelerationScheme(){}

      /// Destructor.
      virtual ~ConstAverageAccelerationScheme(){}
      
      DEMIntegrationScheme* CloneRaw() const override {
            DEMIntegrationScheme* cloned_scheme(new ConstAverageAccelerationScheme(*this));
            return cloned_scheme;
      }
      
      DEMIntegrationScheme::Pointer CloneShared() const override {
            DEMIntegrationScheme::Pointer cloned_scheme(new ConstAverageAccelerationScheme(*this));
            return cloned_scheme;
        }
      
      void SetIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose = true) const override {
        if(verbose) std::cout << "\nAssigning ConstAverageAccelerationScheme to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }
            
     /// Its the same to do a loop`in nodes or element??? Need to be compared.  
     /// Need to check if the velocity or the dispalcement are the degree of freedon. Talk to M. Celigueta
     /*void Calculate(ModelPart& model_part, bool TRotationOption, int StepFlag = -1)  override
     {
        KRATOS_TRY
        
	ProcessInfo& r_process_info  = model_part.GetProcessInfo();
	NodesArrayType& pNodes           = model_part.Nodes(); 
        
	//double aux          = 0;
        array_1d<double, 3 >  new_accel;
        array_1d<double, 3 >  prev_accel;
	double delta_t      =  r_process_info[DELTA_TIME];

        vector<unsigned int> node_partition;
	//NodesArrayType::iterator it_begin = pNodes.ptr_begin();
	//NodesArrayType::iterator it_end   = pNodes.ptr_end();
	int number_of_threads             = 1; //OpenMPUtils::GetNumThreads();
	OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);
	
	
        #pragma omp parallel for 
	for(int k=0; k<number_of_threads; k++)
	{
	  NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
	  NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
	  for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)      
	  {
               
	     array_1d<double, 3 > & vel             = i->FastGetSolutionStepValue(VELOCITY);
	     array_1d<double, 3 > & displ           = i->FastGetSolutionStepValue(DISPLACEMENT);
             array_1d<double, 3 > & delta_displ     = i->FastGetSolutionStepValue(DELTA_DISPLACEMENT); // (DISPLACEMENT,1) can not be used as its the same as DISPLACEMENT
             array_1d<double, 3 > & coor            = i->Coordinates();
  	     array_1d<double, 3 > & initial_coor    = i->GetInitialPosition();
  	     array_1d<double, 3 > & force           = i->FastGetSolutionStepValue(RHS);
	     array_1d<double, 3 > & prev_force      = i->FastGetSolutionStepValue(RHS,1); // (RHS,1) is different from RHS which has been calculated in the previous step.
        double mass                            = i->FastGetSolutionStepValue(NODAL_MASS);

             
	     //             aux = delta_t / mass;
            
	     new_accel = force / mass;
             prev_accel = prev_force / mass;

             //velocidad = c1->GetVelocidad() + 0.5 * dt * (c1->GetAceleracion() + accel);
	     //desplazamiento = dt * c1->GetVelocidad() + 0.5 * dt * dt * (c1->GetAceleracion() + accel);
	     //if( i->pGetDof(VELOCITY_X)->IsFixed() == false )
             if( i->IsNot(DEMFlags::FIXED_VEL_X))
             {

                 delta_displ[0] = delta_t * vel[0] + 0.25 * delta_t * delta_t * (prev_accel[0] + new_accel[0]);

                 displ[0]  += delta_displ[0];

	         vel[0]    = vel[0] + 0.5 * delta_t * (prev_accel[0] + new_accel[0]);
                 
	         coor[0]   = initial_coor[0] + displ[0];
                 

		 //prev_accel[0] = new_accel[0];
                 
             }
                             
             else
             {
               
                 delta_displ[0] = delta_t * vel[0];
                 
                 displ[0]  += delta_displ[0];
                 
                 coor[0]   = initial_coor[0] + displ[0];
                 

             }
	     
	     //if( i->pGetDof(VELOCITY_Y)->IsFixed() == false )
             if( i->IsNot(DEMFlags::FIXED_VEL_Y))
             {

                 delta_displ[1] = delta_t * vel[1] + 0.25 * delta_t * delta_t * (prev_accel[1] + new_accel[1]);

                 displ[1]  += delta_displ[1];

	         vel[1]    = vel[1] + 0.5 * delta_t * (prev_accel[1] + new_accel[1]);

	         coor[1]   = initial_coor[1] + displ[1];
                 

		 //prev_accel[1] = new_accel[1];
	     }
              else
             {
                 delta_displ[1] = delta_t * vel[1];

                 displ[1]  += delta_displ[1];

                 coor[1]   = initial_coor[1] + displ[1];
                 
             }
	     
             //if( i->pGetDof(VELOCITY_Z)->IsFixed() == false )
             if( i->IsNot(DEMFlags::FIXED_VEL_Z))
	     {
                 
	         displ[2]  += delta_t * vel[2] + 0.25 * delta_t * delta_t * (prev_accel[2] + new_accel[2]);
                 
	         vel[2]    = vel[2] + 0.5 * delta_t * (prev_accel[2] + new_accel[2]);
                
	         coor[2]   = initial_coor[2] + displ[2];
                 
                 
		 //prev_accel[2] = new_accel[2];
                
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
     }  */
    
    
     
      /// Turn back information as a string.
      virtual std::string Info() const override
      {
	std::stringstream buffer;
        buffer << "ConstAverageAccelerationScheme" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const override{rOStream << "ConstAverageAccelerationScheme";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const override{}
      
            
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
     ConstAverageAccelerationScheme& operator=(ConstAverageAccelerationScheme const& rOther)
     {
       return *this;
     }

      /// Copy constructor.
      ConstAverageAccelerationScheme(ConstAverageAccelerationScheme const& rOther)
      {
	*this = rOther;
      }

        
      ///@}    
        
    }; // Class ConstAverageAccelerationScheme

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    ConstAverageAccelerationScheme& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const ConstAverageAccelerationScheme& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_CONSTANT_ACERAGE_ACCELERATION_SCHEME_H_INCLUDED  defined 

