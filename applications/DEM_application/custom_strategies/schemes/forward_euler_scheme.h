//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson Lafontaine $
//   Date:                $Date: 2012-04-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//



#if !defined(KRATOS_FORWARD_EULER_SCHEME_H_INCLUDED )
#define  KRATOS_FORWARD_EULER_SCHEME_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cfloat>

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
  
  
  class ForwardEulerScheme :  public IntegrationScheme
    {
    public:
      ///@name Type Definitions
      ///@{
      
      typedef ModelPart::NodesContainerType NodesArrayType;
    
    
      /// Pointer definition of ForwardEulerScheme
      KRATOS_CLASS_POINTER_DEFINITION(ForwardEulerScheme);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      ForwardEulerScheme(){}

      /// Destructor.
      virtual ~ForwardEulerScheme(){}
      
      virtual NodesArrayType& GetNodes(ModelPart& model_part)
      {
          return model_part.GetCommunicator().LocalMesh().Nodes(); 
      }
      
      virtual NodesArrayType& GetGhostNodes(ModelPart& model_part)
      {
          return model_part.GetCommunicator().GhostMesh().Nodes(); 
      }
      
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

          ProcessInfo& rCurrentProcessInfo  = model_part.GetProcessInfo();
          NodesArrayType& pNodes            = GetNodes(model_part);

          double delta_t              = rCurrentProcessInfo[DELTA_TIME];

          double virtual_mass_coeff   = rCurrentProcessInfo[NODAL_MASS_COEFF];
          bool if_virtual_mass_option = (bool) rCurrentProcessInfo[VIRTUAL_MASS_OPTION];

          vector<unsigned int> node_partition;
          //NodesArrayType::iterator it_begin = pNodes.ptr_begin();
          //NodesArrayType::iterator it_end   = pNodes.ptr_end();

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
                  array_1d<double, 3 > & vel             = i->FastGetSolutionStepValue(VELOCITY);
                  array_1d<double, 3 > & displ           = i->FastGetSolutionStepValue(DISPLACEMENT);
                  array_1d<double, 3 > & delta_displ     = i->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                  array_1d<double, 3 > & coor            = i->Coordinates();
                  array_1d<double, 3 > & initial_coor    = i->GetInitialPosition();
                  array_1d<double, 3 > & force           = i->FastGetSolutionStepValue(TOTAL_FORCES);

                  double mass                            = i->FastGetSolutionStepValue(SQRT_OF_MASS);
                  mass                                  *= mass;
                  
                  double aux = delta_t / mass;     
                  
                  //i->FastGetSolutionStepValue(OLD_COORDINATES) = coor; //saving the coordinates in order to optimize some functions (specially de previous step coordinates)  

                  if (if_virtual_mass_option)
                  {
                      aux = (1 - virtual_mass_coeff)* (delta_t / mass);
                      if (aux<0.0) KRATOS_ERROR(std::runtime_error,"The coefficient assigned for vitual mass is larger than one, virtual_mass_coeff= ",virtual_mass_coeff)
                  }
                  
                  //unsigned int pos = i->FastGetSolutionStepValue(VELOCITY_X_DOF_POS);
                  //if( i->GetDof(VELOCITY_X, pos).IsFixed() == false ) // equivalently:  i->IsFixed(VELOCITY_X) == false
                  if( i->IsNot(DEMFlags::FIXED_VEL_X))
                  {    
                      vel[0] += aux * force[0];

                      displ[0] +=  delta_displ[0];

                      delta_displ[0] = delta_t * vel[0];

                      coor[0] = initial_coor[0] + displ[0];


                  }
                  else
                  {
                      
                      displ[0] += delta_displ[0];
                      
                      delta_displ[0] = delta_t * vel[0];

                      coor[0] = initial_coor[0] + displ[0];
                 
                  }
                  //pos = i->FastGetSolutionStepValue(VELOCITY_Y_DOF_POS);
                  //if( i->GetDof(VELOCITY_Y, pos).IsFixed() == false ) 
                  if( i->IsNot(DEMFlags::FIXED_VEL_Y))
                  {    
                      vel[1] += aux * force[1];

                      displ[1] +=  delta_displ[1];

                      delta_displ[1] = delta_t * vel[1];

                      coor[1] = initial_coor[1] + displ[1];
                      

                  }
                  else
                  {
                      
                      displ[1] += delta_displ[1];
                      
                      delta_displ[1] = delta_t * vel[1];

                      coor[1] = initial_coor[1] + displ[1];
                 
                  }
                  
                    //pos = i->FastGetSolutionStepValue(VELOCITY_Z_DOF_POS);
                    //if( i->GetDof(VELOCITY_Z, pos).IsFixed() == false ) 
                    if( i->IsNot(DEMFlags::FIXED_VEL_Z))
                  {    
                      vel[2] += aux * force[2];

                      displ[2] +=  delta_displ[2];

                      delta_displ[2] = delta_t * vel[2];

                      coor[2] = initial_coor[2] + displ[2];


                  }
                  else
                  {
                      
                      displ[2] += delta_displ[2];
                      
                      delta_displ[2] = delta_t * vel[2];

                      coor[2] = initial_coor[2] + displ[2];
                 
                  }                    
              }
          }
          
          NodesArrayType& pGNodes            = GetGhostNodes(model_part);
      
    
          //NodesArrayType::iterator it_begin = pNodes.ptr_begin();
          //NodesArrayType::iterator it_end   = pNodes.ptr_end();
          
          OpenMPUtils::CreatePartition(number_of_threads, pGNodes.size(), node_partition);
                            
          #pragma omp parallel for shared(delta_t) 
          for(int k=0; k<number_of_threads; k++)
          {
              NodesArrayType::iterator i_begin=pGNodes.ptr_begin()+node_partition[k];
              NodesArrayType::iterator i_end=pGNodes.ptr_begin()+node_partition[k+1];
             
              for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)      
              {   
                
                  array_1d<double, 3 > & vel             = i->FastGetSolutionStepValue(VELOCITY);
                  array_1d<double, 3 > & displ           = i->FastGetSolutionStepValue(DISPLACEMENT);
                  array_1d<double, 3 > & delta_displ     = i->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                  array_1d<double, 3 > & coor            = i->Coordinates();
                  array_1d<double, 3 > & initial_coor    = i->GetInitialPosition();
                  array_1d<double, 3 > & force           = i->FastGetSolutionStepValue(TOTAL_FORCES);
                  
                  double mass                            = i->FastGetSolutionStepValue(SQRT_OF_MASS);
                  mass                                  *= mass;
                  double aux = delta_t / mass;     

                  //i->FastGetSolutionStepValue(OLD_COORDINATES) = coor; //saving the coordinates in order to optimize some functions (specially de previous step coordinates)  s
                  
                  //unsigned int pos = i->FastGetSolutionStepValue(VELOCITY_X_DOF_POS);
                  //if( i->GetDof(VELOCITY_X, pos).IsFixed() == false ) // equivalently:  i->IsFixed(VELOCITY_X) == false
                  if( i->IsNot(DEMFlags::FIXED_VEL_X))
                  {    
                      vel[0] += aux * force[0];

                      displ[0] +=  delta_displ[0];

                      delta_displ[0] = delta_t * vel[0];

                      coor[0] = initial_coor[0] + displ[0];


                  }
                  else
                  {
                      
                      displ[0] += delta_displ[0];
                      
                      delta_displ[0] = delta_t * vel[0];

                      coor[0] = initial_coor[0] + displ[0];
                 
                  }
                  //pos = i->FastGetSolutionStepValue(VELOCITY_Y_DOF_POS);
                  //if( i->GetDof(VELOCITY_Y, pos).IsFixed() == false ) 
                  if( i->IsNot(DEMFlags::FIXED_VEL_Y))
                  {    
                      vel[1] += aux * force[1];

                      displ[1] +=  delta_displ[1];

                      delta_displ[1] = delta_t * vel[1];

                      coor[1] = initial_coor[1] + displ[1];


                  }
                  else
                  {
                      
                      displ[1] += delta_displ[1];
                      
                      delta_displ[1] = delta_t * vel[1];

                      coor[1] = initial_coor[1] + displ[1];
                 
                  }
                  
                    //pos = i->FastGetSolutionStepValue(VELOCITY_Z_DOF_POS);
                    //if( i->GetDof(VELOCITY_Z, pos).IsFixed() == false ) 
                    if( i->IsNot(DEMFlags::FIXED_VEL_Z))    
                  {    
                      vel[2] += aux * force[2];

                      displ[2] +=  delta_displ[2];

                      delta_displ[2] = delta_t * vel[2];

                      coor[2] = initial_coor[2] + displ[2];


                  }
                  else
                  {
                      
                      displ[2] += delta_displ[2];
                      
                      delta_displ[2] = delta_t * vel[2];

                      coor[2] = initial_coor[2] + displ[2];
                 
                  }                    
              }
          }


          KRATOS_CATCH(" ")
      }


      void CalculateRotationalMotion(ModelPart& model_part)
      {
          KRATOS_TRY   

          ProcessInfo& rCurrentProcessInfo  = model_part.GetProcessInfo();
          NodesArrayType& pNodes            = GetNodes(model_part);
          double delta_t =  rCurrentProcessInfo[DELTA_TIME];
          bool if_virtual_mass_option = (bool) rCurrentProcessInfo[VIRTUAL_MASS_OPTION];
          vector<unsigned int> node_partition;
          bool if_trihedron_option = (bool) rCurrentProcessInfo[TRIHEDRON_OPTION];
          double coeff            = rCurrentProcessInfo[NODAL_MASS_COEFF];

    #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
    #else
        int number_of_threads = 1;
    #endif
    OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

    #pragma omp parallel for
    for(int k=0; k<number_of_threads; k++)
    {
            NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
            NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];
            for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
            {

                //double PMass            = i->FastGetSolutionStepValue(NODAL_MASS);

                double PMomentOfInertia = i->FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);                

                array_1d<double, 3 > & AngularVel             = i->FastGetSolutionStepValue(ANGULAR_VELOCITY);
                array_1d<double, 3 > & RotaMoment             = i->FastGetSolutionStepValue(PARTICLE_MOMENT);
                //array_1d<double, 3 > & delta_rotation_displ   = i->FastGetSolutionStepValue(DELTA_ROTA_DISPLACEMENT);
                array_1d<double, 3 > & Rota_Displace          = i->FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
                
                array_1d<double, 3 > delta_rotation_displ;  
                
                //double & Orientation_real                     = i->FastGetSolutionStepValue(ORIENTATION_REAL); 
                double Orientation_real;
                //array_1d<double, 3 > & Orientation_imag       = i->FastGetSolutionStepValue(ORIENTATION_IMAG);                
                array_1d<double, 3 >  Orientation_imag;
                
                bool If_Fix_Rotation[3] = {false, false, false};

                //unsigned int pos = i->FastGetSolutionStepValue(ANGULAR_VELOCITY_X_DOF_POS);                  
                //If_Fix_Rotation[0] = i->GetDof(ANGULAR_VELOCITY_X, pos).IsFixed();
                //pos = i->FastGetSolutionStepValue(ANGULAR_VELOCITY_Y_DOF_POS);
                //If_Fix_Rotation[1] = i->GetDof(ANGULAR_VELOCITY_Y, pos).IsFixed();
                //pos = i->FastGetSolutionStepValue(ANGULAR_VELOCITY_Z_DOF_POS);
                //If_Fix_Rotation[2] = i->GetDof(ANGULAR_VELOCITY_Z, pos).IsFixed(); 
                
                If_Fix_Rotation[0] = i->Is(DEMFlags::FIXED_ANG_VEL_X);
                If_Fix_Rotation[1] = i->Is(DEMFlags::FIXED_ANG_VEL_Y);
                If_Fix_Rotation[2] = i->Is(DEMFlags::FIXED_ANG_VEL_Z);
                
                for(std::size_t iterator = 0 ; iterator < 3; iterator++)
                {
                    if(If_Fix_Rotation[iterator] == false)
                    {
                         double RotaAcc = 0.0;
                         
                         RotaAcc = (RotaMoment[iterator]) / (PMomentOfInertia);

                         if(if_virtual_mass_option)
                         {
                                  RotaAcc = RotaAcc * ( 1 - coeff );
                         }
                       
                         delta_rotation_displ[iterator] = AngularVel[iterator] * delta_t;

                         AngularVel[iterator] += RotaAcc * delta_t;
                         
                         Rota_Displace[iterator] +=  delta_rotation_displ[iterator];                         
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
          
                if(if_trihedron_option)
                {
                    double theta[3] = {0.0};
                    
                    theta[0] = Rota_Displace[0] * 0.5;
                    theta[1] = Rota_Displace[1] * 0.5;
                    theta[2] = Rota_Displace[2] * 0.5;                

                    double thetaMag = sqrt(theta[0] * theta[0] + theta[1] * theta[1] + theta[2] * theta[2]);
                    if(thetaMag * thetaMag * thetaMag * thetaMag / 24.0 < DBL_EPSILON)  //Taylor: low angle
                    {
                        Orientation_real = 1 + thetaMag * thetaMag / 2;
                        Orientation_imag[0] = theta[0] * ( 1 - thetaMag * thetaMag / 6 );
                        Orientation_imag[1] = theta[1] * ( 1 - thetaMag * thetaMag / 6 );
                        Orientation_imag[2] = theta[2] * ( 1 - thetaMag * thetaMag / 6 );
                    }
                    
                    else
                    {
                        Orientation_real = cos (thetaMag);
                        Orientation_imag[0] = (theta[0] / thetaMag) * sin (thetaMag);
                        Orientation_imag[1] = (theta[1] / thetaMag) * sin (thetaMag);
                        Orientation_imag[2] = (theta[2] / thetaMag) * sin (thetaMag);
                    }
                    
                    array_1d<double,3>& EulerAngles       = i->FastGetSolutionStepValue(EULER_ANGLES);    
                
                    double test = Orientation_imag[0] * Orientation_imag[1] + Orientation_imag[2] * Orientation_real;
                    
                    if ( test > 0.49999999 )               // singularity at north pole
                    {
                        EulerAngles[0] = 2 * atan2 ( Orientation_imag[0] , Orientation_real );
                        EulerAngles[1] = pi;
                        EulerAngles[2] = 0.0;
                    }
                    
                    else if (test < -0.49999999 )          // singularity at south pole
                    {                 
                        EulerAngles[0] = -2 * atan2 ( Orientation_imag[0] , Orientation_real );
                        EulerAngles[1] = -pi;
                        EulerAngles[2] = 0.0;
                    }

                    else                
                    {
                        EulerAngles[0] = atan2( 2 * Orientation_real * Orientation_imag[0] + 2 * Orientation_imag[1] * Orientation_imag[2] , 1 - 2 * Orientation_imag[0] * Orientation_imag[0] - 2 * Orientation_imag[1] * Orientation_imag[1] );
                        EulerAngles[1] = asin ( 2 * Orientation_real * Orientation_imag[1] - 2 * Orientation_imag[2] * Orientation_imag[0] );
                        EulerAngles[2] = -atan2( 2 * Orientation_real * Orientation_imag[2] + 2 * Orientation_imag[0] * Orientation_imag[1] , 1 - 2 * Orientation_imag[1] * Orientation_imag[1] - 2 * Orientation_imag[2] * Orientation_imag[2] );
                    }
                    
                  }// Trihedron Option
                
                }//for Node
                
            }//for OMP
        
        KRATOS_CATCH(" ")
        
     }//rotational_motion
    

      /// Turn back information as a string.
      virtual std::string Info() const
      {
    std::stringstream buffer;
        buffer << "ForwardEulerScheme" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "ForwardEulerScheme";}

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
     ForwardEulerScheme& operator=(ForwardEulerScheme const& rOther)
     {
       return *this;
     }

      /// Copy constructor.
      ForwardEulerScheme(ForwardEulerScheme const& rOther)
      {
    *this = rOther;
      }

        
      ///@}    
        
    }; // Class ForwardEulerScheme 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
                    ForwardEulerScheme& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
                    const ForwardEulerScheme& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_FORWARD_EULER_SCHEME_H_INCLUDED  defined 
