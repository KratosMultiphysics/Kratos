#ifndef DEM_FEM_UTILITES_H
#define DEM_FEM_UTILITES_H

// /* External includes */

// System includes

// Project includes
#include "utilities/timer.h"
#include "includes/variables.h"
#include "DEM_application.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>
#include <iostream>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "GeometryFunctions.h"

namespace Kratos
{

class DEMFEMUtilities
    {
     public:

     typedef ModelPart::ElementsContainerType                          ElementsArrayType;
     typedef ModelPart::NodesContainerType                             NodesArrayType;
     typedef WeakPointerVector<Element>                                ParticleWeakVectorType;
     typedef WeakPointerVector<Element >::iterator                     ParticleWeakIteratorType;

     KRATOS_CLASS_POINTER_DEFINITION(DEMFEMUtilities);

      /// Default constructor.

      DEMFEMUtilities()
      {
          
      }

      /// Destructor.

      virtual ~DEMFEMUtilities(){}
      
      void CrossProduct( const array_1d<double,3>& u, const array_1d<double,3>& v, array_1d<double,3>& ReturnVector)
      {
            ReturnVector[0] = u[1]*v[2] - u[2]*v[1];
            ReturnVector[1] = v[0]*u[2] - u[0]*v[2];
            ReturnVector[2] = u[0]*v[1] - u[1]*v[0];
      }
      
      void RotateRightHandedBasisAroundAxis(const array_1d<double, 3 >& e1,  const array_1d<double, 3 >& e2,  const array_1d<double, 3 >& axis, const double ang, 
                                                  array_1d<double, 3 >& new_axes1, array_1d<double, 3 >& new_axes2, array_1d<double, 3 >& new_axes3){
          
            double cang         = cos(ang);
            double sang         = sin(ang);
            
            new_axes1[0] = axis[0] * (axis[0] * e1[0] + axis[1] * e1[1] + axis[2] * e1[2]) * (1 - cang) + e1[0] * cang + (- axis[2] * e1[1] + axis[1] * e1[2]) * sang;
            new_axes1[1] = axis[1] * (axis[0] * e1[0] + axis[1] * e1[1] + axis[2] * e1[2]) * (1 - cang) + e1[1] * cang +   (axis[2] * e1[0] - axis[0] * e1[2]) * sang;
            new_axes1[2] = axis[2] * (axis[0] * e1[0] + axis[1] * e1[1] + axis[2] * e1[2]) * (1 - cang) + e1[2] * cang + (- axis[1] * e1[0] + axis[0] * e1[1]) * sang;
             
            new_axes2[0] = axis[0] * (axis[0] * e2[0] + axis[1] * e2[1] + axis[2] * e2[2]) * (1 - cang) + e2[0] * cang + (- axis[2] * e2[1] + axis[1] * e2[2]) * sang;
            new_axes2[1] = axis[1] * (axis[0] * e2[0] + axis[1] * e2[1] + axis[2] * e2[2]) * (1 - cang) + e2[1] * cang +   (axis[2] * e2[0] - axis[0] * e2[2]) * sang;
            new_axes2[2] = axis[2] * (axis[0] * e2[0] + axis[1] * e2[1] + axis[2] * e2[2]) * (1 - cang) + e2[2] * cang + (- axis[1] * e2[0] + axis[0] * e2[1]) * sang  ;  
            
            CrossProduct(new_axes1, new_axes2, new_axes3);
                      
      }
      

      void MoveAllMeshes(ModelPart& r_model_part, double time)
      {
          
          if ( r_model_part.NumberOfMeshes() > 1 ) {

            for (unsigned int mesh_number = 1; mesh_number < r_model_part.NumberOfMeshes(); mesh_number++){
                
                NodesArrayType& pNodes         = r_model_part.GetMesh(mesh_number).Nodes();

                array_1d<double, 3 >& linear_velocity  = r_model_part.GetMesh(mesh_number)[VELOCITY];
                double                linear_period    = r_model_part.GetMesh(mesh_number)[VELOCITY_PERIOD];
                array_1d<double, 3 >& angular_velocity = r_model_part.GetMesh(mesh_number)[ANGULAR_VELOCITY];
                double                angular_period   = r_model_part.GetMesh(mesh_number)[ANGULAR_VELOCITY_PERIOD];
                array_1d<double, 3 >& initial_center   = r_model_part.GetMesh(mesh_number)[ROTATION_CENTER];
                bool                  fixed_mesh       = r_model_part.GetMesh(mesh_number)[FIXED_MESH_OPTION];

                array_1d<double, 3 > center_position;
                array_1d<double, 3 > linear_velocity_changed;
                array_1d<double, 3 > angular_velocity_changed;
                array_1d<double, 3 > angle;

                if (linear_period > 0.0){
                    double linear_omega = 2.0 * KRATOS_M_PI / linear_period;
                    double inv_linear_omega = 1.0/linear_omega;
                    center_position = initial_center + linear_velocity * sin(linear_omega * time)* inv_linear_omega;
                    linear_velocity_changed = linear_velocity * cos(linear_omega * time);
                }
                else {
                    center_position = initial_center + time * linear_velocity;
                    linear_velocity_changed = linear_velocity;
                }

                if (angular_period > 0.0){
                    double angular_omega = 2.0 * KRATOS_M_PI / angular_period;
                    double inv_angular_omega = 1.0/angular_omega;
                    angle = angular_velocity * sin(angular_omega * time) * inv_angular_omega;
                    angular_velocity_changed = angular_velocity * cos(angular_omega * time);
                }
                else {
                    angle = angular_velocity * time;
                    angular_velocity_changed = angular_velocity;
                }


                double mod_angular_velocity = MathUtils<double>::Norm3(angular_velocity);
                array_1d<double, 3 > relative_position; relative_position[0] = 0.0; relative_position[1] = 0.0; relative_position[2] = 0.0;
                
                array_1d<double, 3 > new_axes1;
                new_axes1[0] = 1.0;
                new_axes1[1] = 0.0;
                new_axes1[2] = 0.0;

                array_1d<double, 3 > new_axes2;
                new_axes2[0] = 0.0;
                new_axes2[1] = 1.0;
                new_axes2[2] = 0.0;

                array_1d<double, 3 > new_axes3;
                new_axes3[0] = 0.0;
                new_axes3[1] = 0.0;
                new_axes3[2] = 1.0;                

                if (mod_angular_velocity > 0.0){
                    double ang = MathUtils<double>::Norm3(angle);
                    array_1d<double, 3 > rotation_axis = angular_velocity/MathUtils<double>::Norm3(angular_velocity);
                    array_1d<double, 3 > e1;
                    e1[0] = 1.0;
                    e1[1] = 0.0;
                    e1[2] = 0.0;

                    array_1d<double, 3 > e2;
                    e2[0] = 0.0;
                    e2[1] = 1.0;
                    e2[2] = 0.0;

                    RotateRightHandedBasisAroundAxis(e1, e2, rotation_axis, ang, new_axes1, new_axes2, new_axes3)  ;
                }                

                if (mod_angular_velocity>0.0 || MathUtils<double>::Norm3(linear_velocity)>0.0) {
                    
                    
                    vector<unsigned int> node_partition;  
                    
                    #ifdef _OPENMP
                    int number_of_threads = omp_get_max_threads();
                    #else
                    int number_of_threads = 1;
                    #endif
                    OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

                    #pragma omp parallel for
                    for(int k=0; k<number_of_threads; k++) {
                    
                        array_1d<double, 3 >  local_coordinates; local_coordinates[0] = 0.0; local_coordinates[1] = 0.0;  local_coordinates[2] = 0.0; 
                        array_1d<double, 3 >  relative_position; relative_position[0] = 0.0; relative_position[1] = 0.0;  relative_position[2] = 0.0; 
                        
                        NodesArrayType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
                        NodesArrayType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

                        for(ModelPart::NodeIterator node=i_begin; node!= i_end; ++node)  {    
                            
                            local_coordinates = node->GetInitialPosition().Coordinates() - initial_center;
                            relative_position = new_axes1 * local_coordinates[0] + new_axes2 * local_coordinates[1] + new_axes3 * local_coordinates[2];                            
                            
                            array_1d<double, 3 > displacement;
                            if(!fixed_mesh){
                                // NEW POSITION
                                node->Coordinates() = center_position + relative_position;

                                // DISPLACEMENT
                                displacement = node->Coordinates() - node->GetInitialPosition().Coordinates();                            
                            }
                            else{
                                displacement[0]=0.0; displacement[1]=0.0; displacement[2]=0.0;
                            }

                            array_1d<double, 3 > velocity_due_to_rotation;  
                            CrossProduct(angular_velocity_changed , relative_position, velocity_due_to_rotation);
                            
                            // NEW VELOCITY                    
                            array_1d<double, 3 > vel;
                            vel[0] = linear_velocity_changed[0] + velocity_due_to_rotation[0];
                            vel[1] = linear_velocity_changed[1] + velocity_due_to_rotation[1];
                            vel[2] = linear_velocity_changed[2] + velocity_due_to_rotation[2];

                            //update VELOCITY
                            node->FastGetSolutionStepValue(VELOCITY) = vel;

                            //update DISPLACEMENT                            
                             node->FastGetSolutionStepValue(DISPLACEMENT) = displacement;
                            
                                                                                      
                        }
                        
                    }
                        
                }
            } //for (unsigned int mesh_number = 1; mesh_number < r_model_part.NumberOfMeshes(); mesh_number++)
            
          } //if ( r_model_part.NumberOfMeshes() > 1 )

          
          
          
          

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

        /// Turn back information as a stemplate<class T, std::size_t dim> tring.

        virtual std::string Info() const
        {
            return "";
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
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
        ///@name Protected static Member rVariables
        ///@{


        ///@}
        ///@name Protected member rVariables
        ///@{ template<class T, std::size_t dim>


        ///@}
        ///@name Protected Operators
        ///@{


        ///@}
        ///@name Protected Operations
        ///@{


        ///@}
        ///@name Protected  Access
        ///@{
        vector<unsigned int>                mElementPartition;

        ///@}
        ///@name Protected Inquiry
        ///@{


        ///@}
        ///@name Protected LifeCycle
        ///@{


        ///@}

    private:


        ///@name Static Member rVariables
        ///@{


        ///@}
        ///@name Member rVariables
        ///@{

        array_1d<double, 3> mInitialCenterOfMassAndMass;
        double mInitialMass;


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
        DEMFEMUtilities & operator=(DEMFEMUtilities const& rOther);


        ///@}

    }; // Class DEMFEMUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{




/// output stream function
// 	template<std::size_t TDim>
// 	inline std::ostream& operator << (std::ostream& rOStream)
// 	{
// 		rThis.PrintInfo(rOStream);
// 		rOStream << std::endl;
// 		rThis.PrintData(rOStream);
//
// 		return rOStream;
// 	}
///@}


} // namespace Kratos.

#endif // DEM_FEM_UTILITES_H
