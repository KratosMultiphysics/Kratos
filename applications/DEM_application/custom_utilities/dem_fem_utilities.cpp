// /* External includes */

// System includes

// Project includes
#include "utilities/timer.h"
#include "includes/variables.h"
#include "DEM_application.h"

// System includes
#include <limits>
#include <iostream>
#include <iomanip>

// External includes
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "GeometryFunctions.h"

#include "dem_fem_utilities.h"

namespace Kratos
{

    //DEMFEMUtilities::DEMFEMUtilities() {}

    /// Destructor

    DEMFEMUtilities::~DEMFEMUtilities(){}
      
    void DEMFEMUtilities::CrossProduct( const array_1d<double,3>& u, const array_1d<double,3>& v, array_1d<double,3>& ReturnVector) {

        ReturnVector[0] = u[1]*v[2] - u[2]*v[1];
        ReturnVector[1] = u[2]*v[0] - u[0]*v[2];
        ReturnVector[2] = u[0]*v[1] - u[1]*v[0];
    }
      
    void DEMFEMUtilities::RotateRightHandedBasisAroundAxis(const array_1d<double, 3 >& e1,  const array_1d<double, 3 >& e2,  const array_1d<double, 3 >& axis,
                                                           const double ang, array_1d<double, 3 >& new_axes1, array_1d<double, 3 >& new_axes2,
                                                           array_1d<double, 3 >& new_axes3) {
        double cang = cos(ang);
        double sang = sin(ang);
            
        new_axes1[0] = axis[0] * (axis[0] * e1[0] + axis[1] * e1[1] + axis[2] * e1[2]) * (1 - cang) + e1[0] * cang + (- axis[2] * e1[1] + axis[1] * e1[2]) * sang;
        new_axes1[1] = axis[1] * (axis[0] * e1[0] + axis[1] * e1[1] + axis[2] * e1[2]) * (1 - cang) + e1[1] * cang +   (axis[2] * e1[0] - axis[0] * e1[2]) * sang;
        new_axes1[2] = axis[2] * (axis[0] * e1[0] + axis[1] * e1[1] + axis[2] * e1[2]) * (1 - cang) + e1[2] * cang + (- axis[1] * e1[0] + axis[0] * e1[1]) * sang;
          
        new_axes2[0] = axis[0] * (axis[0] * e2[0] + axis[1] * e2[1] + axis[2] * e2[2]) * (1 - cang) + e2[0] * cang + (- axis[2] * e2[1] + axis[1] * e2[2]) * sang;
        new_axes2[1] = axis[1] * (axis[0] * e2[0] + axis[1] * e2[1] + axis[2] * e2[2]) * (1 - cang) + e2[1] * cang +   (axis[2] * e2[0] - axis[0] * e2[2]) * sang;
        new_axes2[2] = axis[2] * (axis[0] * e2[0] + axis[1] * e2[1] + axis[2] * e2[2]) * (1 - cang) + e2[2] * cang + (- axis[1] * e2[0] + axis[0] * e2[1]) * sang;  
            
        CrossProduct(new_axes1, new_axes2, new_axes3);
    }
      
    void DEMFEMUtilities::MoveAllMeshes(ModelPart& r_model_part, double time, double dt) {

        if (r_model_part.NumberOfMeshes() > 1) {

            for (unsigned int mesh_number = 1; mesh_number < r_model_part.NumberOfMeshes(); mesh_number++) {

                NodesArrayType& pNodes = r_model_part.GetMesh(mesh_number).Nodes();

                const double velocity_start_time = r_model_part.GetMesh(mesh_number)[VELOCITY_START_TIME];
                const double velocity_stop_time = r_model_part.GetMesh(mesh_number)[VELOCITY_STOP_TIME];
                const double angular_velocity_start_time = r_model_part.GetMesh(mesh_number)[ANGULAR_VELOCITY_START_TIME];
                const double angular_velocity_stop_time = r_model_part.GetMesh(mesh_number)[ANGULAR_VELOCITY_STOP_TIME];

                array_1d<double, 3>& previous_displ = r_model_part.GetMesh(mesh_number)[DISPLACEMENT];
                const array_1d<double, 3>& linear_velocity = r_model_part.GetMesh(mesh_number)[VELOCITY];
                const double linear_period = r_model_part.GetMesh(mesh_number)[VELOCITY_PERIOD];
                const array_1d<double, 3>& angular_velocity = r_model_part.GetMesh(mesh_number)[ANGULAR_VELOCITY];
                const double angular_period = r_model_part.GetMesh(mesh_number)[ANGULAR_VELOCITY_PERIOD];
                const array_1d<double, 3>& initial_center = r_model_part.GetMesh(mesh_number)[ROTATION_CENTER];
                const bool fixed_mesh = r_model_part.GetMesh(mesh_number)[FIXED_MESH_OPTION];
                array_1d<double, 3> center_position;
                array_1d<double, 3> linear_velocity_changed;
                array_1d<double, 3> angular_velocity_changed;
                array_1d<double, 3> angle = ZeroVector(3);
                double sign_angle = 1.0;
                array_1d<double, 3> final_angle = ZeroVector(3);

                if (time < velocity_start_time || time > velocity_stop_time) {
                    center_position[0] = initial_center[0] + previous_displ[0];
                    center_position[1] = initial_center[1] + previous_displ[1];
                    center_position[2] = initial_center[2] + previous_displ[2];
                    linear_velocity_changed = ZeroVector(3);
                } else {
                    if (linear_period > 0.0) {
                        double linear_omega = 2.0 * KRATOS_M_PI / linear_period;
                        double inv_linear_omega = 1.0 / linear_omega;
                        noalias(center_position) = initial_center + linear_velocity * sin(linear_omega * time) * inv_linear_omega;
                        noalias(linear_velocity_changed) = linear_velocity * cos(linear_omega * time);
                        noalias(previous_displ) = center_position - initial_center;
                    } else {
                        center_position[0] = initial_center[0] + previous_displ[0] + dt * linear_velocity[0];
                        center_position[1] = initial_center[1] + previous_displ[1] + dt * linear_velocity[1];
                        center_position[2] = initial_center[2] + previous_displ[2] + dt * linear_velocity[2];
                        previous_displ[0] += dt * linear_velocity[0];
                        previous_displ[1] += dt * linear_velocity[1];
                        previous_displ[2] += dt * linear_velocity[2];
                        linear_velocity_changed = linear_velocity;
                    }
                }

                if (time < angular_velocity_start_time) angular_velocity_changed = ZeroVector(3);

                else if (((time - angular_velocity_start_time) > 0.0) && ((time - angular_velocity_stop_time) < 0.0)) {

                if (angular_period > 0.0) {
                    double angular_omega = 2.0 * KRATOS_M_PI / angular_period;
                    double inv_angular_omega = 1.0 / angular_omega;
                    noalias(angle) = angular_velocity * sin(angular_omega * (time - angular_velocity_start_time)) * inv_angular_omega;
                    sign_angle = sin(angular_omega * (time - angular_velocity_start_time)) / fabs(sin(angular_omega * (time - angular_velocity_start_time)));
                    noalias(angular_velocity_changed) = angular_velocity * cos(angular_omega * (time - angular_velocity_start_time));
                    final_angle = angle;
                } else {
                    noalias(angle) = angular_velocity * (time - angular_velocity_start_time);
                    angular_velocity_changed = angular_velocity;
                    }
                }
                  else { //if ((time - angular_velocity_stop_time) > 0.0) {
                    angular_velocity_changed = ZeroVector(3);

                    if (angular_period > 0.0) {
                        angle = final_angle;
                    } else {
                        noalias(angle) = angular_velocity * (angular_velocity_stop_time - angular_velocity_start_time);
                        }
                    }

                double mod_angular_velocity = MathUtils<double>::Norm3(angular_velocity);

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

                if (mod_angular_velocity > 0.0) {

                    double ang = sign_angle * MathUtils<double>::Norm3(angle);
                    array_1d<double, 3 > rotation_axis;
                    noalias(rotation_axis) = angular_velocity / MathUtils<double>::Norm3(angular_velocity);
                    array_1d<double, 3 > e1;
                    e1[0] = 1.0;
                    e1[1] = 0.0;
                    e1[2] = 0.0;

                    array_1d<double, 3 > e2;
                    e2[0] = 0.0;
                    e2[1] = 1.0;
                    e2[2] = 0.0;

                    RotateRightHandedBasisAroundAxis(e1, e2, rotation_axis, ang, new_axes1, new_axes2, new_axes3);
                }

                if (mod_angular_velocity > 0.0 || MathUtils<double>::Norm3(linear_velocity) > 0.0) {

                    vector<unsigned int> node_partition;

                    #ifdef _OPENMP
                    int number_of_threads = omp_get_max_threads();
                    #else
                    int number_of_threads = 1;
                    #endif
                    OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

                    #pragma omp parallel for
                    for (int k = 0; k < number_of_threads; k++) {

                        array_1d<double, 3 > local_coordinates;
                        local_coordinates[0] = 0.0;
                        local_coordinates[1] = 0.0;
                        local_coordinates[2] = 0.0;
                        array_1d<double, 3 > relative_position;
                        relative_position[0] = 0.0;
                        relative_position[1] = 0.0;
                        relative_position[2] = 0.0;

                        NodesArrayType::iterator i_begin = pNodes.ptr_begin() + node_partition[k];
                        NodesArrayType::iterator i_end = pNodes.ptr_begin() + node_partition[k + 1];

                        for (ModelPart::NodeIterator node = i_begin; node != i_end; ++node) {

                            noalias(local_coordinates) = node->GetInitialPosition().Coordinates() - initial_center;
                            noalias(relative_position) = new_axes1 * local_coordinates[0] + new_axes2 * local_coordinates[1] + new_axes3 * local_coordinates[2];

                            array_1d<double, 3 > displacement;
                            array_1d<double, 3 > old_coordinates;
                            noalias(old_coordinates) = node->Coordinates();
                                                                
                            array_1d<double, 3 > velocity_due_to_rotation;
                            CrossProduct(angular_velocity_changed, relative_position, velocity_due_to_rotation);
                               
                            array_1d<double, 3 >& velocity = node->FastGetSolutionStepValue(VELOCITY);
                            noalias(velocity) = linear_velocity_changed + velocity_due_to_rotation;  
                                                               
                            if (!fixed_mesh) {
                                // NEW POSITION
                                noalias(node->Coordinates()) = center_position + relative_position;

                                // DISPLACEMENT
                                noalias( node->FastGetSolutionStepValue(DISPLACEMENT) ) = node->Coordinates() - node->GetInitialPosition().Coordinates();
                                noalias( node->FastGetSolutionStepValue(DELTA_DISPLACEMENT) ) = node->Coordinates() - old_coordinates;
                            } else {
                                (node->FastGetSolutionStepValue(DISPLACEMENT)).clear(); //Set values to zero
                                noalias( node->FastGetSolutionStepValue(DELTA_DISPLACEMENT) ) = velocity * dt; //But still there must be some delta_displacement (or motion won't be detected by the spheres!)
                            }
                        }
                    }
                }
            } //for (unsigned int mesh_number = 1; mesh_number < r_model_part.NumberOfMeshes(); mesh_number++)
        } //if ( r_model_part.NumberOfMeshes() > 1 )
    }        


    void DEMFEMUtilities::CreateRigidFacesFromAllElements(ModelPart& r_model_part, PropertiesType::Pointer pProps) {
            
        ElementsArrayType&  all_elements = r_model_part.Elements();            
        
        for (unsigned int i = 0; i < all_elements.size(); i++) {
            
            ConditionType::Pointer pCondition;
            ElementsArrayType::iterator pElement = all_elements.ptr_begin()+i;
            pCondition = ConditionType::Pointer(new RigidFace3D( pElement->Id(), pElement->pGetGeometry(), pProps));        
            r_model_part.Conditions().push_back(pCondition); 
        }
    }
        
        
        
        /// Turn back information as a string.
    std::string DEMFEMUtilities::Info() const {
            return "";
    }

        /// Print information about this object.
    void DEMFEMUtilities::PrintInfo(std::ostream& rOStream) const {}

        /// Print object's data.
    void DEMFEMUtilities::PrintData(std::ostream& rOStream) const {}


   

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


} // namespace Kratos
