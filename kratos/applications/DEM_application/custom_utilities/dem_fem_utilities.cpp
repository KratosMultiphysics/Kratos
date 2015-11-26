//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Salva Latorre$
//   Date:                $Date: 2015-10-26 09:56:42 $
//   Revision:            $Revision: 1.5 $
//

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
    //using namespace GeometryFunctions;
    
    // Constructor
    //DEMFEMUtilities::DEMFEMUtilities() {}

    // Destructor

    DEMFEMUtilities::~DEMFEMUtilities() {}
      
    void DEMFEMUtilities::MoveAllMeshes(ModelPart& r_model_part, double time, double dt) {

        if (r_model_part.NumberOfMeshes() > 1) {

            for (unsigned int mesh_number = 1; mesh_number < r_model_part.NumberOfMeshes(); mesh_number++) {

                NodesArrayType& pNodes = r_model_part.GetMesh(mesh_number).Nodes();
                
                array_1d<double, 3>& previous_displ = r_model_part.GetMesh(mesh_number)[DISPLACEMENT];
                
                const array_1d<double, 3>& linear_velocity = r_model_part.GetMesh(mesh_number)[VELOCITY];
                const double velocity_start_time = r_model_part.GetMesh(mesh_number)[VELOCITY_START_TIME];
                const double velocity_stop_time = r_model_part.GetMesh(mesh_number)[VELOCITY_STOP_TIME];
                const double linear_period = r_model_part.GetMesh(mesh_number)[VELOCITY_PERIOD];
                const bool fixed_mesh = r_model_part.GetMesh(mesh_number)[FIXED_MESH_OPTION];
                
                const array_1d<double, 3>& angular_velocity = r_model_part.GetMesh(mesh_number)[ANGULAR_VELOCITY];
                const double angular_velocity_start_time = r_model_part.GetMesh(mesh_number)[ANGULAR_VELOCITY_START_TIME];
                const double angular_velocity_stop_time = r_model_part.GetMesh(mesh_number)[ANGULAR_VELOCITY_STOP_TIME];
                const double angular_period = r_model_part.GetMesh(mesh_number)[ANGULAR_VELOCITY_PERIOD];
                const array_1d<double, 3>& initial_center = r_model_part.GetMesh(mesh_number)[ROTATION_CENTER];
                                
                array_1d<double, 3> center_position;
                array_1d<double, 3> linear_velocity_changed;
                array_1d<double, 3> angular_velocity_changed;
                double mod_angular_velocity = MathUtils<double>::Norm3(angular_velocity);
                array_1d<double, 3> new_axes1;
                array_1d<double, 3> new_axes2;
                array_1d<double, 3> new_axes3;
                
                GeometryFunctions::TranslateGridOfNodes(time, velocity_start_time, velocity_stop_time, center_position, initial_center, previous_displ,
                                                        linear_velocity_changed, linear_period, dt, linear_velocity);

                GeometryFunctions::RotateGridOfNodes(time, angular_velocity_start_time, angular_velocity_stop_time, angular_velocity_changed,
                                                     angular_period, mod_angular_velocity, angular_velocity, new_axes1, new_axes2, new_axes3);

                GeometryFunctions::UpdateKinematicVariablesOfAGridOfNodes(mod_angular_velocity, linear_velocity, initial_center, new_axes1,
                                                                          new_axes2, new_axes3, angular_velocity_changed, linear_velocity_changed, center_position,
                                                                          fixed_mesh, dt, pNodes);
                
            } //for (unsigned int mesh_number = 1; mesh_number < r_model_part.NumberOfMeshes(); mesh_number++)
        } //if ( r_model_part.NumberOfMeshes() > 1 )
    }
    
    void DEMFEMUtilities::CreateRigidFacesFromAllElements(ModelPart& r_model_part, PropertiesType::Pointer pProps) {
            
        ElementsArrayType& all_elements = r_model_part.Elements();            
        
        for (unsigned int i = 0; i < all_elements.size(); i++) {
            
            ConditionType::Pointer pCondition;
            ElementsArrayType::iterator pElement = all_elements.ptr_begin() + i;
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
