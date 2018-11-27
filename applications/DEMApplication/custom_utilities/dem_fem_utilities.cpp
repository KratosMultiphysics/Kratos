//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Salva Latorre$
//   Date:                $Date: 2015-10-26 09:56:42 $
//   Revision:            $Revision: 1.5 $
//

// Project includes

// System includes
#include <limits>
#include <iostream>
#include <iomanip>

// External includes
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "dem_fem_utilities.h"
#include "custom_conditions/RigidFace.h"

namespace Kratos
{
    // Constructor
    //DEMFEMUtilities::DEMFEMUtilities() {}

    // Destructor

    DEMFEMUtilities::~DEMFEMUtilities() {}

    void DEMFEMUtilities::MoveAllMeshes(ModelPart& r_model_part, double time, double dt) {

        if (r_model_part.NumberOfSubModelParts()) {
            for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = r_model_part.SubModelPartsBegin();
                                                                 sub_model_part != r_model_part.SubModelPartsEnd(); ++sub_model_part) {

                ModelPart& submp = *sub_model_part;

                const bool rigid_body_motion = submp[RIGID_BODY_MOTION];
                if (!rigid_body_motion) continue;

                NodesArrayType& rNodes = submp.Nodes();
                array_1d<double, 3>& previous_displ = submp[DISPLACEMENT];
                const array_1d<double, 3>& linear_velocity = submp[LINEAR_VELOCITY];
                const double velocity_start_time = submp[VELOCITY_START_TIME];
                const double velocity_stop_time = submp[VELOCITY_STOP_TIME];
                const double linear_period = submp[VELOCITY_PERIOD];
                const bool fixed_mesh = submp[FIXED_MESH_OPTION];

                const array_1d<double, 3>& angular_velocity = submp[ANGULAR_VELOCITY];
                const double angular_velocity_start_time = submp[ANGULAR_VELOCITY_START_TIME];
                const double angular_velocity_stop_time = submp[ANGULAR_VELOCITY_STOP_TIME];
                const double angular_period = submp[ANGULAR_VELOCITY_PERIOD];
                const array_1d<double, 3>& initial_center = submp[ROTATION_CENTER];

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
                                                                          fixed_mesh, dt, rNodes);
            } //for ModelPart::SubModelPartsContainerType::iterator
        } //if r_model_part.NumberOfMeshes() > 1
    }

    // void DEMFEMUtilities::MoveAllMeshesUsingATable(ModelPart& r_model_part, double time, double dt) {

    //     if (r_model_part.NumberOfSubModelParts()) {

    //         for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = r_model_part.SubModelPartsBegin();
    //                                                              sub_model_part != r_model_part.SubModelPartsEnd(); ++sub_model_part) {
    //             ModelPart& submp = *sub_model_part;
    //             NodesArrayType& rNodes = submp.Nodes();
    //             const int table_number = submp[TABLE_NUMBER];
    //             if (!table_number) continue;
    //             const int velocity_component = submp[TABLE_VELOCITY_COMPONENT];

    //             #pragma omp parallel for
    //             for (int k = 0; k < (int)rNodes.size(); k++) {

    //                 array_1d<double, 3> old_coordinates = ZeroVector(3);
    //                 ModelPart::NodeIterator node = rNodes.begin() + k;

    //                 old_coordinates[0] = node->Coordinates()[0];
    //                 old_coordinates[1] = node->Coordinates()[1];
    //                 old_coordinates[2] = node->Coordinates()[2];

    //                 array_1d<double, 3> velocity = ZeroVector(3);
    //                 velocity[0] = 0.0;
    //                 velocity[1] = 0.0;
    //                 velocity[2] = 0.0;
    //                 velocity[velocity_component] = r_model_part.GetTable(table_number).GetValue(time);
    //                 noalias(node->FastGetSolutionStepValue(VELOCITY)) = velocity;

    //                 node->Coordinates()[0] = old_coordinates[0] + velocity[0] * dt;
    //                 node->Coordinates()[1] = old_coordinates[1] + velocity[1] * dt;
    //                 node->Coordinates()[2] = old_coordinates[2] + velocity[2] * dt;

    //                 array_1d<double, 3> displacement = ZeroVector(3);
    //                 displacement[0] = node->Coordinates()[0] - node->GetInitialPosition().Coordinates()[0];
    //                 displacement[1] = node->Coordinates()[1] - node->GetInitialPosition().Coordinates()[1];
    //                 displacement[2] = node->Coordinates()[2] - node->GetInitialPosition().Coordinates()[2];
    //                 noalias(node->FastGetSolutionStepValue(DISPLACEMENT)) = displacement;
    //             } // for loop
    //         } // for ModelPart::SubModelPartsContainerType::iterator
    //     } // if r_model_part.NumberOfMeshes()
    // } // function

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
