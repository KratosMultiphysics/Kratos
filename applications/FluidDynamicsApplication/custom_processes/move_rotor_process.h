//
//  Main authors:    Miguel Ángel Celigueta
//
//

#ifndef KRATOS_MOVE_ROTOR_PROCESS_H
#define KRATOS_MOVE_ROTOR_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"

// Application includes
namespace Kratos
{

class MoveRotorProcess : public Process
{
public:

    /// Pointer definition of MoveRotorProcess
    KRATOS_CLASS_POINTER_DEFINITION(MoveRotorProcess);

    typedef Node<3>                     NodeType;
    typedef Geometry<NodeType>      GeometryType;

    /// Constructor.
    MoveRotorProcess(ModelPart& rModelPart,
                    const double angular_velocity_with_respect_to_stator_center,
                    const double coordinates_of_stator_center_x,
                    const double coordinates_of_stator_center_y,
                    const double initial_coordinates_of_rotor_center_x,
                    const double initial_coordinates_of_rotor_center_y,
                    const unsigned int number_of_rotor_lobules): Process(), mrModelPart(rModelPart)
    {
        mW1[0] = 0.0;
        mW1[1] = 0.0;
        mW1[2] = angular_velocity_with_respect_to_stator_center;
        mCoordinatesOfStatorCenter[0] = coordinates_of_stator_center_x;
        mCoordinatesOfStatorCenter[1] = coordinates_of_stator_center_y;
        mCoordinatesOfStatorCenter[2] = 0.0;
        mInitialCoordinatesOfRotorCenter[0] = initial_coordinates_of_rotor_center_x;
        mInitialCoordinatesOfRotorCenter[1] = initial_coordinates_of_rotor_center_y;
        mInitialCoordinatesOfRotorCenter[2] = 0.0;
        const array_1d<double,3> centers_distance = mInitialCoordinatesOfRotorCenter - mCoordinatesOfStatorCenter;
        mEccentricity = sqrt(centers_distance[0]*centers_distance[0] + centers_distance[1]*centers_distance[1]);
        mW2[0] = 0.0;
        mW2[1] = 0.0;
        mW2[2] = -mW1[2] / number_of_rotor_lobules;
    }

    MoveRotorProcess(ModelPart& rModelPart,
                    Parameters rParameters): Process(), mrModelPart(rModelPart)
    {
        Parameters default_parameters( R"(
            {
                "model_part_name":"rotor",
                "angular_velocity_with_respect_to_stator_center":0.0,
                "coordinates_of_stator_center_x":0.0,
                "coordinates_of_stator_center_y":0.0,
                "initial_coordinates_of_rotor_center_x":0.0,
                "initial_coordinates_of_rotor_center_y":0.0,
                "number_of_rotor_lobules": 0
            }  )" );

        rParameters.ValidateAndAssignDefaults(default_parameters);

        mW1[0] = 0.0;
        mW1[1] = 0.0;
        mW1[2] = rParameters["angular_velocity_with_respect_to_stator_center"].GetDouble();
        mCoordinatesOfStatorCenter[0] = rParameters["coordinates_of_stator_center_x"].GetDouble();
        mCoordinatesOfStatorCenter[1] = rParameters["coordinates_of_stator_center_y"].GetDouble();
        mCoordinatesOfStatorCenter[2] = 0.0;
        mInitialCoordinatesOfRotorCenter[0] = rParameters["initial_coordinates_of_rotor_center_x"].GetDouble();
        mInitialCoordinatesOfRotorCenter[1] = rParameters["initial_coordinates_of_rotor_center_y"].GetDouble();
        mInitialCoordinatesOfRotorCenter[2] = 0.0;
        const array_1d<double,3> centers_distance = mInitialCoordinatesOfRotorCenter - mCoordinatesOfStatorCenter;
        mEccentricity = sqrt(centers_distance[0]*centers_distance[0] + centers_distance[1]*centers_distance[1]);
        mW2[0] = 0.0;
        mW2[1] = 0.0;
        mW2[2] = -mW1[2] / rParameters["number_of_rotor_lobules"].GetInt();

    }

    /// Destructor.
    ~MoveRotorProcess() override{}


    void ExecuteBeforeSolutionLoop() override
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }


    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
	    const double& rCurrentTime = rCurrentProcessInfo[TIME];

        if (!mrModelPart.NodesBegin()->SolutionStepsDataHas(DISPLACEMENT))
            KRATOS_ERROR << "DISPLACEMENT variable is not in move rotor submodelpart.";

        if (!mrModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITY))
            KRATOS_ERROR << "VELOCITY variable is not in move rotor submodelpart.";

        //UPDATE POSITION OF ROTOR AXIS:
        const double initial_angle = atan2( (mInitialCoordinatesOfRotorCenter[1] - mCoordinatesOfStatorCenter[1]), (mInitialCoordinatesOfRotorCenter[0] - mCoordinatesOfStatorCenter[0]) );
        const double rotated_angle1 = mW1[2] * rCurrentTime;
        const double current_angle1 = initial_angle + rotated_angle1;
        array_1d<double,3> vector_to_rotor_center;
        vector_to_rotor_center[0] = mEccentricity*cos(current_angle1);
        vector_to_rotor_center[1] = mEccentricity*sin(current_angle1);
        vector_to_rotor_center[2] = 0.0;
        const array_1d<double,3> current_rotor_position = mCoordinatesOfStatorCenter + vector_to_rotor_center;

        //UPDATE VELOCITY OF ROTOR AXIS:
        array_1d<double, 3> rotor_velocity;
        MathUtils<double>::CrossProduct(rotor_velocity, mW1, vector_to_rotor_center);

        //UPDATE LOCAL AXES (ROTATE THEM AROUND (0,0,1) THE ROTATED ANGLE )
        const double rotated_angle2 = mW2[2] * rCurrentTime;
        array_1d<double,3> current_local_axis_1;
        array_1d<double,3> current_local_axis_2;
        array_1d<double,3> current_local_axis_3;

        array_1d<double,3> vertical_vector; vertical_vector[0] = 0.0; vertical_vector[1] = 0.0; vertical_vector[2] = 1.0;
        array_1d<double,3> initial_local_axis_1; initial_local_axis_1[0] = 1.0; initial_local_axis_1[1] = 0.0; initial_local_axis_1[2] = 0.0; //(local axes are assumed oriented as global axes at the beginning)
        array_1d<double,3> initial_local_axis_2; initial_local_axis_2[0] = 0.0; initial_local_axis_2[1] = 1.0; initial_local_axis_2[2] = 0.0; //(local axes are assumed oriented as global axes at the beginning)

        RotateAVectorAGivenAngleAroundAUnitaryVector(initial_local_axis_1, vertical_vector, rotated_angle2, current_local_axis_1);
        RotateAVectorAGivenAngleAroundAUnitaryVector(initial_local_axis_2, vertical_vector, rotated_angle2, current_local_axis_2);

        //UPDATE POSITION AND VELOCITY OF ALL NODES
        for (ModelPart::NodesContainerType::iterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); node_i++) {
            //Get local coordinates at the beginning (local axes are assumed oriented as global axes at the beginning)
            array_1d<double,3> local_coordinates;
            local_coordinates[0] = node_i->X0() - mInitialCoordinatesOfRotorCenter[0];
            local_coordinates[1] = node_i->Y0() - mInitialCoordinatesOfRotorCenter[1];

            //Use local coordinates with the updated local axes
            array_1d<double,3> from_rotor_center_to_node;
            noalias(from_rotor_center_to_node) = local_coordinates[0] * current_local_axis_1 + local_coordinates[1] * current_local_axis_2;

            array_1d<double,3>& current_node_position = node_i->Coordinates();
            current_node_position[0] = current_rotor_position[0] + local_coordinates[0] * current_local_axis_1[0] + local_coordinates[1] * current_local_axis_2[0];
            current_node_position[1] = current_rotor_position[1] + local_coordinates[0] * current_local_axis_1[1] + local_coordinates[1] * current_local_axis_2[1];

            noalias(node_i->FastGetSolutionStepValue(DISPLACEMENT)) = node_i->Coordinates() - node_i->GetInitialPosition();

            array_1d<double,3>& current_node_velocity = node_i->FastGetSolutionStepValue(VELOCITY);

            array_1d<double, 3> aux_velocity;
            MathUtils<double>::CrossProduct(aux_velocity, mW2, from_rotor_center_to_node);
            noalias(current_node_velocity) = rotor_velocity + aux_velocity;

        }//end of loop over nodes
        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "MoveRotorProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "MoveRotorProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ModelPart&                                 mrModelPart;
    array_1d<double,3>                                 mW1;
    array_1d<double,3>                                 mW2;
    double                                   mEccentricity;
    array_1d<double,3>                         mLocalAxis1;
    array_1d<double,3>                         mLocalAxis2;
    array_1d<double,3>                         mLocalAxis3;
    array_1d<double,3>    mInitialCoordinatesOfRotorCenter;
    array_1d<double,3>          mCoordinatesOfStatorCenter;

    void RotateAVectorAGivenAngleAroundAUnitaryVector(const array_1d<double, 3>& old_vec, const array_1d<double, 3>& axis,
                                                                const double ang, array_1d<double, 3>& new_vec) {
        double cang = cos(ang);
        double sang = sin(ang);

        new_vec[0] = axis[0] * (axis[0] * old_vec[0] + axis[1] * old_vec[1] + axis[2] * old_vec[2]) * (1 - cang) + old_vec[0] * cang + (-axis[2] * old_vec[1] + axis[1] * old_vec[2]) * sang;
        new_vec[1] = axis[1] * (axis[0] * old_vec[0] + axis[1] * old_vec[1] + axis[2] * old_vec[2]) * (1 - cang) + old_vec[1] * cang + ( axis[2] * old_vec[0] - axis[0] * old_vec[2]) * sang;
        new_vec[2] = axis[2] * (axis[0] * old_vec[0] + axis[1] * old_vec[1] + axis[2] * old_vec[2]) * (1 - cang) + old_vec[2] * cang + (-axis[1] * old_vec[0] + axis[0] * old_vec[1]) * sang;
    }



private:

    /// Assignment operator.
    MoveRotorProcess& operator=(MoveRotorProcess const& rOther){return *this;}

}; // Class MoveRotorProcess

};  // namespace Kratos.

#endif // KRATOS_MOVE_ROTOR_PROCESS_H
