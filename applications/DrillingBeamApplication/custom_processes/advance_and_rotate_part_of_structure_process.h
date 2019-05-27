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

class AdvanceAndRotatePartOfStructureProcess : public Process
{
public:

    /// Pointer definition of AdvanceAndRotatePartOfStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(AdvanceAndRotatePartOfStructureProcess);

    typedef Node<3>                     NodeType;
    typedef Geometry<NodeType>      GeometryType;

    /// Constructor.
    AdvanceAndRotatePartOfStructureProcess(ModelPart& rModelPart,
                    const double angular_velocity_x,
                    const double angular_velocity_y,
                    const double angular_velocity_z,
                    const double coordinates_rotation_center_x,
                    const double coordinates_rotation_center_y,
                    const double coordinates_rotation_center_z,
                    const double advance_velocity): Process(), mrModelPart(rModelPart)
    {
        mAngularVelocity[0] = angular_velocity_x;
        mAngularVelocity[1] = angular_velocity_y;
        mAngularVelocity[2] = angular_velocity_z;
        mCoordinatesOfCenter[0] = coordinates_rotation_center_x;
        mCoordinatesOfCenter[1] = coordinates_rotation_center_y;
        mCoordinatesOfCenter[2] = coordinates_rotation_center_z;
        mAdvanceVelocity = advance_velocity;
    }

    AdvanceAndRotatePartOfStructureProcess(ModelPart& rModelPart,
                    Parameters rParameters): Process(), mrModelPart(rModelPart)
    {
        Parameters default_parameters( R"(
            {
                "model_part_name":"rotor",
                "angular_velocity_x":0.0,
                "angular_velocity_y":0.0,
                "angular_velocity_z":0.0,
                "coordinates_rotation_center_x":0.0,
                "coordinates_rotation_center_y":0.0,
                "coordinates_rotation_center_z":0.0,
                "advance_velocity": 0.0
            }  )" );

        rParameters.ValidateAndAssignDefaults(default_parameters);

        mAngularVelocity[0] = rParameters["angular_velocity_x"].GetDouble();
        mAngularVelocity[1] = rParameters["angular_velocity_y"].GetDouble();
        mAngularVelocity[2] = rParameters["angular_velocity_z"].GetDouble();
        mCoordinatesOfCenter[0] = rParameters["coordinates_rotation_center_x"].GetDouble();
        mCoordinatesOfCenter[1] = rParameters["coordinates_rotation_center_y"].GetDouble();
        mCoordinatesOfCenter[2] = rParameters["coordinates_rotation_center_z"].GetDouble();;
        mAdvanceVelocity = rParameters["advance_velocity"].GetDouble();;
    }

    /// Destructor.
    ~AdvanceAndRotatePartOfStructureProcess() override{}


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
            KRATOS_ERROR << "DISPLACEMENT variable is not in the submodelpart to be rotated.";

        if (!mrModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITY))
            KRATOS_ERROR << "VELOCITY variable is not in the submodelpart to be rotated.";

        //UPDATE LOCAL AXES (ROTATE THEM AROUND mAngularVelocity THE ROTATED ANGLE )
        const double norm_of_angular_velocity = MathUtils<double>::Norm3(mAngularVelocity);
        const double rotated_angle = norm_of_angular_velocity * rCurrentTime;
        array_1d<double,3> unitary_angular_velocity = mAngularVelocity * (1.0 / norm_of_angular_velocity);
        array_1d<double,3> initial_local_axis_1; initial_local_axis_1[0] = 1.0; initial_local_axis_1[1] = 0.0; initial_local_axis_1[2] = 0.0;
        array_1d<double,3> initial_local_axis_2; initial_local_axis_2[0] = 0.0; initial_local_axis_2[1] = 1.0; initial_local_axis_2[2] = 0.0;
        array_1d<double,3> initial_local_axis_3; initial_local_axis_3[0] = 0.0; initial_local_axis_3[1] = 0.0; initial_local_axis_3[2] = 1.0;
        array_1d<double,3> current_local_axis_1;
        array_1d<double,3> current_local_axis_2;
        array_1d<double,3> current_local_axis_3;

        array_1d<double,3> current_coordinates_of_center;
        current_coordinates_of_center[0] = mCoordinatesOfCenter[0] + mAdvanceVelocity * rCurrentTime;
        current_coordinates_of_center[1] = 0.0;
        current_coordinates_of_center[2] = 0.0;
        array_1d<double,3> current_velocity_of_center;
        current_velocity_of_center[0] = mAdvanceVelocity;
        current_velocity_of_center[1] = 0.0;
        current_velocity_of_center[2] = 0.0;

        RotateAVectorAGivenAngleAroundAUnitaryVector(initial_local_axis_1, unitary_angular_velocity, rotated_angle, current_local_axis_1);
        RotateAVectorAGivenAngleAroundAUnitaryVector(initial_local_axis_2, unitary_angular_velocity, rotated_angle, current_local_axis_2);
        RotateAVectorAGivenAngleAroundAUnitaryVector(initial_local_axis_3, unitary_angular_velocity, rotated_angle, current_local_axis_3);

        //UPDATE POSITION AND VELOCITY OF ALL NODES
        for (ModelPart::NodesContainerType::iterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); node_i++) {

            //Get local coordinates at the beginning (local axes are assumed oriented as global axes at the beginning)
            array_1d<double,3> local_coordinates;
            local_coordinates[0] = node_i->X0() - mCoordinatesOfCenter[0];
            local_coordinates[1] = node_i->Y0() - mCoordinatesOfCenter[1];
            local_coordinates[2] = node_i->Z0() - mCoordinatesOfCenter[2];

            //Use local coordinates with the updated local axes
            array_1d<double,3> from_center_to_node;
            noalias(from_center_to_node) = local_coordinates[0] * current_local_axis_1 + local_coordinates[1] * current_local_axis_2 + local_coordinates[2] * current_local_axis_3;

            array_1d<double,3>& current_node_position = node_i->Coordinates();
            noalias(current_node_position) = current_coordinates_of_center + from_center_to_node;

            noalias(node_i->FastGetSolutionStepValue(DISPLACEMENT)) = node_i->Coordinates() - node_i->GetInitialPosition();

            array_1d<double,3>& current_node_velocity = node_i->FastGetSolutionStepValue(VELOCITY);

            array_1d<double, 3> aux_velocity;
            MathUtils<double>::CrossProduct(aux_velocity, mAngularVelocity, from_center_to_node);
            noalias(current_node_velocity) = aux_velocity + current_velocity_of_center;

        }//end of loop over nodes
        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "AdvanceAndRotatePartOfStructureProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "AdvanceAndRotatePartOfStructureProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ModelPart&                  mrModelPart;
    array_1d<double,3>     mAngularVelocity;
    array_1d<double,3> mCoordinatesOfCenter;
    double mAdvanceVelocity;

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
    AdvanceAndRotatePartOfStructureProcess& operator=(AdvanceAndRotatePartOfStructureProcess const& rOther){return *this;}

}; // Class AdvanceAndRotatePartOfStructureProcess

};  // namespace Kratos.

#endif // KRATOS_MOVE_ROTOR_PROCESS_H
