//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             MACeligueta $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//

#ifndef KRATOS_CONSTANT_ROTATION_PROCESS_H_INCLUDED
#define KRATOS_CONSTANT_ROTATION_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"

// Application includes
namespace Kratos
{

class ConstantRotationProcess : public Process
{
public:

    /// Pointer definition of ConstantRotationProcess
    KRATOS_CLASS_POINTER_DEFINITION(ConstantRotationProcess);

    typedef Node<3>                     NodeType;
    typedef Geometry<NodeType>      GeometryType;

    /// Constructor.
    ConstantRotationProcess(ModelPart& rModelPart,
                    const double angular_velocity_x,
                    const double angular_velocity_y,
                    const double angular_velocity_z,
                    const double center_x,
                    const double center_y,
                    const double center_z): Process(), mrModelPart(rModelPart)
    {
        mW[0] = angular_velocity_x;
        mW[1] = angular_velocity_y;
        mW[2] = angular_velocity_z;
        mCenter[0] = center_x;
        mCenter[1] = center_y;
        mCenter[2] = center_z;
    }

    ConstantRotationProcess(ModelPart& rModelPart,
                    Parameters rParameters): Process(), mrModelPart(rModelPart)
    {
        Parameters default_parameters( R"(
            {
                "model_part_name":"rotating_part",
                "angular_velocity_x":0.0,
                "angular_velocity_y":0.0,
                "angular_velocity_z":1.0,
                "center_x":0.0,
                "center_y":0.0,
                "center_z":0.0,
                "interval": [0,"End"]
            }  )" );

        rParameters.ValidateAndAssignDefaults(default_parameters);

        mW[0] = rParameters["angular_velocity_x"].GetDouble();
        mW[1] = rParameters["angular_velocity_y"].GetDouble();
        mW[2] = rParameters["angular_velocity_z"].GetDouble();
        mCenter[0] = rParameters["center_x"].GetDouble();
        mCenter[1] = rParameters["center_y"].GetDouble();
        mCenter[2] = rParameters["center_z"].GetDouble();

        if(rParameters.Has("interval")) {
            if(rParameters["interval"][1].IsString()) {
                if(rParameters["interval"][1].GetString() == "End") rParameters["interval"][1].SetDouble(1e30);
                else KRATOS_THROW_ERROR(std::runtime_error, "The second value of interval can be \"End\" or a number", 0);
            }
            mInitialTime = rParameters["interval"][0].GetDouble();
            mFinalTime   = rParameters["interval"][1].GetDouble();
        }
    }

    /// Destructor.
    virtual ~ConstantRotationProcess(){}


    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
	const double& rCurrentTime = rCurrentProcessInfo[TIME];

        if(rCurrentTime < mInitialTime || rCurrentTime > mFinalTime) {
            SetAllNodesVelocityToZero();
            return;
        }

        array_1d<double,3> current_local_axis_1;
        array_1d<double,3> current_local_axis_2;
        array_1d<double,3> current_local_axis_3;

        const double modulus_omega = sqrt(mW[0]*mW[0] + mW[1]*mW[1] + mW[2]*mW[2]);
        const double rotated_angle = modulus_omega * rCurrentTime;

        array_1d<double,3> unitary_omega;
        noalias(unitary_omega) = mW / modulus_omega;

        array_1d<double,3> initial_local_axis_1; initial_local_axis_1[0] = 1.0; initial_local_axis_1[1] = 0.0; initial_local_axis_1[2] = 0.0; //(local axes are assumed oriented as global axes at the beginning)
        array_1d<double,3> initial_local_axis_2; initial_local_axis_2[0] = 0.0; initial_local_axis_2[1] = 1.0; initial_local_axis_2[2] = 0.0; //(local axes are assumed oriented as global axes at the beginning)
        array_1d<double,3> initial_local_axis_3; initial_local_axis_2[0] = 0.0; initial_local_axis_3[1] = 0.0; initial_local_axis_3[2] = 1.0; //(local axes are assumed oriented as global axes at the beginning)

        RotateAVectorAGivenAngleAroundAUnitaryVector(initial_local_axis_1, unitary_omega, rotated_angle, current_local_axis_1);
        RotateAVectorAGivenAngleAroundAUnitaryVector(initial_local_axis_2, unitary_omega, rotated_angle, current_local_axis_2);
        RotateAVectorAGivenAngleAroundAUnitaryVector(initial_local_axis_3, unitary_omega, rotated_angle, current_local_axis_3);

        //UPDATE POSITION AND VELOCITY OF ALL NODES
        for (ModelPart::NodesContainerType::iterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); ++node_i) {
            //Get local coordinates at the beginning (local axes are assumed oriented as global axes at the beginning)
            array_1d<double,3> local_coordinates;
            local_coordinates[0] = node_i->X0() - mCenter[0];
            local_coordinates[1] = node_i->Y0() - mCenter[1];
            local_coordinates[2] = node_i->Z0() - mCenter[2];

            //Use local coordinates with the updated local axes
            array_1d<double,3> new_from_center_to_node;
            new_from_center_to_node = local_coordinates[0] * current_local_axis_1 + local_coordinates[1] * current_local_axis_2 + local_coordinates[2] * current_local_axis_3;

            array_1d<double,3>& current_node_position = node_i->Coordinates();
            noalias(current_node_position) = mCenter + new_from_center_to_node;

            node_i->pGetDof(VELOCITY_X)->FixDof();
            node_i->pGetDof(VELOCITY_Y)->FixDof();
            node_i->pGetDof(VELOCITY_Z)->FixDof();
            array_1d<double,3>& current_node_velocity = node_i->FastGetSolutionStepValue(VELOCITY);
            //noalias(current_node_velocity) = MathUtils<double>::CrossProduct(new_from_center_to_node, mW);
            MathUtils<double>::CrossProduct(current_node_velocity, mW, new_from_center_to_node);
        }//end of loop over nodes
        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ConstantRotationProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "ConstantRotationProcess";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}


protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ModelPart&                                 mrModelPart;
    array_1d<double,3>                                  mW;
    array_1d<double,3>                             mCenter;
    double                                    mInitialTime;
    double                                      mFinalTime;

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
    ConstantRotationProcess& operator=(ConstantRotationProcess const& rOther){return *this;}

    void SetAllNodesVelocityToZero(){
        KRATOS_TRY;
        for (ModelPart::NodesContainerType::iterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); ++node_i) {
            node_i->pGetDof(VELOCITY_X)->FixDof();
            node_i->pGetDof(VELOCITY_Y)->FixDof();
            node_i->pGetDof(VELOCITY_Z)->FixDof();
            array_1d<double,3>& current_node_velocity = node_i->FastGetSolutionStepValue(VELOCITY);
            current_node_velocity[0] = 0.0;
            current_node_velocity[1] = 0.0;
            current_node_velocity[2] = 0.0;
        }//end of loop over nodes
        KRATOS_CATCH("");

    }

}; // Class ConstantRotationProcess

};  // namespace Kratos.

#endif // KRATOS_CONSTANT_ROTATION_PROCESS_H_INCLUDED
