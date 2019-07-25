//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
// ==============================================================================
//

#ifndef ROTATE_REGION_PROCESS_H
#define ROTATE_REGION_PROCESS_H

// System includes
#include <string>
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "utilities/quaternion.h"

// Application includes

namespace Kratos
{

class KRATOS_API(CHIMERA_APPLICATION) RotateRegionProcess : public Process
{
public:
    /// Pointer definition of MoveRotorProcess
    KRATOS_CLASS_POINTER_DEFINITION(RotateRegionProcess);

    typedef ModelPart::NodeIterator NodeIteratorType;
    typedef ProcessInfo ProcessInfoType;
    typedef ProcessInfo::Pointer ProcessInfoPointerType;
    typedef Matrix MatrixType;
    typedef Vector VectorType;

    /// Constructor.
    RotateRegionProcess(ModelPart &rModelPart,
                        Parameters rParameters) : Process(Flags()), mrModelPart(rModelPart), mParameters(rParameters)
    {

        Parameters default_parameters(R"(
            {
                "model_part_name":"SPECIFY_MODELPART_NAME",
                "center_of_rotation":[],
                "angular_velocity_radians":0.0,
                "axis_of_rotation":[],
                "is_ale" : false
            }  )");

        mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mAngularVelocityRadians = mParameters["angular_velocity_radians"].GetDouble();
        mCenterOfRotation = mParameters["center_of_rotation"].GetVector();
        auto axis_of_rotation_raw = mParameters["axis_of_rotation"].GetVector();
        double norm = norm_2(axis_of_rotation_raw);
        mAxisOfRotationVector = axis_of_rotation_raw/norm;
        mTheta = 0.0;
    }

    /// Destructor.
    virtual ~RotateRegionProcess()
    {
    }

    void ExecuteBeforeSolutionLoop() override
    {

    }

    void SetAngularVelocity(const double NewAngularVelocity)
    {
        mAngularVelocityRadians = NewAngularVelocity;
    }

    void ExecuteFinalizeSolutionStep() override
    {
        KRATOS_TRY;
        const auto &r_process_info = mrModelPart.GetProcessInfo();
        const double dt = r_process_info[DELTA_TIME];
        const int domain_size = r_process_info[DOMAIN_SIZE];
        mTheta += mAngularVelocityRadians * dt;

        const int num_nodes = mrModelPart.NumberOfNodes();
        const NodeIteratorType it_slave_node_begin = mrModelPart.NodesBegin();

#pragma omp parallel for schedule(guided, 512)
        for (int i_node = 0; i_node < num_nodes; ++i_node)
        {
            NodeIteratorType it_node = it_slave_node_begin;
            std::advance(it_node, i_node);

            /// Calculating the displacement of the current node
            array_1d<double, 3> transformed_coordinates;
            TransformNode(it_node->GetInitialPosition().Coordinates(), transformed_coordinates);

            it_node->X() = transformed_coordinates[0];
            it_node->Y() = transformed_coordinates[1];
            if (domain_size > 2)
                it_node->Z() = transformed_coordinates[2];

            // Computing the linear velocity at this it_node
            DenseVector<double> radius(3);
            DenseVector<double> linearVelocity(3);
            radius[0] = it_node->X() - mCenterOfRotation[0];
            radius[1] = it_node->Y() - mCenterOfRotation[1];
            radius[2] = it_node->Z() - mCenterOfRotation[2];
            CalculateLinearVelocity(mAxisOfRotationVector, radius, linearVelocity);
            if (mParameters["is_ale"].GetBool())
            {
                it_node->FastGetSolutionStepValue(MESH_VELOCITY_X, 0) = mAngularVelocityRadians * linearVelocity[0];
                it_node->FastGetSolutionStepValue(MESH_VELOCITY_Y, 0) = mAngularVelocityRadians * linearVelocity[1];
                if (domain_size > 2)
                    it_node->FastGetSolutionStepValue(MESH_VELOCITY_Z, 0) = mAngularVelocityRadians * linearVelocity[2];

                if (it_node->IsFixed(VELOCITY_X))
                    it_node->FastGetSolutionStepValue(VELOCITY_X, 0) = it_node->FastGetSolutionStepValue(MESH_VELOCITY_X, 0);

                if (it_node->IsFixed(VELOCITY_Y))
                    it_node->FastGetSolutionStepValue(VELOCITY_Y, 0) = it_node->FastGetSolutionStepValue(MESH_VELOCITY_Y, 0);

                if (domain_size > 2)
                    if (it_node->IsFixed(VELOCITY_Z))
                        it_node->FastGetSolutionStepValue(VELOCITY_Z, 0) = it_node->FastGetSolutionStepValue(MESH_VELOCITY_Z, 0);
            }
        }

        KRATOS_CATCH("");
    }

    void ExecuteInitializeSolutionStep() override
    {
    }

    void ExecuteAfterOutputStep() override
    {
    }

    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "RotateRegionProcess";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override { rOStream << "RotateRegionProcess"; }

    /// Print object's data.
    void PrintData()
    {
    }

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}

private:
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ModelPart &mrModelPart;
    Parameters mParameters;
    std::string mSubModelPartName;
    double mAngularVelocityRadians;
    DenseVector<double> mAxisOfRotationVector;
    DenseVector<double> mCenterOfRotation;
    double mTheta;

    ///@}

    /// Assignment operator.
    RotateRegionProcess &operator=(RotateRegionProcess const &rOther) { return *this; }

    /*
     * @brief Calculates the linear velocity v = r x w
     * @param rAxisOfRotationVector the axis of rotation vector
     * @param rRadius the radius vector of the point for which v is to be calculated.
     * @out   rLinearVelocity the calculated linear velocity.
     */
    void CalculateLinearVelocity(const DenseVector<double> &rAxisOfRotationVector,
                                 const DenseVector<double> &rRadius,
                                 DenseVector<double> &rLinearVelocity)
    {
        assert(rAxisOfRotationVector.size() == rRadius.size());
        rLinearVelocity[0] = rAxisOfRotationVector[1] * rRadius[2] - rAxisOfRotationVector[2] * rRadius[1];
        rLinearVelocity[1] = rAxisOfRotationVector[2] * rRadius[0] - rAxisOfRotationVector[0] * rRadius[2];
        rLinearVelocity[2] = rAxisOfRotationVector[0] * rRadius[1] - rAxisOfRotationVector[1] * rRadius[0];
    }

    /*
     * @brief Rotates the given node by mTheta around the mAxisOfRotationVector and mCenterOfRotation
     *          This function uses Quaternion for rotations.
     * @param rCoordinates The nodal coordinates which are to be rotated.
     * @out rTransformedCoordinates the rotated nodal coordinates.
     */
    void TransformNode(const array_1d<double, 3> &rCoordinates, array_1d<double, 3> &rTransformedCoordinates) const
    {
        // Changing the origin to the center of rotation
        auto new_coordinates = rCoordinates - mCenterOfRotation;
        Quaternion<double> mQuaternion = Quaternion<double>::FromAxisAngle(mAxisOfRotationVector(0),
                                                                           mAxisOfRotationVector(1),
                                                                           mAxisOfRotationVector(2), mTheta);
        mQuaternion.RotateVector3(new_coordinates, rTransformedCoordinates);
        // Bringing back to the original coordinate system.
        rTransformedCoordinates = rTransformedCoordinates + mCenterOfRotation;
    }
}; // Class MoveRotorProcess

}; // namespace Kratos.

#endif // KRATOS_MOVE_ROTOR_PROCESS_H
