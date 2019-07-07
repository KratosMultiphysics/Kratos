//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//  The implementation follows : http://paulbourke.net/geometry/rotate/

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
                "center_of_rotation":[],
                "angular_velocity_radians":0.0,
                "axis_of_rotation":[],
                "is_ale" : false
            }  )");

        mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mAngularVelocityRadians = mParameters["angular_velocity_radians"].GetDouble();
        mCenterOfRotation = mParameters["center_of_rotation"].GetVector();
        mAxisOfRotationVector = mParameters["axis_of_rotation"].GetVector();
        mAxisOfRoationVectorNormalized.resize(3);

        mTransformationMatrix.resize(4, 4, false);

        // normalizing the axis of roatation
        double norm = 0.0;
        for (std::size_t d = 0; d < 3; ++d)
            norm += mAxisOfRotationVector[d] * mAxisOfRotationVector[d];
        norm = sqrt(norm);
        for (std::size_t d = 0; d < 3; ++d)
            mAxisOfRoationVectorNormalized[d] = (mAxisOfRotationVector[d] / norm);

        CalculateRotationMatrix(-1 * mAngularVelocityRadians, mTransformationMatrix);
        mTheta = 0.0;
    }

    /// Destructor.
    virtual ~RotateRegionProcess()
    {
    }

    void ExecuteBeforeSolutionLoop() override
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

    void SetAngularVelocity(const double NewAngularVelocity)
    {
        mAngularVelocityRadians = NewAngularVelocity;
    }

    void ExecuteInitializeSolutionStep() override
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
            TransformNode(it_node->Coordinates(), transformed_coordinates);
            //TransformNodeWithQuaternion(it_node->Coordinates(), transformed_coordinates);

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
                it_node->FastGetSolutionStepValue(MESH_VELOCITY_X, 0) += mAngularVelocityRadians * linearVelocity[0];
                it_node->FastGetSolutionStepValue(MESH_VELOCITY_Y, 0) += mAngularVelocityRadians * linearVelocity[1];
                if (domain_size > 2)
                    it_node->FastGetSolutionStepValue(MESH_VELOCITY_Z, 0) += mAngularVelocityRadians * linearVelocity[2];

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

    void ExecuteFinalizeSolutionStep() override
    {
        const auto &r_process_info = mrModelPart.GetProcessInfo();

        int domain_size = r_process_info[DOMAIN_SIZE];
        const int num_nodes = mrModelPart.NumberOfNodes();
        const NodeIteratorType it_slave_node_begin = mrModelPart.NodesBegin();

#pragma omp parallel for schedule(guided, 512)
        for (int i_node = 0; i_node < num_nodes; ++i_node)
        {
            NodeIteratorType it_node = it_slave_node_begin;
            std::advance(it_node, i_node);
            if (mParameters["is_ale"].GetBool())
            {
                it_node->FastGetSolutionStepValue(MESH_VELOCITY_X, 0) = 0.0;
                it_node->FastGetSolutionStepValue(MESH_VELOCITY_Y, 0) = 0.0;
                if (domain_size > 2)
                    it_node->FastGetSolutionStepValue(MESH_VELOCITY_Z, 0) = 0.0;
            }
        }
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
    MatrixType mTransformationMatrix;
    DenseVector<double> mAxisOfRotationVector;
    DenseVector<double> mCenterOfRotation;
    DenseVector<double> mAxisOfRoationVectorNormalized;
    double mTheta;

    ///@}

    /// Assignment operator.
    RotateRegionProcess &operator=(RotateRegionProcess const &rOther) { return *this; }

    /*
     * Calculates the cross product of two vectors
     *  c = axb
     */
    void CalculateLinearVelocity(const DenseVector<double> &mAxisOfRotationVector,
                                 const DenseVector<double> &rRadius,
                                 DenseVector<double> &rLinearVelocity)
    {
        assert(mAxisOfRotationVector.size() == rRadius.size());
        rLinearVelocity[0] = mAxisOfRotationVector[1] * rRadius[2] - mAxisOfRotationVector[2] * rRadius[1];
        rLinearVelocity[1] = mAxisOfRotationVector[2] * rRadius[0] - mAxisOfRotationVector[0] * rRadius[2];
        rLinearVelocity[2] = mAxisOfRotationVector[0] * rRadius[1] - mAxisOfRotationVector[1] * rRadius[0];
    }

    void TransformNode(const array_1d<double, 3> &rCoordinates, array_1d<double, 3> &rTransformedCoordinates) const
    {
        DenseVector<double> original_node(4, 0.0);
        DenseVector<double> transformed_node(4, 0.0);

        original_node[0] = rCoordinates(0);
        original_node[1] = rCoordinates(1);
        original_node[2] = rCoordinates(2);
        original_node[3] = 1.0;
        // Multiplying the point to get the rotated point
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                transformed_node[i] += mTransformationMatrix(i, j) * original_node[j];
            }
        }
        rTransformedCoordinates(0) = transformed_node[0];
        rTransformedCoordinates(1) = transformed_node[1];
        rTransformedCoordinates(2) = transformed_node[2];
    }

    void TransformNodeWithQuaternion(const array_1d<double, 3> &rCoordinates, array_1d<double, 3> &rTransformedCoordinates) const
    {
        Quaternion<double> mQuaternion = Quaternion<double>::FromAxisAngle(mAxisOfRotationVector(1), mAxisOfRotationVector(2), mAxisOfRotationVector(3), mTheta);
        mQuaternion.RotateVector3(rCoordinates, rTransformedCoordinates);
    }

    void CalculateRotationMatrix(const double Theta, MatrixType &rMatrix)
    {
        DenseVector<double> U(3); // normalized axis of rotation
        // normalizing the axis of rotation
        double norm = 0.0;
        for (IndexType d = 0; d < 3; ++d)
            norm += mAxisOfRotationVector[d] * mAxisOfRotationVector[d];
        norm = sqrt(norm);
        KRATOS_ERROR_IF(norm < std::numeric_limits<double>::epsilon()) << "Norm of the provided axis of rotation is Zero !" << std::endl;
        for (IndexType d = 0; d < 3; ++d)
            U[d] = mAxisOfRotationVector[d] / norm;

        // Constructing the transformation matrix
        const double x1 = mCenterOfRotation[0];
        const double y1 = mCenterOfRotation[1];
        const double z1 = mCenterOfRotation[2];

        const double a = U[0];
        const double b = U[1];
        const double c = U[2];

        const double t2 = std::cos(Theta);
        const double t3 = std::sin(Theta);
        const double t4 = a * a;
        const double t5 = b * b;
        const double t6 = c * c;
        const double t7 = a * b;
        const double t8 = t5 + t6;
        const double t9 = std::abs(t8) < std::numeric_limits<double>::epsilon() ? 1.0e8 : 1.0 / t8;
        const double t10 = a * c;
        const double t11 = b * t3;
        const double t12 = a * t3 * t5;
        const double t13 = a * t3 * t6;
        const double t14 = b * c * t2;
        rMatrix(0, 0) = t4 + t2 * t8;
        rMatrix(0, 1) = t7 - c * t3 - a * b * t2;
        rMatrix(0, 2) = t10 + t11 - a * c * t2;
        rMatrix(0, 3) = x1 - t4 * x1 - a * b * y1 - a * c * z1 - b * t3 * z1 + c * t3 * y1 - t2 * t5 * x1 - t2 * t6 * x1 + a * b * t2 * y1 + a * c * t2 * z1;
        rMatrix(1, 0) = t7 + c * t3 - a * b * t2;
        rMatrix(1, 1) = t9 * (t2 * t6 + t5 * t8 + t2 * t4 * t5);
        rMatrix(1, 2) = -t9 * (t12 + t13 + t14 - b * c * t8 - b * c * t2 * t4);
        rMatrix(1, 3) = -t9 * (-t8 * y1 + t2 * t6 * y1 + t5 * t8 * y1 + a * b * t8 * x1 - b * c * t2 * z1 + b * c * t8 * z1 - a * t3 * t5 * z1 - a * t3 * t6 * z1 + c * t3 * t8 * x1 + t2 * t4 * t5 * y1 - a * b * t2 * t8 * x1 + b * c * t2 * t4 * z1);
        rMatrix(2, 0) = t10 - t11 - a * c * t2;
        rMatrix(2, 1) = t9 * (t12 + t13 - t14 + b * c * t8 + b * c * t2 * t4);
        rMatrix(2, 2) = t9 * (t2 * t5 + t6 * t8 + t2 * t4 * t6);
        rMatrix(2, 3) = -t9 * (-t8 * z1 + t2 * t5 * z1 + t6 * t8 * z1 + a * c * t8 * x1 - b * c * t2 * y1 + b * c * t8 * y1 + a * t3 * t5 * y1 + a * t3 * t6 * y1 - b * t3 * t8 * x1 + t2 * t4 * t6 * z1 - a * c * t2 * t8 * x1 + b * c * t2 * t4 * y1);
        rMatrix(3, 0) = 0.0;
        rMatrix(3, 1) = 0.0;
        rMatrix(3, 2) = 0.0;
        rMatrix(3, 3) = 1.0;
    }
}; // Class MoveRotorProcess

}; // namespace Kratos.

#endif // KRATOS_MOVE_ROTOR_PROCESS_H
