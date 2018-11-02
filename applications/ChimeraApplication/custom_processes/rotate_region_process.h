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

// Application includes

namespace Kratos
{

class RotateRegionProcess : public Process
{
  public:
    /// Pointer definition of MoveRotorProcess
    KRATOS_CLASS_POINTER_DEFINITION(RotateRegionProcess);

    typedef ModelPart::NodeIterator NodeIteratorType;
    typedef ProcessInfo ProcessInfoType;
    typedef ProcessInfo::Pointer ProcessInfoPointerType;

    /// Constructor.
    RotateRegionProcess(ModelPart &model_part,
                        Parameters rParameters) : Process(Flags()), mr_model_part(model_part), m_parameters(rParameters)
    {

        Parameters default_parameters(R"(
            {
                "movement_name":"default",
                "sub_model_part_name":"default",
                "sub_model_part_boundary_name":"default",
                "center_of_rotation":[],
                "angular_velocity_radians":0.0,
                "axis_of_rotation":[],
                "is_ale" : false
            }  )");

        m_sub_model_part_name = m_parameters["sub_model_part_name"].GetString();
        m_sub_model_part_boundary_name = m_parameters["sub_model_part_boundary_name"].GetString();
        m_angular_velocity_radians = m_parameters["angular_velocity_radians"].GetDouble();

        mCenterOfRotation.push_back(m_parameters["center_of_rotation"][0].GetDouble());
        mCenterOfRotation.push_back(m_parameters["center_of_rotation"][1].GetDouble());
        mCenterOfRotation.push_back(m_parameters["center_of_rotation"][2].GetDouble());

        mAxisOfRoationVector.push_back(m_parameters["axis_of_rotation"][0].GetDouble());
        mAxisOfRoationVector.push_back(m_parameters["axis_of_rotation"][1].GetDouble());
        mAxisOfRoationVector.push_back(m_parameters["axis_of_rotation"][2].GetDouble());

        // normalizing the axis of roatation
        double norm = 0.0;
        for (std::size_t d = 0; d < 3; ++d)
            norm += mAxisOfRoationVector[d] * mAxisOfRoationVector[d];
        norm = sqrt(norm);
        for (std::size_t d = 0; d < 3; ++d)
            mAxisOfRoationVectorNormalized.push_back(mAxisOfRoationVector[d] / norm);

        mTheta = 0.0;
    }

    /**
		Activates the constraint set or deactivates
		@arg isActive true/false
		*/
    void SetActive(bool isActive = true)
    {
        mIsActive = isActive;
    }

    /**
		Sets the name of the constraint set
		@arg name
		*/
    void SetName(std::string name)
    {
        this->mName = name;
    }

    void SetCentreOfRotation(double x, double y, double z)
    {

        mCenterOfRotation[0] = x;
        mCenterOfRotation[1] = y;
        mCenterOfRotation[2] = z;
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
    
    void ChangeAngularVelocity(double ang_velocity) 
    {
        m_angular_velocity_radians = ang_velocity;

    }

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        ModelPart &sub_model_part = mr_model_part.GetSubModelPart(m_sub_model_part_name);
        double delta_t;
        int domain_size;

        ProcessInfoPointerType info = mr_model_part.pGetProcessInfo();
        delta_t = (*info)[DELTA_TIME];
        domain_size = (*info)[DOMAIN_SIZE];
        mTheta += m_angular_velocity_radians * delta_t;

#pragma omp parallel for
        for (NodeIteratorType node = sub_model_part.NodesBegin(); node < sub_model_part.NodesEnd(); node++)
        {

            /// Calculating the displacement of the current node
            std::vector<double> node_cords = {node->X0(), node->Y0(), node->Z0()};
            std::vector<double> rotatedNode = RotateNode2(node_cords, mTheta);

            node->X() = rotatedNode[0];
            node->Y() = rotatedNode[1];
            if (domain_size > 2)
                node->Z() = rotatedNode[2];

            // This is done for visualization
            /*node->FastGetSolutionStepValue(MESH_DISPLACEMENT_X) = node->X() - node->X0();
            node->FastGetSolutionStepValue(MESH_DISPLACEMENT_Y) = node->Y() - node->Y0();
            if (domain_size > 2)
                node->FastGetSolutionStepValue(MESH_DISPLACEMENT_Z) = node->Z() - node->Z0();*/

            // Computing the linear velocity at this node
            std::vector<double> r(3);
            std::vector<double> linearVelocity(3);
            r[0] = node->X() - mCenterOfRotation[0];
            r[1] = node->Y() - mCenterOfRotation[1];
            r[2] = node->Z() - mCenterOfRotation[2];
            CrossProduct(mAxisOfRoationVector, r, linearVelocity);

            if (m_parameters["is_ale"].GetBool())
            {
                node->FastGetSolutionStepValue(MESH_VELOCITY_X, 0) += m_angular_velocity_radians * linearVelocity[0];
                node->FastGetSolutionStepValue(MESH_VELOCITY_Y, 0) += m_angular_velocity_radians * linearVelocity[1];
                if (domain_size > 2)
                    node->FastGetSolutionStepValue(MESH_VELOCITY_Z, 0) += m_angular_velocity_radians * linearVelocity[2];
            }

            if (node->IsFixed(VELOCITY_X))
                node->FastGetSolutionStepValue(VELOCITY_X, 0) = node->FastGetSolutionStepValue(MESH_VELOCITY_X, 0);

            if (node->IsFixed(VELOCITY_Y))
                node->FastGetSolutionStepValue(VELOCITY_Y, 0) = node->FastGetSolutionStepValue(MESH_VELOCITY_Y, 0);

            if (domain_size > 2)
                if (node->IsFixed(VELOCITY_Z))
                    node->FastGetSolutionStepValue(VELOCITY_Z, 0) = node->FastGetSolutionStepValue(MESH_VELOCITY_Z, 0);
        }

        /*if (m_parameters["is_ale"].GetBool())
        {

            ModelPart &sub_model_part_boundary = mr_model_part.GetSubModelPart(m_sub_model_part_boundary_name);
#pragma omp parallel for
            for (NodeIteratorType node = sub_model_part_boundary.NodesBegin(); node < sub_model_part_boundary.NodesEnd(); node++)
            {
                node->FastGetSolutionStepValue(MESH_VELOCITY_X, 0) = 0;
                node->FastGetSolutionStepValue(MESH_VELOCITY_Y, 0) = 0;
                if (domain_size > 2)
                    node->FastGetSolutionStepValue(MESH_VELOCITY_Z, 0) = 0;
            }
        }*/

        KRATOS_CATCH("");
    }

    void ExecuteFinalizeSolutionStep() override
    {
        ProcessInfoPointerType info = mr_model_part.pGetProcessInfo();

        int domain_size = (*info)[DOMAIN_SIZE];

        ModelPart &sub_model_part = mr_model_part.GetSubModelPart(m_sub_model_part_name);
#pragma omp parallel for
        for (NodeIteratorType node = sub_model_part.NodesBegin(); node < sub_model_part.NodesEnd(); node++)
        {

            if (m_parameters["is_ale"].GetBool())
            {
                node->FastGetSolutionStepValue(MESH_VELOCITY_X, 0) = 0.0;
                node->FastGetSolutionStepValue(MESH_VELOCITY_Y, 0) = 0.0;
                if (domain_size > 2)
                    node->FastGetSolutionStepValue(MESH_VELOCITY_Z, 0) = 0.0;
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
    ModelPart &mr_model_part;
    Parameters m_parameters;
    std::string m_sub_model_part_name;
    std::string m_sub_model_part_boundary_name;
    double m_angular_velocity_radians;
    std::string mName;
    bool mIsActive;
    std::vector<double> mCenterOfRotation;
    std::vector<double> mAxisOfRoationVector;
    std::vector<double> mAxisOfRoationVectorNormalized;
    double mTheta;

  private:
    /// Assignment operator.
    RotateRegionProcess &operator=(RotateRegionProcess const &rOther) { return *this; }

    /*
     * Calculates the cross product of two vectors
     *  c = axb
     */
    void CrossProduct(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c)
    {
        assert(a.size() == b.size());
        c[0] = a[1] * b[2] - a[2] * b[1];
        c[1] = a[2] * b[0] - a[0] * b[2];
        c[2] = a[0] * b[1] - a[1] * b[0];
    }

    /*
     * Rotates a given point(node_cords) in space around a given mAxisOfRoationVector by an angle thetha
     */
    std::vector<double> RotateNode(std::vector<double> node_cords, double theta)
    {
        std::vector<double> rotatedNode(4);
        std::vector<double> U(3); // normalized axis of rotation
        node_cords.push_back(1.0);
        
        // normalizing the axis of roatation
        double norm = 0.0;
        for (std::size_t d = 0; d < 3; ++d)
            norm += mAxisOfRoationVector[d] * mAxisOfRoationVector[d];
        norm = sqrt(norm);
        for (std::size_t d = 0; d < 3; ++d)
            U[d] = mAxisOfRoationVector[d] / norm;

        // Constructing the transformation matrix
        double A0[4][4];
        double x1 = mCenterOfRotation[0];
        double y1 = mCenterOfRotation[1];
        double z1 = mCenterOfRotation[2];

        double a = U[0];
        double b = U[1];
        double c = U[2];

        double t2 = cos(theta);
        double t3 = sin(theta);
        double t4 = a * a;
        double t5 = b * b;
        double t6 = c * c;
        double t7 = a * b;
        double t8 = t5 + t6;
        double t9 = 1.0 / t8;
        double t10 = a * c;
        double t11 = b * t3;
        double t12 = a * t3 * t5;
        double t13 = a * t3 * t6;
        double t14 = b * c * t2;
        A0[0][0] = t4 + t2 * t8;
        A0[0][1] = t7 - c * t3 - a * b * t2;
        A0[0][2] = t10 + t11 - a * c * t2;
        A0[0][3] = x1 - t4 * x1 - a * b * y1 - a * c * z1 - b * t3 * z1 + c * t3 * y1 - t2 * t5 * x1 - t2 * t6 * x1 + a * b * t2 * y1 + a * c * t2 * z1;
        A0[1][0] = t7 + c * t3 - a * b * t2;
        A0[1][1] = t9 * (t2 * t6 + t5 * t8 + t2 * t4 * t5);
        A0[1][2] = -t9 * (t12 + t13 + t14 - b * c * t8 - b * c * t2 * t4);
        A0[1][3] = -t9 * (-t8 * y1 + t2 * t6 * y1 + t5 * t8 * y1 + a * b * t8 * x1 - b * c * t2 * z1 + b * c * t8 * z1 - a * t3 * t5 * z1 - a * t3 * t6 * z1 + c * t3 * t8 * x1 + t2 * t4 * t5 * y1 - a * b * t2 * t8 * x1 + b * c * t2 * t4 * z1);
        A0[2][0] = t10 - t11 - a * c * t2;
        A0[2][1] = t9 * (t12 + t13 - t14 + b * c * t8 + b * c * t2 * t4);
        A0[2][2] = t9 * (t2 * t5 + t6 * t8 + t2 * t4 * t6);
        A0[2][3] = -t9 * (-t8 * z1 + t2 * t5 * z1 + t6 * t8 * z1 + a * c * t8 * x1 - b * c * t2 * y1 + b * c * t8 * y1 + a * t3 * t5 * y1 + a * t3 * t6 * y1 - b * t3 * t8 * x1 + t2 * t4 * t6 * z1 - a * c * t2 * t8 * x1 + b * c * t2 * t4 * y1);
        A0[3][3] = 1.0;

        // Multiplying the point to get the rotated point
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                rotatedNode[i] += (A0[i][j] * node_cords[j]);
            }
        }

        return rotatedNode;
    }

    /*
   Rotate a point p by angle theta around an arbitrary axis r
   Return the rotated point.
   Positive angles are anticlockwise looking down the axis
   towards the origin.
   Assume right hand coordinate system.
*/
    std::vector<double> RotateNode2(std::vector<double> p, double theta)
    {
        std::vector<double> q(3);
        double costheta, sintheta;

        std::vector<double> r = mAxisOfRoationVectorNormalized;
        costheta = cos(theta);
        sintheta = sin(theta);

        p[0] -= mCenterOfRotation[0];
        p[1] -= mCenterOfRotation[1];
        p[2] -= mCenterOfRotation[2];

        q[0] += (costheta + (1 - costheta) * r[0] * r[0]) * p[0];
        q[0] += ((1 - costheta) * r[0] * r[1] - r[2] * sintheta) * p[1];
        q[0] += ((1 - costheta) * r[0] * r[2] + r[1] * sintheta) * p[2];

        q[1] += ((1 - costheta) * r[0] * r[1] + r[2] * sintheta) * p[0];
        q[1] += (costheta + (1 - costheta) * r[1] * r[1]) * p[1];
        q[1] += ((1 - costheta) * r[1] * r[2] - r[0] * sintheta) * p[2];

        q[2] += ((1 - costheta) * r[0] * r[2] - r[1] * sintheta) * p[0];
        q[2] += ((1 - costheta) * r[1] * r[2] + r[0] * sintheta) * p[1];
        q[2] += (costheta + (1 - costheta) * r[2] * r[2]) * p[2];

        q[0] += mCenterOfRotation[0];
        q[1] += mCenterOfRotation[1];
        q[2] += mCenterOfRotation[2];

        return (q);
    }

}; // Class MoveRotorProcess
}; // namespace Kratos.

#endif // KRATOS_MOVE_ROTOR_PROCESS_H
