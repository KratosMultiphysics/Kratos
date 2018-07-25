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
//
//

#ifndef APPLY_PERIODIC_CONDITION_PROCESS_H
#define APPLY_PERIODIC_CONDITION_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/point.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "utilities/binbased_fast_point_locator.h"
#include "spatial_containers/spatial_containers.h"


namespace Kratos
{

class ApplyPeriodicConditionProcess : public Process
{

  public:
    /// Pointer definition of ApplyPeriodicConditionProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyPeriodicConditionProcess);

    typedef Dof<double> *DofPointerType;
    typedef Dof<double> DofType;

    typedef ModelPart::VariableComponentType VariableComponentType;
    typedef KratosComponents<Variable<array_1d<double, 3>>> VectorVariableType;
    typedef ProcessInfo ProcessInfoType;
    typedef ProcessInfo::Pointer ProcessInfoPointerType;
    typedef unsigned int IndexType;

    typedef ModelPart::DoubleVariableType VariableType;
    typedef ModelPart::NodeIterator NodeIterator;
    typedef Element ElementType;
    typedef Node<3> NodeType;

    /// Constructor.
    ApplyPeriodicConditionProcess(ModelPart &model_part,
                                  Parameters rParameters) : Process(Flags()), mrMainModelPart(model_part), m_parameters(rParameters)
    {
        Parameters default_parameters(R"(
                                            {
                                            "constraint_set_name":"default",
                                            "master_sub_model_part_name":"default_master",
                                            "slave_sub_model_part_name":"default_slave",
                                            "variable_names":[],
                                            "center":[0,0,0],
                                            "axis_of_rotation":[0.0,0.0,0.0],
                                            "angle":0.0,
                                            "flip":false,
                                            "dir_of_translation":[0.0,0.0,0.0],
                                            "magnitude":0.0
                                            }  )");

        // Initializing
        m_parameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        // Use the salve model part and newly created model part to form MPC constraints
        mSlaveSubModelPartName = m_parameters["slave_sub_model_part_name"].GetString();
        mMasterSubModelPartName = m_parameters["master_sub_model_part_name"].GetString();

        mCenterOfRotation.push_back(m_parameters["center"][0].GetDouble());
        mCenterOfRotation.push_back(m_parameters["center"][1].GetDouble());
        mCenterOfRotation.push_back(m_parameters["center"][2].GetDouble());

        mAxisOfRoationVector.push_back(m_parameters["axis_of_rotation"][0].GetDouble());
        mAxisOfRoationVector.push_back(m_parameters["axis_of_rotation"][1].GetDouble());
        mAxisOfRoationVector.push_back(m_parameters["axis_of_rotation"][2].GetDouble());

        mModulus = m_parameters["magnitude"].GetDouble();
        mDirOfTranslation.push_back(m_parameters["dir_of_translation"][0].GetDouble());
        mDirOfTranslation.push_back(m_parameters["dir_of_translation"][1].GetDouble());
        mDirOfTranslation.push_back(m_parameters["dir_of_translation"][2].GetDouble());

        if (m_parameters["flip"].GetBool())
            mSign = -1;
        else
            mSign = 1;

        // normalizing the axis of roatation
        double norm = 0.0;
        for (unsigned int d = 0; d < 3; ++d)
            norm += mAxisOfRoationVector[d] * mAxisOfRoationVector[d];
        norm = sqrt(norm);
        for (unsigned int d = 0; d < 3; ++d)
            mAxisOfRoationVectorNormalized.push_back(mAxisOfRoationVector[d] / norm);

        mTheta = m_parameters["angle"].GetDouble() * 2 * 3.1416 / 360.0;

        mTransformationMatrix.resize(4);
        for (auto &it : mTransformationMatrix)
        {
            it.resize(4, 0.0);
        }

        if (mTheta == 0 && mModulus != 0)
            mType = "translation";
        else if (mTheta != 0.0 && mModulus == 0.0)
            mType = "rotation";
        else if (mTheta == 0.0 && mModulus == 0.0)
            KRATOS_THROW_ERROR(std::runtime_error, "Both angle of rotation and modulus of translation cannot be zero. Please check the input", "");

        CalculateTransformationMatrix();
        mIsInitialized = false;
    }
    ~ApplyPeriodicConditionProcess()
    {
    }

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        if (!mIsInitialized)
        {
            ModelPart &masterModelPart = mrMainModelPart.GetSubModelPart(m_parameters["master_sub_model_part_name"].GetString());

            ProcessInfoPointerType info = mrMainModelPart.pGetProcessInfo();
            int probDim = info->GetValue(DOMAIN_SIZE);
            // Rotate the master so it goes to the slave
            if (probDim == 2)
            {
                GetRotatedMaster<2>(masterModelPart);
                ApplyConstraintsForPeriodicConditions<2>();
            }
            else if (probDim == 3)
            {
                GetRotatedMaster<3>(masterModelPart);
                ApplyConstraintsForPeriodicConditions<3>();
            }
            mIsInitialized = true;
        }

        // Once master and slave nodes are on the same surface we interpolate and appply constraints
        KRATOS_CATCH("");
    }

  private:
    std::vector<std::vector<double>> mTransformationMatrix;
    ModelPart &mrMainModelPart;
    Parameters m_parameters;
    ModelPart mpRotatedMasterModelPart; // This is new modelpart
    std::string mSlaveSubModelPartName;
    std::string mMasterSubModelPartName;
    double mTheta;
    bool mIsInitialized;
    std::vector<double> mCenterOfRotation;
    std::vector<double> mAxisOfRoationVector;
    std::vector<double> mAxisOfRoationVectorNormalized;
    int mSign;
    std::string mType;
    double mModulus;
    std::vector<double> mDirOfTranslation;

    template <int TDim>
    void ApplyConstraintsForPeriodicConditions()
    {
    }



    void CalculateTransformationMatrix()
    {

        if (mType == "translation")
            CalculateTranslationMatrix();
        else if (mType == "rotation")
            CalculateRoatationMatrix();
    }

    void CalculateTranslationMatrix()
    {
        mTransformationMatrix[0][0] = 1;
        mTransformationMatrix[0][1] = 0;
        mTransformationMatrix[0][2] = 0;
        mTransformationMatrix[0][3] = mModulus * mDirOfTranslation[0];
        mTransformationMatrix[1][0] = 0;
        mTransformationMatrix[1][1] = 1.0;
        mTransformationMatrix[1][2] = 0;
        mTransformationMatrix[1][3] = mModulus * mDirOfTranslation[1];
        mTransformationMatrix[2][0] = 0;
        mTransformationMatrix[2][1] = 0;
        mTransformationMatrix[2][2] = 1.0;
        mTransformationMatrix[2][3] = mModulus * mDirOfTranslation[2];
        mTransformationMatrix[3][3] = 1.0;
    }

    void CalculateRoatationMatrix()
    {
        std::vector<double> U(3); // normalized axis of rotation
        // normalizing the axis of roatation
        double norm = 0.0;
        for (unsigned int d = 0; d < 3; ++d)
            norm += mAxisOfRoationVector[d] * mAxisOfRoationVector[d];
        norm = sqrt(norm);
        for (unsigned int d = 0; d < 3; ++d)
            U[d] = mAxisOfRoationVector[d] / norm;

        // Constructing the transformation matrix
        double x1 = mCenterOfRotation[0];
        double y1 = mCenterOfRotation[1];
        double z1 = mCenterOfRotation[2];

        double a = U[0];
        double b = U[1];
        double c = U[2];
        //double eps = 1e-12;

        double t2 = cos(mTheta);
        double t3 = sin(mTheta);
        double t4 = a * a;
        double t5 = b * b;
        double t6 = c * c;
        double t7 = a * b;
        double t8 = t5 + t6;
        double t9 = 1.0 / t8;
        if(isnan(t9))
            t9 = 0.0;
        double t10 = a * c;
        double t11 = b * t3;
        double t12 = a * t3 * t5;
        double t13 = a * t3 * t6;
        double t14 = b * c * t2;
        mTransformationMatrix[0][0] = t4 + t2 * t8;
        mTransformationMatrix[0][1] = t7 - c * t3 - a * b * t2;
        mTransformationMatrix[0][2] = t10 + t11 - a * c * t2;
        mTransformationMatrix[0][3] = x1 - t4 * x1 - a * b * y1 - a * c * z1 - b * t3 * z1 + c * t3 * y1 - t2 * t5 * x1 - t2 * t6 * x1 + a * b * t2 * y1 + a * c * t2 * z1;
        mTransformationMatrix[1][0] = t7 + c * t3 - a * b * t2;
        mTransformationMatrix[1][1] = t9 * (t2 * t6 + t5 * t8 + t2 * t4 * t5);
        mTransformationMatrix[1][2] = -t9 * (t12 + t13 + t14 - b * c * t8 - b * c * t2 * t4);
        mTransformationMatrix[1][3] = -t9 * (-t8 * y1 + t2 * t6 * y1 + t5 * t8 * y1 + a * b * t8 * x1 - b * c * t2 * z1 + b * c * t8 * z1 - a * t3 * t5 * z1 - a * t3 * t6 * z1 + c * t3 * t8 * x1 + t2 * t4 * t5 * y1 - a * b * t2 * t8 * x1 + b * c * t2 * t4 * z1);
        mTransformationMatrix[2][0] = t10 - t11 - a * c * t2;
        mTransformationMatrix[2][1] = t9 * (t12 + t13 - t14 + b * c * t8 + b * c * t2 * t4);
        mTransformationMatrix[2][2] = t9 * (t2 * t5 + t6 * t8 + t2 * t4 * t6);
        mTransformationMatrix[2][3] = -t9 * (-t8 * z1 + t2 * t5 * z1 + t6 * t8 * z1 + a * c * t8 * x1 - b * c * t2 * y1 + b * c * t8 * y1 + a * t3 * t5 * y1 + a * t3 * t6 * y1 - b * t3 * t8 * x1 + t2 * t4 * t6 * z1 - a * c * t2 * t8 * x1 + b * c * t2 * t4 * y1);
        mTransformationMatrix[3][3] = 1.0;
    }



};

}