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
    typedef Matrix MatrixType;
    typedef Vector VectorType;

    ApplyPeriodicConditionProcess(ModelPart &model_part,
                                  Parameters rParameters) : Process(Flags()), mrMainModelPart(model_part), mParameters(rParameters)
    {
        Parameters default_parameters(R"(
                                            {
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
        mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        // Use the salve model part and newly created model part to form MPC constraints
        mSlaveSubModelPartName = mParameters["slave_sub_model_part_name"].GetString();
        mMasterSubModelPartName = mParameters["master_sub_model_part_name"].GetString();

        mCenterOfRotation.push_back(mParameters["center"][0].GetDouble());
        mCenterOfRotation.push_back(mParameters["center"][1].GetDouble());
        mCenterOfRotation.push_back(mParameters["center"][2].GetDouble());

        mAxisOfRoationVector.push_back(mParameters["axis_of_rotation"][0].GetDouble());
        mAxisOfRoationVector.push_back(mParameters["axis_of_rotation"][1].GetDouble());
        mAxisOfRoationVector.push_back(mParameters["axis_of_rotation"][2].GetDouble());

        mModulus = mParameters["magnitude"].GetDouble();
        mDirOfTranslation.push_back(mParameters["dir_of_translation"][0].GetDouble());
        mDirOfTranslation.push_back(mParameters["dir_of_translation"][1].GetDouble());
        mDirOfTranslation.push_back(mParameters["dir_of_translation"][2].GetDouble());

        if (mParameters["flip"].GetBool())
            mSign = -1;
        else
            mSign = 1;

        mTheta = mParameters["angle"].GetDouble() * 2 * 3.1416 / 360.0;

        mTransformationMatrix.resize(4,4);

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
            ProcessInfoPointerType info = mrMainModelPart.pGetProcessInfo();
            int probDim = info->GetValue(DOMAIN_SIZE);
            // Rotate the master so it goes to the slave
            if (probDim == 2)
            {
                ApplyConstraintsForPeriodicConditions<2>();
            }
            else if (probDim == 3)
            {
                ApplyConstraintsForPeriodicConditions<3>();
            }
            mIsInitialized = true;
        }

        // Once master and slave nodes are on the same surface we interpolate and appply constraints
        KRATOS_CATCH("");
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream <<"ApplyPeriodicConditionProcess Process "<<std::endl;
    }

  private:
    MatrixType mTransformationMatrix; // This can be either for rotating or for translating the master or slave
    ModelPart &mrMainModelPart;       // the MainModelPart to which the master-slave constraints are added.
    Parameters mParameters;          // parameters
    ModelPart mpRotatedMasterModelPart; // This is new modelpart
    std::string mSlaveSubModelPartName;
    std::string mMasterSubModelPartName;
    double mTheta;
    bool mIsInitialized;
    std::vector<double> mCenterOfRotation;
    std::vector<double> mAxisOfRoationVector;
    int mSign;
    std::string mType;
    double mModulus;
    std::vector<double> mDirOfTranslation;

    /**
     * @brief   The function to figure out how the master and slave model parts relate together and add master-slave constraints 
     *          to the rModelPart
     */
    template <int TDim>
    void ApplyConstraintsForPeriodicConditions()
    {

        ModelPart &slave_model_part = mrMainModelPart.GetSubModelPart(mParameters["slave_sub_model_part_name"].GetString());
        ModelPart &master_model_part = mrMainModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
        int num_vars = mParameters["variable_names"].size();
        BinBasedFastPointLocator<TDim> bin_based_point_locator(master_model_part);
        bin_based_point_locator.UpdateSearchDatabase();

        // for bin based point locator
		VectorType shape_function_values;
		const int max_results = 1000;
		typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        Element::Pointer p_host_elem;
        typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

        for (int j = 0; j < num_vars; j++)
        {
            std::string var_name = mParameters["variable_names"][j].GetString();
            if (KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Has(var_name + "_X"))
            {   // Checking if the variable is a vector variable
                VariableComponentType r_var_x = KratosComponents<VariableComponentType>::Get(var_name + std::string("_X"));
                VariableComponentType r_var_y = KratosComponents<VariableComponentType>::Get(var_name + std::string("_Y"));
                VariableComponentType r_var_z = KratosComponents<VariableComponentType>::Get(var_name + std::string("_Z"));

                for (auto& slave_node : slave_model_part.Nodes())
                {
                    // Finding the host element for this node
                	bool is_found = false;
				    is_found = bin_based_point_locator.FindPointOnMesh(slave_node.Coordinates(), shape_function_values, p_host_elem, result_begin, max_results);
                    if(is_found)
                    {
                        auto& geometry = p_host_elem->GetGeometry();
                        IndexType master_index = 0;
                        for (auto& master_node : geometry)
                        {
                            int current_num_constraint = mrMainModelPart.NumberOfMasterSlaveConstraints();

                            double master_weight = mSign * shape_function_values(master_index);
                            double constant_x(0.0), constant_y(0.0), constant_z(0.0);

                            constant_x = master_weight * mTransformationMatrix(0,3);
                            constant_y = master_weight * mTransformationMatrix(1,3);
                            constant_z = master_weight * mTransformationMatrix(2,3);

                            mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_x, slave_node, r_var_x, master_weight * mTransformationMatrix(0,0), constant_x);
                            mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_y, slave_node, r_var_x, master_weight * mTransformationMatrix(0,1), constant_x);
                            mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_z, slave_node, r_var_x, master_weight * mTransformationMatrix(0,2), constant_x);

                            mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_x, slave_node, r_var_y, master_weight * mTransformationMatrix(1,0), constant_y);
                            mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_y, slave_node, r_var_y, master_weight * mTransformationMatrix(1,1), constant_y);
                            mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_z, slave_node, r_var_y, master_weight * mTransformationMatrix(1,2), constant_y);

                            if (TDim == 3)
                            {
                                mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_x, slave_node, r_var_z, master_weight * mTransformationMatrix(2,0), constant_z);
                                mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_y, slave_node, r_var_z, master_weight * mTransformationMatrix(2,1), constant_z);
                                mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_z, slave_node, r_var_z, master_weight * mTransformationMatrix(2,2), constant_z);
                            }

                            master_index++;
                        }
                    }
                }
            }
        }

    }

    /**
     * @brief   Fuctions calculates the transformation matrix to account for the moving the two periodic condition modelparts together.
     */
    void CalculateTransformationMatrix()
    {
        if (mType == "translation")
            CalculateTranslationMatrix();
        else if (mType == "rotation")
            CalculateRoatationMatrix();
    }

    /**
     * @brief   Calculates the transformation matrix which translates the given vector alone mDirOfTranslation by mModulus
     */
    void CalculateTranslationMatrix()
    {
        mTransformationMatrix(0,0) = 1;
        mTransformationMatrix(0,1) = 0;
        mTransformationMatrix(0,2) = 0;
        mTransformationMatrix(0,3) = mModulus * mDirOfTranslation[0];
        mTransformationMatrix(1,0) = 0;
        mTransformationMatrix(1,1) = 1.0;
        mTransformationMatrix(1,2) = 0;
        mTransformationMatrix(1,3) = mModulus * mDirOfTranslation[1];
        mTransformationMatrix(2,0) = 0;
        mTransformationMatrix(2,1) = 0;
        mTransformationMatrix(2,2) = 1.0;
        mTransformationMatrix(2,3) = mModulus * mDirOfTranslation[2];
        mTransformationMatrix(3,3) = 1.0;
    }

    /**
     * @brief   Calculates the transformation matrix which rotates the given vector around mAxisOfRoationVector and mCenterOfRotation
     *          by mTheta
     */
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
        mTransformationMatrix(0,0) = t4 + t2 * t8;
        mTransformationMatrix(0,1) = t7 - c * t3 - a * b * t2;
        mTransformationMatrix(0,2) = t10 + t11 - a * c * t2;
        mTransformationMatrix(0,3) = x1 - t4 * x1 - a * b * y1 - a * c * z1 - b * t3 * z1 + c * t3 * y1 - t2 * t5 * x1 - t2 * t6 * x1 + a * b * t2 * y1 + a * c * t2 * z1;
        mTransformationMatrix(1,0) = t7 + c * t3 - a * b * t2;
        mTransformationMatrix(1,1) = t9 * (t2 * t6 + t5 * t8 + t2 * t4 * t5);
        mTransformationMatrix(1,2) = -t9 * (t12 + t13 + t14 - b * c * t8 - b * c * t2 * t4);
        mTransformationMatrix(1,3) = -t9 * (-t8 * y1 + t2 * t6 * y1 + t5 * t8 * y1 + a * b * t8 * x1 - b * c * t2 * z1 + b * c * t8 * z1 - a * t3 * t5 * z1 - a * t3 * t6 * z1 + c * t3 * t8 * x1 + t2 * t4 * t5 * y1 - a * b * t2 * t8 * x1 + b * c * t2 * t4 * z1);
        mTransformationMatrix(2,0) = t10 - t11 - a * c * t2;
        mTransformationMatrix(2,1) = t9 * (t12 + t13 - t14 + b * c * t8 + b * c * t2 * t4);
        mTransformationMatrix(2,2) = t9 * (t2 * t5 + t6 * t8 + t2 * t4 * t6);
        mTransformationMatrix(2,3) = -t9 * (-t8 * z1 + t2 * t5 * z1 + t6 * t8 * z1 + a * c * t8 * x1 - b * c * t2 * y1 + b * c * t8 * y1 + a * t3 * t5 * y1 + a * t3 * t6 * y1 - b * t3 * t8 * x1 + t2 * t4 * t6 * z1 - a * c * t2 * t8 * x1 + b * c * t2 * t4 * y1);
        mTransformationMatrix(3,3) = 1.0;
    }

};

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyPeriodicConditionProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyPeriodicConditionProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_APPLY_CONSTANT_VECTORVALUE_PROCESS_H_INCLUDED  defined