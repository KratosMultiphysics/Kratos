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

// System includes
#include <string>
#include <iostream>
#include <chrono>

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "processes/apply_periodic_boundary_condition_process.h"
#include "utilities/binbased_fast_point_locator_conditions.h"

namespace Kratos
{

    ApplyPeriodicConditionProcess::ApplyPeriodicConditionProcess(ModelPart &rMasterModelPart, ModelPart &rSlaveModelPart,
                                  Parameters rParameters) : Process(Flags()), mrMasterModelPart(rMasterModelPart),
                                   mrSlaveModelPart(rSlaveModelPart), mParameters(rParameters)
    {
        Parameters default_parameters(R"(
                                            {
                                                "variable_names":[],
                                                "transformation_settings":{
                                                    "rotation_settings":{
                                                        "center":[0,0,0],
                                                        "axis_of_rotation":[0.0,0.0,0.0],
                                                        "angle_degree":0.0
                                                    },
                                                    "translation_settings":{
                                                        "dir_of_translation":[0.0,0.0,0.0],
                                                        "magnitude":0.0
                                                    }
                                                },
                                                "search_settings":{
                                                    "max_results":100000,
                                                    "tolerance": 1E-6
                                                }
                                            }  )");

        // Initializing
        mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mCenterOfRotation = mParameters["transformation_settings"]["rotation_settings"]["center"].GetVector();
        mAxisOfRotationVector = mParameters["transformation_settings"]["rotation_settings"]["axis_of_rotation"].GetVector();
        mDirOfTranslation = mParameters["transformation_settings"]["translation_settings"]["dir_of_translation"].GetVector();

        mMagnitude = mParameters["transformation_settings"]["translation_settings"]["magnitude"].GetDouble();
        mAngleOfRotation = mParameters["transformation_settings"]["rotation_settings"]["angle_degree"].GetDouble() * 2 * Globals::Pi / 360.0;

        mTransformationMatrix.resize(4,4);
        mTransformationMatrixVariable.resize(4,4);

        mSearchMaxResults = mParameters["search_settings"]["max_results"].GetInt();
        mSearchTolerance = mParameters["search_settings"]["tolerance"].GetDouble();

        RemoveCommonNodesFromSlaveModelPart();

        if (mAngleOfRotation == 0 && mMagnitude != 0)
            mType = "translation";
        else if (mAngleOfRotation != 0.0 && mMagnitude == 0.0)
            mType = "rotation";
        else
            KRATOS_ERROR_IF(mAngleOfRotation == 0.0 && mMagnitude == 0.0)<<"Both angle of rotation and modulus of translation cannot be zero. Please check the input"<<std::endl;
        else
            KRATOS_ERROR_IF(mAngleOfRotation != 0.0 && mMagnitude != 0.0)<<"Both angle of rotation and modulus of translation cannot be specified at the same time. Please check the input"<<std::endl;

        CalculateTransformationMatrix();
        mIsInitialized = false;
    }

    ApplyPeriodicConditionProcess::~ApplyPeriodicConditionProcess()
    {
    }

    /**
     * @brief  Function to remove the common nodes of slave and master modelparts from the slave modelpart
     */
    void ApplyPeriodicConditionProcess::RemoveCommonNodesFromSlaveModelPart()
    {
        const int num_slave_nodes = mrSlaveModelPart.NumberOfNodes();
        const NodeIteratorType it_slave_node_begin = mrSlaveModelPart.NodesBegin();
        for(int i_node = 0; i_node<num_slave_nodes; ++i_node)
        {
            NodeIteratorType it_slave_node = it_slave_node_begin;
            std::advance(it_slave_node, i_node);
            auto iterator = mrMasterModelPart.Nodes().find(it_slave_node->Id());
            if(iterator != mrMasterModelPart.NodesEnd()){
                it_slave_node->Set(TO_ERASE);
            }
        }
        mrSlaveModelPart.RemoveNodes(TO_ERASE);
    }

    /**
     * @brief Function initializes the process
     */
    void ApplyPeriodicConditionProcess::ExecuteInitialize()
    {
        KRATOS_TRY;
        if (!mIsInitialized)
        {
            const ProcessInfoPointerType info = mrMasterModelPart.pGetProcessInfo();
            const int prob_dim = info->GetValue(DOMAIN_SIZE);
            // Rotate the master so it goes to the slave
            if (prob_dim == 2){
                ApplyConstraintsForPeriodicConditions<2>();
            }
            else if (prob_dim == 3){
                ApplyConstraintsForPeriodicConditions<3>();
            }
            mIsInitialized = true;
        }
        std::cout<<"Initialization of Periodic conditions finished .. !"<<std::endl;

        KRATOS_CATCH("");
    }

    /**
     * @brief Function initializes the solution step
     */
    void ApplyPeriodicConditionProcess::ExecuteInitializeSolutionStep()
    {
    }

    /**
     * @brief Function to print the information about this current process
     */
    void ApplyPeriodicConditionProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream <<"ApplyPeriodicConditionProcess Process "<<std::endl;
    }

    template <int TDim>
    void ApplyPeriodicConditionProcess::ApplyConstraintsForPeriodicConditions()
    {
        const double start_apply = OpenMPUtils::GetCurrentTime();
        const int num_vars = mParameters["variable_names"].size();
        BinBasedFastPointLocatorConditions<TDim> bin_based_point_locator(mrMasterModelPart);
        bin_based_point_locator.UpdateSearchDatabase();

        const int num_slave_nodes = mrSlaveModelPart.NumberOfNodes();
        const NodeIteratorType it_slave_node_begin = mrSlaveModelPart.NodesBegin();

        unsigned int num_slaves_found = 0;

        #pragma omp parallel for schedule(guided, 512) reduction( + : num_slaves_found )
        for(int i_node = 0; i_node<num_slave_nodes; ++i_node)
        {
            Condition::Pointer p_host_cond;
            VectorType shape_function_values;
            NodeIteratorType it_slave_node = it_slave_node_begin;
            std::advance(it_slave_node, i_node);
            array_1d<double, 3 > transformed_slave_coordinates;
            TransformNode(it_slave_node->Coordinates(), transformed_slave_coordinates);

            // Finding the host element for this node
            const bool is_found = bin_based_point_locator.FindPointOnMeshSimplified(transformed_slave_coordinates, shape_function_values, p_host_cond, mSearchMaxResults, mSearchTolerance);
            if(is_found)
            {
                ++num_slaves_found;
                for (int j = 0; j < num_vars; j++)
                {
                    const std::string var_name = mParameters["variable_names"][j].GetString();
                    // Checking if the variable is a vector variable
                    if (KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Has(var_name + "_X"))
                    {   // TODO: Look for a better alternative to do this.
                        ConstraintSlaveNodeWithConditionForVectorVariable<TDim>(*it_slave_node, p_host_cond->GetGeometry() , shape_function_values, var_name);
                    } else if (KratosComponents<VariableType>::Has(var_name))
                    {
                        ConstraintSlaveNodeWithConditionForScalarVariable<TDim>(*it_slave_node, p_host_cond->GetGeometry() , shape_function_values, var_name);
                    }
                }
            }
        }
        KRATOS_WARNING_IF("",num_slaves_found != mrSlaveModelPart.NumberOfNodes())<<"Periodic condition cannot be applied for all the nodes."<<std::endl;
        const double end_apply = OpenMPUtils::GetCurrentTime();
        std::cout<<"Applying periodic boundary conditions took : "<<end_apply - start_apply<<" seconds." <<std::endl;
    }

    template <int TDim>
    void ApplyPeriodicConditionProcess::ConstraintSlaveNodeWithConditionForVectorVariable(NodeType& rSlaveNode, const GeometryType& rHostedGeometry, const VectorType& rWeights,
                                                                                            const std::string& rVarName )
    {
        const VariableComponentType& r_var_x = KratosComponents<VariableComponentType>::Get(rVarName + std::string("_X"));
        const VariableComponentType& r_var_y = KratosComponents<VariableComponentType>::Get(rVarName + std::string("_Y"));
        const VariableComponentType& r_var_z = KratosComponents<VariableComponentType>::Get(rVarName + std::string("_Z"));

        // Reference constraint
        const auto& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");

        IndexType master_index = 0;
        for (auto& master_node : rHostedGeometry)
        {
            const double master_weight = rWeights(master_index);

            const double constant_x = master_weight * mTransformationMatrixVariable(0,3);
            const double constant_y = master_weight * mTransformationMatrixVariable(1,3);
            const double constant_z = master_weight * mTransformationMatrixVariable(2,3);

            #pragma omp critical
            {
                int current_num_constraint = mrMasterModelPart.GetRootModelPart().NumberOfMasterSlaveConstraints();
                auto constraint1 = r_clone_constraint.Create(++current_num_constraint, master_node, r_var_x, rSlaveNode, r_var_x, master_weight * mTransformationMatrixVariable(0,0), constant_x);
                auto constraint2 = r_clone_constraint.Create(++current_num_constraint, master_node, r_var_y, rSlaveNode, r_var_x, master_weight * mTransformationMatrixVariable(0,1), constant_x);
                auto constraint3 = r_clone_constraint.Create(++current_num_constraint, master_node, r_var_z, rSlaveNode, r_var_x, master_weight * mTransformationMatrixVariable(0,2), constant_x);

                auto constraint4 = r_clone_constraint.Create(++current_num_constraint, master_node, r_var_x, rSlaveNode, r_var_y, master_weight * mTransformationMatrixVariable(1,0), constant_y);
                auto constraint5 = r_clone_constraint.Create(++current_num_constraint, master_node, r_var_y, rSlaveNode, r_var_y, master_weight * mTransformationMatrixVariable(1,1), constant_y);
                auto constraint6 = r_clone_constraint.Create(++current_num_constraint, master_node, r_var_z, rSlaveNode, r_var_y, master_weight * mTransformationMatrixVariable(1,2), constant_y);

                mrMasterModelPart.AddMasterSlaveConstraint(constraint1);
                mrMasterModelPart.AddMasterSlaveConstraint(constraint2);
                mrMasterModelPart.AddMasterSlaveConstraint(constraint3);
                mrMasterModelPart.AddMasterSlaveConstraint(constraint4);
                mrMasterModelPart.AddMasterSlaveConstraint(constraint5);
                mrMasterModelPart.AddMasterSlaveConstraint(constraint6);


                if (TDim == 3)
                {
                    auto constraint7 = r_clone_constraint.Create(++current_num_constraint, master_node, r_var_x, rSlaveNode, r_var_z, master_weight * mTransformationMatrixVariable(2,0), constant_z);
                    auto constraint8 = r_clone_constraint.Create(++current_num_constraint, master_node, r_var_y, rSlaveNode, r_var_z, master_weight * mTransformationMatrixVariable(2,1), constant_z);
                    auto constraint9 = r_clone_constraint.Create(++current_num_constraint, master_node, r_var_z, rSlaveNode, r_var_z, master_weight * mTransformationMatrixVariable(2,2), constant_z);

                    mrMasterModelPart.AddMasterSlaveConstraint(constraint7);
                    mrMasterModelPart.AddMasterSlaveConstraint(constraint8);
                    mrMasterModelPart.AddMasterSlaveConstraint(constraint9);
                }
            }

            master_index++;
        }
    }

    template <int TDim>
    void ApplyPeriodicConditionProcess::ConstraintSlaveNodeWithConditionForScalarVariable(NodeType& rSlaveNode, const GeometryType& rHostedGeometry, const VectorType& rWeights,
                                                                                            const std::string& rVarName )
    {
        const VariableType r_var = KratosComponents<VariableType>::Get(rVarName);

        // Reference constraint
        const auto& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");

        IndexType master_index = 0;
        for (auto& master_node : rHostedGeometry)
        {
            const double master_weight = rWeights(master_index);
            #pragma omp critical
            {
                int current_num_constraint = mrMasterModelPart.GetRootModelPart().NumberOfMasterSlaveConstraints();
                auto constraint = r_clone_constraint.Create(++current_num_constraint,master_node, r_var, rSlaveNode, r_var, master_weight, 0.0);
                mrMasterModelPart.AddMasterSlaveConstraint(constraint);
            }
            master_index++;
        }
    }

    void ApplyPeriodicConditionProcess::CalculateTransformationMatrix()
    {
        if (mType == "translation"){
            CalculateTranslationMatrix(-1*mMagnitude, mTransformationMatrix);
            CalculateTranslationMatrix(0.0, mTransformationMatrixVariable);
        }
        else if (mType == "rotation"){
            CalculateRotationMatrix(-1*mAngleOfRotation, mTransformationMatrix);
            CalculateRotationMatrix(mAngleOfRotation, mTransformationMatrixVariable);
        }
    }

    void ApplyPeriodicConditionProcess::CalculateTranslationMatrix(const double Modulus, MatrixType& rMatrix)
    {
        rMatrix(0,0) = 1.0;
        rMatrix(0,1) = 0.0;
        rMatrix(0,2) = 0.0;
        rMatrix(0,3) = Modulus *mDirOfTranslation[0];
        rMatrix(1,0) = 0.0;
        rMatrix(1,1) = 1.0;
        rMatrix(1,2) = 0.0;
        rMatrix(1,3) = Modulus *mDirOfTranslation[1];
        rMatrix(2,0) = 0.0;
        rMatrix(2,1) = 0.0;
        rMatrix(2,2) = 1.0;
        rMatrix(2,3) = Modulus *mDirOfTranslation[2];
        rMatrix(3,0) = 0.0;
        rMatrix(3,1) = 0.0;
        rMatrix(3,2) = 0.0;
        rMatrix(3,3) = 1.0;
    }

    void ApplyPeriodicConditionProcess::CalculateRotationMatrix(const double Theta, MatrixType& rMatrix )
    {
        std::vector<double> U(3); // normalized axis of rotation
        // normalizing the axis of rotation
        double norm = 0.0;
        for (unsigned int d = 0; d < 3; ++d)
            norm += mAxisOfRotationVector[d] * mAxisOfRotationVector[d];
        norm = sqrt(norm);
        for (unsigned int d = 0; d < 3; ++d)
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
        rMatrix(0,0) = t4 + t2 * t8;
        rMatrix(0,1) = t7 - c * t3 - a * b * t2;
        rMatrix(0,2) = t10 + t11 - a * c * t2;
        rMatrix(0,3) = x1 - t4 * x1 - a * b * y1 - a * c * z1 - b * t3 * z1 + c * t3 * y1 - t2 * t5 * x1 - t2 * t6 * x1 + a * b * t2 * y1 + a * c * t2 * z1;
        rMatrix(1,0) = t7 + c * t3 - a * b * t2;
        rMatrix(1,1) = t9 * (t2 * t6 + t5 * t8 + t2 * t4 * t5);
        rMatrix(1,2) = -t9 * (t12 + t13 + t14 - b * c * t8 - b * c * t2 * t4);
        rMatrix(1,3) = -t9 * (-t8 * y1 + t2 * t6 * y1 + t5 * t8 * y1 + a * b * t8 * x1 - b * c * t2 * z1 + b * c * t8 * z1 - a * t3 * t5 * z1 - a * t3 * t6 * z1 + c * t3 * t8 * x1 + t2 * t4 * t5 * y1 - a * b * t2 * t8 * x1 + b * c * t2 * t4 * z1);
        rMatrix(2,0) = t10 - t11 - a * c * t2;
        rMatrix(2,1) = t9 * (t12 + t13 - t14 + b * c * t8 + b * c * t2 * t4);
        rMatrix(2,2) = t9 * (t2 * t5 + t6 * t8 + t2 * t4 * t6);
        rMatrix(2,3) = -t9 * (-t8 * z1 + t2 * t5 * z1 + t6 * t8 * z1 + a * c * t8 * x1 - b * c * t2 * y1 + b * c * t8 * y1 + a * t3 * t5 * y1 + a * t3 * t6 * y1 - b * t3 * t8 * x1 + t2 * t4 * t6 * z1 - a * c * t2 * t8 * x1 + b * c * t2 * t4 * y1);
        rMatrix(3,0) = 0.0;
        rMatrix(3,1) = 0.0;
        rMatrix(3,2) = 0.0;
        rMatrix(3,3) = 1.0;
    }

    void ApplyPeriodicConditionProcess::TransformNode(array_1d<double, 3 >& rCoordinates, array_1d<double, 3 >& rTransformedCoordinates)
    {
        std::vector<double> original_node(4, 0.0f);
        std::vector<double> transformed_node(4, 0.0f);

        original_node[0] = rCoordinates(0); original_node[1] = rCoordinates(1); original_node[2] = rCoordinates(2); original_node[3] = 1.0;
        // Multiplying the point to get the rotated point
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                transformed_node[i] += mTransformationMatrix(i,j) * original_node[j];
            }
        }

        rTransformedCoordinates(0) = transformed_node[0]; rTransformedCoordinates(1) = transformed_node[1]; rTransformedCoordinates(2) = transformed_node[2];
    }

}  // namespace Kratos.
