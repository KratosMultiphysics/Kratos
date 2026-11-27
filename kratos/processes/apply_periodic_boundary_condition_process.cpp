//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "processes/apply_periodic_boundary_condition_process.h"
#include "utilities/binbased_fast_point_locator_conditions.h"
#include "utilities/geometrical_transformation_utilities.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "constraints/linear_master_slave_constraint.h"

namespace Kratos
{

ApplyPeriodicConditionProcess::ApplyPeriodicConditionProcess(ModelPart &rMasterModelPart, ModelPart &rSlaveModelPart,
                                Parameters Settings) : Process(Flags()), mrMasterModelPart(rMasterModelPart),
                                mrSlaveModelPart(rSlaveModelPart), mParameters(Settings)
{
    // Initializing
    const Parameters default_parameters = this->GetDefaultParameters();
    mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mCenterOfRotation = mParameters["transformation_settings"]["rotation_settings"]["center"].GetVector();
    mAxisOfRotationVector = mParameters["transformation_settings"]["rotation_settings"]["axis_of_rotation"].GetVector();
    mDirOfTranslation = mParameters["transformation_settings"]["translation_settings"]["dir_of_translation"].GetVector();

    mDistance = mParameters["transformation_settings"]["translation_settings"]["magnitude"].GetDouble();
    mAngleOfRotation = mParameters["transformation_settings"]["rotation_settings"]["angle_degree"].GetDouble() * 2 * Globals::Pi / 360.0;

    mTransformationMatrix.resize(4,4, false);
    mTransformationMatrixVariable.resize(4,4, false);

    mSearchMaxResults = mParameters["search_settings"]["max_results"].GetInt();
    mSearchTolerance = mParameters["search_settings"]["tolerance"].GetDouble();

    const double eps = std::numeric_limits<double>::epsilon();

    RemoveCommonNodesFromSlaveModelPart();

    if (std::abs(mAngleOfRotation) < eps && std::abs(mDistance) > eps){
        mTransformationType = ApplyPeriodicConditionProcess::TransformationType::TRANSLATION;
    }else if (std::abs(mAngleOfRotation) > eps && std::abs(mDistance) < eps){
        mTransformationType = ApplyPeriodicConditionProcess::TransformationType::ROTATION;
    }else{
        KRATOS_ERROR_IF(std::abs(mAngleOfRotation) < eps && std::abs(mDistance) < eps)<<"Both angle of rotation and modulus of translation cannot be zero. Please check the input"<<std::endl;
        KRATOS_ERROR_IF(std::abs(mAngleOfRotation) > eps && std::abs(mDistance) > eps)<<"Both angle of rotation and modulus of translation cannot be specified at the same time. Please check the input"<<std::endl;
    }

    CalculateTransformationMatrix();
}

ApplyPeriodicConditionProcess::~ApplyPeriodicConditionProcess()
{
}

Process::Pointer ApplyPeriodicConditionProcess::Create(Model &rModel, Parameters ThisParameters)
{
    const std::string master_model_part_name = ThisParameters["master_model_part_name"].GetString();
    const std::string slave_model_part_name = ThisParameters["slave_model_part_name"].GetString();
    return Kratos::make_shared<ApplyPeriodicConditionProcess>(rModel.GetModelPart(master_model_part_name), rModel.GetModelPart(slave_model_part_name), ThisParameters);
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
    const int domain_size = mrMasterModelPart.GetProcessInfo()[DOMAIN_SIZE];
    // Rotate the master so it goes to the slave
    if (domain_size == 2){
        ApplyConstraintsForPeriodicConditions<2>();
    }
    else if (domain_size == 3){
        ApplyConstraintsForPeriodicConditions<3>();
    } else {
        KRATOS_ERROR <<"Periodic conditions are designed onyl for 2 and 3 Dimensional cases ! "<<std::endl;
    }

    KRATOS_CATCH("");
}

/**
    * @brief Function initializes the solution step
    */
void ApplyPeriodicConditionProcess::ExecuteInitializeSolutionStep()
{
}

const Parameters ApplyPeriodicConditionProcess::GetDefaultParameters() const
{
    const Parameters default_parameters(R"(
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
    return default_parameters;
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
    const auto timer = BuiltinTimer();
    const unsigned int num_vars = mParameters["variable_names"].size();
    BinBasedFastPointLocatorConditions<TDim> bin_based_point_locator(mrMasterModelPart);
    bin_based_point_locator.UpdateSearchDatabase();

    // Define auxiliary functions
    using MyFunction = std::function<void(NodeType&, const GeometryType&, const VectorType&, const std::vector<const VariableType*>&)>;
    const MyFunction function_for_vector_variable =
    [this](NodeType& rSlaveNode, const GeometryType& rHostedGeometry, const VectorType& rWeights, const std::vector<const VariableType*>& rVars) {ConstraintSlaveNodeWithConditionForVectorVariable<TDim>(rSlaveNode, rHostedGeometry, rWeights, rVars);};
    const MyFunction function_for_scalar_variable =
    [this](NodeType& rSlaveNode, const GeometryType& rHostedGeometry, const VectorType& rWeights, const std::vector<const VariableType*>& rVars) {ConstraintSlaveNodeWithConditionForScalarVariable<TDim>(rSlaveNode, rHostedGeometry, rWeights, rVars);};

    // Fill an auxiliary vector with the functions
    std::vector<const MyFunction*> functions_required(num_vars, nullptr);
    std::vector<std::vector<const VariableType*>> variables_vector(num_vars, std::vector<const VariableType*>());
    for (unsigned int j = 0; j < num_vars; j++) {
        const std::string& r_var_name = mParameters["variable_names"][j].GetString();
        if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_var_name)) {
            functions_required[j] = &function_for_vector_variable;
            variables_vector[j].push_back(&KratosComponents<VariableType>::Get(r_var_name + "_X"));
            variables_vector[j].push_back(&KratosComponents<VariableType>::Get(r_var_name + "_Y"));
            variables_vector[j].push_back(&KratosComponents<VariableType>::Get(r_var_name + "_Z"));
        } else if (KratosComponents<VariableType>::Has(r_var_name)) {
            functions_required[j] = &function_for_scalar_variable;
            variables_vector[j].push_back(&KratosComponents<VariableType>::Get(r_var_name));
        } else {
            KRATOS_ERROR << "This only works with scalar and 3 components vector variables. Variable considered: " << r_var_name << std::endl;
        }
    }

    struct TLS_container {IndexType counter = 0; Condition::Pointer p_host_cond; VectorType shape_function_values; array_1d<double, 3 > transformed_slave_coordinates;};
    const IndexType num_slaves_found = block_for_each<SumReduction<IndexType>>(mrSlaveModelPart.Nodes(), TLS_container(), [&](Node& rNode, TLS_container& tls_container){
        tls_container.counter = 0;
        TransformNode(rNode.Coordinates(), tls_container.transformed_slave_coordinates);

        // Finding the host element for this node
        const bool is_found = bin_based_point_locator.FindPointOnMeshSimplified(tls_container.transformed_slave_coordinates, tls_container.shape_function_values, tls_container.p_host_cond, mSearchMaxResults, mSearchTolerance);
        if(is_found) {
            ++tls_container.counter;
            for (unsigned int j = 0; j < num_vars; j++) {
                (*functions_required[j])(rNode, tls_container.p_host_cond->GetGeometry(), tls_container.shape_function_values, variables_vector[j]);
            }
        }
        return tls_container.counter;
    });

    KRATOS_WARNING_IF("ApplyPeriodicConditionProcess",num_slaves_found != mrSlaveModelPart.NumberOfNodes())<<"Periodic condition cannot be applied for all the nodes."<<std::endl;
    KRATOS_INFO("ApplyPeriodicConditionProcess")<<"Applying periodic boundary conditions took : "<< timer.ElapsedSeconds() <<" seconds." <<std::endl;
}

template <int TDim>
void ApplyPeriodicConditionProcess::ConstraintSlaveNodeWithConditionForVectorVariable(
    NodeType& rSlaveNode, 
    const GeometryType& rHostedGeometry, 
    const VectorType& rWeights,
    const std::vector<const VariableType*>& rVars
    )
{
    // Reference constraint
    const auto& r_clone_constraint = LinearMasterSlaveConstraint();

    // Constant values
    array_1d<double, TDim> constants;

    IndexType master_index = 0;
    for (auto& r_master_node : rHostedGeometry) {
        const double master_weight = rWeights(master_index);

        for (unsigned int i = 0; i < TDim; ++i) {
            constants[i] = master_weight * mTransformationMatrixVariable(i,3);
        }

        #pragma omp critical
        {
            int current_num_constraint = mrMasterModelPart.GetRootModelPart().NumberOfMasterSlaveConstraints();
            for (unsigned int i = 0; i < 3; ++i) {
                for (unsigned int j = 0; j < TDim; ++j) {
                    auto p_constraint = r_clone_constraint.Create(++current_num_constraint, r_master_node, (*rVars[i]), rSlaveNode, (*rVars[j]), master_weight * mTransformationMatrixVariable(j,i), constants[j]);
                    mrMasterModelPart.AddMasterSlaveConstraint(p_constraint);
                }
            }
        }

        master_index++;
    }
}

template <int TDim>
void ApplyPeriodicConditionProcess::ConstraintSlaveNodeWithConditionForScalarVariable(
    NodeType& rSlaveNode, 
    const GeometryType& rHostedGeometry, 
    const VectorType& rWeights,
    const std::vector<const VariableType*>& rVars 
    )
{
    const VariableType& r_var = (*rVars[0]);

    // Reference constraint
    const auto& r_clone_constraint = LinearMasterSlaveConstraint();

    IndexType master_index = 0;
    for (auto& r_master_node : rHostedGeometry)
    {
        const double master_weight = rWeights(master_index);
        #pragma omp critical
        {
            int current_num_constraint = mrMasterModelPart.GetRootModelPart().NumberOfMasterSlaveConstraints();
            auto constraint = r_clone_constraint.Create(++current_num_constraint,r_master_node, r_var, rSlaveNode, r_var, master_weight, 0.0);
            mrMasterModelPart.AddMasterSlaveConstraint(constraint);
        }
        master_index++;
    }
}

void ApplyPeriodicConditionProcess::CalculateTransformationMatrix()
{
    if (mTransformationType == ApplyPeriodicConditionProcess::TransformationType::TRANSLATION){
        GeometricalTransformationUtilities::CalculateTranslationMatrix(-1*mDistance, mTransformationMatrix, mDirOfTranslation);
        GeometricalTransformationUtilities::CalculateTranslationMatrix(0.0, mTransformationMatrixVariable, mDirOfTranslation);
    }
    else if (mTransformationType == ApplyPeriodicConditionProcess::TransformationType::ROTATION){
        GeometricalTransformationUtilities::CalculateRotationMatrix(-1*mAngleOfRotation, mTransformationMatrix, mAxisOfRotationVector, mCenterOfRotation);
        GeometricalTransformationUtilities::CalculateRotationMatrix(mAngleOfRotation, mTransformationMatrixVariable, mAxisOfRotationVector, mCenterOfRotation);
    }
}

void ApplyPeriodicConditionProcess::TransformNode(const array_1d<double, 3 >& rCoordinates, array_1d<double, 3 >& rTransformedCoordinates) const
{
    DenseVector<double> original_node(4, 0.0);
    DenseVector<double> transformed_node(4, 0.0);

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
