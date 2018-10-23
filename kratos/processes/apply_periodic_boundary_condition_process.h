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
#include <chrono>

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/point.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "utilities/binbased_fast_point_locator_conditions.h"



namespace Kratos
{

class ApplyPeriodicConditionProcess : public Process
{

  public:
    /// Pointer definition of ApplyPeriodicConditionProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyPeriodicConditionProcess);

    typedef Dof<double>*                                    DofPointerType;
    typedef Dof<double>                                     DofType;
    typedef ModelPart::VariableComponentType                VariableComponentType;
    typedef KratosComponents<Variable<array_1d<double, 3>>> VectorVariableType;
    typedef ProcessInfo                                     ProcessInfoType;
    typedef ProcessInfo::Pointer                            ProcessInfoPointerType;
    typedef unsigned int                                    IndexType;
    typedef ModelPart::DoubleVariableType                   VariableType;
    typedef ModelPart::NodeIterator                         NodeIterator;
    typedef Element                                         ElementType;
    typedef Node<3>                                         NodeType;
    typedef Matrix                                          MatrixType;
    typedef Vector                                          VectorType;
    typedef Geometry<NodeType>                              GeometryType;

    /**
     * @brief Constructor of the process to apply periodic boundary condition
     * @param rModelPart The main model part where the boundary conditions are applied.
     * @param rParameters parameters for the periodic condition to be applied
     */
    ApplyPeriodicConditionProcess(ModelPart &rMasterModelPart, ModelPart &rSlaveModelPart,
                                  Parameters rParameters) : Process(Flags()), mrMasterModelPart(rMasterModelPart),
                                   mrSlaveModelPart(rSlaveModelPart), mParameters(rParameters)
    {
        Parameters default_parameters(R"(
                                            {
                                            "variable_names":[],
                                            "center":[0,0,0],
                                            "axis_of_rotation":[0.0,0.0,0.0],
                                            "angle":0.0,
                                            "dir_of_translation":[0.0,0.0,0.0],
                                            "magnitude":0.0
                                            }  )");

        // Initializing
        mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        mCenterOfRotation.push_back(mParameters["center"][0].GetDouble());
        mCenterOfRotation.push_back(mParameters["center"][1].GetDouble());
        mCenterOfRotation.push_back(mParameters["center"][2].GetDouble());

        mAxisOfRotationVector.push_back(mParameters["axis_of_rotation"][0].GetDouble());
        mAxisOfRotationVector.push_back(mParameters["axis_of_rotation"][1].GetDouble());
        mAxisOfRotationVector.push_back(mParameters["axis_of_rotation"][2].GetDouble());

        mModulus = mParameters["magnitude"].GetDouble();
        mDirOfTranslation.push_back(mParameters["dir_of_translation"][0].GetDouble());
        mDirOfTranslation.push_back(mParameters["dir_of_translation"][1].GetDouble());
        mDirOfTranslation.push_back(mParameters["dir_of_translation"][2].GetDouble());

        mTheta = mParameters["angle"].GetDouble() * 2 * 3.1416 / 360.0;

        mTransformationMatrix.resize(4,4);
        mTransformationMatrixVariable.resize(4,4);

        RemoveCommonNodesFromSlaveModelPart();

        if (mTheta == 0 && mModulus != 0)
            mType = "translation";
        else if (mTheta != 0.0 && mModulus == 0.0)
            mType = "rotation";
        else
            KRATOS_ERROR_IF(mTheta == 0.0 && mModulus == 0.0)<<"Both angle of rotation and modulus of translation cannot be zero. Please check the input"<<std::endl;
        else
            KRATOS_ERROR_IF(mTheta != 0.0 && mModulus != 0.0)<<"Both angle of rotation and modulus of translation cannot be specified at the same time. Please check the input"<<std::endl;

        CalculateTransformationMatrix();
        mIsInitialized = false;
        std::cout<<"########### New Implementation :: "<<std::endl;
    }

    /**
     * @brief Destructor of the process class
     */
    ~ApplyPeriodicConditionProcess()
    {
    }

    /**
     * @brief  Function to remove the common nodes of slave and master modelparts from the slave modelpart
     */
    void RemoveCommonNodesFromSlaveModelPart()
    {
        const int num_slave_nodes = mrSlaveModelPart.NumberOfNodes();
        const ModelPart::NodeIterator it_slave_node_begin = mrSlaveModelPart.NodesBegin();
        for(int i_node = 0; i_node<num_slave_nodes; ++i_node)
        {
            ModelPart::NodeIterator it_slave_node = it_slave_node_begin;
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
    virtual void ExecuteInitialize() override
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
    void ExecuteInitializeSolutionStep() override
    {
    }

    /**
     * @brief Function to print the information about this current process
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream <<"ApplyPeriodicConditionProcess Process "<<std::endl;
    }

  private:
    MatrixType mTransformationMatrix; // This can be either for rotating or for translating the slave geometry to Master geometry
    MatrixType mTransformationMatrixVariable; // This can be either for rotating or for translating the master variable to slave geometry
    ModelPart &mrMasterModelPart;       // the master modelpart to which the master-slave constraints are added.
    ModelPart &mrSlaveModelPart;
    Parameters mParameters;          // parameters
    double mTheta;
    bool mIsInitialized;
    std::vector<double> mCenterOfRotation;
    std::vector<double> mAxisOfRotationVector;
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
        const double start_apply = OpenMPUtils::GetCurrentTime();
        const int num_vars = mParameters["variable_names"].size();
        BinBasedFastPointLocatorConditions<TDim> bin_based_point_locator(mrMasterModelPart);
        bin_based_point_locator.UpdateSearchDatabase();

        // for bin based point locator
		VectorType shape_function_values;
		const int max_results = 10000;
		typename BinBasedFastPointLocatorConditions<TDim>::ResultContainerType results(max_results);
        typename BinBasedFastPointLocatorConditions<TDim>::ResultIteratorType result_begin = results.begin();

        const int num_slave_nodes = mrSlaveModelPart.NumberOfNodes();
        const ModelPart::NodeIterator it_slave_node_begin = mrSlaveModelPart.NodesBegin();

        unsigned int num_slaves_found = 0;
#pragma omp parallel for schedule(guided, 512)
        for(int i_node = 0; i_node<num_slave_nodes; ++i_node)
        {
            Condition::Pointer p_host_cond;
            ModelPart::NodeIterator it_slave_node = it_slave_node_begin;
            std::advance(it_slave_node, i_node);
            array_1d<double, 3 > transformed_slave_coordinates;
            TransformNode(it_slave_node->Coordinates(), transformed_slave_coordinates);

            // Finding the host element for this node
            const bool is_found = bin_based_point_locator.FindPointOnMesh(transformed_slave_coordinates, shape_function_values, p_host_cond, result_begin, max_results);
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

    /**
     * @brief   The function add the linear master-slave constraint to rModelPart. This function is specifically for applying
     *          periodic conditions for vector variable. This distinction is because, for vector variables the variable should also be
     *          transformed according to the transfromation specified in the settings.
     * @param rSalveNode The slave node which is to be connected to the rHostGeometry.
     * @param rHostedGeometry the Host geometry which has the rSlaveNode.
     * @param rWeights The weights with which the rSlaveNode is connected to the rHostedGeometry's nodes.
     * @param rVarName The name of the vector variable on which periodic boundary condition can be applied.
     */
    template <int TDim>
    void ConstraintSlaveNodeWithConditionForVectorVariable(NodeType& rSalveNode, const GeometryType& rHostedGeometry, const VectorType& rWeights, const std::string& rVarName )
    {
        const VariableComponentType r_var_x = KratosComponents<VariableComponentType>::Get(rVarName + std::string("_X"));
        const VariableComponentType r_var_y = KratosComponents<VariableComponentType>::Get(rVarName + std::string("_Y"));
        const VariableComponentType r_var_z = KratosComponents<VariableComponentType>::Get(rVarName + std::string("_Z"));

        IndexType master_index = 0;
        for (auto& master_node : rHostedGeometry)
        {
                int current_num_constraint = mrMasterModelPart.GetRootModelPart().NumberOfMasterSlaveConstraints();

                const double master_weight = rWeights(master_index);

                const double constant_x = master_weight * mTransformationMatrixVariable(0,3);
                const double constant_y = master_weight * mTransformationMatrixVariable(1,3);
                const double constant_z = master_weight * mTransformationMatrixVariable(2,3);

                #pragma omp critical
                {
                    mrMasterModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_x, rSalveNode, r_var_x, master_weight * mTransformationMatrixVariable(0,0), constant_x);
                    mrMasterModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_y, rSalveNode, r_var_x, master_weight * mTransformationMatrixVariable(0,1), constant_x);
                    mrMasterModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_z, rSalveNode, r_var_x, master_weight * mTransformationMatrixVariable(0,2), constant_x);

                    mrMasterModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_x, rSalveNode, r_var_y, master_weight * mTransformationMatrixVariable(1,0), constant_y);
                    mrMasterModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_y, rSalveNode, r_var_y, master_weight * mTransformationMatrixVariable(1,1), constant_y);
                    mrMasterModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_z, rSalveNode, r_var_y, master_weight * mTransformationMatrixVariable(1,2), constant_y);

                    if (TDim == 3)
                    {
                        mrMasterModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_x, rSalveNode, r_var_z, master_weight * mTransformationMatrixVariable(2,0), constant_z);
                        mrMasterModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_y, rSalveNode, r_var_z, master_weight * mTransformationMatrixVariable(2,1), constant_z);
                        mrMasterModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var_z, rSalveNode, r_var_z, master_weight * mTransformationMatrixVariable(2,2), constant_z);
                    }
                }

            master_index++;
        }
    }

    /**
     * @brief   The function add the linear master-slave constraint to rModelPart. This function is specifically for applying
     *          periodic conditions for scalar variable. This distinction is because, for scalar variables the variable need NOT be
     *          transformed according to the transfromation specified in the settings.
     * @param rSalveNode The slave node which is to be connected to the rHostGeometry.
     * @param rHostedGeometry the Host geometry which has the rSlaveNode.
     * @param rWeights The weights with which the rSlaveNode is connected to the rHostedGeometry's nodes.
     * @param rVarName The name of the scalar variable on which periodic boundary condition can be applied.
     */
    template <int TDim>
    void ConstraintSlaveNodeWithConditionForScalarVariable(NodeType& rSalveNode, const GeometryType& rHostedGeometry, const VectorType& rWeights, const std::string& rVarName )
    {
        const VariableType r_var = KratosComponents<VariableType>::Get(rVarName);

        IndexType master_index = 0;
        for (auto& master_node : rHostedGeometry)
        {
            int current_num_constraint = mrMasterModelPart.NumberOfMasterSlaveConstraints();

            const double master_weight = rWeights(master_index);
            #pragma omp critical
            {
                mrMasterModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", current_num_constraint++, master_node, r_var, rSalveNode, r_var, master_weight, 0.0);
            }
            master_index++;
        }
    }

    /**
     * @brief   Fuctions calculates the transformation matrix to account for the moving the two periodic condition modelparts together.
     */
    void CalculateTransformationMatrix()
    {
        if (mType == "translation"){
            CalculateTranslationMatrix(-1*mModulus, mTransformationMatrix);
            CalculateTranslationMatrix(0.0, mTransformationMatrixVariable);
        }
        else if (mType == "rotation"){
            CalculateRotationMatrix(-1*mTheta, mTransformationMatrix);
            CalculateRotationMatrix(mTheta, mTransformationMatrixVariable);
        }
    }

    /**
     * @brief   Calculates the transformation matrix which translates the given vector alone mDirOfTranslation by mModulus
     * @param   Modulus is the magnitude by which the translation should happen in the direction of mDirOfTranslation.
     * @param   rMatrix is the transformation matrix which will be calculated in this function. This should be of correct size (4x4).
     */
    void CalculateTranslationMatrix(const double Modulus, MatrixType& rMatrix)
    {
        const double tolerance = std::numeric_limits<double>::epsilon();
        rMatrix(0,0) = 1.0;
        rMatrix(0,1) = 0.0;
        rMatrix(0,2) = 0.0;
        rMatrix(0,3) = (abs(Modulus)<tolerance) ?  1.0 : mModulus *mDirOfTranslation[0];
        rMatrix(1,0) = 0.0;
        rMatrix(1,1) = 1.0;
        rMatrix(1,2) = 0.0;
        rMatrix(1,3) = (abs(Modulus)<tolerance) ?  1.0 : mModulus *mDirOfTranslation[1];
        rMatrix(2,0) = 0.0;
        rMatrix(2,1) = 0.0;
        rMatrix(2,2) = 1.0;
        rMatrix(2,3) = (abs(Modulus)<tolerance) ?  1.0 : mModulus *mDirOfTranslation[2];
        rMatrix(3,0) = 0.0;
        rMatrix(3,1) = 0.0;
        rMatrix(3,2) = 0.0;
        rMatrix(3,3) = 1.0;
    }

    /**
     * @brief   Calculates the transformation matrix which rotates the given vector around mAxisOfRotationVector and mCenterOfRotation
     *          by provided Theta and stores the result in rMatrix The following code is generated from MATLAB and is adapted here.
     * @param   Theta is the angle of rotation about mAxisOfRotationVector and mCenterOfRotation.
     * @param   rMatrix is the transformation matrix which will be calculated in this function. This should be of correct size (4x4).
     */
    void CalculateRotationMatrix(const double Theta, MatrixType& rMatrix )
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

        const double t2 = cos(Theta);
        const double t3 = sin(Theta);
        const double t4 = a * a;
        const double t5 = b * b;
        const double t6 = c * c;
        const double t7 = a * b;
        const double t8 = t5 + t6;
        double t9 = 1.0 / t8;
        if(isnan(t9))
            t9 = 0.0;
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

    /*
     * @brief Rotates a given point(node_cords) in space around a given mAxisOfRotationVector by an angle thetha
     * @param rCoordinates The original coordinates which have to be transformed
     * @param rTransformedCoordinates The new coordinates which are transformed with rTransformationMatrix.
     */
    void TransformNode(array_1d<double, 3 >& rCoordinates, array_1d<double, 3 >& rTransformedCoordinates)
    {
        std::vector<double> original_node(4, 0.0f);
        std::vector<double> transformed_node(4, 0.0f);

        original_node[0] = rCoordinates(0); original_node[1] = rCoordinates(1); original_node[2] = rCoordinates(2); original_node[3] = 1.0;
        // Multiplying the point to get the rotated point
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                transformed_node[i] += mTransformationMatrix(i,j) * original_node[j];
            }
        }

        rTransformedCoordinates(0) = transformed_node[0]; rTransformedCoordinates(1) = transformed_node[1]; rTransformedCoordinates(2) = transformed_node[2];
    }

}; // Class


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