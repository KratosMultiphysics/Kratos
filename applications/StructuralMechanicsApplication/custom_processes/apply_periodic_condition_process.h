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
#include "includes/element.h"
#include "includes/node.h"
#include "geometries/point.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "utilities/binbased_fast_point_locator.h"
#include "spatial_containers/spatial_containers.h"
#include "includes/linear_master_slave_constraint.h"

namespace Kratos
{
/**
 * @class ApplyPeriodicConditionProcess
 *
 * @ingroup StructuralMechanicsApplication
 *
 * @brief This method computes and assigns the periodic boundary conditions on
 * specified boundary conditions process.
 *
 * @author Aditya Ghantasala
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ApplyPeriodicConditionProcess
    : public Process
{
  public:
    /// Pointer definition of ApplyPeriodicConditionProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyPeriodicConditionProcess);

    typedef Element ElementType;
    typedef Dof<double> *DofPointerType;
    typedef Dof<double> DofType;
    typedef Geometry<ModelPart::NodeType> GeometryType;

    typedef typename ModelPart::VariableComponentType VariableComponentType;
    typedef KratosComponents<Variable<array_1d<double, 3>>> VectorVariableType;
    typedef ProcessInfo ProcessInfoType;
    typedef ProcessInfo::Pointer ProcessInfoPointerType;
    typedef std::size_t IndexType;

    typedef typename ModelPart::DoubleVariableType VariableType;

    /// Constructor.
    ApplyPeriodicConditionProcess(ModelPart &rModelPart,
                                  Parameters rParameters) : Process(Flags()), mrMainModelPart(rModelPart), mParameters(rParameters)
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
        mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        // using the given data first rotate the master surface to match with slave surface. obtain a new model part.

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

        // normalizing the axis of roatation
        double norm = 0.0;
        for (unsigned int d = 0; d < 3; ++d)
            norm += mAxisOfRoationVector[d] * mAxisOfRoationVector[d];
        norm = sqrt(norm);
        for (unsigned int d = 0; d < 3; ++d)
            mAxisOfRoationVectorNormalized.push_back(mAxisOfRoationVector[d] / norm);

        mTheta = mParameters["angle"].GetDouble() * 2 * 3.1416 / 360.0;
        mpRotatedMasterModelPart = ModelPart::Pointer(new ModelPart("rotatedMaster"));

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
        // Multiplying the point to get the rotated point
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                std::cout << mTransformationMatrix[i][j] << " :: ";
            }
            std::cout << std::endl;
        }

        mSearchRadius = 0.0;
    }
    ~ApplyPeriodicConditionProcess()
    {
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
       rOStream <<" ApplyPeriodicConditionProcess "<<std::endl;
    }

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        mConstraintId = 0;
        if (!mIsInitialized)
        {
            ModelPart &masterModelPart = mrMainModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());

            ProcessInfoPointerType info = mrMainModelPart.pGetProcessInfo();
            int prob_dim = info->GetValue(DOMAIN_SIZE);
            // Rotate the master so it goes to the slave
            if (prob_dim == 2)
            {
                GetRotatedMaster<2>(masterModelPart);
                ApplyConstraintsForPeriodicConditions<2>();
            }
            else if (prob_dim == 3)
            {
                GetRotatedMaster<3>(masterModelPart);
                std::cout << "applying periodic conditions >> START " << std::endl;
                ApplyConstraintsForPeriodicConditions<3>();
                std::cout << "applying periodic conditions >> END " << std::endl;
            }

            mIsInitialized = true;
        }

        std::cout<<"##################################### "<<std::endl;
        std::cout<<"##################################### :: "<<mrMainModelPart.MasterSlaveConstraints().size()<<std::endl;
        std::cout<<"##################################### "<<std::endl;


        KRATOS_CATCH("");
    }

  private:
    IndexType mConstraintId;
    std::vector<std::vector<double>> mTransformationMatrix;
    ModelPart &mrMainModelPart;
    Parameters mParameters;
    ModelPart::Pointer mpRotatedMasterModelPart; // This is new modelpart
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
    double mSearchRadius;

  private:
    class ConditionCenterPoint : public Point
    {

      public:
        KRATOS_CLASS_POINTER_DEFINITION(ConditionCenterPoint);

        ConditionCenterPoint() : Point(), mrCondition(  *(Kratos::make_shared<Condition>()) )
        {
        }

        ConditionCenterPoint(Point& rPoint) : Point(rPoint), mrCondition(  *(Kratos::make_shared<Condition>()) )
        {
        }

        ConditionCenterPoint(Condition& rCondition) : Point(), mrCondition(rCondition)
        {
            double xC(0.0), yC(0.0), zC(0.0);
            int number_points = rCondition.GetGeometry().size();
            auto geometry = rCondition.GetGeometry();

            for (int i = 0; i < number_points; i++)
            {
                xC += geometry[i].X();
                yC += geometry[i].Y();
                zC += geometry[i].Z();
            }

            this->X() = xC / number_points;
            this->Y() = yC / number_points;
            this->Z() = zC / number_points;
        }

        Condition& GetCondition()
        {
            return mrCondition;
        }

        ~ConditionCenterPoint()
        {
        }

      private:
        Condition& mrCondition;
    };

    typedef ConditionCenterPoint ElementCenterPointType;
    typedef ConditionCenterPoint::Pointer ElementCenterPointTypePointer;
    typedef Point::Pointer PointPointerType;
    typedef std::vector<ElementCenterPointTypePointer> ElementCenterPointerVector;
    typedef std::vector<ElementCenterPointTypePointer>::iterator CenterIterator;
    typedef std::vector<double>::iterator DoubleVectorIterator;

    template <int TDim>
    void ApplyConstraintsForPeriodicConditions()
    {
        ModelPart &slave_model_part = mrMainModelPart.GetSubModelPart(mParameters["slave_sub_model_part_name"].GetString());
        ModelPart &master_model_part = (*mpRotatedMasterModelPart);
        int num_vars = mParameters["variable_names"].size();

        // Type definitions for tree-search
        typedef Bucket<TDim, ElementCenterPointType, ElementCenterPointerVector, ElementCenterPointTypePointer, CenterIterator, DoubleVectorIterator> BucketType;
        typedef Tree<KDTreePartition<BucketType>> KDTree;
        std::size_t bucket_size = 100;
        unsigned int max_number_of_neighbors = 3;
        ElementCenterPointerVector centers_vector(master_model_part.NumberOfConditions());
        this->CalculateElementCenterVectors(master_model_part, centers_vector);

        typename KDTree::Pointer p_search_tree = Kratos::make_shared<KDTree>(centers_vector.begin(), centers_vector.end(), bucket_size);

        for (int var_index = 0; var_index < num_vars; var_index++)
        {
            std::string var_name = mParameters["variable_names"][var_index].GetString();
            if (KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Has(var_name + "_X"))
            { // its a variable with components, we can apply periodic condition with vector variables

                // Loop over the slave nodes
                for (auto &node : slave_model_part.Nodes())
                {
                    ElementCenterPointerVector neighbor_nodes(max_number_of_neighbors);
                    std::vector<double> resulting_squared_distances(max_number_of_neighbors);
                    // Marking the node as a slave
                    node.Set(SLAVE);
                    ConditionCenterPoint center_point(node);
                    std::size_t number_of_neighbors = p_search_tree->SearchInRadius(center_point,
                                                                                    mSearchRadius,
                                                                                    neighbor_nodes.begin(),
                                                                                    resulting_squared_distances.begin(),
                                                                                    max_number_of_neighbors);

                    // Now we know which 10 elements are close now we project the slave point on to these elements
                    for (unsigned int nei_index = 0; nei_index < number_of_neighbors; nei_index++)
                    {
                        auto &cond = neighbor_nodes[nei_index]->GetCondition();
                        typename Point::CoordinatesArrayType local_coordinates;
                        Vector shape_function_values;
                        bool is_found = cond.GetGeometry().IsInside(node.Coordinates(), local_coordinates);
                        cond.GetGeometry().ShapeFunctionsValues(shape_function_values, local_coordinates);
                        // If found apply the constraint and break
                        if (is_found)
                        {
                            ApplyCondition<TDim>(var_name, node, cond, shape_function_values);
                            break;
                        }
                    }
                }
            }
            else
            {
                KRATOS_ERROR << "Provided variable :: " << var_name << ". Periodic boundary conditions are only applied for variables with components. " << std::endl;
            }
        }
        // void SearchNearestPoint(PointType const& rThisPoint, PointerType& rResult, CoordinateType& rResultDistance )
    }

    void CalculateElementCenterVectors(ModelPart &rModelPart, ElementCenterPointerVector &rCentersVector)
    {
        std::size_t index = 0;
        for (auto &cond : rModelPart.Conditions())
        {
            rCentersVector[index] = Kratos::make_shared<ElementCenterPointType>(cond);
            auto elem_length = cond.GetGeometry().Length();
            mSearchRadius = mSearchRadius < elem_length ? elem_length : mSearchRadius;
            index++;
        }
    }

    template <int TDim>
    void ApplyCondition(std::string &rVarName, Node<3> &rSlaveNode, Condition &rCondition, Vector &rShapeFunctionValues)
    {

        const VariableComponentType rVarX = KratosComponents<VariableComponentType>::Get(rVarName + std::string("_X"));
        const VariableComponentType rVarY = KratosComponents<VariableComponentType>::Get(rVarName + std::string("_Y"));
        const VariableComponentType rVarZ = KratosComponents<VariableComponentType>::Get(rVarName + std::string("_Z"));

        std::size_t slave_node_id = rSlaveNode.Id();
        GeometryType &geom = rCondition.GetGeometry();
        std::size_t num_master_nodes = geom.size();


        double constant_cumulative_x(0.0), constant_cumulative_y(0.0), constant_cumulative_z(0.0);

        for (unsigned int i = 0; i < num_master_nodes; i++)
        {
            std::size_t master_node_id = geom[i].Id();
            // This can happen for example nodes on the axis
            if (slave_node_id == master_node_id)
            {
                return;
            }
        }

        for (unsigned int i = 0; i < num_master_nodes; i++)
        {
            double master_weight = mSign * rShapeFunctionValues[i];

            constant_cumulative_x += master_weight * mTransformationMatrix[0][3];
            constant_cumulative_y += master_weight * mTransformationMatrix[1][3];
            constant_cumulative_z += master_weight * mTransformationMatrix[2][3];

            double constant_x(0.0), constant_y(0.0), constant_z(0.0);
            if (i == num_master_nodes - 1)
            {
                constant_x = constant_cumulative_x;
                constant_y = constant_cumulative_y;
                constant_z = constant_cumulative_z;
            }

            if (! geom[i].HasDofFor(rVarX) )
            {
                KRATOS_ERROR << "Error while applying periodic boundary condition. Master Node " <<geom[i].Id()<<" does not have DOF "<<rVarX<<std::endl;
            }
            if (! rSlaveNode.HasDofFor(rVarX) )
            {
                KRATOS_ERROR << "Error while applying periodic boundary condition. Slave Node " <<rSlaveNode.Id()<<" does not have DOF "<<rVarX<<std::endl;
            }

            mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", mConstraintId++, geom[i], rVarX, rSlaveNode, rVarX, master_weight * mTransformationMatrix[0][0], constant_x);
            mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", mConstraintId++, geom[i], rVarY, rSlaveNode, rVarX, master_weight * mTransformationMatrix[0][1], constant_x);
            mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", mConstraintId++, geom[i], rVarZ, rSlaveNode, rVarX, master_weight * mTransformationMatrix[0][2], constant_x);

            mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", mConstraintId++, geom[i], rVarX, rSlaveNode, rVarY, master_weight * mTransformationMatrix[1][0], constant_y);
            mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", mConstraintId++, geom[i], rVarY, rSlaveNode, rVarY, master_weight * mTransformationMatrix[1][1], constant_y);
            mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", mConstraintId++, geom[i], rVarZ, rSlaveNode, rVarY, master_weight * mTransformationMatrix[1][2], constant_y);

            if (TDim == 3)
            {
                mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", mConstraintId++, geom[i], rVarX, rSlaveNode, rVarZ, master_weight * mTransformationMatrix[2][0], constant_z);
                mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", mConstraintId++, geom[i], rVarY, rSlaveNode, rVarZ, master_weight * mTransformationMatrix[2][1], constant_z);
                mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", mConstraintId++, geom[i], rVarZ, rSlaveNode, rVarZ, master_weight * mTransformationMatrix[2][2], constant_z);
            }
        }
    }

    // Functions which use two variable components
    template <int TDim>
    void GetRotatedMaster(ModelPart &master_model_part)
    {

        long int n_master_nodes = master_model_part.Nodes().size();
         for (int i = 0; i < n_master_nodes; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = (master_model_part).NodesBegin() + i;
            Node<3>::Pointer pnode = *(iparticle.base());
            auto new_node = pnode->Clone();
            (*mpRotatedMasterModelPart).AddNode(new_node);
        }

        long int n_master_elem = master_model_part.Elements().size();
        for (int i = 0; i < n_master_elem; i++)
        {
            Element::NodesArrayType new_nodes_array;
            ModelPart::ElementsContainerType::iterator iparticle = (master_model_part).ElementsBegin() + i;
            for (unsigned int ii = 0; ii < (*iparticle).GetGeometry().size(); ii++)
            {
                new_nodes_array.push_back((*mpRotatedMasterModelPart).pGetNode((*iparticle).GetGeometry()[ii].Id()));
            }
            Element::Pointer pElem = (*iparticle).Clone((*iparticle).Id(), new_nodes_array);
            (*mpRotatedMasterModelPart).AddElement(pElem);
        }

        long int n_master_cond = master_model_part.Conditions().size();
        for (int i = 0; i < n_master_cond; i++)
        {
            Condition::NodesArrayType new_nodes_array;
            ModelPart::ConditionsContainerType::iterator iparticle = (master_model_part).ConditionsBegin() + i;
            for (unsigned int ii = 0; ii < (*iparticle).GetGeometry().size(); ii++)
            {
                new_nodes_array.push_back((*mpRotatedMasterModelPart).pGetNode((*iparticle).GetGeometry()[ii].Id()));
            }
            Condition::Pointer pCond = (*iparticle).Clone((*iparticle).Id(), new_nodes_array);
            (*mpRotatedMasterModelPart).AddCondition(pCond);
        }

        // Rotating the nodes
        for (int i = 0; i < n_master_nodes; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = (*mpRotatedMasterModelPart).NodesBegin() + i;
            Node<3>::Pointer p_master_node = *(iparticle.base());
            std::vector<double> masterNode = {p_master_node->X(), p_master_node->Y(), p_master_node->Z(), 1};

            std::vector<double> rotatedMasterNode = TransformNode(i, masterNode, mTheta);

            p_master_node->X() = rotatedMasterNode[0];
            p_master_node->Y() = rotatedMasterNode[1];
            if (TDim > 2)
                p_master_node->Z() = rotatedMasterNode[2];

            p_master_node->X0() = rotatedMasterNode[0];
            p_master_node->Y0() = rotatedMasterNode[1];
            if (TDim > 2)
            {
                p_master_node->Z0() = rotatedMasterNode[2];
            }
        }
    }

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

    /*
     * Rotates a given point(node_cords) in space around a given mAxisOfRoationVector by an angle thetha
     */
    std::vector<double> TransformNode(long int index, std::vector<double> node_cords, double theta)
    {
        std::vector<double> rotatedNode(4, 0.0f);

        // Multiplying the point to get the rotated point
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                rotatedNode[i] += (mTransformationMatrix[i][j] * node_cords[j]);
            }
        }

        return rotatedNode;
    }
};
}

#endif