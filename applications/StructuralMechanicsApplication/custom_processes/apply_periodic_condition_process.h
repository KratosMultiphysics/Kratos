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
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "utilities/binbased_fast_point_locator.h"
#include "processes/apply_multi_point_constraints_process.h"

#include "../../MappingApplication/custom_utilities/mapper.h"
#include "../../MappingApplication/custom_utilities/mapper_communicator.h"


namespace Kratos
{

class ApplyPeriodicConditionProcess : public Process
{
  public:
    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyPeriodicConditionProcess);

    typedef Dof<double> *DofPointerType;
    typedef Dof<double> DofType;

    typedef ApplyMultipointConstraintsProcess::VariableComponentType VariableComponentType;
    typedef KratosComponents<Variable<array_1d<double, 3>>> VectorVariableType;
    typedef ProcessInfo ProcessInfoType;
    typedef ProcessInfo::Pointer ProcessInfoPointerType;
    typedef unsigned int IndexType;

    typedef ApplyMultipointConstraintsProcess::VariableType VariableType;
    typedef ModelPart::NodeIterator NodeIterator;
    typedef ApplyMultipointConstraintsProcess::Pointer MpcProcessType;

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
                                            "magnitude":0.0,
                                            "interpolation_type":"nodes",                           
                                            "interpolation_settings":{
                                                            "search_radius" : 0.5,
                                                            "search_iterations":10,
                                                            "approximation_tolerance":1e-3,
                                                            "echo_level":3,
                                                            "mapper_type":"nearest_element",
                                                            "interface_submodel_part_origin":"default_master",
                                                            "interface_submodel_part_destination":"default_slave"
                                                             }                                                                                                            
                                            }  )");

        // Initializing
        m_parameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        std::cout << "Parameters afeter validating are :: " << m_parameters << std::endl;

        // using the given data first rotate the master surface to match with slave surface. obtain a new model part.

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

        mpMpcProcess = ApplyMultipointConstraintsProcess::Pointer(new ApplyMultipointConstraintsProcess(mrMainModelPart, m_parameters["constraint_set_name"].GetString()));
        CalculateTransformationMatrix();
        mIsInitialized = false;
        //mpMpcProcess->SetWeak();
        // Multiplying the point to get the rotated point
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                std::cout << mTransformationMatrix[i][j] << " :: ";
            }
            std::cout << std::endl;
        }
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
    ModelPart::Pointer mpRotatedMasterModelPart; // This is new modelpart
    std::string mSlaveSubModelPartName;
    std::string mMasterSubModelPartName;
    double mTheta;
    bool mIsInitialized;
    std::vector<double> mCenterOfRotation;
    std::vector<double> mAxisOfRoationVector;
    std::vector<double> mAxisOfRoationVectorNormalized;
    MpcProcessType mpMpcProcess;
    int mSign;
    std::string mType;
    double mModulus;
    std::vector<double> mDirOfTranslation;

    /*
    *   Structrue which contain the slave master information. This is used in conjection with the mapper communicator
    *   
    */
    struct MasterSlaveRelation
    {
        std::vector<int> MastersNodeIds;
        std::vector<double> MastersNodeWeights;
        std::vector<double> MasterConstants;
    };

    template <int TDim>
    void ApplyConstraintsForPeriodicConditions()
    {
        ModelPart &slaveModelPart = mrMainModelPart.GetSubModelPart(m_parameters["slave_sub_model_part_name"].GetString());
        ModelPart &master_model_part = (*mpRotatedMasterModelPart);
        int numVars = m_parameters["variable_names"].size();
        Parameters mapper_parameters = m_parameters["interpolation_settings"]; // TODO: find out how to give these settings
        MapperCommunicator::Pointer mpMapperCommunicator = MapperCommunicator::Pointer(
            new MapperCommunicator(master_model_part,
                                   slaveModelPart,
                                   mapper_parameters));
        if (m_parameters["interpolation_type"].GetString() == "elements")
        {
            mpMapperCommunicator->InitializeOrigin(MapperUtilities::Condition_Center);
            mpMapperCommunicator->InitializeDestination(MapperUtilities::Node_Coords);
            mpMapperCommunicator->Initialize();
        }
        else
        {
            mpMapperCommunicator->InitializeOrigin(MapperUtilities::Node_Coords);
            mpMapperCommunicator->InitializeDestination(MapperUtilities::Node_Coords);
            mpMapperCommunicator->Initialize();
        }

        for (int j = 0; j < numVars; j++)
        {
            std::string varName = m_parameters["variable_names"][j].GetString();
            if (KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Has(varName + "_X"))
            {
                if (m_parameters["interpolation_type"].GetString() == "elements")
                {
                    auto function_pointer_origin = std::bind(&GetMasterRelationInformationFromElementVectorVariable,
                                                             std::placeholders::_1,
                                                             std::placeholders::_2);
                    auto function_pointer_destination = std::bind(&SetMpcDataAtNodeVectorVariable<void *>,
                                                                  std::placeholders::_1,
                                                                  std::placeholders::_2,
                                                                  mpMpcProcess,
                                                                  mTransformationMatrix,
                                                                  varName, TDim, mSign);
                    mpMapperCommunicator->TransferVariableData(function_pointer_origin, function_pointer_destination);
                }
                else
                {
                    auto function_pointer_origin = std::bind(&GetMasterRelationInformationFromNodeVectorVariable,
                                                             std::placeholders::_1,
                                                             std::placeholders::_2);
                    auto function_pointer_destination = std::bind(&SetMpcDataAtNodeVectorVariable<void *>,
                                                                  std::placeholders::_1,
                                                                  std::placeholders::_2,
                                                                  mpMpcProcess,
                                                                  mTransformationMatrix,
                                                                  varName, TDim, mSign);
                    mpMapperCommunicator->TransferVariableData(function_pointer_origin, function_pointer_destination);
                }
            }
        }
    }

    /*
    * Function to be used in realation with the nearest node mapper. Slave side 
    */
    template <typename T>
    static void SetMpcDataAtNodeVectorVariable(InterfaceObject *pInterfaceObject, T rValue, MpcProcessType pMpcProcess, const std::vector<std::vector<double>> &transformationMatrix, std::string varName, int TDim, int sign)
    {
        VariableComponentType rVarX = KratosComponents<VariableComponentType>::Get(varName + std::string("_X"));
        VariableComponentType rVarY = KratosComponents<VariableComponentType>::Get(varName + std::string("_Y"));
        VariableComponentType rVarZ = KratosComponents<VariableComponentType>::Get(varName + std::string("_Z"));

        Node<3> *p_base_node = static_cast<InterfaceNode *>(pInterfaceObject)->pGetBase();
        unsigned int slaveNodeId = p_base_node->Id();

        MasterSlaveRelation *pMasterSlaveRelation = static_cast<MasterSlaveRelation *>(rValue);
        KRATOS_ERROR_IF_NOT(p_base_node) << "Base Pointer is nullptr!!!" << std::endl;
        // Marking the node as a slave
        p_base_node->Set(SLAVE);

        double constantCumulativeX(0.0), constantCumulativeY(0.0), constantCumulativeZ(0.0);

        for (unsigned int i = 0; i < pMasterSlaveRelation->MastersNodeIds.size(); i++)
        {
            int masterNodeNumber = pMasterSlaveRelation->MastersNodeIds[i];
            double masterweight = sign * pMasterSlaveRelation->MastersNodeWeights[i];

            constantCumulativeX += masterweight * transformationMatrix[0][3];
            constantCumulativeY += masterweight * transformationMatrix[1][3];
            constantCumulativeZ += masterweight * transformationMatrix[2][3];

            double constantX(0.0), constantY(0.0), constantZ(0.0);
            if (i == pMasterSlaveRelation->MastersNodeIds.size() - 1)
            {
                constantX = constantCumulativeX;
                constantY = constantCumulativeY;
                constantZ = constantCumulativeZ;
            }
            ///std::cout<<"rVarX :: "<<rVarX<<" :: rVarY ::"<<rVarY<<" :: rVarZ :: "<<rVarZ<<std::endl;

            pMpcProcess->AddMasterSlaveRelationWithNodeIdsAndVariableComponents(masterNodeNumber, rVarX, slaveNodeId, rVarX, masterweight * transformationMatrix[0][0], constantX);
            pMpcProcess->AddMasterSlaveRelationWithNodeIdsAndVariableComponents(masterNodeNumber, rVarY, slaveNodeId, rVarX, masterweight * transformationMatrix[0][1], constantX);
            pMpcProcess->AddMasterSlaveRelationWithNodeIdsAndVariableComponents(masterNodeNumber, rVarZ, slaveNodeId, rVarX, masterweight * transformationMatrix[0][2], constantX);

            pMpcProcess->AddMasterSlaveRelationWithNodeIdsAndVariableComponents(masterNodeNumber, rVarX, slaveNodeId, rVarY, masterweight * transformationMatrix[1][0], constantY);
            pMpcProcess->AddMasterSlaveRelationWithNodeIdsAndVariableComponents(masterNodeNumber, rVarY, slaveNodeId, rVarY, masterweight * transformationMatrix[1][1], constantY);
            pMpcProcess->AddMasterSlaveRelationWithNodeIdsAndVariableComponents(masterNodeNumber, rVarZ, slaveNodeId, rVarY, masterweight * transformationMatrix[1][2], constantY);

            if (TDim == 3)
            {
                pMpcProcess->AddMasterSlaveRelationWithNodeIdsAndVariableComponents(masterNodeNumber, rVarX, slaveNodeId, rVarZ, masterweight * transformationMatrix[2][0], constantZ);
                pMpcProcess->AddMasterSlaveRelationWithNodeIdsAndVariableComponents(masterNodeNumber, rVarY, slaveNodeId, rVarZ, masterweight * transformationMatrix[2][1], constantZ);
                pMpcProcess->AddMasterSlaveRelationWithNodeIdsAndVariableComponents(masterNodeNumber, rVarZ, slaveNodeId, rVarZ, masterweight * transformationMatrix[2][2], constantZ);
            }
        }

        delete pMasterSlaveRelation;
    }

    /*
    * Function to be used in realation with the nearest node mapper. Master side 
    */
    static MasterSlaveRelation *GetMasterRelationInformationFromElementVectorVariable(InterfaceObject *pInterfaceObject,
                                                                                      const std::vector<double> &rShapeFunctionValues)
    {
        MasterSlaveRelation *pMasterSlaveRelation = new MasterSlaveRelation();
        Geometry<Node<3>> *p_base_geometry = static_cast<InterfaceGeometryObject *>(pInterfaceObject)->pGetBase();
        KRATOS_ERROR_IF_NOT(p_base_geometry) << "Base Pointer is nullptr!!!" << std::endl;
        double constant = 0.0;
        for (std::size_t i = 0; i < p_base_geometry->PointsNumber(); ++i)
        {
            unsigned int nodeId = p_base_geometry->GetPoint(i).Id();
            pMasterSlaveRelation->MastersNodeIds.push_back(nodeId);

            pMasterSlaveRelation->MastersNodeWeights.push_back(rShapeFunctionValues[i]);
            pMasterSlaveRelation->MasterConstants.push_back(constant);
        }

        return pMasterSlaveRelation;
    }

    static MasterSlaveRelation *GetMasterRelationInformationFromNodeVectorVariable(InterfaceObject *pInterfaceObject,
                                                                                   const std::vector<double> &rShapeFunctionValues)
    {
        MasterSlaveRelation *pMasterSlaveRelation = new MasterSlaveRelation();
        Node<3> *p_base_node = static_cast<InterfaceNode *>(pInterfaceObject)->pGetBase();
        KRATOS_ERROR_IF_NOT(p_base_node) << "Base Pointer is nullptr!!!" << std::endl;
        double constant = 0.0;
        unsigned int nodeId = p_base_node->Id();
        pMasterSlaveRelation->MastersNodeIds.push_back(nodeId);

        pMasterSlaveRelation->MastersNodeWeights.push_back(1.0);
        pMasterSlaveRelation->MasterConstants.push_back(constant);

        return pMasterSlaveRelation;
    }

    // Functions which use two variable components
    template <int TDim>
    void GetRotatedMaster(ModelPart &master_model_part)
    {
        // iterating over slave nodes to find thecorresponding masters
        long int n_master_nodes = master_model_part.Nodes().size();
        for (int i = 0; i < n_master_nodes; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = (master_model_part).NodesBegin() + i;
            Node<3>::Pointer pnode = *(iparticle.base());
            (*mpRotatedMasterModelPart).CreateNewNode(pnode->Id(), *pnode);
        }

        // iterating over slave nodes to find thecorresponding masters
        long int n_master_elem = master_model_part.Elements().size();
        for (int i = 0; i < n_master_elem; i++)
        {
            Element::NodesArrayType newNodesArray;
            ModelPart::ElementsContainerType::iterator iparticle = (master_model_part).ElementsBegin() + i;
            for (unsigned int ii = 0; ii < (*iparticle).GetGeometry().size(); ii++)
            {
                newNodesArray.push_back((*mpRotatedMasterModelPart).pGetNode((*iparticle).GetGeometry()[ii].Id()));
            }
            Element::Pointer pElem = (*iparticle).Clone((*iparticle).Id(), newNodesArray);
            (*mpRotatedMasterModelPart).AddElement(pElem);
        }

        long int n_master_cond = master_model_part.Conditions().size();
        for (int i = 0; i < n_master_cond; i++)
        {
            Condition::NodesArrayType newNodesArray;
            ModelPart::ConditionsContainerType::iterator iparticle = (master_model_part).ConditionsBegin() + i;
            for (unsigned int ii = 0; ii < (*iparticle).GetGeometry().size(); ii++)
            {
                newNodesArray.push_back((*mpRotatedMasterModelPart).pGetNode((*iparticle).GetGeometry()[ii].Id()));
            }
            Condition::Pointer pCond = (*iparticle).Clone((*iparticle).Id(), newNodesArray);
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
        double eps = 1e-12;

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