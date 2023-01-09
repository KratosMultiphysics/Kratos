

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "custom_utilities/geometrical/symmetry_utility.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"
#include "optimization_application.h"

#define PI 3.14159265

// ==============================================================================

namespace Kratos
{

SymmetryUtility::SymmetryUtility(std::string Name, ModelPart& rModelPart, Parameters SymmetrySettings )
    : mUtilName(Name), mrModelPart( rModelPart ), mSymmetrySettings( SymmetrySettings )
{

    if(mSymmetrySettings["type"].GetString()=="plane_symmetry"){
        mPlaneSymmetry = true;
        mPlaneSymmetryData.Point = mSymmetrySettings["settings"]["point"].GetVector();
        array_3d normal = mSymmetrySettings["settings"]["normal"].GetVector();
        const double norm_normal = norm_2(normal);
        if (norm_normal < std::numeric_limits<double>::epsilon())
            KRATOS_ERROR<<"SymmetryUtility: norm of normal is close to zero"<<std::endl;
            
        mPlaneSymmetryData.Normal = normal/norm_normal;
        mPlaneSymmetryData.ReflectionMatrix = IdentityMatrix(3) - (2*outer_prod(mPlaneSymmetryData.Normal, mPlaneSymmetryData.Normal));
    }
    else if(mSymmetrySettings["type"].GetString()=="axis_symmetry"){
        mAxisSymmetry = true;
        mAxisSymmetryData.Point = mSymmetrySettings["settings"]["point"].GetVector();
        array_3d axis = mSymmetrySettings["settings"]["axis"].GetVector();
        const double norm_axis = norm_2(axis);
        if (norm_axis < std::numeric_limits<double>::epsilon())
            KRATOS_ERROR<<"SymmetryUtility: norm of axis is close to zero"<<std::endl;
            
        mAxisSymmetryData.Axis = axis/norm_axis;
        mAxisSymmetryData.Angle = mSymmetrySettings["settings"]["angle"].GetDouble();        
    }
    else
        KRATOS_ERROR<<"SymmetryUtility: type should be either plane or axis"<<std::endl;

}



void SymmetryUtility::Initialize()
{

    NodeVector mListOfNodesInModelPart;
    mListOfNodesInModelPart.resize(mrModelPart.Nodes().size());
    int counter = 0;
    for (ModelPart::NodesContainerType::iterator node_it = mrModelPart.NodesBegin(); node_it != mrModelPart.NodesEnd(); ++node_it)
    {
        NodeTypePointer pnode = *(node_it.base());
        mListOfNodesInModelPart[counter++] = pnode;
    }

    const size_t bucket_size = 100;
    KDTree search_tree(mListOfNodesInModelPart.begin(), mListOfNodesInModelPart.end(), bucket_size);

    if(mAxisSymmetry){
        mAxisSymmetryData.Map.resize(mrModelPart.Nodes().size());
        int counter = 0;
        for (auto& r_node : mrModelPart.Nodes()){
            int num_rotations = (360.0/mAxisSymmetryData.Angle);
            NodeVector p_neighbor_rot_nodes;
            p_neighbor_rot_nodes.resize(num_rotations-1);
            for(int r_i=1;r_i<num_rotations;r_i++){
               NodeTypePointer p_rot_node = GetRotatedNode(r_node,r_i*mAxisSymmetryData.Angle);
               double distance;
               NodeTypePointer p_neighbor_rot_node = search_tree.SearchNearestPoint(*p_rot_node, distance);
               p_neighbor_rot_nodes[r_i-1] = p_neighbor_rot_node;
            }

            mAxisSymmetryData.Map[counter] = p_neighbor_rot_nodes;
            counter++;
        }
    }
    if(mPlaneSymmetry){
        mPlaneSymmetryData.Map.resize(mrModelPart.Nodes().size());
        int counter = 0;
        // #pragma omp parallel for reduction(+:counter)
        for (auto& r_node : mrModelPart.Nodes()){
            NodeTypePointer p_ref_node = GetReflectedNode(r_node);
            double distance;
            NodeTypePointer p_neighbor_ref_node = search_tree.SearchNearestPoint(*p_ref_node, distance);
            mPlaneSymmetryData.Map[counter] =  p_neighbor_ref_node;
            counter ++;          
        }
    }

}

SymmetryUtility::NodeTypePointer SymmetryUtility::GetRotatedNode(NodeType& rNode, double angle){
    NodeTypePointer p_new_node = Kratos::make_intrusive<NodeType>(rNode.Id(), rNode[0], rNode[1], rNode[2]);

    Matrix rot_mat = ZeroMatrix(3,3);
    const array_3d v1 = rNode.Coordinates()-mAxisSymmetryData.Point;

    const double c=cos(angle*PI/180);
    const double s=sin(angle*PI/180);
    const double t = 1-c;
    const array_3d& k = mAxisSymmetryData.Axis;

    rot_mat(0,0) = t*k[0]*k[0] + c;         rot_mat(0,1) = t*k[0]*k[1] - k[2]*s;    rot_mat(0,2) = t*k[0]*k[2] + k[1]*s;
    rot_mat(1,0) = t*k[0]*k[1] + k[2]*s;    rot_mat(1,1) = t*k[1]*k[1] + c;         rot_mat(1,2) = t*k[1]*k[2] - k[0]*s;
    rot_mat(2,0) = t*k[0]*k[2] - k[1]*s;    rot_mat(2,1) = t*k[1]*k[2] + k[0]*s;    rot_mat(2,2) = t*k[2]*k[2] + c;

    p_new_node->Coordinates() = prod(rot_mat, v1) + mAxisSymmetryData.Point;

    return p_new_node;
}

SymmetryUtility::NodeTypePointer SymmetryUtility::GetReflectedNode(NodeType& rNode){
    NodeTypePointer p_new_node = Kratos::make_intrusive<NodeType>(rNode.Id(), rNode[0], rNode[1], rNode[2]);
    const array_3d v = rNode.Coordinates() - mPlaneSymmetryData.Point;
    p_new_node->Coordinates() = prod(mPlaneSymmetryData.ReflectionMatrix, v) + mPlaneSymmetryData.Point;
    return p_new_node;
}

void SymmetryUtility::Update()
{
    // for(auto& node_i : mrModelPartToDamp.Nodes())
    // {
    //     node_i.SetValue(DAMPING_FACTOR_X,1.0);
    //     node_i.SetValue(DAMPING_FACTOR_Y,1.0);
    //     node_i.SetValue(DAMPING_FACTOR_Z,1.0);
    // }
}

// void SymmetryUtility::SetDampingFactorsForAllDampingRegions()
// {
//     KRATOS_INFO("") << std::endl;

//     BuiltinTimer timer;
//     KRATOS_INFO("ShapeOpt") << "Starting to prepare damping..." << std::endl;

//     // Loop over all regions for which damping is to be applied
//     for (const auto& r_region_parameters : mSymmetrySettings["damping_regions"])
//     {
//         const std::string& sub_model_part_name = r_region_parameters["sub_model_part_name"].GetString();
//         ModelPart& dampingRegion = mrModelPartToDamp.GetRootModelPart().GetSubModelPart(sub_model_part_name);

//         const bool dampX = r_region_parameters["damp_X"].GetBool();
//         const bool dampY = r_region_parameters["damp_Y"].GetBool();
//         const bool dampZ = r_region_parameters["damp_Z"].GetBool();
//         const auto& dampingFunctionType = r_region_parameters["damping_function_type"].GetString();
//         const double dampingRadius = r_region_parameters["damping_radius"].GetDouble();

//         const auto p_damping_function = CreateDampingFunction(dampingFunctionType, dampingRadius );

//         // Loop over all nodes in specified damping sub-model part
//         block_for_each(dampingRegion.Nodes(), [&](const ModelPart::NodeType& rNode) {
//             NodeVector neighbor_nodes( mMaxNeighborNodes );
//             const unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( rNode,
//                                                                                 dampingRadius,
//                                                                                 neighbor_nodes.begin(),
//                                                                                 mMaxNeighborNodes );

//             ThrowWarningIfNodeNeighborsExceedLimit( rNode, number_of_neighbors );

//             // Loop over all nodes in radius (including node on damping region itself)
//             for(unsigned int j_itr = 0 ; j_itr<number_of_neighbors ; j_itr++)
//             {
//                 ModelPart::NodeType& neighbor_node = *neighbor_nodes[j_itr];
//                 const double dampingFactor = 1.0 - p_damping_function->ComputeWeight( rNode.Coordinates(), neighbor_node.Coordinates());

//                 // For every specified damping direction we check if new damping factor is smaller than the assigned one for current node.
//                 // In case yes, we overwrite the value. This ensures that the damping factor of a node is computed by its closest distance to the damping region
//                 array_3d& damping_factor_variable = neighbor_node.GetValue(DAMPING_FACTOR);
//                 neighbor_node.SetLock();
//                 if(dampX && dampingFactor < damping_factor_variable[0])
//                     damping_factor_variable[0] = dampingFactor;
//                 if(dampY && dampingFactor < damping_factor_variable[1])
//                     damping_factor_variable[1] = dampingFactor;
//                 if(dampZ && dampingFactor < damping_factor_variable[2])
//                     damping_factor_variable[2] = dampingFactor;
//                 neighbor_node.UnSetLock();
//             }
//         });
//     }

//     KRATOS_INFO("ShapeOpt") << "Finished preparation of damping in: " << timer.ElapsedSeconds() << " s" << std::endl;
// }

// FilterFunction::Pointer SymmetryUtility::CreateDampingFunction( std::string damping_type, double damping_radius ) const
// {
//     return Kratos::make_unique<FilterFunction>(damping_type, damping_radius);
// }

// void SymmetryUtility::ThrowWarningIfNodeNeighborsExceedLimit( const ModelPart::NodeType& given_node, const unsigned int number_of_neighbors ) const
// {
//     if(number_of_neighbors >= mMaxNeighborNodes)
//         KRATOS_WARNING("ShapeOpt::SymmetryUtility") << "For node " << given_node.Id() << " and specified damping radius, maximum number of neighbor nodes (=" << mMaxNeighborNodes << " nodes) reached!" << std::endl;
// }

void SymmetryUtility::ApplyOnVectorField( const Variable<array_3d> &rNodalVariable )
{
    // block_for_each(mrModelPartToDamp.Nodes(), [&](ModelPart::NodeType& rNode) {
    //     const array_3d& damping_factor = rNode.GetValue(DAMPING_FACTOR);
    //     array_3d& nodalVariable = rNode.FastGetSolutionStepValue(rNodalVariable);

    //     nodalVariable[0] *= damping_factor[0];
    //     nodalVariable[1] *= damping_factor[1];
    //     nodalVariable[2] *= damping_factor[2];
    // });
}

void SymmetryUtility::ApplyOnScalarField( const Variable<double> &rNodalVariable )
{
    Vector aux_field = ZeroVector(mrModelPart.Nodes().size());

    if(mAxisSymmetry){        
        int counter = 0;
        for (auto& r_node : mrModelPart.Nodes()){
            for(int p_rot_c= 0; p_rot_c<mAxisSymmetryData.Map[counter].size();p_rot_c++)
                aux_field[counter] += mAxisSymmetryData.Map[counter][p_rot_c]->FastGetSolutionStepValue(rNodalVariable);                
;    
            counter ++;          
        }
        counter = 0;
        for (auto& r_node : mrModelPart.Nodes()){
            r_node.FastGetSolutionStepValue(rNodalVariable) += aux_field[counter];
            r_node.FastGetSolutionStepValue(rNodalVariable) /= (mAxisSymmetryData.Map[counter].size()+1);
            counter ++;          
        }
    }


    if(mPlaneSymmetry){        
        int counter = 0;
        for (auto& r_node : mrModelPart.Nodes()){
            aux_field[counter] = mPlaneSymmetryData.Map[counter]->FastGetSolutionStepValue(rNodalVariable);
            counter ++;          
        }
        counter = 0;
        for (auto& r_node : mrModelPart.Nodes()){
            r_node.FastGetSolutionStepValue(rNodalVariable) += aux_field[counter];
            r_node.FastGetSolutionStepValue(rNodalVariable) /= 2.0;
            counter ++;          
        }
    }

}

} // namespace Kratos
