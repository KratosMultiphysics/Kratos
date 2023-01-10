

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

    BuiltinTimer timer;
    KRATOS_INFO("SymmetryUtility::Initialize: ") << "Starting initialization of symmetry utility "<<mUtilName<<" ..." << std::endl;

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

        #pragma omp declare reduction (merge : std::vector<std::pair <NodeTypePointer,NodeVector>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))        
        std::vector<std::pair <NodeTypePointer,NodeVector>>& rMap = mAxisSymmetryData.Map;
        #pragma omp parallel for reduction(merge: rMap)
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
            std::pair <NodeTypePointer,NodeVector> nodes_pair;
            nodes_pair = std::make_pair(mrModelPart.pGetNode(r_node.Id()),p_neighbor_rot_nodes);
            rMap.push_back(nodes_pair); 
        }
    }

    if(mPlaneSymmetry){
        #pragma omp declare reduction (merge : std::vector<std::pair <NodeTypePointer,NodeTypePointer>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))        
        std::vector<std::pair <NodeTypePointer,NodeTypePointer>>& rMap = mPlaneSymmetryData.Map;
        #pragma omp parallel for reduction(merge: rMap)
        for (auto& r_node : mrModelPart.Nodes()){
            NodeTypePointer p_ref_node = GetReflectedNode(r_node);
            double distance;
            NodeTypePointer p_neighbor_ref_node = search_tree.SearchNearestPoint(*p_ref_node, distance);
            std::pair <NodeTypePointer,NodeTypePointer> nodes_pair;
            nodes_pair = std::make_pair(mrModelPart.pGetNode(r_node.Id()),p_neighbor_ref_node);
            rMap.push_back(nodes_pair);          
        }
    }

    KRATOS_INFO("SymmetryUtility::Initialize: ") << " Finished initialization of symmetry utility "<<mUtilName<<" in " << timer.ElapsedSeconds() << " s." << std::endl;    
}

SymmetryUtility::NodeTypePointer SymmetryUtility::GetRotatedNode(NodeType& rNode, double angle){
    NodeTypePointer p_new_node = Kratos::make_intrusive<NodeType>(rNode.Id(), rNode[0], rNode[1], rNode[2]);

    const array_3d v1 = rNode.Coordinates()-mAxisSymmetryData.Point;
    Matrix rot_mat;
    GetRotationMatrix(angle,rot_mat);
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

}

void SymmetryUtility::GetRotationMatrix(double angle, Matrix& rRotMat){
    rRotMat.resize(3,3);
    rRotMat = ZeroMatrix(3,3);

    const double c=cos(angle*PI/180);
    const double s=sin(angle*PI/180);
    const double t = 1-c;
    const array_3d& k = mAxisSymmetryData.Axis;

    rRotMat(0,0) = t*k[0]*k[0] + c;         rRotMat(0,1) = t*k[0]*k[1] - k[2]*s;    rRotMat(0,2) = t*k[0]*k[2] + k[1]*s;
    rRotMat(1,0) = t*k[0]*k[1] + k[2]*s;    rRotMat(1,1) = t*k[1]*k[1] + c;         rRotMat(1,2) = t*k[1]*k[2] - k[0]*s;
    rRotMat(2,0) = t*k[0]*k[2] - k[1]*s;    rRotMat(2,1) = t*k[1]*k[2] + k[0]*s;    rRotMat(2,2) = t*k[2]*k[2] + c;
}

void SymmetryUtility::ApplyOnVectorField( const Variable<array_3d> &rNodalVariable )
{
    if(mAxisSymmetry){        
        std::vector<Vector> aux;
        #pragma omp parallel
        {
            std::vector<Vector> aux_private;
            #pragma omp for nowait schedule(static)
            for(auto& r_pair : mAxisSymmetryData.Map) {
                Vector vec_val = r_pair.first->FastGetSolutionStepValue(rNodalVariable);
                Matrix rot_mat;
                for(int n_c = 0; n_c<r_pair.second.size();n_c++){
                    GetRotationMatrix((n_c+1)*mAxisSymmetryData.Angle,rot_mat);
                    vec_val += prod(trans(rot_mat),r_pair.second[n_c]->FastGetSolutionStepValue(rNodalVariable));
                }                
                    
                vec_val /= (r_pair.second.size()+1);
                aux_private.push_back(vec_val);
            }
            #pragma omp for schedule(static) ordered
            for(int i=0; i<omp_get_num_threads(); i++) {
                #pragma omp ordered
                aux.insert(aux.end(), aux_private.begin(), aux_private.end());
            }
        }

        #pragma omp parallel for 
        for(int i=0;i<aux.size();i++)
            mAxisSymmetryData.Map[i].first->FastGetSolutionStepValue(rNodalVariable) = aux[i];
    }

    if(mPlaneSymmetry){  

        std::vector<Vector> aux;
        #pragma omp parallel
        {
            std::vector<Vector> aux_private;
            #pragma omp for nowait schedule(static)
            for(auto& r_pair : mPlaneSymmetryData.Map) {
                Vector val = r_pair.first->FastGetSolutionStepValue(rNodalVariable);
                Vector ref_node_val = r_pair.second->FastGetSolutionStepValue(rNodalVariable);
                ref_node_val = prod(trans(mPlaneSymmetryData.ReflectionMatrix),ref_node_val);
                val += ref_node_val;
                val /= 2.0;
                aux_private.push_back(val);
            }
            #pragma omp for schedule(static) ordered
            for(int i=0; i<omp_get_num_threads(); i++) {
                #pragma omp ordered
                aux.insert(aux.end(), aux_private.begin(), aux_private.end());
            }
        }

        #pragma omp parallel for 
        for(int i=0;i<aux.size();i++)
            mPlaneSymmetryData.Map[i].first->FastGetSolutionStepValue(rNodalVariable) = aux[i];
        
    }
}

void SymmetryUtility::ApplyOnScalarField( const Variable<double> &rNodalVariable )
{    

    if(mAxisSymmetry){        
        std::vector<double> aux;
        #pragma omp parallel
        {
            std::vector<double> aux_private;
            #pragma omp for nowait schedule(static)
            for(auto& r_pair : mAxisSymmetryData.Map) {
                double val = r_pair.first->FastGetSolutionStepValue(rNodalVariable);
                for(int n_c = 0; n_c<r_pair.second.size();n_c++)                
                    val += r_pair.second[n_c]->FastGetSolutionStepValue(rNodalVariable);
                val /= (r_pair.second.size()+1);
                aux_private.push_back(val);
            }
            #pragma omp for schedule(static) ordered
            for(int i=0; i<omp_get_num_threads(); i++) {
                #pragma omp ordered
                aux.insert(aux.end(), aux_private.begin(), aux_private.end());
            }
        }

        #pragma omp parallel for 
        for(int i=0;i<aux.size();i++)
            mAxisSymmetryData.Map[i].first->FastGetSolutionStepValue(rNodalVariable) = aux[i];
    }

    if(mPlaneSymmetry){  

        std::vector<double> aux;
        #pragma omp parallel
        {
            std::vector<double> aux_private;
            #pragma omp for nowait schedule(static)
            for(auto& r_pair : mPlaneSymmetryData.Map) {
                double val = r_pair.first->FastGetSolutionStepValue(rNodalVariable);
                val += r_pair.second->FastGetSolutionStepValue(rNodalVariable);
                val /= 2.0;
                aux_private.push_back(val);
            }
            #pragma omp for schedule(static) ordered
            for(int i=0; i<omp_get_num_threads(); i++) {
                #pragma omp ordered
                aux.insert(aux.end(), aux_private.begin(), aux_private.end());
            }
        }

        #pragma omp parallel for 
        for(int i=0;i<aux.size();i++)
            mPlaneSymmetryData.Map[i].first->FastGetSolutionStepValue(rNodalVariable) = aux[i];
        
    }

}

} // namespace Kratos
