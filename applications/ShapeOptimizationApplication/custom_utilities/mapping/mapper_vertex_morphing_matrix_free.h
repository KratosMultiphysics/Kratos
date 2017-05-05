// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_MATRIX_FREE_H
#define MAPPER_VERTEX_MORPHING_MATRIX_FREE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "../../kratos/includes/define.h"
#include "../../kratos/processes/process.h"
#include "../../kratos/includes/node.h"
#include "../../kratos/includes/element.h"
#include "../../kratos/includes/model_part.h"
#include "../../kratos/includes/kratos_flags.h"
#include "../../kratos/spatial_containers/spatial_containers.h"
#include "../../kratos/utilities/timer.h"
#include "../../kratos/processes/node_erase_process.h"
#include "../../kratos/utilities/binbased_fast_point_locator.h"
#include "../../kratos/utilities/normal_calculation_utils.h"
#include "../../kratos/spaces/ublas_space.h"
#include "shape_optimization_application.h"
#include "filter_function.h"
//#include "damping_function.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    typedef boost::python::extract<ModelPart&> extractModelPart;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

*/

class MapperVertexMorphingMatrixFree
{
public:
    ///@name Type Definitions
    ///@{

    // Type definitions for better reading later
    typedef array_1d<double,3> array_3d;
    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double> DoubleVector;
    typedef std::vector<double>::iterator DoubleVectorIterator;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    // Type definitions for linear algebra including sparse systems
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef SparseSpaceType::MatrixType SparseMatrixType;
    typedef SparseSpaceType::VectorType VectorType;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;    

    /// Pointer definition of MapperVertexMorphingMatrixFree
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingMatrixFree);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingMatrixFree( ModelPart& designSurface, boost::python::dict dampingRegions, Parameters& optimizationSettings )
        : mrDesignSurface( designSurface ),
          mNumberOfDesignVariables( designSurface.Nodes().size() ),
          mFilterRadius( optimizationSettings["design_variables"]["filter"]["filter_radius"].GetDouble() ),
          mFilterType( optimizationSettings["design_variables"]["filter"]["filter_function_type"].GetString() )
          //mPerformDamping( optimizationSettings["design_variables"]["damping"]["perform_damping"].GetBool() )
    {
        boost::timer timer;        
        std::cout << "\n> Creating search tree..." << std::endl;
        CreateListOfNodesOfDesignSurface();
        CreateSearchTreeWithAllNodesOnDesignSurface();
        std::cout << "> Search tree created in: " << timer.elapsed() << " s" << std::endl;

        CreateFilterFunction();
        AssignMappingIds();

        //InitalizeDampingFactorsToHaveNoInfluence(); 
        //if(mPerformDamping)   
        //    SetDampingFactorsForAllDampingRegions( dampingRegions, optimizationSettings["design_variables"]["damping"]["damping_regions"] );         
    }

    /// Destructor.
    virtual ~MapperVertexMorphingMatrixFree()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void CreateListOfNodesOfDesignSurface()
    {
        for (ModelPart::NodesContainerType::iterator node_it = mrDesignSurface.NodesBegin(); node_it != mrDesignSurface.NodesEnd(); ++node_it)
        {
            NodeTypePointer pnode = *(node_it.base());
            mListOfNodesOfDesignSurface.push_back(pnode);
        }
    }

    // --------------------------------------------------------------------------
    void CreateSearchTreeWithAllNodesOnDesignSurface()
    {
        mpSearchTree = boost::shared_ptr<KDTree>(new KDTree(mListOfNodesOfDesignSurface.begin(), mListOfNodesOfDesignSurface.end(), mBucketSize));
    }   

    // --------------------------------------------------------------------------
    void CreateFilterFunction()
    {
        mpFilterFunction = boost::shared_ptr<FilterFunction>(new FilterFunction(mFilterType, mFilterRadius));
    }     

    // --------------------------------------------------------------------------
    void AssignMappingIds()
    {
        unsigned int i = 0;
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
            node_i->SetValue(MAPPING_ID,i++);
    }

    // --------------------------------------------------------------------------
/*    void initalizeDampingFactorsToHaveNoInfluence()


    {
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            node_i->SetValue(DAMPING_FACTOR_X,1.0);    
            node_i->SetValue(DAMPING_FACTOR_Y,1.0);  
            node_i->SetValue(DAMPING_FACTOR_Z,1.0);  
        } 
    }
*/
    // --------------------------------------------------------------------------
/*    void setDampingFactorsForAllDampingRegions( boost::python::dict dampingRegions, Parameters dampingSettings )
    {
        std::cout << "\n> Starting to prepare damping..." << std::endl;

        UpdateSearchTreeIfGeometryHasChanged();        

        // Loop over all regions for which damping is to be applied
        for (unsigned int regionNumber = 0; regionNumber < len(dampingRegions); regionNumber++)
        {
            // Extract sub-model part for damping
            std::string dampingRegionSubModelPartName = dampingSettings[regionNumber]["sub_model_part_name"].GetString();
            ModelPart& dampingRegion = extractModelPart( dampingRegions[dampingRegionSubModelPartName] );

            // Read settings for current edge
            bool dampX = dampingSettings[regionNumber]["damp_X"].GetBool();
            bool dampY = dampingSettings[regionNumber]["damp_Y"].GetBool();
            bool dampZ = dampingSettings[regionNumber]["damp_Z"].GetBool();
            std::string dampingFunctionType = dampingSettings[regionNumber]["damping_function_type"].GetString();
            double dampingRadius = dampingSettings[regionNumber]["damping_radius"].GetDouble();
            DampingFunction::Pointer mpDampingFunction = CreateDampingFunction( dampingFunctionType, dampingRadius );

            // Loop over all nodes in specified damping sub-model part 
            for ( ModelPart::NodeIterator node_itr = dampingRegion.NodesBegin(); node_itr != dampingRegion.NodesEnd(); ++node_itr )
            {
                ModelPart::NodeType& currentDampingNode = *node_itr;
                NodeVector neighbor_nodes( mMaxNeighborNodes );
                DoubleVector resulting_squared_distances( mMaxNeighborNodes,0.0 );
                unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( currentDampingNode,
                                                                                 mFilterRadius, 
                                                                                 neighbor_nodes.begin(),
                                                                                 resulting_squared_distances.begin(), 
                                                                                 mMaxNeighborNodes );

                ThrowWarningIfMaxNodeNeighborsReached( currentDampingNode, number_of_neighbors );                                                                               

                // Loop over all nodes in radius (including node on damping region itself)
                for(unsigned int j_itr = 0 ; j_itr<number_of_neighbors ; j_itr++)
                {
                    ModelPart::NodeType& neighbor_node = *neighbor_nodes[j_itr];
                    double dampingFactor = mpDampingFunction->compute_damping_factor( currentDampingNode.Coordinates(), neighbor_node.Coordinates());

                    // For every specified damping direction we check if new damping factor is smaller than the assigned one for current node. 
                    // In case yes, we overwrite the value. This ensures that the damping factor of a node is computed by its closest distance to the damping region
                    auto& damping_factor_variable = neighbor_node.GetValue(DAMPING_FACTOR);
                    if(dampX && dampingFactor < damping_factor_variable[0])
                        damping_factor_variable[0] = dampingFactor;
                    if(dampY && dampingFactor < damping_factor_variable[1])       
                        damping_factor_variable[1] = dampingFactor;   
                    if(dampZ && dampingFactor < damping_factor_variable[2])       
                        damping_factor_variable[2] = dampingFactor;                            
                }
            }
        }

        std::cout << "> Finished preparation of damping." << std::endl;
    }
*/
    // --------------------------------------------------------------------------
/*    DampingFunction::Pointer CreateDampingFunction( std::string damping_type, double damping_radius )
    {
        return boost::shared_ptr<DampingFunction>(new DampingFunction(damping_type, damping_radius));
    }    
*/
    // --------------------------------------------------------------------------
/*    void MapToDesignSpace( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInDesignSpace )
    {
        //if(mPerformDamping)
        //    dampNodalVariable( rNodalVariable );
        PerformMappingToDesignSpace( rNodalVariable, rNodalVariableInDesignSpace );
    }
*/
    // --------------------------------------------------------------------------
/*    void MapToGeometrySpace( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInGeometrySpace )
    {
        PerformMappingToGeometrySpace( rNodalVariable, rNodalVariableInGeometrySpace );
        //if(mPerformDamping)
        //    dampNodalVariable( rNodalVariable );
    }
*/
    // --------------------------------------------------------------------------
    void MapToDesignSpace( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInDesignSpace )
    {
        boost::timer mapping_time;
        std::cout << "\n> Starting to map " << rNodalVariable.Name() << " to design space..." << std::endl;

        VectorType x_variables_in_design_space, y_variables_in_design_space, z_variables_in_design_space;
        x_variables_in_design_space.resize(mNumberOfDesignVariables);
        y_variables_in_design_space.resize(mNumberOfDesignVariables);
        z_variables_in_design_space.resize(mNumberOfDesignVariables);       
        
        UpdateSearchTreeIfGeometryHasChanged();

        MapComponentwiseVariablesToDesignSpace( rNodalVariable,
                                                x_variables_in_design_space, 
                                                y_variables_in_design_space, 
                                                z_variables_in_design_space );

        AssignComponentwiseMappedVariablesToNodes( x_variables_in_design_space, 
                                                   y_variables_in_design_space, 
                                                   z_variables_in_design_space,
                                                   rNodalVariableInDesignSpace );

        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void MapToGeometrySpace( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInGeometrySpace )
    {
        boost::timer mapping_time;
        std::cout << "\n> Starting to map " << rNodalVariable.Name() << " to geometry space..." << std::endl;

        VectorType x_variables_in_geometry_space, y_variables_in_geometry_space, z_variables_in_geometry_space;
        x_variables_in_geometry_space.resize(mNumberOfDesignVariables);
        y_variables_in_geometry_space.resize(mNumberOfDesignVariables);
        z_variables_in_geometry_space.resize(mNumberOfDesignVariables);

        UpdateSearchTreeIfGeometryHasChanged();

        MapComponentwiseVariablesToGeometrySpace( rNodalVariable,
                                                  x_variables_in_geometry_space, 
                                                  y_variables_in_geometry_space, 
                                                  z_variables_in_geometry_space );

        AssignComponentwiseMappedVariablesToNodes( x_variables_in_geometry_space, 
                                                   y_variables_in_geometry_space, 
                                                   z_variables_in_geometry_space,
                                                   rNodalVariableInGeometrySpace );

        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
/*    void dampNodalVariable( const Variable<array_3d> &rNodalVariable )
    {
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {   
            auto& damping_factor = node_i->GetValue(DAMPING_FACTOR);
            auto& nodalVariable = node_i->FastGetSolutionStepValue(rNodalVariable);

            nodalVariable[0] *= damping_factor[0];
            nodalVariable[1] *= damping_factor[1];
            nodalVariable[2] *= damping_factor[2];  
        }  
    } 
*/
    // --------------------------------------------------------------------------
    void UpdateSearchTreeIfGeometryHasChanged()
    {
        if(HasGeometryChanged())
        {
            DeleteSearchTree();
            CreateSearchTreeWithAllNodesOnDesignSurface();
        }
    }

    // --------------------------------------------------------------------------
    bool HasGeometryChanged()
    {
        return true;
    }

    // --------------------------------------------------------------------------
    void DeleteSearchTree()
    {
        mpSearchTree.reset();
    }     

    // --------------------------------------------------------------------------
    void MapComponentwiseVariablesToDesignSpace( const Variable<array_3d>& rNodalVariable,
                                                 VectorType& x_variables_in_design_space, 
                                                 VectorType& y_variables_in_design_space, 
                                                 VectorType& z_variables_in_design_space )
    {
        for (ModelPart::NodeIterator node_itr = mrDesignSurface.NodesBegin(); node_itr != mrDesignSurface.NodesEnd(); ++node_itr)
        {
            ModelPart::NodeType& design_node = *node_itr;

            NodeVector neighbor_nodes( mMaxNeighborNodes );
            DoubleVector resulting_squared_distances( mMaxNeighborNodes );
            unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( design_node,
                                                                             mFilterRadius, 
                                                                             neighbor_nodes.begin(),
                                                                             resulting_squared_distances.begin(), 
                                                                             mMaxNeighborNodes );

            ThrowWarningIfMaxNodeNeighborsReached( design_node, number_of_neighbors );                                                                               

            DoubleVector list_of_weights( number_of_neighbors, 0.0 );
            double sum_of_weights = 0.0;
            ComputeWeightForAllNeighbors( design_node, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
            
            // Local transpose mapping including post-scaling
            array_3d& nodal_variable = design_node.FastGetSolutionStepValue(rNodalVariable);
            for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
            {
                ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
                int neighbor_node_mapping_id = neighbor_node.GetValue(MAPPING_ID);

                double weight = list_of_weights[neighbor_itr] / sum_of_weights;

                x_variables_in_design_space[neighbor_node_mapping_id] += weight*nodal_variable[0];
                y_variables_in_design_space[neighbor_node_mapping_id] += weight*nodal_variable[1];
                z_variables_in_design_space[neighbor_node_mapping_id] += weight*nodal_variable[2];
            }                          
        } 
    }   

    // --------------------------------------------------------------------------
    void MapComponentwiseVariablesToGeometrySpace( const Variable<array_3d>& rNodalVariable,
                                                   VectorType& x_variables_in_geometry_space, 
                                                   VectorType& y_variables_in_geometry_space, 
                                                   VectorType& z_variables_in_geometry_space )
    {
        for (ModelPart::NodeIterator node_itr = mrDesignSurface.NodesBegin(); node_itr != mrDesignSurface.NodesEnd(); ++node_itr)
        {
            ModelPart::NodeType& design_node = *node_itr;

            NodeVector neighbor_nodes(mMaxNeighborNodes);
            DoubleVector resulting_squared_distances(mMaxNeighborNodes);
            unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( design_node,
                                                                             mFilterRadius, 
                                                                             neighbor_nodes.begin(),
                                                                             resulting_squared_distances.begin(), 
                                                                             mMaxNeighborNodes );

            ThrowWarningIfMaxNodeNeighborsReached( design_node, number_of_neighbors );                                                                               
            
            DoubleVector list_of_weights( number_of_neighbors, 0.0 );
            double sum_of_weights = 0.0;
            ComputeWeightForAllNeighbors( design_node, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );

            // Local mapping including post-scaling
            int design_node_mapping_id = design_node.GetValue(MAPPING_ID);
            for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
            {
                double weight = list_of_weights[neighbor_itr] / sum_of_weights;

                ModelPart::NodeType& node_j = *neighbor_nodes[neighbor_itr];
                array_3d& nodal_variable = node_j.FastGetSolutionStepValue(rNodalVariable);

                x_variables_in_geometry_space[design_node_mapping_id] += weight*nodal_variable[0];
                y_variables_in_geometry_space[design_node_mapping_id] += weight*nodal_variable[1];
                z_variables_in_geometry_space[design_node_mapping_id] += weight*nodal_variable[2];
            }                          
        }        
    }

    // --------------------------------------------------------------------------
    void AssignComponentwiseMappedVariablesToNodes( VectorType& vector_field_x, 
                                                    VectorType& vector_field_y, 
                                                    VectorType& vector_field_z,
                                                    const Variable<array_3d> &rNodalVariable )
    {
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            int i = node_i->GetValue(MAPPING_ID);

            VectorType node_vector = ZeroVector(3);
            node_vector(0) = vector_field_x[i];
            node_vector(1) = vector_field_y[i];
            node_vector(2) = vector_field_z[i];
            node_i->FastGetSolutionStepValue(rNodalVariable) = node_vector;
        }
    }    

    // --------------------------------------------------------------------------
    void ThrowWarningIfMaxNodeNeighborsReached( ModelPart::NodeType& given_node, unsigned int number_of_neighbors )
    {
        if(number_of_neighbors >= mMaxNeighborNodes)
            std::cout << "\n> WARNING!!!!! For node " << given_node.Id() << " and specified filter radius, maximum number of neighbor nodes (=" << mMaxNeighborNodes << " nodes) reached!" << std::endl;
    }

    // --------------------------------------------------------------------------
    void ComputeWeightForAllNeighbors(  ModelPart::NodeType& design_node, 
                                        NodeVector& neighbor_nodes, 
                                        unsigned int number_of_neighbors,
                                        DoubleVector& list_of_weights, 
                                        double& sum_of_weights )
    {
        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            double weight = mpFilterFunction->compute_weight( design_node.Coordinates(), neighbor_node.Coordinates() );

            list_of_weights[neighbor_itr] = weight;
            sum_of_weights += weight;
        }
    }

    // ==============================================================================

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "MapperVertexMorphingMatrixFree";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MapperVertexMorphingMatrixFree";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    ModelPart& mrDesignSurface;
    const unsigned int mNumberOfDesignVariables;
    double mFilterRadius;
    std::string mFilterType;
    FilterFunction::Pointer mpFilterFunction;
    //bool mPerformDamping;

    // ==============================================================================
    // Variables for spatial search
    // ==============================================================================
    unsigned int mBucketSize = 100;
    unsigned int mMaxNeighborNodes = 10000;            
    NodeVector mListOfNodesOfDesignSurface;
    KDTree::Pointer mpSearchTree;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
//      MapperVertexMorphingMatrixFree& operator=(MapperVertexMorphingMatrixFree const& rOther);

    /// Copy constructor.
//      MapperVertexMorphingMatrixFree(MapperVertexMorphingMatrixFree const& rOther);


    ///@}

}; // Class MapperVertexMorphingMatrixFree

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_MATRIX_FREE_H
