// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_H
#define MAPPER_VERTEX_MORPHING_H

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

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

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

class MapperVertexMorphing
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

    /// Pointer definition of MapperVertexMorphing
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphing);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphing( ModelPart& designSurface, Parameters& optimizationSettings )
        : mrDesignSurface( designSurface ),
          mNumberOfDesignVariables(designSurface.Nodes().size()),
          mFilterRadius( optimizationSettings["design_variables"]["filter"]["filter_radius"].GetDouble() ),
          mFilterType( optimizationSettings["design_variables"]["filter"]["filter_function_type"].GetString() )
    {
        CreateListOfNodesOfDesignSurface();
        InitializeMappingMatrix();                
        ComputeMappingMatrix();
        CreateFilterFunction();
        AssignMappingIds();
    }

    /// Destructor.
    virtual ~MapperVertexMorphing()
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
    void InitializeMappingMatrix()
    {
        mMappingMatrix.resize(mNumberOfDesignVariables,mNumberOfDesignVariables);
        mMappingMatrix.clear();    
    }    

    // --------------------------------------------------------------------------
    void ComputeMappingMatrix()
    {
        boost::timer timer;        
        std::cout << "> Creating mapping matrix to perform mapping..." << std::endl;        
        CreateSearchTreeWithAllNodesOnDesignSurface();
        ComputeEntriesOfMappingMatrix(); 
        std::cout << "> Mapping matrix created in: " << timer.elapsed() << " s" << std::endl;
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
    void ComputeEntriesOfMappingMatrix()
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

            DoubleVector list_of_weights( number_of_neighbors, 0.0 );
            double sum_of_weights = 0.0;
            
            ThrowWarningIfMaxNodeNeighborsReached( design_node, number_of_neighbors );                                                                               
            ComputeWeightForAllNeighbors( design_node, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
            FillMappingMatrixWithWeights( design_node, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
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

    // --------------------------------------------------------------------------
    void FillMappingMatrixWithWeights(  ModelPart::NodeType& design_node, 
                                        NodeVector& neighbor_nodes, 
                                        unsigned int number_of_neighbors,
                                        DoubleVector& list_of_weights, 
                                        double& sum_of_weights )
    {
        unsigned int row_id = design_node.GetValue(MAPPING_ID);
        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            int collumn_id = neighbor_node.GetValue(MAPPING_ID);

            double weight = list_of_weights[neighbor_itr] / sum_of_weights;
            mMappingMatrix.push_back(row_id,collumn_id,weight);
        }        
    }

    // --------------------------------------------------------------------------
    void MapToDesignSpace( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInDesignSpace )
    {
        boost::timer mapping_time;
        std::cout << "\n> Starting to map " << rNodalVariable.Name() << " to design space..." << std::endl;

        RecomputeMappingMatrixIfGeometryHasChanged();

    //     // First we apply edge damping to the sensitivities if specified
    //     if(mPerformDamping)
    //         for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
    //         {
    //             node_i->FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY_X) *= node_i->GetValue(DAMPING_FACTOR_X);
    //             node_i->FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY_Y) *= node_i->GetValue(DAMPING_FACTOR_Y);
    //             node_i->FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY_Z) *= node_i->GetValue(DAMPING_FACTOR_Z);
    //             if(constraint_given)
    //             {
    //                 node_i->FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY_X) *= node_i->GetValue(DAMPING_FACTOR_X);
    //                 node_i->FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY_Y) *= node_i->GetValue(DAMPING_FACTOR_Y);
    //                 node_i->FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY_Z) *= node_i->GetValue(DAMPING_FACTOR_Z);                    
    //             }                
    //         }

    //     // Map objective sensitivities

    //     // Assign nodal sensitivities to vectors used for mapping
    //     VectorType dJdX, dJdY, dJdZ;
    //     dJdX.resize(mNumberOfDesignVariables);
    //     dJdY.resize(mNumberOfDesignVariables);
    //     dJdZ.resize(mNumberOfDesignVariables);
    //     for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
    //     {
    //         int i = node_i->GetValue(MAPPING_ID);
    //         array_3d& node_sens = node_i->FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY);
    //         dJdX[i] = node_sens[0];
    //         dJdY[i] = node_sens[1];
    //         dJdZ[i] = node_sens[2];
    //     }
    //     VectorType dCdX, dCdY, dCdZ;
    //     if(constraint_given)
    //     {
    //         dCdX.resize(mNumberOfDesignVariables);
    //         dCdY.resize(mNumberOfDesignVariables);
    //         dCdZ.resize(mNumberOfDesignVariables);
    //         for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
    //         {
    //             int i = node_i->GetValue(MAPPING_ID);
    //             array_3d& node_sens = node_i->FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY);
    //             dCdX[i] = node_sens[0];
    //             dCdY[i] = node_sens[1];
    //             dCdZ[i] = node_sens[2];
    //         }
    //     }

    //     // Perform mapping of objective sensitivities
    //     VectorType dJdsx, dJdsy, dJdsz;
    //     dJdsx.resize(mNumberOfDesignVariables);
    //     dJdsy.resize(mNumberOfDesignVariables);
    //     dJdsz.resize(mNumberOfDesignVariables);
    //     SparseSpaceType::TransposeMult(mMappingMatrix,dJdX,dJdsx);
    //     SparseSpaceType::TransposeMult(mMappingMatrix,dJdY,dJdsy);
    //     SparseSpaceType::TransposeMult(mMappingMatrix,dJdZ,dJdsz);

    //     // Assign results to nodal variables
    //     for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
    //     {
    //         int i = node_i->GetValue(MAPPING_ID);
    //         VectorType dJds_i = ZeroVector(3);
    //         dJds_i(0) = dJdsx[i];
    //         dJds_i(1) = dJdsy[i];
    //         dJds_i(2) = dJdsz[i];
    //         node_i->FastGetSolutionStepValue(MAPPED_OBJECTIVE_SENSITIVITY) = dJds_i;
    //     }

    //     // Repeat mapping to map constraint sensitivities
    //     if(constraint_given)
    //     {
    //         // Perform mapping
    //         VectorType dCdsx, dCdsy, dCdsz;
    //         dCdsx.resize(mNumberOfDesignVariables);
    //         dCdsy.resize(mNumberOfDesignVariables);
    //         dCdsz.resize(mNumberOfDesignVariables);
    //         SparseSpaceType::TransposeMult(mMappingMatrix,dCdX,dCdsx);
    //         SparseSpaceType::TransposeMult(mMappingMatrix,dCdY,dCdsy);
    //         SparseSpaceType::TransposeMult(mMappingMatrix,dCdZ,dCdsz);

    //         // Assign results to nodal variables
    //         for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
    //         {
    //             int i = node_i->GetValue(MAPPING_ID);
    //             VectorType dCds_i = ZeroVector(3);
    //             dCds_i(0) = dCdsx[i];
    //             dCds_i(1) = dCdsy[i];
    //             dCds_i(2) = dCdsz[i];
    //             node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY) = dCds_i;
    //         }
    //     }

        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void MapToGeometrySpace( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInGeometrySpace )
    {
        boost::timer mapping_time;
        std::cout << "\n> Starting to map " << rNodalVariable.Name() << " to geometry space..." << std::endl;

        RecomputeMappingMatrixIfGeometryHasChanged();        

    //     // Assign design update to vector which shall be mapped (depending specified mapping matrix )
    //     VectorType dsx, dsy, dsz;
    //     dsx.resize(mNumberOfDesignVariables);
    //     dsy.resize(mNumberOfDesignVariables);
    //     dsz.resize(mNumberOfDesignVariables);
    //     for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
    //     {
    //         int i = node_i->GetValue(MAPPING_ID);
    //         dsx[i] = node_i->FastGetSolutionStepValue(DESIGN_UPDATE_X);
    //         dsy[i] = node_i->FastGetSolutionStepValue(DESIGN_UPDATE_Y);
    //         dsz[i] = node_i->FastGetSolutionStepValue(DESIGN_UPDATE_Z);
    //     }

    //     // Perform mapping to compute shape update
    //     VectorType dx, dy, dz;
    //     dx.resize(mNumberOfDesignVariables);
    //     dy.resize(mNumberOfDesignVariables);
    //     dz.resize(mNumberOfDesignVariables);
    //     noalias(dx) = prod(mMappingMatrix,dsx);
    //     noalias(dy) = prod(mMappingMatrix,dsy);
    //     noalias(dz) = prod(mMappingMatrix,dsz);

    //     // Assign dx as nodal shape updates and update coordinates accordingly
    //     for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
    //     {
    //         // If shape update deactivated, set it to zero
    //         if(node_i->FastGetSolutionStepValue(SHAPE_UPDATES_DEACTIVATED))
    //         {
    //             array_3d zero_array(3,0.0);
    //             noalias(node_i->FastGetSolutionStepValue(SHAPE_UPDATE)) = zero_array;
    //         }
    //         // In case it is not deactivated, it is checked if it is on a specified boundary beyond which no update is wanted
    //         else
    //         {
    //             // Read shape update from solution vector
    //             int i = node_i->GetValue(MAPPING_ID);
    //             array_3d shape_update;
    //             shape_update[0] = dx[i];
    //             shape_update[1] = dy[i];
    //             shape_update[2] = dz[i];

    //             // If node is on specified boundary, project shape update on specified boundary plane
    //             // I.e. remove component that is normal to the boundary plane
    //             if(node_i->FastGetSolutionStepValue(IS_ON_BOUNDARY))
    //             {
    //                 array_3d boundary_plane = node_i->FastGetSolutionStepValue(BOUNDARY_PLANE);
    //                 shape_update = shape_update - (inner_prod(shape_update,boundary_plane))*boundary_plane/norm_2(boundary_plane);
    //             }

    //             // Assign shape update to nodal variable
    //             noalias(node_i->FastGetSolutionStepValue(SHAPE_UPDATE)) = shape_update;
    //         }
    //     }

    //     // Apply edge damping if specified
    //     if(mPerformDamping)
    //         for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
    //         {   
    //             node_i->FastGetSolutionStepValue(SHAPE_UPDATE_X) *= node_i->GetValue(DAMPING_FACTOR_X);
    //             node_i->FastGetSolutionStepValue(SHAPE_UPDATE_Y) *= node_i->GetValue(DAMPING_FACTOR_Y);
    //             node_i->FastGetSolutionStepValue(SHAPE_UPDATE_Z) *= node_i->GetValue(DAMPING_FACTOR_Z);
    //         }

    //     // Update model part and record absolute shape change 
    //     for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
    //     {
    //         array_3d& shape_update = node_i->FastGetSolutionStepValue(SHAPE_UPDATE);
                
    //         // Update coordinates
    //         node_i->X() += shape_update[0];
    //         node_i->Y() += shape_update[1];
    //         node_i->Z() += shape_update[2];

    //         // Add shape update to previous updates
    //         noalias(node_i->FastGetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE)) += shape_update;
    //     }

        std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
    }

    // --------------------------------------------------------------------------
    void RecomputeMappingMatrixIfGeometryHasChanged()
    {
        if(HasGeometryChanged())
        {
            ResetComputationOfMappingMatrix();
            ComputeMappingMatrix();
        }
    }

    // --------------------------------------------------------------------------
    bool HasGeometryChanged()
    {
        return true;
    }   

    // --------------------------------------------------------------------------
    void ResetComputationOfMappingMatrix()
    {
        mpSearchTree.reset();
        mMappingMatrix.clear();
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
        return "MapperVertexMorphing";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MapperVertexMorphing";
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

    // // ==============================================================================
    // // Initialized by class constructor
    // // ==============================================================================
    ModelPart& mrDesignSurface;
    const unsigned int mNumberOfDesignVariables;
    double mFilterRadius;
    std::string mFilterType;
    FilterFunction::Pointer mpFilterFunction;

    // ==============================================================================
    // Variables for spatial search
    // ==============================================================================
    unsigned int mBucketSize = 100;
    unsigned int mMaxNeighborNodes = 10000;            
    NodeVector mListOfNodesOfDesignSurface;
    KDTree::Pointer mpSearchTree;

    // // ==============================================================================
    // // General working arrays
    // // ==============================================================================
    SparseMatrixType mMappingMatrix;

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
//      MapperVertexMorphing& operator=(MapperVertexMorphing const& rOther);

    /// Copy constructor.
//      MapperVertexMorphing(MapperVertexMorphing const& rOther);


    ///@}

}; // Class MapperVertexMorphing

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_H
