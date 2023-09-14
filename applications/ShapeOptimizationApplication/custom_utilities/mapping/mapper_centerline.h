// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:
//
// ==============================================================================

#ifndef MAPPER_CENTERLINE_H
#define MAPPER_CENTERLINE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "mapper_base.h"
#include "custom_utilities/filter_function.h"
#include "custom_utilities/heat_method_utilities.h"
#include "linear_solvers/linear_solver.h"

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

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) MapperCenterline : public Mapper
{
public:
    ///@name Type Definitions
    ///@{

    // Type definitions for better reading later
    typedef Node NodeType;
    typedef Node::Pointer NodeTypePointer;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double>::iterator DoubleVectorIterator;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;


    // Type definitions for linear algebra including sparse systems
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef SparseSpaceType::MatrixType SparseMatrixType;
    typedef SparseSpaceType::VectorType VectorType;
    typedef UblasSpace<double, Matrix, Vector> DenseSpace;
    typedef UblasSpace<double, SparseMatrix, Vector> SparseSpace;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// Pointer definition of MapperCenterline
    KRATOS_CLASS_POINTER_DEFINITION(MapperCenterline);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperCenterline( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters MapperSettings)
        : mrOriginModelPart( rOriginModelPart ),
          mrDestinationModelPart( rDestinationModelPart ),
          mMapperSettings( MapperSettings ),
          mFilterRadius( MapperSettings["filter_radius"].GetDouble() ),
          mMaxNumberOfNeighbors( MapperSettings["max_nodes_in_filter_radius"].GetInt())
    {
    }

    /// Destructor.
    virtual ~MapperCenterline()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void Initialize() override;

    // --------------------------------------------------------------------------
    void Map( const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable ) override;

    // --------------------------------------------------------------------------
    void Map( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable ) override;

    // --------------------------------------------------------------------------
    void InverseMap( const Variable<array_3d> &rDestinationVariable, const Variable<array_3d> &rOriginVariable ) override;

    // --------------------------------------------------------------------------
    void InverseMap( const Variable<double> &rDestinationVariable, const Variable<double> &rOriginVariable ) override;

    // --------------------------------------------------------------------------
    void Update(LinearSolver<SparseSpaceType, LocalSpaceType>& rSolver, ModelPart& rOriginModelPart);

    // --------------------------------------------------------------------------

    ///@}
    ///@name Access
    ///@{

    FilterFunction::UniquePointer mpFilterFunction;

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "MapperCenterline";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MapperCenterline";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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

    // Initialized by class constructor
    ModelPart& mrOriginModelPart;
    ModelPart& mrDestinationModelPart;
    Parameters mMapperSettings;
    double mFilterRadius;
    unsigned int mMaxNumberOfNeighbors;

    // Variables for spatial search
    unsigned int mBucketSize = 100;
    NodeVector mListOfNodesInOriginModelPart;
    KDTree::Pointer mpSearchTree;

    // Variables for mapping
    SparseMatrixType mMappingMatrix;
    std::vector<Vector> mValuesOrigin;
    std::vector<Vector> mValuesDestination;
    bool mIsMappingInitialized = false;

    // Variables for heat method
    HeatMethodUtilities::UniquePointer mpHeatMethod;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // --------------------------------------------------------------------------
    void CreateListOfNodesInOriginModelPart();

    // --------------------------------------------------------------------------
    void CreateFilterFunction();

    // --------------------------------------------------------------------------
    void InitializeMappingVariables();

    // --------------------------------------------------------------------------
    void AssignMappingIds();
    // --------------------------------------------------------------------------
    void CreateSearchTreeWithAllNodesInOriginModelPart();

    // --------------------------------------------------------------------------
    void ThrowWarningIfNumberOfNeighborsExceedsLimit(ModelPart::NodeType& given_node, unsigned int number_of_neighbors);

    // --------------------------------------------------------------------------
    virtual void ComputeWeightForAllNeighbors(  const ModelPart::NodeType& design_node,
                                        const NodeVector& neighbor_nodes,
                                        const unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights );

    // --------------------------------------------------------------------------
    double GetVertexMorphingRadius(const NodeType& rNode) const override
    {
        return mFilterRadius;
    }

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
//      MapperCenterline& operator=(MapperCenterline const& rOther);

    /// Copy constructor.
//      MapperCenterline(MapperCenterline const& rOther);


    ///@}

}; // Class MapperCenterline

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_CENTERLINE_H
