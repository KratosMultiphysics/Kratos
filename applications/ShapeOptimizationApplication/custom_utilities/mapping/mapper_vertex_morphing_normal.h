// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_NORMAL_H
#define MAPPER_VERTEX_MORPHING_NORMAL_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "mapper_vertex_morphing.h"
#include "custom_utilities/geometry_utilities.h"

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

class MapperVertexMorphingNormal : public MapperVertexMorphing
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperVertexMorphingNormal
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingNormal);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingNormal( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters MapperSettings )
        : MapperVertexMorphing(rOriginModelPart, rDestinationModelPart, MapperSettings)
    {
    }

    /// Destructor.
    virtual ~MapperVertexMorphingNormal()
    {
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void Map( const Variable<double> &rOriginVariable, const Variable<array_3d> &rDestinationVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        BuiltinTimer timer;
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("ShapeOpt") << "Starting mapping of " << rOriginVariable.Name() << "..." << std::endl;

        // Prepare vectors for mapping
        mOriginVector.clear();
        mDestinationVector.clear();

        for(auto& node_i : mrOriginModelPart.Nodes())
        {
            const int i = node_i.GetValue(MAPPING_ID);
            mOriginVector[i] = node_i.FastGetSolutionStepValue(rOriginVariable);
        }

        // Perform mapping
        noalias(mDestinationVector) = prod(mMappingMatrix, mOriginVector);

        // Assign results to nodal variable
        for(auto& node_i : mrDestinationModelPart.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);

            array_3d& r_node_vector = node_i.FastGetSolutionStepValue(rDestinationVariable);
            r_node_vector(0) = mDestinationVector[i*3+0];
            r_node_vector(1) = mDestinationVector[i*3+1];
            r_node_vector(2) = mDestinationVector[i*3+2];
        }

        KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void Map( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable) override
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    // --------------------------------------------------------------------------
    void Map( const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable) override
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    // --------------------------------------------------------------------------
    void InverseMap( const Variable<array_3d> &rDestinationVariable, const Variable<double> &rOriginVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        BuiltinTimer timer;
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("ShapeOpt") << "Starting inverse mapping of " << rDestinationVariable.Name() << "..." << std::endl;

        // Prepare vectors for mapping
        mOriginVector.clear();
        mDestinationVector.clear();

        for(auto& node_i : mrDestinationModelPart.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);
            array_3d& r_nodal_variable = node_i.FastGetSolutionStepValue(rDestinationVariable);
            mDestinationVector[i*3+0] = r_nodal_variable[0];
            mDestinationVector[i*3+1] = r_nodal_variable[1];
            mDestinationVector[i*3+2] = r_nodal_variable[2];
        }

        // Perform mapping
        if(mMapperSettings["consistent_mapping"].GetBool())
        {
            KRATOS_ERROR << "Not Implemented";
        }
        else
        {
            SparseSpaceType::TransposeMult(mMappingMatrix,mDestinationVector,mOriginVector);
        }

        // Assign results to nodal variable
        for(auto& node_i : mrOriginModelPart.Nodes())
        {
            int i = node_i.GetValue(MAPPING_ID);
            node_i.FastGetSolutionStepValue(rOriginVariable) = mOriginVector[i];
        }

        KRATOS_INFO("ShapeOpt") << "Finished mapping in " << timer.ElapsedSeconds() << " s." << std::endl;
    }

    // --------------------------------------------------------------------------
    void InverseMap(const Variable<array_3d> &rDestinationVariable, const Variable<array_3d> &rOriginVariable) override
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    // --------------------------------------------------------------------------
    void InverseMap(const Variable<double> &rDestinationVariable, const Variable<double> &rOriginVariable) override
    {
        KRATOS_ERROR << "Not implemented!" << std::endl;
    }

    // --------------------------------------------------------------------------

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
    virtual std::string Info() const override
    {
        return "MapperVertexMorphingNormal";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MapperVertexMorphingNormal";
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

    // --------------------------------------------------------------------------
    virtual void ComputeMappingMatrix()
    {
        GeometryUtilities(mrOriginModelPart).ComputeUnitSurfaceNormals();
        MapperVertexMorphing::ComputeMappingMatrix();
    }

    // --------------------------------------------------------------------------
    void FillMappingMatrixWithWeights(  ModelPart::NodeType& destination_node,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights )
    {
        unsigned int row_id = destination_node.GetValue(MAPPING_ID);
        for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            int collumn_id = neighbor_node.GetValue(MAPPING_ID);
            const array_3d& normal = neighbor_node.FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);

            array_3d weight = normal * list_of_weights[neighbor_itr] / sum_of_weights;
            mMappingMatrix.insert_element(row_id*3+0,collumn_id,weight[0]);
            mMappingMatrix.insert_element(row_id*3+1,collumn_id,weight[1]);
            mMappingMatrix.insert_element(row_id*3+2,collumn_id,weight[2]);
        }
    }

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

    // Variables for mapping
    Vector mOriginVector;
    Vector mDestinationVector;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // --------------------------------------------------------------------------
    void InitializeMappingVariables()
    {
        const unsigned int origin_node_number = mrOriginModelPart.Nodes().size();
        const unsigned int destination_node_number = mrDestinationModelPart.Nodes().size();

        mMappingMatrix.resize(destination_node_number*3,origin_node_number,false);
        mMappingMatrix.clear();

        mOriginVector = ZeroVector(origin_node_number);
        mDestinationVector = ZeroVector(destination_node_number*3);
    }

    // --------------------------------------------------------------------------

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
//      MapperVertexMorphingNormal& operator=(MapperVertexMorphingNormal const& rOther);

    /// Copy constructor.
//      MapperVertexMorphingNormal(MapperVertexMorphingNormal const& rOther);


    ///@}

}; // Class MapperVertexMorphingNormal

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_NORMAL_H
