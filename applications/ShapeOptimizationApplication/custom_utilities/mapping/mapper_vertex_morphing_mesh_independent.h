// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

#ifndef MAPPER_VERTEX_MORPHING_MESH_INDEPENDENT_H
#define MAPPER_VERTEX_MORPHING_MESH_INDEPENDENT_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "mapper_vertex_morphing.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "linear_solvers/amgcl_solver.h"
#include "processes/find_conditions_neighbours_process.h"
#include "utilities/math_utils.h"
#include "custom_utilities/optimization_utilities.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
using ReordererType = Reorderer<SparseSpaceType, LocalSpaceType>;
using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
using SkylineLUFactorizationSolverType = SkylineLUFactorizationSolver<SparseSpaceType, LocalSpaceType, ReordererType>;
using AMGCLSolverType = AMGCLSolver<SparseSpaceType, LocalSpaceType, ReordererType>;

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

class MapperVertexMorphingMeshIndependent : public MapperVertexMorphing
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperVertexMorphingMeshIndependent
    KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphingMeshIndependent);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperVertexMorphingMeshIndependent( ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, Parameters MapperSettings )
        : MapperVertexMorphing(rOriginModelPart, rDestinationModelPart, MapperSettings)
    {
        mWeighting = MapperSettings["mesh_independent_weighting"].GetBool();
        mWeightingFullMassOrigin = MapperSettings["mesh_independent_weighting_full_mass_origin"].GetBool();
        mWeightingFullMassDestination = MapperSettings["mesh_independent_weighting_full_mass_destination"].GetBool();
        if (mWeightingFullMassOrigin || mWeightingFullMassDestination) {
            mWeightingFullMass = true;
        }
    }

    /// Destructor.
    virtual ~MapperVertexMorphingMeshIndependent()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        if (mIsMappingInitialized == false)
        {
            FindNeighbourConditions(mrOriginModelPart);
            FindNeighbourConditions(mrDestinationModelPart);
        }

        MapperVertexMorphing::Initialize();
    }

    void ComputeMappingMatrix() override
    {
        MapperVertexMorphing::ComputeMappingMatrix();

        KRATOS_WATCH("diagonal")

        // compute the scaling for the mesh indpendency
        Vector diagonal_mass_matrix(mrOriginModelPart.Nodes().size());
        SparseSpaceType::TransposeMult(mMappingMatrix, mDestinationNodalAreas, diagonal_mass_matrix);

        for(auto& node_i : mrOriginModelPart.Nodes())
        {
            const int mapping_id = node_i.GetValue(MAPPING_ID);
            node_i.SetValue(VARIABLE_SCALING_FACTOR, std::sqrt(1.0/diagonal_mass_matrix[mapping_id]));
        }
    }

    void Map( const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        if (!mWeighting)
        {
            ScaleOriginValues(rOriginVariable);
        }
        MapperVertexMorphing::Map(rOriginVariable, rDestinationVariable);
    }

    void Map( const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        if (!mWeighting)
        {
            ScaleOriginValues(rOriginVariable);
        }
        MapperVertexMorphing::Map(rOriginVariable, rDestinationVariable);
    }

    void InverseMap( const Variable<array_3d> &rDestinationVariable, const Variable<array_3d> &rOriginVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        if (mWeightingFullMassDestination) {
            WeightDestinationValues(rDestinationVariable);
        }

        MapperVertexMorphing::InverseMap(rDestinationVariable, rOriginVariable);

        if (mWeightingFullMassOrigin) {
            WeightOriginValues(rOriginVariable);
        }

        if (!mWeightingFullMass) {
            ScaleOriginValues(rOriginVariable);
            if (mWeighting)
            {
                ScaleOriginValues(rOriginVariable);
                KRATOS_WARNING("ShapeOpt::MapperVertexMorphing") << "Sensitivity Weighting!" << std::endl;
            }
        }

    }

    void InverseMap(const Variable<double> &rDestinationVariable, const Variable<double> &rOriginVariable) override
    {
        if (mIsMappingInitialized == false)
            Initialize();

        // if (mWeightingFullMassDestination) {
        //     WeightDestinationValues(rDestinationVariable);
        // }

        MapperVertexMorphing::InverseMap(rDestinationVariable, rOriginVariable);

        // if (mWeightingFullMassOrigin) {
        //     WeightOriginValues(rOriginVariable);
        // }

        if (!mWeightingFullMass) {
            ScaleOriginValues(rOriginVariable);
            if (mWeighting)
            {
                ScaleOriginValues(rOriginVariable);
                KRATOS_WARNING("ShapeOpt::MapperVertexMorphing") << "Sensitivity Weighting!" << std::endl;
            }
        }

    }


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
    std::string Info() const override
    {
        return "MapperVertexMorphingMeshIndependent";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MapperVertexMorphingMeshIndependent";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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
    void ComputeWeightForAllNeighbors(  ModelPart::NodeType& node_i,
                                        NodeVector& neighbor_nodes,
                                        unsigned int number_of_neighbors,
                                        std::vector<double>& list_of_weights,
                                        double& sum_of_weights ) override
    {
        for(unsigned int j_itr = 0 ; j_itr<number_of_neighbors ; j_itr++)
        {
            // Get node information
            ModelPart::NodeType& node_j = *neighbor_nodes[j_itr];

            // Get all neighbour conditions
            double Aij = mpFilterFunction->ComputeWeight(node_j.Coordinates(),node_i.Coordinates(), GetVertexMorphingRadius(node_i));
            Aij *= mOriginNodalAreas[node_j.GetValue(MAPPING_ID)];

            // Add values to list
            list_of_weights[j_itr] += Aij;

            // Computed for integration of weighting function later using post-scaling
            sum_of_weights += Aij;
        }
    }

    void InitializeComputationOfMappingMatrix() override
    {
        // from base class
        MapperVertexMorphing::InitializeComputationOfMappingMatrix();

        ComputeNodalAreas(mrOriginModelPart, mOriginNodalAreas);
        ComputeNodalAreas(mrDestinationModelPart, mDestinationNodalAreas);

        const unsigned int origin_node_number = mrOriginModelPart.Nodes().size();
        const unsigned int destination_node_number = mrDestinationModelPart.Nodes().size();

        if (mWeightingFullMass) {
            KRATOS_WARNING("ShapeOpt::MapperVertexMorphingMeshIndependent") << "Computing full mass matrix!" << std::endl;
            mOriginMassMatrix.resize(destination_node_number,origin_node_number,false);
            mOriginMassMatrix.clear();

            mDestinationMassMatrix.resize(destination_node_number,origin_node_number,false);
            mDestinationMassMatrix.clear();

            ComputeMassMatrix(mrOriginModelPart, mOriginMassMatrix);
            ComputeMassMatrix(mrDestinationModelPart, mDestinationMassMatrix);

            for (unsigned int i = 0; i < 10; i++) {
                for (unsigned int j = 0; j < 10; j++) {
                    std::cout << "mOriginMassMatrix(" << i << "," << j << "): " << mOriginMassMatrix(i,j) << " " << std::endl;
                }
            }
            for (unsigned int i = 0; i < 10; i++) {
                for (unsigned int j = 0; j < 10; j++) {
                    std::cout << "mDestinationMassMatrix(" << i << "," << j << "): " << mOriginMassMatrix(i,j) << " " << std::endl;
                }
            }
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

    Vector mOriginNodalAreas;
    Vector mDestinationNodalAreas;
    bool mWeighting = false;

    SparseMatrixType mOriginMassMatrix;
    SparseMatrixType mDestinationMassMatrix;
    bool mWeightingFullMass = false;
    bool mWeightingFullMassOrigin = false;
    bool mWeightingFullMassDestination = false;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void FindNeighbourConditions(ModelPart& rModelPart)
    {
        FindConditionsNeighboursProcess find_conditions_neighbours_process(rModelPart, rModelPart.GetProcessInfo()[DOMAIN_SIZE]);
        find_conditions_neighbours_process.Execute();
    }

    void ComputeNodalAreas(ModelPart& rModelPart, Vector& rNodalAreas){
        rNodalAreas.resize(rModelPart.Nodes().size(), false);
        noalias(rNodalAreas) = ZeroVector(rNodalAreas.size());
        for(auto& node_i : rModelPart.Nodes())
        {
            const int& i = node_i.GetValue(MAPPING_ID);

            // Get all neighbour conditions
            const GlobalPointersVector<Condition>& rConditions = node_i.GetValue(NEIGHBOUR_CONDITIONS);

            // loop conditions
            for(unsigned int c_itr=0; c_itr<rConditions.size(); c_itr++)
            {
                // Get geometry of current condition
                Condition rCondition = rConditions[c_itr];
                Condition::GeometryType& geom_i = rCondition.GetGeometry();
                rNodalAreas[i] += geom_i.DomainSize() / geom_i.size();
            }
        }
    }

    void ComputeMassMatrix(ModelPart& rModelPart, SparseMatrixType& rMassMatrix)
    {
        for (const auto& rCondition : rModelPart.Conditions())
        {
            const auto& rGeometry = rCondition.GetGeometry();
            const double element_area = rGeometry.DomainSize();

            for (const auto& rNodeI : rGeometry)
            {
                const int i = rNodeI.GetValue(MAPPING_ID);

                for (const auto& rNodeJ : rGeometry)
                {
                    const int j = rNodeJ.GetValue(MAPPING_ID);

                    if (j < i)
                        continue; // only compute upper triangle, the matrix is symmetric

                    if (i == j)
                    {
                        rMassMatrix(i, j) += (2.0 / 12.0) * element_area;
                    }
                    else
                    {
                        rMassMatrix(i, j) += (1.0 / 12.0) * element_area;
                        rMassMatrix(j, i) += (1.0 / 12.0) * element_area;
                    }
                }
            }
        }
    }


    template <typename T>
    void ScaleOriginValues(const T& rOriginVariable) {
        for (auto& node_i : mrOriginModelPart.Nodes())
        {
            node_i.FastGetSolutionStepValue(rOriginVariable) *= node_i.GetValue(VARIABLE_SCALING_FACTOR);
        }
    }

    // template <typename T>
    void WeightOriginValues(const Variable<array_3d>& rOriginVariable) {

        Parameters empty_parameters =  Parameters(R"({})");
        LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new AMGCLSolverType(empty_parameters) );

        Vector raw_sens = ZeroVector(mrOriginModelPart.Nodes().size()*3);
        Vector weighted_sens = ZeroVector(mrOriginModelPart.Nodes().size()*3);
        OptimizationUtilities::AssembleVector(mrOriginModelPart, raw_sens, rOriginVariable);

        for (int dim = 0; dim < 3; dim++) {

            Vector raw_sens_dim = ZeroVector(mrOriginModelPart.Nodes().size());
            for (int i = 0; i < mrOriginModelPart.Nodes().size(); i++) {
                raw_sens_dim[i] = raw_sens[i*3 + dim];
                if (i < 10) {
                    std::cout << "raw_sens_dim[" << i << "]: " << raw_sens_dim[i] << " " << std::endl;
                }
            }

            Vector weighted_sens_dim = ZeroVector(mrOriginModelPart.Nodes().size());
            psolver->Solve(mOriginMassMatrix, weighted_sens_dim, raw_sens_dim);

            for (int i = 0; i < mrOriginModelPart.Nodes().size(); i++) {
                weighted_sens[i*3 + dim] = weighted_sens_dim[i];
                if (i < 10) {
                    std::cout << "weighted_sens_dim[" << i << "]: " << weighted_sens_dim[i] << " " << std::endl;
                }
            }
        }
        OptimizationUtilities::AssignVectorToVariable(mrOriginModelPart, weighted_sens, rOriginVariable);
    }

    // template <typename T>
    void WeightDestinationValues(const Variable<array_3d>& rDestinationVariable) {
        Parameters empty_parameters =  Parameters(R"({})");
        LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new AMGCLSolverType(empty_parameters) );

        Vector raw_sens = ZeroVector(mrDestinationModelPart.Nodes().size()*3);
        Vector weighted_sens = ZeroVector(mrDestinationModelPart.Nodes().size()*3);
        OptimizationUtilities::AssembleVector(mrDestinationModelPart, raw_sens, rDestinationVariable);

        for (int dim = 0; dim < 3; dim++) {

            Vector raw_sens_dim = ZeroVector(mrDestinationModelPart.Nodes().size());
            for (int i = 0; i < mrDestinationModelPart.Nodes().size(); i++) {
                raw_sens_dim[i] = raw_sens[i*3 + dim];
                if (i < 10) {
                    std::cout << "raw_sens_dim[" << i << "]: " << raw_sens_dim[i] << " " << std::endl;
                }
            }

            Vector weighted_sens_dim = ZeroVector(mrDestinationModelPart.Nodes().size());
            psolver->Solve(mDestinationMassMatrix, weighted_sens_dim, raw_sens_dim);

            for (int i = 0; i < mrDestinationModelPart.Nodes().size(); i++) {
                weighted_sens[i*3 + dim] = weighted_sens_dim[i];
                if (i < 10) {
                    std::cout << "weighted_sens_dim[" << i << "]: " << weighted_sens_dim[i] << " " << std::endl;
                }
            }
        }
        OptimizationUtilities::AssignVectorToVariable(mrDestinationModelPart, weighted_sens, rDestinationVariable);
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

    ///@}

}; // Class MapperVertexMorphingMeshIndependent

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MAPPER_VERTEX_MORPHING_MESH_INDEPENDENT_H
